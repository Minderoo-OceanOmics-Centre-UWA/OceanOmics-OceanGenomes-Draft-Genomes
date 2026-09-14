#!/usr/bin/env python3
import argparse
import configparser
import csv
import glob
import os
import re
import sys

import psycopg2

# Sibling import: Nextflow bind-mounts bin/ into the task container and puts it
# on PATH, so this resolves via sys.path[0]. Kept byte-identical to the copy in
# the mitogenome pipeline.
from taxdump_lineage import TaxdumpLineage

def load_db_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)

    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port')
    }

def get_species_info(db_params, og_id):
    """
    Return (nominal_species_id, taxon_id, tax_class) for an og_id,
    or None if not found.
    """
    query = """
    WITH sample_q AS (
        SELECT
            s.og_id,
            s.nominal_species_id,
            trim(s.nominal_species_id) AS nominal_name,
            split_part(trim(s.nominal_species_id), ' ', 1) AS nominal_genus
        FROM sample s
        WHERE s.og_id = %s
    )
    SELECT
        s.og_id,
        s.nominal_species_id,
        m.ncbi_taxon_id,
        m.class,
        m.species_matched,
        m.sim,
        m.match_level
    FROM sample_q s
    LEFT JOIN LATERAL (
        SELECT *
        FROM (
            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                1 AS priority,
                1.0 AS sim,
                'species_exact' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND lower(sp.species) = lower(s.nominal_name)

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                2 AS priority,
                1.0 AS sim,
                'genus_exact' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND lower(sp.genus) = lower(s.nominal_genus)

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                3 AS priority,
                1.0 AS sim,
                'family_exact' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND lower(sp.family) = lower(s.nominal_name)

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                4 AS priority,
                1.0 AS sim,
                'order_exact' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND lower(sp.ordr) = lower(s.nominal_name)

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                5 AS priority,
                similarity(sp.species, s.nominal_name) AS sim,
                'species_fuzzy_same_genus' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND lower(sp.genus) = lower(s.nominal_genus)
              AND sp.species %% s.nominal_name

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                6 AS priority,
                similarity(sp.family, s.nominal_name) AS sim,
                'family_fuzzy' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND sp.family %% s.nominal_name
              AND similarity(sp.family, s.nominal_name) >= 0.65

            UNION ALL

            SELECT
                sp.ncbi_taxon_id,
                sp.class,
                sp.species AS species_matched,
                7 AS priority,
                similarity(sp.ordr, s.nominal_name) AS sim,
                'order_fuzzy' AS match_level
            FROM species sp
            WHERE sp.ncbi_taxon_id IS NOT NULL
              AND sp.ordr %% s.nominal_name
              AND similarity(sp.ordr, s.nominal_name) >= 0.65
        ) ranked_matches
        ORDER BY
            priority,
            sim DESC,
            ncbi_taxon_id
        LIMIT 1
    ) m ON TRUE
    ORDER BY s.og_id;
    """
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            # IMPORTANT: pass a tuple, not a bare string
            cur.execute(query, (og_id,))
            result = cur.fetchone()
            if result:
                (
                    ogid,
                    nominal_species_id,
                    taxon_id,
                    tax_class,
                    species_matched,
                    sim,
                    match_level
                ) = result

                # Deliberately quiet on a miss: the caller falls back to the
                # NCBI taxdump, so warning here would report a failure that
                # mostly does not happen. Per-sample provenance goes to the
                # resolution report, and a real miss is reported at the end.
                return (nominal_species_id, taxon_id, tax_class)
            else:
                print(f"[WARN] No sample row found for OG ID: {og_id}")
                return None
    except Exception as e:
        print(f"[ERROR] Database query failed for {og_id}: {e}")
        return None
    finally:
        if conn:
            conn.close()

def parse_run_id(run_id: str):
    """
    Assumes format like NOVA_251031_TWAD.
    Returns {'date': '251031'} (or '' if not present).
    """
    parts = run_id.split('_')
    return {"date": parts[1] if len(parts) > 1 else ""}

def detect_read_role(filename: str):
    """
    Try to detect whether this is R1 or R2 from the filename.
    Returns 'R1', 'R2', or None if unknown.
    """
    base = os.path.basename(filename)

    # Common patterns like sample.R1.fq.gz or sample.R2.fastq.gz
    if ".R1." in base or "_R1." in base or base.endswith("_R1") or base.endswith(".R1"):
        return "R1"
    if ".R2." in base or "_R2." in base or base.endswith("_R2") or base.endswith(".R2"):
        return "R2"

    # Fallback: look for 'R1' or 'R2' as separate tokens
    tokens = base.replace('.', '_').split('_')
    if "R1" in tokens:
        return "R1"
    if "R2" in tokens:
        return "R2"

    return None

def discover_pairs_from_dir(fastq_dir: str):
    """
    Scan a directory for FASTQ files and return a list of
    (og_id, [fastq_1, fastq_2]) tuples.

    og_id is taken as the first part of the filename before the first '.'.
    """
    patterns = [
        os.path.join(fastq_dir, "*.fastq.gz"),
        os.path.join(fastq_dir, "*.fq.gz"),
        os.path.join(fastq_dir, "*.fastq"),
        os.path.join(fastq_dir, "*.fq"),
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat))

    if not files:
        print(f"[ERROR] No FASTQ files found in directory: {fastq_dir}")
        sys.exit(1)

    # Map og_id -> {'R1': path, 'R2': path}
    by_og = {}

    for path in sorted(files):
        base = os.path.basename(path)
        og_id = base.split('.')[0]  # first part before first '.'

        read_role = detect_read_role(base)
        if read_role not in ("R1", "R2"):
            print(f"[WARN] Could not determine R1/R2 for file: {base}; skipping.")
            continue

        if og_id not in by_og:
            by_og[og_id] = {}

        if read_role in by_og[og_id]:
            print(
                f"[WARN] Duplicate {read_role} for OG {og_id}: {base}; "
                f"existing: {os.path.basename(by_og[og_id][read_role])}"
            )
        # Keep the first one we saw
        by_og[og_id].setdefault(read_role, path)

    pairs = []
    for og_id, reads in by_og.items():
        r1 = reads.get("R1")
        r2 = reads.get("R2")
        if not r1 or not r2:
            print(
                f"[WARN] Missing pair for OG {og_id}: "
                f"R1={bool(r1)}, R2={bool(r2)}; skipping."
            )
            continue
        pairs.append((og_id, [r1, r2]))

    if not pairs:
        print("[ERROR] No complete (R1,R2) pairs found in directory.")
        sys.exit(1)

    return pairs

def create_samplesheet(rows, output_file):
    """
    Write a multi-row samplesheet CSV.

    Each element in rows is a dict with keys:
      'sample','run','date','prefix','nom_species_id','taxon_id','class','fastq_1','fastq_2'
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Header
        writer.writerow([
            'sample',
            'run',
            'date',
            'prefix',
            'nom_species_id',
            'taxon_id',
            'class',
            'fastq_1',
            'fastq_2'
        ])

        # Rows
        for row in rows:
            writer.writerow([
                row['sample'],
                row['run'],
                row['date'],
                row['prefix'],
                row['nom_species_id'],
                row['taxon_id'],
                row['class'],
                row['fastq_1'],
                row['fastq_2']
            ])


def clean_species_name(name):
    """
    Strip the qualifier noise that makes a nominal_species_id unusable as an
    NCBI taxonomy query (e.g. 'Alepes vari (TBC)', 'Astronesthes spp.',
    "Nesogobius sp. `groove cheek`"). Returns the best plain query string we can
    salvage: 'Genus species' if a binomial survives, otherwise the bare genus.
    Used only as a fallback when the species table yields no match.

    Kept in step with the copy in the mitogenome pipeline's create_samplesheet.py.
    """
    if not name:
        return ""
    s = str(name)
    s = re.sub(r"\([^)]*\)", " ", s)        # drop parentheticals: (TBC), (Class)
    s = re.sub(r"`[^`]*`", " ", s)           # drop backtick descriptors: `groove cheek`
    s = re.sub(r"[`'\"]", " ", s)            # stray quotes
    # drop open-nomenclature qualifiers and undescribed-species markers
    s = re.sub(r"\b(spp?|cf|aff|nr|sp|TBC)\.?\b", " ", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+", " ", s).strip()
    # Strip leading/trailing punctuation from each token and drop any token that
    # has no letters. This removes the stray '.' left behind by e.g. 'spp.'
    # ('Blachea spp.' -> 'Blachea', not 'Blachea .'), which NCBI rejects.
    tokens = [re.sub(r"^[^A-Za-z]+|[^A-Za-z]+$", "", t) for t in s.split()]
    tokens = [t for t in tokens if t]
    if not tokens:
        return ""
    # keep at most a binomial (Genus species); a lone genus is a valid query too
    return " ".join(tokens[:2])


def build_resolver(taxdump_dir):
    """
    One TaxdumpLineage for the whole run, or None when no usable dump was given.

    Parsing nodes.dmp/names.dmp costs ~6s and ~1 GB, so it must not happen per
    sample. index_all_ranks is on because our nominal names sit at whatever rank
    the collection could identify -- 'Caridea' is an infraorder, which the
    default five-rank index does not hold at all. class_falls_back_to_phylum is
    on to match what scripts/taxonomy/load_taxonomy.py writes into
    `species.class` for lineages where NCBI has no class rank.
    """
    if not taxdump_dir:
        return None
    resolver = TaxdumpLineage(
        taxdump_dir,
        index_all_ranks=True,
        class_falls_back_to_phylum=True,
    )
    if not resolver.available:
        print(f"[WARN] --taxdump-dir '{taxdump_dir}' has no nodes.dmp/names.dmp; "
              f"taxonomy fallback disabled")
        return None
    return resolver


def resolve_species_info(db_params, og_id, resolver=None):
    """
    (nominal_species_id, taxon_id, tax_class, source) for a sample.

    The curated `species` table stays authoritative and is tried first. It only
    holds taxa someone has loaded, though, and the invertebrate runs draw from
    most of Metazoa, so a miss there is the normal case rather than the
    exception. Anything the table leaves blank is filled from the NCBI taxdump.

    `source` is 'db', 'db+taxdump', 'taxdump' or 'unresolved'.
    """
    species_info = get_species_info(db_params, og_id)
    if not species_info:
        return None, None, None, 'unresolved'

    nominal_species_id, taxon_id, tax_class = species_info

    db_had_taxon = taxon_id not in (None, "")
    db_had_class = bool(tax_class)
    if db_had_taxon and db_had_class:
        return nominal_species_id, taxon_id, tax_class, 'db'

    lineage = {}
    query_name = clean_species_name(nominal_species_id)
    if resolver is not None and query_name:
        try:
            lineage = resolver.lineage_for_name(query_name)
        except Exception as exc:
            print(f"[WARN] taxdump lookup failed for {og_id} "
                  f"('{query_name}'): {exc}")

    if lineage:
        if not db_had_taxon:
            # The taxid of whatever node matched, at whatever rank. FCS-GX takes
            # a --tax-id at any rank, so a class-rank id beats no id at all.
            taxon_id = lineage.get('matched_taxid') or taxon_id
        if not db_had_class:
            tax_class = lineage.get('class') or tax_class

    resolved = taxon_id not in (None, "") and bool(tax_class)
    if not resolved:
        source = 'unresolved'
    elif db_had_taxon or db_had_class:
        source = 'db+taxdump' if lineage else 'db'
    else:
        source = 'taxdump'

    return nominal_species_id, taxon_id, tax_class, source


# What an unresolved field is written as. Matches the `^(unknown|None)$` pattern
# assets/schema_input.json already allows for taxon_id, so these rows parse and
# reach validateTaxonomy() rather than failing schema validation first.
UNKNOWN = 'unknown'

RESOLUTION_COLUMNS = ('sample', 'nominal_species_id', 'taxon_id', 'class', 'source')


def write_resolution_report(path, rows):
    """
    Per-sample taxonomy provenance, written next to the samplesheet.

    Without it a taxdump-derived class is indistinguishable from a curated one,
    and there is no way to audit which samples leaned on the fallback.
    """
    if not path:
        return
    with open(path, 'w', newline='') as handle:
        # lineterminator: csv.writer defaults to CRLF, which leaves a stray \r on
        # the last column of every row and breaks awk/cut/grep on this file.
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(RESOLUTION_COLUMNS)
        for row in rows:
            writer.writerow([row.get(column, '') for column in RESOLUTION_COLUMNS])


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a draft-genomes samplesheet from a directory of FASTQ files.")
    parser.add_argument("config_file", help="Path to the postgres config file.")
    parser.add_argument("run_id", help="Sequencing run ID, e.g. NOVA_260909_LA.")
    parser.add_argument("fastq_dir", help="Directory of pooled FASTQ files.")
    parser.add_argument("--taxdump-dir", default=None,
                        help="NCBI taxdump directory (nodes.dmp/names.dmp). Used to "
                             "resolve taxon_id/class when the species table has no "
                             "match. Without it, unmatched samples fail the run.")
    parser.add_argument("--resolution-report", default=None,
                        help="Where to write the per-sample taxonomy provenance table "
                             "(default: taxonomy_resolution.tsv next to the samplesheet).")
    return parser.parse_args()


def main():
    args = parse_args()

    db_params = load_db_config(args.config_file)

    # Discover (og_id, [R1,R2]) pairs from directory
    pairs = discover_pairs_from_dir(args.fastq_dir)

    run_info = parse_run_id(args.run_id)
    date = run_info["date"]

    output_file = f"{args.run_id}_samplesheet.csv"
    resolution_report = args.resolution_report or os.path.join(
        os.path.dirname(os.path.abspath(output_file)), "taxonomy_resolution.tsv")

    resolver = build_resolver(args.taxdump_dir)
    if resolver is None:
        print("[WARN] No taxdump available: taxon_id/class come from the species "
              "table only, so any sample it does not carry will fail below.")

    rows = []
    resolution_rows = []
    unresolved = []
    for og_id, files in pairs:
        nominal_species_id, taxon_id, tax_class, source = resolve_species_info(
            db_params, og_id, resolver)

        resolution_rows.append({
            'sample': og_id,
            'nominal_species_id': nominal_species_id or '',
            'taxon_id': taxon_id if taxon_id is not None else '',
            'class': tax_class or '',
            'source': source,
        })

        missing = []
        if not nominal_species_id:
            missing.append("nominal_species_id")
        if taxon_id in (None, ""):
            missing.append("taxon_id")
        if not tax_class:
            missing.append("class")
        if missing:
            # Emit the row anyway, marked UNKNOWN, and keep going. The samplesheet
            # is only published if this task succeeds (publishDir does not run on a
            # failed task), and an operator cannot fix rows they cannot see.
            #
            # This does NOT put an unknown taxon back in front of FCS-GX, which is
            # what these placeholders used to cause: `--tax-id unknown` failed only
            # after MEGAHIT had burned the SUs. The run is stopped instead by
            # validateTaxonomy() in subworkflows/local/prepare_samplesheet, which
            # fires after this sheet is published and before any assembly task is
            # submitted. Keep the two in step: that guard is what makes writing
            # these rows safe.
            unresolved.append((og_id, nominal_species_id,
                               f"missing {', '.join(missing)}"))

        rows.append({
            'sample': og_id,
            'run': args.run_id,
            'date': date,
            'prefix': f"{og_id}.ilmn.{date}",
            'nom_species_id': nominal_species_id or UNKNOWN,
            'taxon_id': taxon_id if taxon_id not in (None, "") else UNKNOWN,
            'class': tax_class or UNKNOWN,
            'fastq_1': str(files[0]),
            'fastq_2': str(files[1]),
        })

    # Written even on the failure path: it is the only record of what the taxdump
    # did and did not rescue, which is what the operator needs to fix the rest.
    write_resolution_report(resolution_report, resolution_rows)
    print(f"[INFO] Wrote taxonomy provenance to: {resolution_report}")

    if not rows:
        # A different failure from unresolved taxonomy, and there is no sheet to
        # write: no FASTQ pair was usable at all.
        print("[ERROR] No rows to write (all entries skipped).")
        return 1

    create_samplesheet(rows, output_file)
    print(f"[INFO] Wrote samplesheet to: {output_file}")

    if unresolved:
        print(f"[ERROR] Taxonomy could not be resolved for {len(unresolved)} of "
              f"{len(rows)} samples. They are written to the samplesheet with "
              f"'{UNKNOWN}' in the taxon_id and class columns:")
        for og_id, nominal_species_id, reason in unresolved:
            print(f"[ERROR]   {og_id}: {reason} "
                  f"(nominal_species_id='{nominal_species_id or ''}')")
        print("[ERROR] The pipeline will stop before any assembly work rather than "
              "run these, so no compute is spent until they are fixed.")
        print("[FIX] These names were not found in the species table OR in the NCBI "
              "taxonomy, which usually means the nominal_species_id is misspelled "
              "(e.g. 'Actinaria' for 'Actiniaria') or is not a taxon name at all "
              "(e.g. 'Larval fish'). Fix it either way round: correct the "
              "nominal_species_id in the sample table and delete the published "
              "samplesheet so it is regenerated, or edit the taxon_id and class "
              "columns in the published samplesheet and re-run. If the name is right "
              "and simply uncurated, load it with scripts/taxonomy/load_taxonomy.py "
              "(e.g. --phylum Mollusca,Echinodermata).")

    taxdump_rescued = sum(1 for r in resolution_rows if r['source'] in ('taxdump', 'db+taxdump'))
    if taxdump_rescued:
        print(f"[INFO] {taxdump_rescued} of {len(resolution_rows)} samples had their "
              f"taxonomy completed from the NCBI taxdump rather than the species table.")


if __name__ == "__main__":
    sys.exit(main() or 0)
