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

# Tiers of a master_species match, best first. The SQL below returns the first
# that hits, and resolve_species_info() decides how far to trust it.
#
# Only an exact species or synonym hit names the sample's own taxon. Every tier
# after that returns SOME row of the right genus/family/order, so its class is
# right but its ncbi_taxon_id belongs to another species -- which is why the
# taxdump is asked before those tiers' taxid is used.
EXACT_SPECIES_LEVELS = ('species_exact', 'synonym_exact')
FUZZY_LEVELS = ('family_fuzzy', 'order_fuzzy')

# master_species, not species: it is the maintained table (OceanGenom_Web's), and
# only read here. It carries no invertebrates; those resolve from the taxdump.
MATCH_QUERY = """
SELECT
    m.ncbi_taxon_id,
    m.class,
    m.species_matched,
    m.match_level
FROM (
    SELECT sp.ncbi_taxon_id, sp.class, sp.species AS species_matched,
           1 AS priority, 1.0 AS sim, 'species_exact' AS match_level
    FROM master_species sp
    WHERE lower(sp.species) IN (lower(%(raw)s), lower(%(clean)s))

    UNION ALL

    -- synonym is a single name, or occasionally a comma-separated list
    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           2, 1.0, 'synonym_exact'
    FROM master_species sp
    WHERE EXISTS (
        SELECT 1 FROM unnest(string_to_array(sp.synonym, ',')) AS syn(name)
        WHERE lower(btrim(syn.name)) IN (lower(%(raw)s), lower(%(clean)s))
    )

    UNION ALL

    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           3, 1.0, 'genus_exact'
    FROM master_species sp
    WHERE lower(sp.genus) = lower(%(genus)s)

    UNION ALL

    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           4, 1.0, 'family_exact'
    FROM master_species sp
    WHERE lower(sp.family) IN (lower(%(raw)s), lower(%(clean)s))

    UNION ALL

    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           5, 1.0, 'order_exact'
    FROM master_species sp
    WHERE lower(sp.ordr) IN (lower(%(raw)s), lower(%(clean)s))

    UNION ALL

    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           6, similarity(sp.family, %(clean)s), 'family_fuzzy'
    FROM master_species sp
    WHERE sp.family %% %(clean)s
      AND similarity(sp.family, %(clean)s) >= 0.65

    UNION ALL

    SELECT sp.ncbi_taxon_id, sp.class, sp.species,
           7, similarity(sp.ordr, %(clean)s), 'order_fuzzy'
    FROM master_species sp
    WHERE sp.ordr %% %(clean)s
      AND similarity(sp.ordr, %(clean)s) >= 0.65
) m
ORDER BY m.priority, m.sim DESC, m.ncbi_taxon_id NULLS LAST, m.species_matched
LIMIT 1
"""


def match_params(nominal_species_id):
    """
    Query parameters for MATCH_QUERY. The cleaned name is matched alongside the
    raw one so qualifier noise ('Gobiidae spp.', 'Cynoglossidae (Family)') hits
    an exact tier instead of only the fuzzy ones.
    """
    raw = (nominal_species_id or "").strip()
    clean = clean_species_name(raw) or raw
    return {'raw': raw, 'clean': clean, 'genus': clean.split(' ')[0] if clean else ''}


def get_species_info(db_params, og_id):
    """
    Return (nominal_species_id, taxon_id, tax_class, match_level) for an og_id,
    or None if the sample row is missing. match_level is None when
    master_species has no match at all; otherwise one of the tiers in MATCH_QUERY.
    """
    conn = None
    try:
        conn = psycopg2.connect(**db_params)
        with conn.cursor() as cur:
            # IMPORTANT: pass a tuple, not a bare string
            cur.execute("SELECT nominal_species_id FROM sample WHERE og_id = %s", (og_id,))
            sample_row = cur.fetchone()
            if not sample_row:
                print(f"[WARN] No sample row found for OG ID: {og_id}")
                return None
            nominal_species_id = sample_row[0]
            if not (nominal_species_id or "").strip():
                return (nominal_species_id, None, None, None)

            cur.execute(MATCH_QUERY, match_params(nominal_species_id))
            match = cur.fetchone()
            # Deliberately quiet on a miss: the caller falls back to the
            # NCBI taxdump, so warning here would report a failure that
            # mostly does not happen. Per-sample provenance goes to the
            # resolution report, and a real miss is reported at the end.
            if not match:
                return (nominal_species_id, None, None, None)
            taxon_id, tax_class, _species_matched, match_level = match
            return (nominal_species_id, taxon_id, tax_class, match_level)
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
    Matched alongside the raw name in MATCH_QUERY, and used as the taxdump query.

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
    (nominal_species_id, taxon_id, tax_class, source, match_level, taxdump_rank)
    for a sample.

    master_species is tried first, but a hit there is only trusted as far as its
    tier deserves (see EXACT_SPECIES_LEVELS):

    - exact species/synonym: class and taxid from the table; a missing taxid
      (40% of master_species rows have none) is filled from the taxdump.
    - genus/family/order exact: the table's class is kept -- it is the curated
      spelling BUSCO routing and the resource profiles key on, and congeners share
      it -- but the taxid comes from the taxdump, which can name the sample's own
      species, or at least the real genus/family node. The table's taxid is some
      other species in that group and is used only when there is no taxdump.
    - fuzzy family/order: a last resort, used only when the taxdump has nothing.
    - no match: everything from the taxdump. This is the normal case for
      invertebrates, which master_species does not carry.

    `source` is 'db', 'db+taxdump', 'taxdump' or 'unresolved'. match_level is
    the master_species tier ('' for none) and taxdump_rank the rank of the
    taxdump node that was used ('' if it was not), so the resolution report shows
    every sample that leaned on a congener or a fuzzy match.
    """
    species_info = get_species_info(db_params, og_id)
    if not species_info:
        return None, None, None, 'unresolved', '', ''

    nominal_species_id, db_taxon, db_class, match_level = species_info
    match_level = match_level or ''
    db_had_taxon = db_taxon not in (None, "")
    db_had_class = bool(db_class)

    if match_level in EXACT_SPECIES_LEVELS and db_had_taxon and db_had_class:
        return nominal_species_id, db_taxon, db_class, 'db', match_level, ''

    lineage = {}
    query_name = clean_species_name(nominal_species_id)
    if resolver is not None and query_name:
        try:
            lineage = resolver.lineage_for_name(query_name)
        except Exception as exc:
            print(f"[WARN] taxdump lookup failed for {og_id} "
                  f"('{query_name}'): {exc}")

    # The taxid of whatever node matched, at whatever rank. FCS-GX takes a
    # --tax-id at any rank, so a class-rank id beats no id at all.
    lineage_taxon = lineage.get('matched_taxid') if lineage else None
    lineage_class = lineage.get('class') if lineage else None

    # A real NCBI name beats a near-miss spelling of a family or order.
    trust_table = bool(match_level) and not (match_level in FUZZY_LEVELS and lineage)

    used_taxdump = False
    if trust_table:
        tax_class = db_class
        if not db_had_class and lineage_class:
            tax_class = lineage_class
            used_taxdump = True
        if match_level in EXACT_SPECIES_LEVELS and db_had_taxon:
            taxon_id = db_taxon
        elif lineage_taxon is not None:
            # Never let a congener's taxid win over a lookup of the sample's own name.
            taxon_id = lineage_taxon
            used_taxdump = True
        else:
            taxon_id = db_taxon
    else:
        taxon_id, tax_class = lineage_taxon, lineage_class
        used_taxdump = bool(lineage)

    resolved = taxon_id not in (None, "") and bool(tax_class)
    if not resolved:
        source = 'unresolved'
    elif trust_table:
        source = 'db+taxdump' if used_taxdump else 'db'
    else:
        source = 'taxdump'

    taxdump_rank = lineage.get('matched_rank', '') if used_taxdump else ''
    return nominal_species_id, taxon_id, tax_class, source, match_level, taxdump_rank


# What an unresolved field is written as. Matches the `^(unknown|None)$` pattern
# assets/schema_input.json already allows for taxon_id, so these rows parse and
# reach validateTaxonomy() rather than failing schema validation first.
UNKNOWN = 'unknown'

RESOLUTION_COLUMNS = ('sample', 'nominal_species_id', 'taxon_id', 'class', 'source',
                      'match_level', 'taxdump_rank')


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
                             "resolve taxon_id/class when master_species has no "
                             "exact match. Without it, unmatched samples fail the run.")
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
        print("[WARN] No taxdump available: taxon_id/class come from master_species "
              "only, so any sample it does not carry will fail below, and a "
              "genus-level match keeps a congener's taxon_id.")

    rows = []
    resolution_rows = []
    unresolved = []
    for og_id, files in pairs:
        (nominal_species_id, taxon_id, tax_class, source,
         match_level, taxdump_rank) = resolve_species_info(db_params, og_id, resolver)

        resolution_rows.append({
            'sample': og_id,
            'nominal_species_id': nominal_species_id or '',
            'taxon_id': taxon_id if taxon_id is not None else '',
            'class': tax_class or '',
            'source': source,
            'match_level': match_level,
            'taxdump_rank': taxdump_rank,
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
        print("[FIX] These names were not found in master_species OR in the NCBI "
              "taxonomy, which usually means the nominal_species_id is misspelled "
              "(e.g. 'Actinaria' for 'Actiniaria') or is not a taxon name at all "
              "(e.g. 'Larval fish'). Fix it either way round: correct the "
              "nominal_species_id in the sample table and delete the published "
              "samplesheet so it is regenerated, or edit the taxon_id and class "
              "columns in the published samplesheet and re-run. If the name is right "
              "and simply missing from NCBI's taxonomy, fill the taxon_id and class in the "
              "published samplesheet by hand.")

    taxdump_rescued = sum(1 for r in resolution_rows if r['source'] in ('taxdump', 'db+taxdump'))
    if taxdump_rescued:
        print(f"[INFO] {taxdump_rescued} of {len(resolution_rows)} samples had their "
              f"taxonomy completed from the NCBI taxdump rather than master_species.")


if __name__ == "__main__":
    sys.exit(main() or 0)
