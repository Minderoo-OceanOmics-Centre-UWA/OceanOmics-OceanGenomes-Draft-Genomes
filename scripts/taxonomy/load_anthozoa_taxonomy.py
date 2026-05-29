#!/usr/bin/env python3
"""
Bulk-load NCBI species-rank taxa into the OceanOmics `species` table.

Default target is class Anthozoa (corals, anemones, sea pens). Source is the
NCBI `new_taxdump` archive (rankedlineage.dmp), which carries the full
class/order/family/genus/species lineage plus the NCBI taxon_id for every
named taxon. Existing rows in the `species` table are never touched
(ON CONFLICT (species) DO NOTHING).

Default mode is dry-run: writes a CSV preview of the rows that would be
inserted. Use --apply <db.cfg> to actually load them.
"""

import argparse
import configparser
import csv
import datetime
import os
import sys
import tarfile
import urllib.request
from pathlib import Path

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
RANKEDLINEAGE_MEMBER = "rankedlineage.dmp"

# rankedlineage.dmp columns, in order:
#   tax_id | name | species | genus | family | order | class | phylum | kingdom | superkingdom
COL_TAX_ID, COL_NAME, COL_SPECIES, COL_GENUS, COL_FAMILY, COL_ORDER, COL_CLASS = range(7)


def load_db_config(config_file):
    config = configparser.ConfigParser()
    config.read(config_file)
    return {
        'dbname': config.get('postgres', 'dbname'),
        'user': config.get('postgres', 'user'),
        'password': config.get('postgres', 'password'),
        'host': config.get('postgres', 'host'),
        'port': config.getint('postgres', 'port'),
    }


def ensure_taxdump(cache_dir: Path, refresh: bool) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    archive = cache_dir / "new_taxdump.tar.gz"
    if archive.exists() and not refresh:
        print(f"[INFO] Using cached taxdump: {archive}", file=sys.stderr)
        return archive
    print(f"[INFO] Downloading {TAXDUMP_URL} -> {archive}", file=sys.stderr)
    tmp = archive.with_suffix(".tar.gz.part")
    with urllib.request.urlopen(TAXDUMP_URL) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.rename(archive)
    return archive


def extract_rankedlineage(archive: Path, cache_dir: Path) -> Path:
    out_path = cache_dir / RANKEDLINEAGE_MEMBER
    if out_path.exists():
        return out_path
    print(f"[INFO] Extracting {RANKEDLINEAGE_MEMBER} from archive", file=sys.stderr)
    with tarfile.open(archive, "r:gz") as tf:
        member = tf.getmember(RANKEDLINEAGE_MEMBER)
        with tf.extractfile(member) as src, open(out_path, "wb") as dst:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                dst.write(chunk)
    return out_path


def parse_rankedlineage_line(line: str):
    # rankedlineage.dmp lines are terminated with "\t|\n"; fields are "\t|\t"-separated.
    line = line.rstrip("\n")
    if line.endswith("\t|"):
        line = line[:-2]
    return [f.strip() for f in line.split("\t|\t")]


def collect_species_rows(rankedlineage_path: Path, target_class: str):
    """
    Yield dicts for species-rank entries belonging to `target_class`.

    In rankedlineage.dmp each lineage column is empty for the rank of the row
    itself and populated for ranks above it. So a species-rank row has an
    empty `species` column and a populated `genus` column; `name` carries the
    binomial. Subspecies/strain rows have both `species` and `genus`
    populated. Genus-rank rows have both empty.
    """
    target_lower = target_class.lower()
    kept = 0
    skipped_noise = 0
    with open(rankedlineage_path, "r", encoding="utf-8") as fh:
        for raw in fh:
            cols = parse_rankedlineage_line(raw)
            if len(cols) < 10:
                continue
            cls = cols[COL_CLASS]
            if cls.lower() != target_lower:
                continue
            species_col = cols[COL_SPECIES]
            genus = cols[COL_GENUS]
            if species_col or not genus:
                # subspecies/strain (species_col populated) or higher rank (genus empty)
                continue
            name = cols[COL_NAME]
            if not name or " " not in name:
                # not a binomial — skip
                continue
            # Drop common "sp."/"cf."/"aff."/hybrid noise that won't match real samples
            lowered = name.lower()
            if (" sp." in lowered or " cf." in lowered or " aff." in lowered
                    or lowered.startswith("unidentified ") or lowered.startswith("uncultured ")
                    or lowered.startswith("unclassified ") or " x " in lowered
                    or "environmental" in lowered):
                skipped_noise += 1
                continue
            try:
                tax_id = int(cols[COL_TAX_ID])
            except ValueError:
                continue
            epithet = name[len(genus):].strip() if name.startswith(genus + " ") else ""
            kept += 1
            yield {
                "species": name,
                "class": cls,
                "ordr": cols[COL_ORDER] or None,
                "family": cols[COL_FAMILY] or None,
                "genus": genus or None,
                "epithet": epithet or None,
                "ncbi_taxon_id": tax_id,
            }
    print(f"[INFO] Matched {kept} species-rank rows in class {target_class}", file=sys.stderr)
    if skipped_noise:
        print(f"[INFO] Skipped {skipped_noise} noisy names (sp./cf./aff./hybrid/unidentified)", file=sys.stderr)


def write_preview_csv(rows, out_path: Path, provenance: str):
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["species", "class", "ordr", "family", "genus", "epithet", "ncbi_taxon_id", "comments"])
        n = 0
        for r in rows:
            w.writerow([
                r["species"], r["class"], r["ordr"], r["family"],
                r["genus"], r["epithet"], r["ncbi_taxon_id"], provenance,
            ])
            n += 1
    return n


def apply_inserts(rows, db_params, provenance: str, batch_size: int = 500):
    try:
        import psycopg2
        from psycopg2.extras import execute_values
    except ImportError:
        print("[ERROR] psycopg2 is required for --apply. Install it or run in dry-run mode.", file=sys.stderr)
        sys.exit(2)

    sql = """
        INSERT INTO species
            (species, class, ordr, family, genus, epithet, ncbi_taxon_id, comments)
        VALUES %s
        ON CONFLICT (species) DO NOTHING
    """
    conn = psycopg2.connect(**db_params)
    inserted_total = 0
    attempted_total = 0
    try:
        with conn:
            with conn.cursor() as cur:
                batch = []
                for r in rows:
                    batch.append((
                        r["species"], r["class"], r["ordr"], r["family"],
                        r["genus"], r["epithet"], r["ncbi_taxon_id"], provenance,
                    ))
                    if len(batch) >= batch_size:
                        execute_values(cur, sql, batch)
                        inserted_total += cur.rowcount
                        attempted_total += len(batch)
                        batch.clear()
                if batch:
                    execute_values(cur, sql, batch)
                    inserted_total += cur.rowcount
                    attempted_total += len(batch)
    finally:
        conn.close()
    return attempted_total, inserted_total


def main():
    p = argparse.ArgumentParser(description="Bulk-load NCBI species-rank taxa for a target class into the OceanOmics species table.")
    p.add_argument("--class", dest="target_class", default="Anthozoa",
                   help="NCBI class to load (default: Anthozoa)")
    p.add_argument("--cache-dir", default=str(Path(__file__).parent / ".taxdump_cache"),
                   help="Where to cache the downloaded NCBI taxdump (default: ./.taxdump_cache next to this script)")
    p.add_argument("--refresh", action="store_true",
                   help="Force re-download even if a cached taxdump exists")
    p.add_argument("--preview-out", default=None,
                   help="Path for the CSV preview written in dry-run mode (default: <class>_preview_<date>.csv in cwd)")
    p.add_argument("--apply", metavar="DB_CFG", default=None,
                   help="Insert rows into the species table using this postgres cfg file. Without this flag, runs in dry-run mode.")
    p.add_argument("--limit", type=int, default=None,
                   help="For testing: only process the first N matching rows")
    args = p.parse_args()

    cache_dir = Path(args.cache_dir)
    archive = ensure_taxdump(cache_dir, refresh=args.refresh)
    rankedlineage_path = extract_rankedlineage(archive, cache_dir)

    today = datetime.date.today().isoformat()
    provenance = f"auto-loaded from NCBI taxdump {today}"

    def row_iter():
        seen = 0
        for r in collect_species_rows(rankedlineage_path, args.target_class):
            if args.limit is not None and seen >= args.limit:
                return
            seen += 1
            yield r

    if args.apply is None:
        out_path = Path(args.preview_out) if args.preview_out else Path(f"{args.target_class.lower()}_preview_{today}.csv")
        n = write_preview_csv(row_iter(), out_path, provenance)
        print(f"[OK] DRY-RUN wrote {n} rows to {out_path}")
        print(f"[NEXT] Inspect the CSV, then re-run with: --apply {args.apply or '<db.cfg>'}")
        return

    db_params = load_db_config(args.apply)
    attempted, inserted = apply_inserts(row_iter(), db_params, provenance)
    skipped = attempted - inserted
    print(f"[OK] Attempted {attempted} rows; inserted {inserted}; skipped (already present) {skipped}")
    print(f"[NEXT] Verify with: SELECT count(*), ordr FROM species WHERE class = '{args.target_class}' GROUP BY ordr ORDER BY count DESC;")


if __name__ == "__main__":
    main()
