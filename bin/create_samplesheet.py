#!/usr/bin/env python3
import psycopg2
import configparser
import sys
import csv
import os
import glob

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
    SELECT
        s.og_id,
        s.nominal_species_id,
        m.ncbi_taxon_id,
        m.class,
        m.species_matched,
        m.sim,
        m.match_level
    FROM sample s
    LEFT JOIN LATERAL (
        SELECT
            sp.ncbi_taxon_id,
            sp.class,
            sp.species AS species_matched,
            -- priority: 1 = species-level match, 2 = only higher-level match
            CASE 
                WHEN sp.species %% s.nominal_species_id THEN 1
                ELSE 2
            END AS priority,
            -- similarity score: max of whichever levels actually match
            GREATEST(
                CASE 
                    WHEN sp.species %% s.nominal_species_id 
                    THEN similarity(sp.species, s.nominal_species_id) 
                    ELSE 0 
                END,
                CASE 
                    WHEN sp.genus %% s.nominal_species_id 
                    THEN similarity(sp.genus, s.nominal_species_id) 
                    ELSE 0 
                END,
                CASE 
                    WHEN sp.family %% s.nominal_species_id 
                    THEN similarity(sp.family, s.nominal_species_id) 
                    ELSE 0 
                END,
                CASE 
                    WHEN sp.ordr %% s.nominal_species_id
                    THEN similarity(sp.ordr, s.nominal_species_id) 
                    ELSE 0 
                END
            ) AS sim,
            CASE 
                WHEN sp.species %% s.nominal_species_id THEN 'species'
                WHEN sp.genus   %% s.nominal_species_id 
                  OR sp.family %% s.nominal_species_id 
                  OR sp.ordr   %% s.nominal_species_id THEN 'higher_taxon'
                ELSE 'unknown'
            END AS match_level
        FROM species sp
        WHERE (
            sp.species %% s.nominal_species_id
            OR sp.genus  %% s.nominal_species_id
            OR sp.family %% s.nominal_species_id
            OR sp.ordr   %% s.nominal_species_id
        )
          AND sp.ncbi_taxon_id IS NOT NULL
        ORDER BY
            priority,   -- species-level matches first
            sim DESC    -- then highest similarity
        LIMIT 1
    ) m ON TRUE
    WHERE s.og_id = %s
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

                # Return exactly what the rest of the script expects
                return (nominal_species_id, taxon_id, tax_class)
            else:
                print(f"[WARN] No matching species data found for OG ID: {og_id}")
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

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage:\n  python create_samplesheet.py <db.cfg> <run_id> <fastq_dir>")
        sys.exit(1)

    config_file = sys.argv[1]
    run_id = sys.argv[2]
    fastq_dir = sys.argv[3]

    db_params = load_db_config(config_file)

    # Discover (og_id, [R1,R2]) pairs from directory
    pairs = discover_pairs_from_dir(fastq_dir)

    run_info = parse_run_id(run_id)
    date = run_info["date"]

    rows = []
    for og_id, files in pairs:
        species_info = get_species_info(db_params, og_id)
        if not species_info:
            # Already warned – emit placeholder values so schema requirements are satisfied
            nominal_species_id = "unknown"
            taxon_id = "unknown"
            tax_class = "unknown"
        else:
            nominal_species_id, taxon_id, tax_class = species_info
            nominal_species_id = nominal_species_id or "unknown"
            taxon_id = taxon_id if taxon_id not in (None, "") else "unknown"
            tax_class = tax_class or "unknown"

        fastq_1 = str(files[0])
        fastq_2 = str(files[1])

        row = {
            'sample': og_id,
            'run': run_id,
            'date': date,
            'prefix': f"{og_id}.ilmn.{date}",
            'nom_species_id': nominal_species_id,
            'taxon_id': taxon_id,
            'class': tax_class,
            'fastq_1': fastq_1,
            'fastq_2': fastq_2
        }
        rows.append(row)

    if not rows:
        print("[ERROR] No rows to write (all entries skipped).")
        sys.exit(1)

    output_file = f"{run_id}_samplesheet.csv"
    create_samplesheet(rows, output_file)
    print(f"[INFO] Wrote samplesheet to: {output_file}")
