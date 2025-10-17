#!/usr/bin/env python3
# Usage (in nf-core module):
#   python 04a_push_gfa_results_to_sqldb.py -c ../configfile.txt -f ${gfastats_sample_tsv}
import sys
import os
import re
import argparse
import psycopg2
import pandas as pd
import numpy as np

# ----------------------------
# Config & helpers
# ----------------------------
def load_kv_config(path: str) -> dict:
    cfg = {}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg

def norm_name(col: str) -> str:
    return re.sub(r"[^a-z0-9_]", "_", col.strip().lower())

def normalise_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [norm_name(c) for c in out.columns]
    return out

def parse_int(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().replace(",", "")
    if s == "" or s.lower() == "nan":
        return None
    try:
        return int(float(s))
    except ValueError:
        return None

def parse_float(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().replace(",", "").replace("%", "")
    if s == "" or s.lower() == "nan":
        return None
    try:
        return float(s)
    except ValueError:
        return None

def to_text_or_none(x):
    return str(x) if (x is not None and not (isinstance(x, float) and np.isnan(x))) else None

def safe_dot_part(s, idx):
    parts = s.split(".")
    return parts[idx] if len(parts) > idx else None

def safe_us_part(s, us_idx, dot_idx=0):
    base = safe_dot_part(s, dot_idx) or ""
    us = base.split("_")
    return us[us_idx] if len(us) > us_idx else None

def ensure_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure og_id, seq_date, stage, haplotype exist.
    If not present but a 'filename' column exists, derive:
      og_id     <- filename split by '_' (first part) from the first dot-part
      seq_date  <- filename split by '_' (second part) from first dot-part, with leading 'v' stripped
      stage     <- integer from the 3rd dot-part
      haplotype <- first '_' part from the 5th dot-part (if present), else None
    """
    df = df.copy()
    cols = set(df.columns)
    if "filename" in cols:
        fname = df["filename"].astype(str)
        if "og_id" not in cols:
            df["og_id"] = fname.apply(lambda x: safe_us_part(x, 0, dot_idx=0))
        if "seq_date" not in cols:
            df["seq_date"] = fname.apply(lambda x: (safe_us_part(x, 1, dot_idx=0) or "").lstrip("v") or None)
        if "stage" not in cols:
            df["stage"] = fname.apply(lambda x: parse_int(safe_dot_part(x, 2)))
        if "haplotype" not in cols:
            df["haplotype"] = fname.apply(lambda x: safe_us_part(x, 0, dot_idx=4))
    return df

# ----------------------------
# SQL
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, stage, haplotype,
    num_contigs, contig_n50, contig_n50_size_mb,
    num_scaffolds, scaffold_n50, scaffold_n50_size_mb,
    largest_scaffold, largest_scaffold_size_mb,
    total_scaffold_length, total_scaffold_length_size_mb,
    gc_content_percent
)
VALUES (
    %(og_id)s, %(seq_date)s,
    %(num_contigs)s, %(contig_n50)s, %(contig_n50_size_mb)s,
    %(num_scaffolds)s, %(scaffold_n50)s, %(scaffold_n50_size_mb)s,
    %(largest_scaffold)s, %(largest_scaffold_size_mb)s,
    %(total_scaffold_length)s, %(total_scaffold_length_size_mb)s,
    %(gc_content_percent)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    num_contigs = EXCLUDED.num_contigs,
    contig_n50 = EXCLUDED.contig_n50,
    contig_n50_size_mb = EXCLUDED.contig_n50_size_mb,
    num_scaffolds = EXCLUDED.num_scaffolds,
    scaffold_n50 = EXCLUDED.scaffold_n50,
    scaffold_n50_size_mb = EXCLUDED.scaffold_n50_size_mb,
    largest_scaffold = EXCLUDED.largest_scaffold,
    largest_scaffold_size_mb = EXCLUDED.largest_scaffold_size_mb,
    total_scaffold_length = EXCLUDED.total_scaffold_length,
    total_scaffold_length_size_mb = EXCLUDED.total_scaffold_length_size_mb,
    gc_content_percent = EXCLUDED.gc_content_percent;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert gfastats per-sample TSV into PostgreSQL (DB creds from key=value config)."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("-f", "--file", required=True,
                    help="Per-sample gfastats TSV")
    args = ap.parse_args()

    # Load DB config
    cfg = load_kv_config(args.config)
    db_params = {
        "dbname":   cfg.get("DB_NAME"),
        "user":     cfg.get("DB_USER"),
        "password": cfg.get("DB_PASSWORD"),
        "host":     cfg.get("DB_HOST"),
        "port":     int(cfg.get("DB_PORT", "5432")),
    }

    # Load TSV
    if not os.path.isfile(args.file):
        sys.exit(f"❌ Input file not found: {args.file}")
    df = pd.read_csv(args.file, sep="\t")
    df = normalise_df(df)
    df = ensure_ids(df)

    # Validate identifiers
    for col in ["og_id", "seq_date"]:
        if col not in df.columns:
            sys.exit(f"❌ Required ID column missing: {col}. Provide it in TSV or include a 'filename' column to derive it.")

    # Coerce/prepare numeric fields (tolerant to missing columns)
    num_int_cols = [
        "num_contigs", "contig_n50", "num_scaffolds", "scaffold_n50",
        "largest_scaffold", "total_scaffold_length"
    ]
    num_float_cols = [
        "contig_n50_size_mb", "scaffold_n50_size_mb",
        "largest_scaffold_size_mb", "total_scaffold_length_size_mb",
        "gc_content_percent"
    ]
    for c in num_int_cols:
        if c in df.columns: df[c] = df[c].apply(parse_int)
    for c in num_float_cols:
        if c in df.columns: df[c] = df[c].apply(parse_float)

    # Connect & upsert
    try:
        conn = psycopg2.connect(**db_params)
    except Exception as e:
        sys.exit(f"❌ Failed to connect to DB: {e}")

    try:
        with conn, conn.cursor() as cur:
            upserted = 0
            for _, row in df.iterrows():
                params = {
                    "og_id":                         to_text_or_none(row.get("og_id")),
                    "seq_date":                      to_text_or_none(row.get("seq_date")),
                    "num_contigs":                   parse_int(row.get("num_contigs")),
                    "contig_n50":                    parse_int(row.get("contig_n50")),
                    "contig_n50_size_mb":            parse_float(row.get("contig_n50_size_mb")),
                    "num_scaffolds":                 parse_int(row.get("num_scaffolds")),
                    "scaffold_n50":                  parse_int(row.get("scaffold_n50")),
                    "scaffold_n50_size_mb":          parse_float(row.get("scaffold_n50_size_mb")),
                    "largest_scaffold":              parse_int(row.get("largest_scaffold")),
                    "largest_scaffold_size_mb":      parse_float(row.get("largest_scaffold_size_mb")),
                    "total_scaffold_length":         parse_int(row.get("total_scaffold_length")),
                    "total_scaffold_length_size_mb": parse_float(row.get("total_scaffold_length_size_mb")),
                    "gc_content_percent":            parse_float(row.get("gc_content_percent")),
                }
                cur.execute(UPSERT_SQL, params)
                upserted += 1

            print(f"✅ Successfully upserted {upserted} row(s) from {args.file}")

    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    main()
