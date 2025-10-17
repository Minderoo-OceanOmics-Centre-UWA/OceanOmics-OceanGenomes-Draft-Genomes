#!/usr/bin/env python3
import sys
import os
import re
import argparse
import psycopg2
import pandas as pd
import numpy as np

# ----------------------------
# Helpers
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

def normalise_name(col: str) -> str:
    # lowercase, trim, replace everything not [a-z0-9_] with "_"
    return re.sub(r"[^a-z0-9_]", "_", col.strip().lower())

def normalise_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [normalise_name(c) for c in df.columns]
    return df

def to_int(x):
    if pd.isna(x) or x == "":
        return 0
    try:
        return int(str(x).replace(",", ""))
    except Exception:
        return 0

def to_float(x):
    if pd.isna(x) or x == "":
        return 0.0
    try:
        val = float(str(x).replace("%", "").replace(",", ""))
        return float(val) if np.isfinite(val) else 0.0
    except Exception:
        return 0.0

# ----------------------------
# SQL
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, passed_filter_reads, low_quality_reads, too_many_n_reads, too_short_reads, too_long_reads,
    raw_total_reads, raw_total_bases, raw_q20_bases, raw_q30_bases, raw_q20_rate, raw_q30_rate, raw_read1_mean_length, raw_read2_mean_length,
    raw_gc_content, total_reads, total_bases, q20_bases, q30_bases, q20_rate, q30_rate, read1_mean_length, read2_mean_length, gc_content, mach, seq_date, initial
)
VALUES (
    %(og_id)s, %(passed_filter_reads)s, %(low_quality_reads)s, %(too_many_n_reads)s, %(too_short_reads)s, %(too_long_reads)s,
    %(raw_total_reads)s, %(raw_total_bases)s, %(raw_q20_bases)s, %(raw_q30_bases)s, %(raw_q20_rate)s, %(raw_q30_rate)s, %(raw_read1_mean_length)s, %(raw_read2_mean_length)s,
    %(raw_gc_content)s, %(total_reads)s, %(total_bases)s, %(q20_bases)s, %(q30_bases)s, %(q20_rate)s, %(q30_rate)s, %(read1_mean_length)s, %(read2_mean_length)s, %(gc_content)s, %(mach)s, %(seq_date)s, %(initial)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    mach = EXCLUDED.mach,
    initial = EXCLUDED.initial,
    passed_filter_reads = EXCLUDED.passed_filter_reads,
    low_quality_reads = EXCLUDED.low_quality_reads,
    too_many_n_reads = EXCLUDED.too_many_n_reads,
    too_short_reads = EXCLUDED.too_short_reads,
    too_long_reads = EXCLUDED.too_long_reads,
    raw_total_reads = EXCLUDED.raw_total_reads,
    raw_total_bases = EXCLUDED.raw_total_bases,
    raw_q20_bases = EXCLUDED.raw_q20_bases,
    raw_q30_bases = EXCLUDED.raw_q30_bases,
    raw_q20_rate = EXCLUDED.raw_q20_rate,
    raw_q30_rate = EXCLUDED.raw_q30_rate,
    raw_read1_mean_length = EXCLUDED.raw_read1_mean_length,
    raw_read2_mean_length = EXCLUDED.raw_read2_mean_length,
    raw_gc_content = EXCLUDED.raw_gc_content,
    total_reads = EXCLUDED.total_reads,
    total_bases = EXCLUDED.total_bases,
    q20_bases = EXCLUDED.q20_bases,
    q30_bases = EXCLUDED.q30_bases,
    q20_rate = EXCLUDED.q20_rate,
    q30_rate = EXCLUDED.q30_rate,
    read1_mean_length = EXCLUDED.read1_mean_length,
    read2_mean_length = EXCLUDED.read2_mean_length,
    gc_content = EXCLUDED.gc_content;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert FastP per-sample TSV into PostgreSQL (DB creds from key=value config)."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("-f", "--file", required=True,
                    help="FastP TSV for a single sample")
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

    # Accept 'sample' as alias for 'og_id'
    if "og_id" not in df.columns and "sample" in df.columns:
        df["og_id"] = df["sample"]

    # Check required identifiers
    missing_ids = [c for c in ["og_id", "seq_date"] if c not in df.columns]
    if missing_ids:
        sys.exit(f"❌ Required column(s) missing in TSV: {missing_ids}. "
                 f"Expected at least og_id (or sample) and seq_date.")

    # Connect DB
    try:
        conn = psycopg2.connect(**db_params)
    except Exception as e:
        sys.exit(f"❌ Failed to connect to DB: {e}")

    try:
        with conn, conn.cursor() as cur:
            upserted = 0
            for _, row in df.iterrows():
                params = {
                    "og_id":              str(row.get("og_id")) if pd.notna(row.get("og_id")) else None,
                    "mach":               str(row.get("mach")) if pd.notna(row.get("mach")) else None,
                    "seq_date":           str(row.get("seq_date")) if pd.notna(row.get("seq_date")) else None,
                    "initial":            str(row.get("initial")) if pd.notna(row.get("initial")) else None,
                    "passed_filter_reads":   to_int(row.get("passed_filter_reads")),
                    "low_quality_reads":     to_int(row.get("low_quality_reads")),
                    "too_many_n_reads":      to_int(row.get("too_many_n_reads")),
                    "too_short_reads":       to_int(row.get("too_short_reads")),
                    "too_long_reads":        to_int(row.get("too_long_reads")),
                    "raw_total_reads":       to_int(row.get("raw_total_reads")),
                    "raw_total_bases":       to_int(row.get("raw_total_bases")),
                    "raw_q20_bases":         to_int(row.get("raw_q20_bases")),
                    "raw_q30_bases":         to_int(row.get("raw_q30_bases")),
                    "raw_q20_rate":          to_float(row.get("raw_q20_rate")),
                    "raw_q30_rate":          to_float(row.get("raw_q30_rate")),
                    "raw_read1_mean_length": to_int(row.get("raw_read1_mean_length")),
                    "raw_read2_mean_length": to_int(row.get("raw_read2_mean_length")),
                    "raw_gc_content":        to_float(row.get("raw_gc_content")),
                    "total_reads":           to_int(row.get("total_reads")),
                    "total_bases":           to_int(row.get("total_bases")),
                    "q20_bases":             to_int(row.get("q20_bases")),
                    "q30_bases":             to_int(row.get("q30_bases")),
                    "q20_rate":              to_float(row.get("q20_rate")),
                    "q30_rate":              to_float(row.get("q30_rate")),
                    "read1_mean_length":     to_int(row.get("read1_mean_length")),
                    "read2_mean_length":     to_int(row.get("read2_mean_length")),
                    "gc_content":            to_float(row.get("gc_content")),
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
