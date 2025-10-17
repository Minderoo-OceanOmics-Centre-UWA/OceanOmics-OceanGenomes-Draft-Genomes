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

def norm_name(col: str) -> str:
    return re.sub(r"[^a-z0-9_]", "_", col.strip().lower())

def normalise_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [norm_name(c) for c in out.columns]
    return out

def to_int_or_none(x):
    if pd.isna(x) or x == "":
        return None
    try:
        return int(str(x).replace(",", ""))
    except Exception:
        return None

def to_float_or_none(x):
    if pd.isna(x) or x == "":
        return None
    try:
        v = float(str(x).replace("%", "").replace(",", ""))
        return float(v) if np.isfinite(v) else None
    except Exception:
        return None

def to_text_or_none(x):
    return str(x) if (x is not None and not (isinstance(x, float) and np.isnan(x))) else None

def ensure_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure og_id and seq_date exist.
    If a 'sample' column like OGID.TECH.SEQDATE exists, split to fill them.
    """
    df = df.copy()
    if "og_id" not in df.columns and "sample" in df.columns:
        split = df["sample"].astype(str).str.split(".", expand=True)
        if split.shape[1] >= 1:
            df["og_id"] = split[0]
        if split.shape[1] >= 3 and "seq_date" not in df.columns:
            df["seq_date"] = split[2]
    return df

# ----------------------------
# SQL
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, complete, single_copy, multi_copy, fragmented,
    missing, n_markers, domain, number_of_scaffolds, number_of_contigs, total_length, percent_gaps, 
    scaffold_n50, contigs_n50
)
VALUES (
    %(og_id)s, %(seq_date)s, %(complete)s, %(single_copy)s, %(multi_copy)s, %(fragmented)s,
    %(missing)s, %(n_markers)s, %(domain)s, %(number_of_scaffolds)s, %(number_of_contigs)s, %(total_length)s,  
    %(percent_gaps)s, %(scaffold_n50)s, %(contigs_n50)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    complete          = EXCLUDED.complete,
    single_copy       = EXCLUDED.single_copy,
    multi_copy        = EXCLUDED.multi_copy,
    fragmented        = EXCLUDED.fragmented,
    missing           = EXCLUDED.missing,
    n_markers         = EXCLUDED.n_markers,
    domain            = EXCLUDED.domain,
    number_of_scaffolds = EXCLUDED.number_of_scaffolds,
    number_of_contigs = EXCLUDED.number_of_contigs,
    total_length      = EXCLUDED.total_length,
    percent_gaps      = EXCLUDED.percent_gaps,
    scaffold_n50      = EXCLUDED.scaffold_n50,
    contigs_n50       = EXCLUDED.contigs_n50;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert BUSCO per-sample TSV into PostgreSQL (DB creds from key=value config)."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("-f", "--file", required=True,
                    help="BUSCO compiled TSV for a single sample")
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

    # Minimal required fields
    required = {"og_id", "seq_date"}
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"❌ Missing required columns in TSV: {missing}. "
                 f"Provide og_id (or sample) and seq_date.")

    # Coerce numeric/text fields safely
    df["complete"]            = df.get("complete").apply(to_float_or_none) if "complete" in df else None
    df["single_copy"]         = df.get("single_copy").apply(to_float_or_none) if "single_copy" in df else None
    df["multi_copy"]          = df.get("multi_copy").apply(to_float_or_none) if "multi_copy" in df else None
    df["fragmented"]          = df.get("fragmented").apply(to_float_or_none) if "fragmented" in df else None
    df["missing"]             = df.get("missing").apply(to_float_or_none) if "missing" in df else None
    df["n_markers"]           = df.get("n_markers").apply(to_int_or_none)   if "n_markers" in df else None
    df["domain"]              = df.get("domain")                             if "domain" in df else None
    df["number_of_scaffolds"] = df.get("number_of_scaffolds").apply(to_int_or_none) if "number_of_scaffolds" in df else None
    df["number_of_contigs"]   = df.get("number_of_contigs").apply(to_int_or_none)   if "number_of_contigs" in df else None
    df["total_length"]        = df.get("total_length").apply(to_int_or_none)        if "total_length" in df else None
    df["percent_gaps"]        = df.get("percent_gaps").apply(to_float_or_none)      if "percent_gaps" in df else None
    df["scaffold_n50"]        = df.get("scaffold_n50").apply(to_float_or_none)      if "scaffold_n50" in df else None
    df["contigs_n50"]         = df.get("contigs_n50").apply(to_int_or_none)         if "contigs_n50" in df else None

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
                    "og_id":              to_text_or_none(row.get("og_id")),
                    "seq_date":           to_text_or_none(row.get("seq_date")),
                    "complete":           row.get("complete"),
                    "single_copy":        row.get("single_copy"),
                    "multi_copy":         row.get("multi_copy"),
                    "fragmented":         row.get("fragmented"),
                    "missing":            row.get("missing"),
                    "n_markers":          row.get("n_markers"),
                    "domain":             to_text_or_none(row.get("domain")),
                    "number_of_scaffolds":row.get("number_of_scaffolds"),
                    "number_of_contigs":  row.get("number_of_contigs"),
                    "total_length":       row.get("total_length"),
                    "percent_gaps":       row.get("percent_gaps"),
                    "scaffold_n50":       row.get("scaffold_n50"),
                    "contigs_n50":        row.get("contigs_n50"),
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
