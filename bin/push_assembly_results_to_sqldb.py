#!/usr/bin/env python3
import sys
import os
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

def normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
        .str.replace(r'[^a-z0-9_]', '_', regex=True)
    )
    return df

def to_int(val):
    if pd.isna(val) or val == "":
        return None
    return int(str(val).replace(",", ""))

def to_float(val):
    if pd.isna(val) or val == "":
        return None
    return float(str(val).replace("%", "").replace(",", ""))

# ----------------------------
# SQL
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, homozygosity, heterozygosity, genomesize, repeatsize, uniquesize, modelfit, errorrate
)
VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    homozygosity = EXCLUDED.homozygosity,
    heterozygosity = EXCLUDED.heterozygosity,
    genomesize = EXCLUDED.genomesize,
    repeatsize = EXCLUDED.repeatsize,
    uniquesize = EXCLUDED.uniquesize,
    modelfit = EXCLUDED.modelfit,
    errorrate = EXCLUDED.errorrate;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Push GenomeScope results for a single sample into PostgreSQL"
    )
    parser.add_argument("-c", "--config", required=True,
                        help="Path to config file with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    parser.add_argument("-f", "--file", required=True,
                        help="TSV file containing results for a single sample")
    args = parser.parse_args()

    # Load config
    cfg = load_kv_config(args.config)
    db_params = {
        "dbname":   cfg.get("DB_NAME"),
        "user":     cfg.get("DB_USER"),
        "password": cfg.get("DB_PASSWORD"),
        "host":     cfg.get("DB_HOST"),
        "port":     int(cfg.get("DB_PORT", "5432")),
    }

    # Load sample TSV
    if not os.path.isfile(args.file):
        sys.exit(f"❌ Input file not found: {args.file}")

    df = pd.read_csv(args.file, sep="\t")
    df = normalise_columns(df)

    if "sample" in df.columns:
        # Split "Sample" into og_id.tech.seq_date
        split = df["sample"].astype(str).str.split(".", expand=True)
        if split.shape[1] >= 3:
            df["og_id"] = split[0]
            df["tech"] = split[1]
            df["seq_date"] = split[2]
        df.drop(columns=["sample", "tech"], errors="ignore", inplace=True)

    # Connect to DB
    try:
        conn = psycopg2.connect(**db_params)
        cur = conn.cursor()

        row_count = 0
        for _, row in df.iterrows():
            params = (
                row.get("og_id"),
                str(row.get("seq_date")) if row.get("seq_date") else None,
                to_float(row.get("homozygosity")),
                to_float(row.get("heterozygosity")),
                to_int(row.get("genomesize")),
                to_int(row.get("repeatsize")),
                to_int(row.get("uniquesize")),
                to_float(row.get("modelfit")),
                to_float(row.get("errorrate")),
            )
            cur.execute(UPSERT_SQL, params)
            row_count += 1

        conn.commit()
        print(f"✅ Successfully upserted {row_count} rows from {args.file}")

    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        if 'cur' in locals():
            cur.close()
        if 'conn' in locals():
            conn.close()

if __name__ == "__main__":
    main()
