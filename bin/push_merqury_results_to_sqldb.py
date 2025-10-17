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
    """Ensure og_id and seq_date exist. If 'sample' like OGID.TECH.SEQDATE is present, split to fill them."""
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
UPSERT_QV_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, unique_k_mers_assembly, k_mers_total, qv, error
)
VALUES (
    %(og_id)s, %(seq_date)s, %(unique_k_mers_assembly)s, %(k_mers_total)s, %(qv)s, %(error)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    unique_k_mers_assembly = EXCLUDED.unique_k_mers_assembly,
    k_mers_total           = EXCLUDED.k_mers_total,
    qv                     = EXCLUDED.qv,
    error                  = EXCLUDED.error;
"""

UPSERT_COMP_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, k_mer_set, solid_k_mers, total_k_mers, completeness
)
VALUES (
    %(og_id)s, %(seq_date)s, %(k_mer_set)s, %(solid_k_mers)s, %(total_k_mers)s, %(completeness)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    k_mer_set    = EXCLUDED.k_mer_set,
    solid_k_mers = EXCLUDED.solid_k_mers,
    total_k_mers = EXCLUDED.total_k_mers,
    completeness = EXCLUDED.completeness;
"""

# ----------------------------
# Loaders
# ----------------------------
def load_qv(tsv_path: str) -> pd.DataFrame:
    """
    Expected columns (case/spacing tolerant): sample (or og_id, seq_date), unique_k_mers_assembly, k_mers_total, qv, error
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df = normalise_df(df)
    df = ensure_ids(df)

    required = {"og_id", "seq_date", "unique_k_mers_assembly", "k_mers_total", "qv", "error"}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{tsv_path}: missing required columns {missing}")

    # Coerce
    df["unique_k_mers_assembly"] = df["unique_k_mers_assembly"].apply(to_int_or_none)
    df["k_mers_total"]           = df["k_mers_total"].apply(to_int_or_none)
    df["qv"]                     = df["qv"].apply(to_float_or_none)
    df["error"]                  = df["error"].apply(to_float_or_none)

    return df

def load_comp(tsv_path: str) -> pd.DataFrame:
    """
    Expected columns (case/spacing tolerant): sample (or og_id, seq_date), k_mer_set, solid_k_mers, total_k_mers, completeness
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df = normalise_df(df)
    df = ensure_ids(df)

    required = {"og_id", "seq_date", "k_mer_set", "solid_k_mers", "total_k_mers", "completeness"}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{tsv_path}: missing required columns {missing}")

    # Coerce
    df["k_mer_set"]     = df["k_mer_set"].apply(to_text_or_none)
    df["solid_k_mers"]  = df["solid_k_mers"].apply(to_int_or_none)
    df["total_k_mers"]  = df["total_k_mers"].apply(to_int_or_none)
    df["completeness"]  = df["completeness"].apply(to_float_or_none)

    return df

# ----------------------------
# Upserters
# ----------------------------
def upsert_qv(conn, df: pd.DataFrame, src_path: str):
    with conn, conn.cursor() as cur:
        n = 0
        for _, row in df.iterrows():
            params = {
                "og_id":                 to_text_or_none(row.get("og_id")),
                "seq_date":              to_text_or_none(row.get("seq_date")),
                "unique_k_mers_assembly":row.get("unique_k_mers_assembly"),
                "k_mers_total":          row.get("k_mers_total"),
                "qv":                    row.get("qv"),
                "error":                 row.get("error"),
            }
            cur.execute(UPSERT_QV_SQL, params)
            n += 1
        print(f"✅ {src_path}: upserted {n} row(s)")

def upsert_comp(conn, df: pd.DataFrame, src_path: str):
    with conn, conn.cursor() as cur:
        n = 0
        for _, row in df.iterrows():
            params = {
                "og_id":        to_text_or_none(row.get("og_id")),
                "seq_date":     to_text_or_none(row.get("seq_date")),
                "k_mer_set":    to_text_or_none(row.get("k_mer_set")),
                "solid_k_mers": row.get("solid_k_mers"),
                "total_k_mers": row.get("total_k_mers"),
                "completeness": row.get("completeness"),
            }
            cur.execute(UPSERT_COMP_SQL, params)
            n += 1
        print(f"✅ {src_path}: upserted {n} row(s)")

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert Merqury per-sample TSVs (QV and/or completeness) into PostgreSQL."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to key=value config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("--qv", help="Merqury QV stats TSV for a single sample")
    ap.add_argument("--comp", help="Merqury completeness stats TSV for a single sample")
    args = ap.parse_args()

    if not args.qv and not args.comp:
        sys.exit("❌ Provide at least one of --qv or --comp")

    cfg = load_kv_config(args.config)
    db_params = {
        "dbname":   cfg.get("DB_NAME"),
        "user":     cfg.get("DB_USER"),
        "password": cfg.get("DB_PASSWORD"),
        "host":     cfg.get("DB_HOST"),
        "port":     int(cfg.get("DB_PORT", "5432")),
    }

    # Connect once
    try:
        conn = psycopg2.connect(**db_params)
    except Exception as e:
        sys.exit(f"❌ Failed to connect to DB: {e}")

    try:
        if args.qv:
            if not os.path.isfile(args.qv):
                raise FileNotFoundError(f"QV TSV not found: {args.qv}")
            qv_df = load_qv(args.qv)
            upsert_qv(conn, qv_df, args.qv)

        if args.comp:
            if not os.path.isfile(args.comp):
                raise FileNotFoundError(f"Completeness TSV not found: {args.comp}")
            comp_df = load_comp(args.comp)
            upsert_comp(conn, comp_df, args.comp)

    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Error: {e}")
    finally:
        conn.close()
        print("Connection Closed")

if __name__ == "__main__":
    main()
