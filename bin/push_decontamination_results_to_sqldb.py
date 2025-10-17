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

def to_text_or_none(x):
    return str(x) if (x is not None and not (isinstance(x, float) and np.isnan(x))) else None

def ensure_ids(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure og_id and seq_date exist. Use 'sample' (OG.TECH.DATE) if needed."""
    df = df.copy()
    if "og_id" not in df.columns and "sample" in df.columns:
        # Split "sample" -> og_id, tech, seq_date when present
        split = df["sample"].astype(str).str.split(".", expand=True)
        if split.shape[1] >= 1:
            df["og_id"] = split[0]
        if split.shape[1] >= 3 and "seq_date" not in df.columns:
            df["seq_date"] = split[2]
    return df

# ----------------------------
# SQL
# ----------------------------
UPSERT_TIARA_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date,
    num_contigs_mitochondrion, num_contigs_plastid, num_contigs_prokarya,
    bp_mitochondrion, bp_plastid, bp_prokarya
)
VALUES (
    %(og_id)s, %(seq_date)s,
    %(num_contigs_mitochondrion)s, %(num_contigs_plastid)s, %(num_contigs_prokarya)s,
    %(bp_mitochondrion)s, %(bp_plastid)s, %(bp_prokarya)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    num_contigs_mitochondrion = EXCLUDED.num_contigs_mitochondrion,
    num_contigs_plastid      = EXCLUDED.num_contigs_plastid,
    num_contigs_prokarya     = EXCLUDED.num_contigs_prokarya,
    bp_mitochondrion         = EXCLUDED.bp_mitochondrion,
    bp_plastid               = EXCLUDED.bp_plastid,
    bp_prokarya              = EXCLUDED.bp_prokarya;
"""

UPSERT_NCBI_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, num_contigs
)
VALUES (
    %(og_id)s, %(seq_date)s, %(num_contigs)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    num_contigs = EXCLUDED.num_contigs;
"""

# ----------------------------
# Tiara processing
# ----------------------------
def load_tiara_pivoted(tsv_path: str) -> pd.DataFrame:
    """
    Input columns expected:
      sample, category, num_contigs, bp
    Produces one row per sample with:
      num_contigs_mitochondrion, num_contigs_plastid, num_contigs_prokarya,
      bp_mitochondrion, bp_plastid, bp_prokarya
    """
    raw = pd.read_csv(tsv_path, sep="\t")
    raw = normalise_df(raw)

    required = {"sample", "category", "num_contigs", "bp"}
    if not required.issubset(raw.columns):
        missing = required - set(raw.columns)
        raise ValueError(f"Tiara TSV missing columns: {sorted(missing)}")

    # Pivot
    piv = raw.pivot(index="sample", columns="category", values=["num_contigs", "bp"])
    piv = piv.fillna(0)

    # Flatten names like ('num_contigs', 'mitochondrion') -> 'num_contigs_mitochondrion'
    piv.columns = [f"{a}_{b}".lower() for a, b in piv.columns.to_flat_index()]
    piv = piv.reset_index()  # restore 'sample' as column
    piv = normalise_df(piv)

    # Ensure IDs
    piv = ensure_ids(piv)

    # Keep only expected tiara output columns if present
    for col in ["num_contigs_mitochondrion", "num_contigs_plastid", "num_contigs_prokarya",
                "bp_mitochondrion", "bp_plastid", "bp_prokarya"]:
        if col not in piv.columns:
            piv[col] = 0  # if a category is absent, treat counts as zero

    # Coerce numeric
    for col in ["num_contigs_mitochondrion", "num_contigs_plastid", "num_contigs_prokarya",
                "bp_mitochondrion", "bp_plastid", "bp_prokarya"]:
        piv[col] = piv[col].apply(to_int_or_none)

    return piv

def upsert_tiara(conn, df: pd.DataFrame, src_path: str):
    df = normalise_df(df)
    df = ensure_ids(df)

    if "og_id" not in df.columns or "seq_date" not in df.columns:
        raise ValueError(f"{src_path}: requires og_id (or sample) and seq_date")

    with conn, conn.cursor() as cur:
        count = 0
        for _, row in df.iterrows():
            params = {
                "og_id": to_text_or_none(row.get("og_id")),
                "seq_date": to_text_or_none(row.get("seq_date")),
                "num_contigs_mitochondrion": to_int_or_none(row.get("num_contigs_mitochondrion")),
                "num_contigs_plastid":      to_int_or_none(row.get("num_contigs_plastid")),
                "num_contigs_prokarya":     to_int_or_none(row.get("num_contigs_prokarya")),
                "bp_mitochondrion":         to_int_or_none(row.get("bp_mitochondrion")),
                "bp_plastid":               to_int_or_none(row.get("bp_plastid")),
                "bp_prokarya":              to_int_or_none(row.get("bp_prokarya")),
            }
            cur.execute(UPSERT_TIARA_SQL, params)
            count += 1
        print(f"✅ {src_path}: upserted {count} row(s)")

# ----------------------------
# NCBI processing
# ----------------------------
def load_ncbi(tsv_path: str) -> pd.DataFrame:
    """
    Expected columns:
      sample (or og_id), num_contigs  [and optionally seq_date]
    If sample is present as OG.TECH.SEQDATE, split to get seq_date.
    """
    df = pd.read_csv(tsv_path, sep="\t")
    df = normalise_df(df)
    df = ensure_ids(df)

    # Validate required fields
    if "og_id" not in df.columns:
        raise ValueError(f"{tsv_path}: requires 'og_id' or 'sample' column")
    if "seq_date" not in df.columns:
        raise ValueError(f"{tsv_path}: requires 'seq_date' (can come from splitting 'sample')")
    if "num_contigs" not in df.columns:
        raise ValueError(f"{tsv_path}: requires 'num_contigs'")

    # Coerce
    df["num_contigs"] = df["num_contigs"].apply(to_int_or_none)
    return df

def upsert_ncbi(conn, df: pd.DataFrame, src_path: str):
    with conn, conn.cursor() as cur:
        count = 0
        for _, row in df.iterrows():
            params = {
                "og_id":    to_text_or_none(row.get("og_id")),
                "seq_date": to_text_or_none(row.get("seq_date")),
                "num_contigs": to_int_or_none(row.get("num_contigs")),
            }
            cur.execute(UPSERT_NCBI_SQL, params)
            count += 1
        print(f"✅ {src_path}: upserted {count} row(s)")

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert Tiara and/or NCBI per-sample TSVs into PostgreSQL."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to key=value config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("--tiara", help="Tiara filter report TSV for a single sample")
    ap.add_argument("--ncbi", help="NCBI >=500bp contig count TSV for a single sample")
    args = ap.parse_args()

    if not args.tiara and not args.ncbi:
        sys.exit("❌ Provide at least one of --tiara or --ncbi")

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
        if args.tiara:
            if not os.path.isfile(args.tiara):
                raise FileNotFoundError(f"Tiara TSV not found: {args.tiara}")
            tiara_df = load_tiara_pivoted(args.tiara)
            upsert_tiara(conn, tiara_df, args.tiara)

        if args.ncbi:
            if not os.path.isfile(args.ncbi):
                raise FileNotFoundError(f"NCBI TSV not found: {args.ncbi}")
            ncbi_df = load_ncbi(args.ncbi)
            upsert_ncbi(conn, ncbi_df, args.ncbi)

    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Error: {e}")
    finally:
        conn.close()
        print("Connection Closed")

if __name__ == "__main__":
    main()
