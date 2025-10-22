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
    if x is None or (isinstance(x, float) and np.isnan(x)) or x == "":
        return None
    try:
        return int(str(x).replace(",", ""))
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return None

def to_text_or_none(x):
    return str(x) if (x is not None and not (isinstance(x, float) and np.isnan(x))) else None

def infer_ids_from_path(path: str):
    """
    Expect filenames like: OGID.TECH.SEQDATE.*  -> return (og_id, seq_date)
    """
    toks = os.path.basename(path).split(".")
    og_id = toks[0] if len(toks) >= 1 else None
    seq_date = toks[2] if len(toks) >= 3 else None
    return og_id, seq_date

def ensure_ids(df: pd.DataFrame, fallback_path: str = None) -> pd.DataFrame:
    """Ensure og_id and seq_date exist. Use 'sample' (OG.TECH.DATE) if needed; else infer from file path."""
    df = df.copy()
    if "og_id" not in df.columns and "sample" in df.columns:
        split = df["sample"].astype(str).str.split(".", expand=True)
        if split.shape[1] >= 1:
            df["og_id"] = split[0]
        if split.shape[1] >= 3 and "seq_date" not in df.columns:
            df["seq_date"] = split[2]
    if (("og_id" not in df.columns) or ("seq_date" not in df.columns)) and fallback_path:
        og_id, seq_date = infer_ids_from_path(fallback_path)
        if "og_id" not in df.columns and og_id:
            df["og_id"] = og_id
        if "seq_date" not in df.columns and seq_date:
            df["seq_date"] = seq_date
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

UPSERT_FILTER_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date,
    num_contigs_exclude, bp_exclude,
    num_contigs_trim,    bp_trim,
    num_contigs_review,  bp_review
)
VALUES (%(og_id)s, %(seq_date)s, %(num_contigs_exclude)s, %(bp_exclude)s,
        %(num_contigs_trim)s, %(bp_trim)s, %(num_contigs_review)s, %(bp_review)s)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    num_contigs_exclude = EXCLUDED.num_contigs_exclude,
    bp_exclude    = EXCLUDED.bp_exclude,
    num_contigs_trim    = EXCLUDED.num_contigs_trim,
    bp_trim       = EXCLUDED.bp_trim,
    num_contigs_review  = EXCLUDED.num_contigs_review,
    bp_review     = EXCLUDED.bp_review;
"""

UPSERT_LT500_SQL = """
INSERT INTO draft_genomes (og_id, seq_date, num_contigs)
VALUES (%(og_id)s, %(seq_date)s, %(num_contigs)s)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    num_contigs = EXCLUDED.num_contigs;
"""

# ----------------------------
# Tiara processing
# ----------------------------
def load_tiara_pivoted(tsv_path: str) -> pd.DataFrame:
    """
    Reads a Tiara summary TSV like:
      Category <tab> num_contigs <tab> bp
      Mitochondrion  5            40908
      Plastid        2            3076
      Prokarya       0            0

    Returns a single-row DataFrame with:
      og_id, seq_date,
      num_contigs_mitochondrion, num_contigs_plastid, num_contigs_prokarya,
      bp_mitochondrion, bp_plastid, bp_prokarya
    """
    # Read + normalise columns
    raw = pd.read_csv(tsv_path, sep="\t", dtype=str)
    raw = normalise_df(raw)

    # Required columns
    required = {"category", "num_contigs", "bp"}
    if not required.issubset(raw.columns):
        missing = required - set(raw.columns)
        raise ValueError(f"Tiara TSV missing columns: {sorted(missing)}")

    # Clean values
    raw["category"] = raw["category"].astype(str).str.strip().str.lower()
    # Allow some aliasing just in case different labels show up
    alias = {
        "mitochondrion": "mitochondrion",
        "mitochondria": "mitochondrion",
        "mito": "mitochondrion",
        "plastid": "plastid",
        "chloroplast": "plastid",
        "prokarya": "prokarya",
        "prokaryote": "prokarya",
        "bacteria": "prokarya",
        "bacterial": "prokarya",
    }
    raw["category"] = raw["category"].map(lambda s: alias.get(s, s))

    # Coerce numerics (tolerant)
    def _to_int(x):
        try:
            return int(str(x).replace(",", "").strip())
        except Exception:
            return 0
    raw["num_contigs"] = raw["num_contigs"].apply(_to_int)
    raw["bp"] = raw["bp"].apply(_to_int)

    # Aggregate
    agg = raw.groupby("category", as_index=False).agg({"num_contigs": "sum", "bp": "sum"})

    # Build output dict with zeros as defaults
    out = {
        "num_contigs_mitochondrion": 0,
        "num_contigs_plastid": 0,
        "num_contigs_prokarya": 0,
        "bp_mitochondrion": 0,
        "bp_plastid": 0,
        "bp_prokarya": 0,
    }
    for _, r in agg.iterrows():
        cat = r["category"]
        if cat == "mitochondrion":
            out["num_contigs_mitochondrion"] = int(r["num_contigs"])
            out["bp_mitochondrion"] = int(r["bp"])
        elif cat == "plastid":
            out["num_contigs_plastid"] = int(r["num_contigs"])
            out["bp_plastid"] = int(r["bp"])
        elif cat == "prokarya":
            out["num_contigs_prokarya"] = int(r["num_contigs"])
            out["bp_prokarya"] = int(r["bp"])

    # IDs from filename OGID.TECH.SEQDATE.*
    og_id, seq_date = infer_ids_from_path(tsv_path)
    if not og_id or not seq_date:
        raise ValueError(f"{tsv_path}: cannot infer og_id/seq_date from filename; expected OGID.TECH.SEQDATE.*")

    row = {"og_id": og_id, "seq_date": seq_date, **out}
    return pd.DataFrame([row], columns=[
        "og_id", "seq_date",
        "num_contigs_mitochondrion", "num_contigs_plastid", "num_contigs_prokarya",
        "bp_mitochondrion", "bp_plastid", "bp_prokarya"
    ])


def upsert_tiara(conn, df: pd.DataFrame, src_path: str):
    if df.empty:
        print(f"ℹ️ {src_path}: no Tiara rows found")
        return
    if df["og_id"].isna().any() or df["seq_date"].isna().any():
        raise ValueError(f"{src_path}: cannot infer og_id/seq_date from filename; expected OGID.TECH.SEQDATE.*")
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
        print(f"✅ {src_path}: upserted {count} Tiara row(s)")

# ----------------------------
# Filter processing (EXCLUDE/TRIM/REVIEW)
# ----------------------------
def parse_filter_report(path: str):
    """
    Lines like:
      EXCLUDE 70 99831
      TRIM 0 0
      REVIEW 0 0
    Returns dict with counts+bp and (og_id, seq_date) from filename.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Filter report not found: {path}")

    og_id, seq_date = infer_ids_from_path(path)
    vals = {
        "num_contigs_exclude": 0, "bp_exclude": 0,
        "num_contigs_trim":    0, "bp_trim":   0,
        "num_contigs_review":  0, "bp_review": 0,
    }
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            m = re.match(r"^(EXCLUDE|TRIM|REVIEW)\s+(\d+)\s+(\d+)$", line, re.IGNORECASE)
            if not m:
                continue
            tag = m.group(1).upper()
            cnt = int(m.group(2))
            bp  = int(m.group(3))
            if tag == "EXCLUDE":
                vals["num_contigs_exclude"] = cnt
                vals["bp_exclude"]    = bp
            elif tag == "TRIM":
                vals["num_contigs_trim"]    = cnt
                vals["bp_trim"]       = bp
            elif tag == "REVIEW":
                vals["num_contigs_review"]  = cnt
                vals["bp_review"]     = bp

    vals["og_id"] = og_id
    vals["seq_date"] = seq_date
    if not og_id or not seq_date:
        raise ValueError(f"{path}: cannot infer og_id/seq_date from filename; expected OGID.TECH.SEQDATE.*")
    return vals

def upsert_filter(conn, rec: dict, src_path: str):
    with conn, conn.cursor() as cur:
        cur.execute(UPSERT_FILTER_SQL, rec)
    print(f"✅ {src_path}: upserted filter stats")

# ----------------------------
# <500 bp contigs processing
# ----------------------------
def parse_lt500(path: str):
    """
    Line like:
      Number of contigs less than 500bp: 44345053
    Returns dict {og_id, seq_date, num_contigs}
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"<500bp file not found: {path}")

    og_id, seq_date = infer_ids_from_path(path)
    n = None
    with open(path, "r") as f:
        for line in f:
            m = re.search(r"(\d+)\s*$", line.strip())
            if m:
                n = int(m.group(1))
                break
    if n is None:
        raise ValueError(f"{path}: could not parse 'contigs <500bp' count")
    if not og_id or not seq_date:
        raise ValueError(f"{path}: cannot infer og_id/seq_date from filename; expected OGID.TECH.SEQDATE.*")
    return {"og_id": og_id, "seq_date": seq_date, "num_contigs": n}

def upsert_lt500(conn, rec: dict, src_path: str):
    with conn, conn.cursor() as cur:
        cur.execute(UPSERT_LT500_SQL, rec)
    print(f"✅ {src_path}: upserted num_contigs={rec['num_contigs']}")

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Upsert decontamination-related results (Tiara/Filter/<500bp) into PostgreSQL."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to key=value config (supports dbname|user|password|host|port OR DB_* variants)")
    ap.add_argument("--tiara", help="Tiara filter summary TSV (Category/num_contigs/bp)")
    ap.add_argument("--filter", help="Filter report text (EXCLUDE/TRIM/REVIEW lines)")
    ap.add_argument("--contigs", help="Text with 'Number of contigs less than 500bp: <N>'")
    args = ap.parse_args()

    if not (args.tiara or args.filter or args.contigs):
        sys.exit("❌ Provide at least one of --tiara, --filter, --contigs")

    cfg = load_kv_config(args.config)
    # support lowercase or DB_* keys
    db_params = {
        "dbname":   cfg.get("dbname")   or cfg.get("DB_NAME"),
        "user":     cfg.get("user")     or cfg.get("DB_USER"),
        "password": cfg.get("password") or cfg.get("DB_PASSWORD"),
        "host":     cfg.get("host")     or cfg.get("DB_HOST"),
        "port":     int(cfg.get("port") or cfg.get("DB_PORT", "5432")),
        "connect_timeout": 10,
    }
    for k in ["dbname", "user", "password", "host"]:
        if not db_params.get(k):
            sys.exit(f"❌ Missing DB config key: {k}")

    if cfg.get("sslmode") or cfg.get("DB_SSLMODE"):
        db_params["sslmode"] = cfg.get("sslmode") or cfg.get("DB_SSLMODE")

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

        if args.filter:
            if not os.path.isfile(args.filter):
                raise FileNotFoundError(f"Filter report not found: {args.filter}")
            filt = parse_filter_report(args.filter)
            upsert_filter(conn, filt, args.filter)

        if args.contigs:
            if not os.path.isfile(args.contigs):
                raise FileNotFoundError(f"<500bp file not found: {args.contigs}")
            rec = parse_lt500(args.contigs)
            upsert_lt500(conn, rec, args.contigs)

    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Error: {e}")
    finally:
        conn.close()
        print("Connection Closed")

if __name__ == "__main__":
    main()
