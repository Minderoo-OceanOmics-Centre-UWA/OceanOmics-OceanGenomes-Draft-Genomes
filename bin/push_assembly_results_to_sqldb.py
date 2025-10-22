#!/usr/bin/env python3
import sys
import os
import re
import argparse
import psycopg2

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

def load_kv_config(path: str) -> dict:
    cfg = {}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            cfg[k.strip()] = v.strip()
    return cfg

def parse_number(s: str) -> float:
    return float(s.replace("bp", "").replace(",", "").strip())

def parse_percent(s: str) -> float:
    return float(s.replace("%", "").strip())

def extract_token_triplet(s: str):
    toks = os.path.basename(s).split(".")
    if len(toks) < 3:
        return (None, None, None)
    og_id, tech, seq_date_raw = toks[0], toks[1], toks[2]
    if not og_id or not seq_date_raw:
        return (None, None, None)
    return (og_id, tech, seq_date_raw)

def clean_seq_date(seq_date_raw: str) -> str:
    # e.g., '250724_genomescope2' -> '250724'
    return seq_date_raw.split("_", 1)[0] if "_" in seq_date_raw else seq_date_raw

def parse_genomescope_report(path: str):
    with open(path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f]

    # Prefer 'name prefix', fall back to 'input file'
    name_prefix = None
    input_file = None
    for ln in lines:
        lower = ln.lower()
        if lower.startswith("name prefix"):
            m = re.search(r"=\s*(.+)$", ln)
            if m: name_prefix = m.group(1).strip()
        elif lower.startswith("input file"):
            m = re.search(r"=\s*(.+)$", ln)
            if m: input_file = m.group(1).strip()

    cand = name_prefix or input_file
    if cand:
        og, tech, seq_date_raw = extract_token_triplet(cand)
        if og and seq_date_raw:
            og_id = og
            seq_date = clean_seq_date(seq_date_raw)
    else:
        og_id = seq_date = None

    if not og_id or not seq_date:
        sys.exit("❌ Could not infer og_id/seq_date from report (checked 'name prefix' then 'input file'). "
                 "Expected something like OGID.TECH.SEQDATE[...].")

    # We’ll split on 2+ spaces to get columns: [label, min, max]
    def split_cols(line: str):
        parts = re.split(r"\s{2,}", line.strip())
        if len(parts) < 3:
            return None
        # Some headers/blank lines may sneak through
        return parts[0], parts[-2], parts[-1]  # label, min, max

    vals = {
        "homozygosity": None,
        "heterozygosity": None,
        "genomesize": None,
        "repeatsize": None,
        "uniquesize": None,
        "modelfit": None,
        "errorrate": None,
    }

    # Label keys we care about (match startswith to allow "(aa)" etc)
    targets = {
        "homozygous": ("homozygosity", parse_percent),
        "heterozygous": ("heterozygosity", parse_percent),
        "genome haploid length": ("genomesize", parse_number),
        "genome repeat length": ("repeatsize", parse_number),
        "genome unique length": ("uniquesize", parse_number),
        "model fit": ("modelfit", parse_percent),
        "read error rate": ("errorrate", parse_percent),
    }

    for ln in lines:
        cols = split_cols(ln)
        if not cols: 
            continue
        label, _min_str, max_str = cols
        key = next((k for k in targets if label.lower().startswith(k)), None)
        if not key:
            continue
        out_key, parser = targets[key]
        try:
            vals[out_key] = parser(max_str)  # <-- use MAX column
        except Exception as e:
            sys.exit(f"❌ Failed to parse line: '{ln}' -> {e}")

    missing = [k for k, v in vals.items() if v is None]
    if missing:
        sys.exit(f"❌ Missing fields in report: {missing}")

    return {"og_id": og_id, "seq_date": seq_date, **vals}

def main():
    ap = argparse.ArgumentParser(description="Push GenomeScope (plaintext) results to PostgreSQL, using MAX column values.")
    ap.add_argument("-c", "--config", required=True, help="Config with dbname, user, password, host, port")
    ap.add_argument("-f", "--file", required=True, help="GenomeScope plaintext report")
    args = ap.parse_args()

    cfg = load_kv_config(args.config)
    for key in ["dbname", "user", "password", "host", "port"]:
        if not cfg.get(key):
            sys.exit(f"❌ Missing DB config key: {key}")

    db_params = {
        "dbname": cfg["dbname"],
        "user": cfg["user"],
        "password": cfg["password"],
        "host": cfg["host"],
        "port": int(cfg["port"]),
        "connect_timeout": 10,
    }
    if cfg.get("sslmode"):
        db_params["sslmode"] = cfg["sslmode"]

    rec = parse_genomescope_report(args.file)

    try:
        conn = psycopg2.connect(**db_params)
        cur = conn.cursor()
        cur.execute(
            UPSERT_SQL,
            (
                rec["og_id"],
                str(rec["seq_date"]),
                float(rec["homozygosity"]),
                float(rec["heterozygosity"]),
                int(round(rec["genomesize"])),
                int(round(rec["repeatsize"])),
                int(round(rec["uniquesize"])),
                float(rec["modelfit"]),
                float(rec["errorrate"]),
            ),
        )
        conn.commit()
        print(f"✅ Upserted og_id={rec['og_id']} seq_date={rec['seq_date']} (MAX values) from {args.file}")
    except Exception as e:
        if 'conn' in locals():
            conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        if 'cur' in locals():
            cur.close()
        if 'conn' in locals():
            conn.close()

if __name__ == "__main__":
    main()
