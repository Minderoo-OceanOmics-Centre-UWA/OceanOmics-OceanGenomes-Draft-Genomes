#!/usr/bin/env python3
import sys
import os
import argparse
import json
import psycopg2

# ----------------------------
# Helpers
# ----------------------------
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

def to_int(x):
    if x is None:
        return 0
    try:
        return int(str(x).replace(",", ""))
    except Exception:
        try:
            return int(float(x))
        except Exception:
            return 0

def to_float(x):
    if x is None:
        return 0.0
    try:
        return float(str(x).replace("%", "").replace(",", ""))
    except Exception:
        return 0.0

def jget(d, *keys, default=None):
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur

def infer_from_filename(path: str):
    """
    Expect filenames like: OGID.<anything>.RUNTOKEN.<rest>.json
    Where RUNTOKEN is formatted as: mach_seqdate_initial
    Returns (og_id, mach, seq_date, initial)
    """
    toks = os.path.basename(path).split(".")
    og_id = toks[0] if len(toks) >= 1 and toks[0] else None
    run_token = toks[2] if len(toks) >= 3 and toks[2] else None

    if not og_id or not run_token:
        return None, None, None, None

    parts = run_token.split("_", 2)  # mach, seq_date, initial (initial may contain underscores -> use maxsplit=2)
    if len(parts) != 3 or not parts[0] or not parts[1] or not parts[2]:
        return og_id, None, None, None

    mach, seq_date, initial = parts[0], parts[1], parts[2]
    return og_id, mach, seq_date, initial

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
    ap = argparse.ArgumentParser(description="Upsert fastp JSON into PostgreSQL.")
    ap.add_argument("-c", "--config", required=True,
                    help="Path to config with DB_NAME, DB_USER, DB_PASSWORD, DB_HOST, DB_PORT")
    ap.add_argument("-f", "--file", required=True,
                    help="Path to fastp JSON (single sample)")
    args = ap.parse_args()

    # Load DB config
    cfg = load_kv_config(args.config)
    db_params = {
        "dbname":   cfg.get('dbname'),
        "user":     cfg.get('user'),
        "password": cfg.get('password'),
        "host":     cfg.get('host'),
        "port":     int(cfg.get('port')),
    }

    # Read JSON
    if not os.path.isfile(args.file):
        sys.exit(f"❌ Input file not found: {args.file}")
    try:
        with open(args.file, "r") as fh:
            data = json.load(fh)
    except Exception as e:
        sys.exit(f"❌ Failed to parse JSON: {e}")

    # Infer identifiers: og_id from token[0]; run token at token[2] as mach_seqdate_initial
    og_id, mach, seq_date, initial = infer_from_filename(args.file)
    if not og_id:
        sys.exit("❌ Could not infer og_id from filename. Expected: OGID.<...>.mach_seqdate_initial.<...>.json")
    if not (mach and seq_date and initial):
        sys.exit("❌ Could not parse run token. Expected 3rd dot-delimited token formatted as 'mach_seqdate_initial' "
                 "(e.g., 'ilmn_2025-10-19_yes').")

    # Extract sections
    fr = jget(data, "filtering_result", default={}) or {}
    bf = jget(data, "summary", "before_filtering", default={}) or {}
    af = jget(data, "summary", "after_filtering", default={}) or {}

    # Normalise key naming (fastp sometimes uses 'too_many_N_reads')
    tmn = fr.get("too_many_n_reads", fr.get("too_many_N_reads", 0))

    # Map to SQL params
    params = {
        "og_id":    str(og_id),
        "mach":     str(mach),
        "seq_date": str(seq_date),
        "initial":  str(initial),

        # filtering_result
        "passed_filter_reads": int(to_int(fr.get("passed_filter_reads"))),
        "low_quality_reads":   int(to_int(fr.get("low_quality_reads"))),
        "too_many_n_reads":    int(to_int(tmn)),
        "too_short_reads":     int(to_int(fr.get("too_short_reads"))),
        "too_long_reads":      int(to_int(fr.get("too_long_reads"))),

        # before_filtering -> raw_*
        "raw_total_reads":       int(to_int(bf.get("total_reads"))),
        "raw_total_bases":       int(to_int(bf.get("total_bases"))),
        "raw_q20_bases":         int(to_int(bf.get("q20_bases"))),
        "raw_q30_bases":         int(to_int(bf.get("q30_bases"))),
        "raw_q20_rate":          float(to_float(bf.get("q20_rate"))),
        "raw_q30_rate":          float(to_float(bf.get("q30_rate"))),
        "raw_read1_mean_length": int(to_int(bf.get("read1_mean_length"))),
        "raw_read2_mean_length": int(to_int(bf.get("read2_mean_length"))),
        "raw_gc_content":        float(to_float(bf.get("gc_content"))),

        # after_filtering -> non-raw
        "total_reads":           int(to_int(af.get("total_reads"))),
        "total_bases":           int(to_int(af.get("total_bases"))),
        "q20_bases":             int(to_int(af.get("q20_bases"))),
        "q30_bases":             int(to_int(af.get("q30_bases"))),
        "q20_rate":              float(to_float(af.get("q20_rate"))),
        "q30_rate":              float(to_float(af.get("q30_rate"))),
        "read1_mean_length":     int(to_int(af.get("read1_mean_length"))),
        "read2_mean_length":     int(to_int(af.get("read2_mean_length"))),
        "gc_content":            float(to_float(af.get("gc_content"))),
    }

    # DB upsert
    try:
        conn = psycopg2.connect(**db_params)
    except Exception as e:
        sys.exit(f"❌ Failed to connect to DB: {e}")

    try:
        with conn, conn.cursor() as cur:
            cur.execute(UPSERT_SQL, params)
            print(f"✅ Upserted og_id={params['og_id']} run(mach_seq_date_initial)={mach}_{seq_date}_{initial} from {args.file}")
    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    main()
