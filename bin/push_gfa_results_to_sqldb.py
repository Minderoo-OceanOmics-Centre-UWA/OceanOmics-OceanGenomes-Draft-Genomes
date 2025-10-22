#!/usr/bin/env python3
# Usage:
#   python push_gfastats_summary_to_sqldb.py -c ../configfile.txt -f /path/to/assembly_summary.txt
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
            cfg[k.strip().lower()] = v.strip()
    return cfg

def parse_int(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip().replace(",", "")
    if s == "" or s.lower() == "nan":
        return None
    try:
        # some values can be floats but represent ints (e.g. "10859.0")
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

def derive_ids_from_filename(path: str):
    base = os.path.basename(path)
    parts = base.split(".")
    og_id = parts[0] if len(parts) > 0 else None
    seq_date = parts[2] if len(parts) > 2 else None
    return og_id, seq_date

# ----------------------------
# Parser for the text summary
# ----------------------------
def parse_gfastats_summary_txt(path: str) -> dict:
    """
    Expects content like:
      '# scaffolds: 125966'
      'Total scaffold length: 482318804'
      'Scaffold N50: 10859'
      'Largest scaffold: 80979'
      '# contigs: 125966'
      'Contig N50: 10859'
      'GC content %: 45.31'  (or 'GC content: 45.31 %')
    """
    # Prepare result container
    res = {
        "num_contigs": None,
        "contig_n50": None,
        "num_scaffolds": None,
        "scaffold_n50": None,
        "largest_scaffold": None,
        "total_scaffold_length": None,
        "gc_content_percent": None,
    }

    # Regexes (case-insensitive, tolerant to spaces)
    RX_INT = r"(-?\d+(?:\.\d+)?)"
    patterns = [
        ("num_scaffolds",          re.compile(r"^\s*#\s*scaffolds\s*:\s*" + RX_INT + r"\s*$", re.I)),
        ("total_scaffold_length",  re.compile(r"^\s*total\s+scaffold\s+length\s*:\s*" + RX_INT + r"\s*$", re.I)),
        ("scaffold_n50",           re.compile(r"^\s*scaffold\s+N50\s*:\s*" + RX_INT + r"\s*$", re.I)),
        ("largest_scaffold",       re.compile(r"^\s*largest\s+scaffold\s*:\s*" + RX_INT + r"\s*$", re.I)),
        ("num_contigs",            re.compile(r"^\s*#\s*contigs\s*:\s*" + RX_INT + r"\s*$", re.I)),
        ("contig_n50",             re.compile(r"^\s*contig\s+N50\s*:\s*" + RX_INT + r"\s*$", re.I)),
        # GC content may appear as 'GC content %: 45.31' or 'GC content: 45.31 %'
        ("gc_content_percent",     re.compile(r"^\s*GC\s*content.*?:\s*([0-9]+(?:\.[0-9]+)?)", re.I)),
    ]

    try:
        with open(path, "r") as f:
            for raw in f:
                line = raw.strip()
                for key, rx in patterns:
                    m = rx.match(line)
                    if m:
                        val = m.group(1)
                        if key == "gc_content_percent":
                            res[key] = parse_float(val)
                        else:
                            # N50/length/count fields first parsed as int if possible
                            ival = parse_int(val)
                            res[key] = ival if ival is not None else parse_float(val)
                        break  # stop checking other regexes for this line
    except Exception as e:
        raise RuntimeError(f"Failed parsing {path}: {e}")

    return res

# ----------------------------
# SQL
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date,
    gfa_num_contigs, gfa_contig_n50, 
    gfa_num_scaffolds, gfa_scaffold_n50, 
    gfa_largest_scaffold, 
    gfa_total_scaffold_length, 
    gfa_gc_content_percent
)
VALUES (
    %(og_id)s, %(seq_date)s,
    %(num_contigs)s, %(contig_n50)s, 
    %(num_scaffolds)s, %(scaffold_n50)s, 
    %(largest_scaffold)s, 
    %(total_scaffold_length)s, 
    %(gc_content_percent)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    gfa_num_contigs = EXCLUDED.gfa_num_contigs,
    gfa_contig_n50 = EXCLUDED.gfa_contig_n50,
    gfa_num_scaffolds = EXCLUDED.gfa_num_scaffolds,
    gfa_scaffold_n50 = EXCLUDED.gfa_scaffold_n50,
    gfa_largest_scaffold = EXCLUDED.gfa_largest_scaffold,
    gfa_total_scaffold_length = EXCLUDED.gfa_total_scaffold_length,
    gfa_gc_content_percent = EXCLUDED.gfa_gc_content_percent;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Parse a gfastats summary text file and upsert key stats into PostgreSQL."
    )
    ap.add_argument("-c", "--config", required=True,
                    help="Path to key=value config with dbname,user,password,host,port")
    ap.add_argument("-f", "--file", required=True,
                    help="Path to gfastats summary text file (e.g., *assembly_summary.txt)")
    args = ap.parse_args()

    if not os.path.isfile(args.file):
        sys.exit(f"❌ Input file not found: {args.file}")

    # Parse the summary text file
    stats = parse_gfastats_summary_txt(args.file)

    # Derive identifiers from basename (OGID.SEQDATE... pattern; we take first dot-part, then underscores)
    og_id, seq_date = derive_ids_from_filename(args.file)
    if not og_id or not seq_date:
        sys.exit("❌ Could not derive og_id/seq_date from filename. Expecting basename like 'OGID_SEQDATE...*.txt'")

    # Build the parameter dict for SQL
    params = {
        "og_id":                         to_text_or_none(og_id),
        "seq_date":                      to_text_or_none(seq_date),
        "num_contigs":                   parse_int(stats.get("num_contigs")),
        "contig_n50":                    parse_int(stats.get("contig_n50")),
        "num_scaffolds":                 parse_int(stats.get("num_scaffolds")),
        "scaffold_n50":                  parse_int(stats.get("scaffold_n50")),
        "largest_scaffold":              parse_int(stats.get("largest_scaffold")),
        "total_scaffold_length":         parse_int(stats.get("total_scaffold_length")),
        "gc_content_percent":            parse_float(stats.get("gc_content_percent")),
    }

    # Load DB config + connect
    cfg = load_kv_config(args.config)
    db_params = {
        "dbname":   cfg.get('dbname'),
        "user":     cfg.get('user'),
        "password": cfg.get('password'),
        "host":     cfg.get('host'),
        "port":     int(cfg.get('port')),
    }

    try:
        conn = psycopg2.connect(**db_params)
    except Exception as e:
        sys.exit(f"❌ Failed to connect to DB: {e}")

    try:
        with conn, conn.cursor() as cur:
            cur.execute(UPSERT_SQL, params)
            print(f"✅ Upserted gfastats summary for {og_id} / {seq_date} from {os.path.basename(args.file)}")
    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        conn.close()
        print("Connection Closed")

if __name__ == "__main__":
    main()
