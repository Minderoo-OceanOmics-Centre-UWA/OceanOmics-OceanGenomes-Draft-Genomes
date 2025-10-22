#!/usr/bin/env python3
import argparse, json, os, re, sys
import psycopg2

# ----------------------------
# Config & helpers
# ----------------------------
def load_kv_config(path):
    cfg = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line: continue
            k, v = line.split("=", 1)
            cfg[k.strip().lower()] = v.strip()
    return cfg

PH = {"ne","n/e","na","n/a","none","null","-",""}  # placeholders → None

def norm_key(s): return re.sub(r"[^a-z0-9_]", "_", s.lower().strip())

def to_int(x):
    if x is None: return None
    s = str(x).strip().replace(",", "")
    if s.lower() in PH: return None
    return int(s) if re.fullmatch(r"[-+]?\d+", s) else None

def to_float(x):
    if x is None: return None
    s = str(x).strip().replace(",", "").replace("%", "")
    if s.lower() in PH: return None
    return float(s) if re.fullmatch(r"[-+]?\d+(\.\d+)?", s) else None

def infer_ids_from_input(inp):
    parts = os.path.basename(inp).split(".")
    og = parts[0] if len(parts) >= 1 else None
    sd = parts[2] if len(parts) >= 3 else None
    return og, sd

def infer_ids_from_path(path):
    base = re.sub(r"\.json$", "", os.path.basename(path), flags=re.I)
    parts = base.split(".")
    og = parts[0] if len(parts) >= 1 else None
    sd = parts[2] if len(parts) >= 3 else None
    return og, sd

# ----------------------------
# Key mapping
# ----------------------------
KEYMAP = {
    "complete":            ["complete_percentage","complete","c"],
    "single_copy":         ["single_copy_percentage","single_percentage","single_copy","single"],
    "multi_copy":          ["multi_copy_percentage","duplicated_percentage","multi_copy","duplicated"],
    "fragmented":          ["fragmented_percentage","fragmented","f"],
    "missing":             ["missing_percentage","missing","m"],
    "n_markers":           ["n_markers","number_of_buscos","total_buscos","lineage_dataset_size"],
    "domain":              ["domain","lineage","lineage_dataset","lineage_name"],
    "number_of_scaffolds": ["number_of_scaffolds","scaffolds"],
    "number_of_contigs":   ["number_of_contigs","contigs"],
    "total_length":        ["total_length","assembly_size","genome_size_bp"],
    "percent_gaps":        ["percent_gaps","gap_percent","n_percent"],
    "scaffold_n50":        ["scaffold_n50","n50_scaffold","scaffold_n50_bp"],
    "contigs_n50":         ["contigs_n50","n50_contig","contig_n50","contig_n50_bp"],
    "internal_stop_codon_count":   ["internal_stop_codon_count"],
    "internal_stop_codon_percent": ["internal_stop_codon_percent"],
}

def pick(norm, keys):
    for k in keys:
        if k in norm and norm[k] is not None:
            return norm[k]
    return None

# ----------------------------
# Parse JSON → row dict
# ----------------------------
def parse_busco_json(path):
    with open(path) as f: data = json.load(f)
    merged = {}
    if isinstance(data.get("metrics"), dict): merged.update(data["metrics"])
    if isinstance(data.get("results"), dict): merged.update(data["results"])
    norm = {norm_key(k): v for k, v in merged.items()}

    row = {}
    for col, keys in KEYMAP.items():
        val = pick(norm, keys)
        if col.endswith("_percent") or col in ("complete","single_copy","multi_copy","fragmented","missing","percent_gaps","scaffold_n50"):
            row[col] = to_float(val)
        elif col.endswith("_count") or col.startswith("n_") or col in ("number_of_scaffolds","number_of_contigs","total_length","contigs_n50"):
            row[col] = to_int(val)
        else:
            row[col] = (val.strip() if isinstance(val,str) and val.strip() else None)

    # infer IDs
    og = sd = None
    inp = (data.get("parameters") or {}).get("in") or (data.get("parameters") or {}).get("out")
    if inp: og, sd = infer_ids_from_input(inp); row["input_file"] = inp
    if not og or not sd:
        og2, sd2 = infer_ids_from_path(path)
        og = og or og2; sd = sd or sd2
    row["og_id"], row["seq_date"] = og, sd
    return row

# ----------------------------
# SQL with new columns
# ----------------------------
UPSERT_SQL = """
INSERT INTO draft_genomes (
    og_id, seq_date, complete, single_copy, multi_copy, fragmented,
    missing, n_markers, domain, number_of_scaffolds, number_of_contigs,
    total_length, percent_gaps, scaffold_n50, contigs_n50,
    internal_stop_codon_count, internal_stop_codon_percent
)
VALUES (
    %(og_id)s, %(seq_date)s, %(complete)s, %(single_copy)s, %(multi_copy)s,
    %(fragmented)s, %(missing)s, %(n_markers)s, %(domain)s,
    %(number_of_scaffolds)s, %(number_of_contigs)s, %(total_length)s,
    %(percent_gaps)s, %(scaffold_n50)s, %(contigs_n50)s,
    %(internal_stop_codon_count)s, %(internal_stop_codon_percent)s
)
ON CONFLICT (og_id, seq_date) DO UPDATE SET
    complete=EXCLUDED.complete, single_copy=EXCLUDED.single_copy,
    multi_copy=EXCLUDED.multi_copy, fragmented=EXCLUDED.fragmented,
    missing=EXCLUDED.missing, n_markers=EXCLUDED.n_markers,
    domain=EXCLUDED.domain, number_of_scaffolds=EXCLUDED.number_of_scaffolds,
    number_of_contigs=EXCLUDED.number_of_contigs, total_length=EXCLUDED.total_length,
    percent_gaps=EXCLUDED.percent_gaps, scaffold_n50=EXCLUDED.scaffold_n50,
    contigs_n50=EXCLUDED.contigs_n50,
    internal_stop_codon_count=EXCLUDED.internal_stop_codon_count,
    internal_stop_codon_percent=EXCLUDED.internal_stop_codon_percent;
"""

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(description="Upsert BUSCO JSON into PostgreSQL.")
    ap.add_argument("-c","--config", required=True)
    ap.add_argument("-f","--file", required=True)
    args = ap.parse_args()
    if not os.path.isfile(args.file): sys.exit(f"❌ No such file: {args.file}")

    row = parse_busco_json(args.file)
    if not row.get("og_id") or not row.get("seq_date"):
        sys.exit("❌ Missing og_id or seq_date")

    print("▶ Upserting:", {k: row.get(k) for k in KEYMAP.keys() if k in row})

    cfg = load_kv_config(args.config)
    try:
        conn = psycopg2.connect(
            dbname=cfg.get("dbname"), user=cfg.get("user"),
            password=cfg.get("password"), host=cfg.get("host"),
            port=int(cfg.get("port", "5432"))
        )
    except Exception as e:
        sys.exit(f"❌ DB connection failed: {e}")

    try:
        with conn, conn.cursor() as cur:
            cur.execute(UPSERT_SQL, row)
            print(f"✅ Upserted ({row['og_id']}, {row['seq_date']}) from {args.file}")
    except Exception as e:
        conn.rollback()
        sys.exit(f"❌ Database error: {e}")
    finally:
        conn.close()

if __name__ == "__main__":
    main()
