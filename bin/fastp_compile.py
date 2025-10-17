#!/usr/bin/env python3
import sys
import os
import json
import pandas as pd
from datetime import datetime

version = "v0.1.0"

def log_write(msg, stream=sys.stdout):
    ts = datetime.now().ctime()
    stream.write(f"[{ts}] {msg}\n")
    stream.flush()

def process_json(json_path: str):
    # Load JSON
    with open(json_path, "r") as f:
        myData = json.load(f)

    # SAMPLE and RUN from basename split by "."
    base = os.path.basename(json_path)
    parts = base.split(".")
    SAMPLE = parts[0] if len(parts) >= 1 else ""
    RUN = parts[2] if len(parts) >= 3 else ""

    log_write(f"Input: {json_path} -> {SAMPLE} {RUN}")

    # filtering_result
    filtering_result = myData.get("filtering_result", {})
    fastp_filt = (
        pd.DataFrame([filtering_result])
        if isinstance(filtering_result, dict)
        else pd.DataFrame(filtering_result)
    )
    if fastp_filt.empty:
        fastp_filt = pd.DataFrame([{}])
    fastp_filt.insert(0, "run", RUN)
    fastp_filt.insert(0, "sample", SAMPLE)

    # before_filtering
    before = myData.get("summary", {}).get("before_filtering", {})
    raw_df = pd.DataFrame([before]) if isinstance(before, dict) else pd.DataFrame(before)
    if raw_df.empty:
        raw_df = pd.DataFrame([{}])
    raw_df = raw_df.add_prefix("raw.")
    raw_df.insert(0, "sample", SAMPLE)

    # after_filtering
    after = myData.get("summary", {}).get("after_filtering", {})
    filt_df = pd.DataFrame([after]) if isinstance(after, dict) else pd.DataFrame(after)
    if filt_df.empty:
        filt_df = pd.DataFrame([{}])
    filt_df.insert(0, "sample", SAMPLE)

    # join
    fastp = fastp_filt.merge(raw_df, on="sample", how="left")
    fastp = fastp.merge(filt_df, on="sample", how="left")

    # write output
    tsv_path = f"{json_path}.tsv"
    fastp.to_csv(tsv_path, sep="\t", index=False)
    log_write(f"Output: {tsv_path}")

def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: fastp_compile.py <file1.json> [file2.json ...]\n")
        sys.exit(1)

    for json_file in sys.argv[1:]:
        if not os.path.isfile(json_file):
            sys.stderr.write(f"File not found: {json_file}\n")
            continue
        try:
            process_json(json_file)
        except Exception as e:
            sys.stderr.write(f"Error processing {json_file}: {e}\n")

if __name__ == "__main__":
    main()
