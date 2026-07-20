#!/usr/bin/env python3
"""Compatibility CLI for uploading one fastp report."""

import argparse
import sys

from draft_genome_stats import parse_fastp, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert fastp JSON into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-f", "--file", required=True)
    args = parser.parse_args()
    try:
        record = parse_fastp(args.file)
        upload_records(args.config, [("fastp", record)])
        print(f"✅ Upserted og_id={record['og_id']} seq_date={record['seq_date']} from {args.file}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
