#!/usr/bin/env python3
"""Compatibility CLI for uploading one gfastats assembly summary."""

import argparse
import sys

from draft_genome_stats import parse_gfastats, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert gfastats results into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-f", "--file", required=True)
    args = parser.parse_args()
    try:
        record = parse_gfastats(args.file)
        upload_records(args.config, [("gfastats", record)])
        print(f"✅ Upserted gfastats summary for {record['og_id']} / {record['seq_date']} from {args.file}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
