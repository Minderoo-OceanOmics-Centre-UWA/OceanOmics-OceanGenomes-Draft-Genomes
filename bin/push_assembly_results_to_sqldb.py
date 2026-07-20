#!/usr/bin/env python3
"""Compatibility CLI for uploading one GenomeScope report."""

import argparse
import sys

from draft_genome_stats import parse_genomescope, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert GenomeScope results into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-f", "--file", required=True)
    args = parser.parse_args()
    try:
        record = parse_genomescope(args.file)
        upload_records(args.config, [("assembly", record)])
        print(f"✅ Upserted og_id={record['og_id']} seq_date={record['seq_date']} from {args.file}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
