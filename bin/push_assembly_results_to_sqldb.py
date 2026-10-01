#!/usr/bin/env python3
"""Compatibility CLI for uploading one GenomeScope report."""

import argparse
import sys

from draft_genome_stats import parse_coverage_summary, parse_genomescope, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert GenomeScope results into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-f", "--file", required=True)
    parser.add_argument("-v", "--coverage-summary",
                        help="published *_coverage_summary.json, for the reliability verdict")
    args = parser.parse_args()
    try:
        records = [("assembly", parse_genomescope(args.file))]
        if args.coverage_summary:
            # A second record for the same (og_id, seq_date); upload_records upserts, so
            # the two merge into one row.
            records.append(("assembly", parse_coverage_summary(args.coverage_summary)))
        upload_records(args.config, records)
        record = records[0][1]
        print(f"✅ Upserted og_id={record['og_id']} seq_date={record['seq_date']} from {args.file}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
