#!/usr/bin/env python3
"""Compatibility CLI for Merqury completeness and QV uploads."""

import argparse
import sys

from draft_genome_stats import _read_merqury, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert Merqury QV and/or completeness into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--qv")
    parser.add_argument("--comp")
    args = parser.parse_args()
    if not (args.qv or args.comp):
        parser.error("provide at least one of --qv or --comp")
    try:
        records = []
        if args.qv:
            records.append(("merqury", _read_merqury(args.qv, "qv")))
        if args.comp:
            records.append(("merqury", _read_merqury(args.comp, "completeness")))
        identities = {(record["og_id"], record["seq_date"]) for _, record in records}
        if len(identities) != 1:
            raise ValueError(f"inputs have mixed identifiers: {sorted(identities)}")
        upload_records(args.config, records)
        og_id, seq_date = next(iter(identities))
        print(f"✅ Upserted Merqury results for {og_id}/{seq_date}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
