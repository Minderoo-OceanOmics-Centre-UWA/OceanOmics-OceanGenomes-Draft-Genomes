#!/usr/bin/env python3
"""Compatibility CLI for decontamination statistics uploads."""

import argparse
import sys

from draft_genome_stats import parse_filter_report, parse_lt500, parse_tiara, upload_records


def main():
    parser = argparse.ArgumentParser(description="Upsert Tiara/filter/<500bp results into PostgreSQL.")
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("--tiara")
    parser.add_argument("--filter")
    parser.add_argument("--contigs")
    args = parser.parse_args()
    if not (args.tiara or args.filter or args.contigs):
        parser.error("provide at least one of --tiara, --filter or --contigs")
    try:
        records = []
        if args.tiara:
            records.append(("decontamination", parse_tiara(args.tiara)))
        if args.filter:
            records.append(("decontamination", parse_filter_report(args.filter)))
        if args.contigs:
            records.append(("decontamination", parse_lt500(args.contigs)))
        identities = {(record["og_id"], record["seq_date"]) for _, record in records}
        if len(identities) != 1:
            raise ValueError(f"inputs have mixed identifiers: {sorted(identities)}")
        upload_records(args.config, records)
        og_id, seq_date = next(iter(identities))
        print(f"✅ Upserted decontamination results for {og_id}/{seq_date}")
    except Exception as exc:
        sys.exit(f"❌ Error: {exc}")


if __name__ == "__main__":
    main()
