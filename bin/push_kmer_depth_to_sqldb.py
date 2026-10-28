#!/usr/bin/env python3
"""Upload read depth and the host partition from *_kmer_depth.json, and nothing else.

Two uses, both after the fact:
  - the depth-only re-run of drafts assembled before MEASURE_KMER_COVERAGE existed, which
    produces a *_kmer_depth.json and no new coverage summary;
  - backfilling busco_read_depth, busco_regions and partition_note for runs whose depth
    was measured but never uploaded.

It writes only DEPTH_COLUMNS, so the GenomeScope numbers, the reliability verdict and every
other family stay as they are. A draft with no draft_genomes row is refused rather than
created: depth without the assembly it describes is not a draft genome.
"""

import argparse
import sys

from draft_genome_stats import DEPTH_COLUMNS, connect, parse_kmer_depth, upsert_record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-c", "--config", required=True)
    parser.add_argument("-k", "--kmer-depth", required=True, nargs="+",
                        help="one or more published *_kmer_depth.json files")
    parser.add_argument("--dry-run", action="store_true",
                        help="parse and check every file, print what would be written, change nothing")
    args = parser.parse_args()

    try:
        records = [parse_kmer_depth(path) for path in args.kmer_depth]
    except Exception as exc:
        sys.exit(f"❌ Could not parse: {exc}")

    connection = connect(args.config)
    try:
        with connection.cursor() as cursor:
            missing = []
            for record in records:
                cursor.execute("SELECT 1 FROM draft_genomes WHERE og_id = %s AND seq_date = %s",
                               (record["og_id"], record["seq_date"]))
                if cursor.fetchone() is None:
                    missing.append(f"{record['og_id']}/{record['seq_date']}")
            if missing:
                raise ValueError(f"no draft_genomes row for: {', '.join(missing)}")

            print("og_id\tseq_date\t" + "\t".join(DEPTH_COLUMNS))
            for record in records:
                print("\t".join([record["og_id"], record["seq_date"]]
                                + ["" if record.get(c) is None else str(record[c]) for c in DEPTH_COLUMNS]))
                if not args.dry_run:
                    upsert_record(cursor, "assembly", record)
        if args.dry_run:
            connection.rollback()
            print(f"--dry-run: {len(records)} checked, nothing written", file=sys.stderr)
        else:
            connection.commit()
            print(f"✅ Upserted depth for {len(records)} draft(s)", file=sys.stderr)
    except Exception as exc:
        connection.rollback()
        sys.exit(f"❌ Error: {exc}")
    finally:
        connection.close()


if __name__ == "__main__":
    main()
