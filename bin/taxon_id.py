#!/usr/bin/env python3
"""
Print one sample's taxonomy as og_id,nominal_species_id,taxon_id,class.

Delegates to create_samplesheet.resolve_species_info so this path resolves a
name exactly as the samplesheet does: master_species first, then the NCBI
taxdump when --taxdump-dir is given.
"""
import argparse
import sys

# Sibling import: Nextflow puts bin/ on PATH, so this resolves via sys.path[0].
from create_samplesheet import build_resolver, load_db_config, resolve_species_info


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("config_file", help="Path to the postgres config file.")
    parser.add_argument("og_id", help="Sample OG ID, e.g. OG1323.")
    parser.add_argument("--taxdump-dir", default=None,
                        help="NCBI taxdump directory (nodes.dmp/names.dmp).")
    args = parser.parse_args()

    db_params = load_db_config(args.config_file)
    resolver = build_resolver(args.taxdump_dir)
    nominal_species_id, taxon_id, tax_class, source, _level, _rank = \
        resolve_species_info(db_params, args.og_id, resolver)

    if source == 'unresolved':
        # The caller treats empty output as 'unknown' and warns.
        print(f"[WARN] No matching species data found for OG ID: {args.og_id}",
              file=sys.stderr)
        return 0
    print("og_id,nominal_species_id,taxon_id,class")
    print(f"{args.og_id},{nominal_species_id},{taxon_id},{tax_class}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
