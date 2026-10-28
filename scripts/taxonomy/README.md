# Taxonomy loader

`load_taxonomy.py` bulk-loads NCBI species-rank taxa into the OceanOmics `species` table from the
NCBI `new_taxdump`, so samples in those clades resolve `taxon_id` and `class` from the curated table.
Target classes (`--class`) and/or phyla (`--phylum`), comma-separated or repeated. It is a dry run by
default (writes a CSV preview); `--apply <db.cfg>` loads the rows. Existing rows are never changed.

```bash
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata --apply ~/postgresql_details/oceanomics.cfg
```

Prefer `--phylum` for invertebrates: NCBI leaves the class empty for many invertebrate lineages, and
those rows get the phylum written into `species.class`.

With `--taxonkit_db_dir` set, the pipeline already falls back to the NCBI taxdump for taxa the table
lacks, so this is for making a name permanent rather than for unblocking a run. See
[usage: where taxonomy comes from](../../docs/usage.md#where-taxonomy-comes-from).
