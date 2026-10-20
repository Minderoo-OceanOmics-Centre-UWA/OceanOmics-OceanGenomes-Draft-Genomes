# Tests

Unit tests for the Python and shell helpers in [`bin/`](../bin). They use only the standard library
and run with the system Python on Setonix:

```bash
python3 -m unittest discover -s tests -p "test_*.py"
```

(`pytest tests` also works where pytest is installed.)

| File | Covers |
| --- | --- |
| `test_create_samplesheet_taxonomy.py` | `create_samplesheet.py` / `taxdump_lineage.py`: taxonomy from the `species` table with the NCBI taxdump fallback, and the `unresolved` abort |
| `test_genomescope_reliability.py` | `genomescope_reliability.py`: the GenomeScope reliability gate, on real NOVA_260909_LA fixtures |
| `test_kmer_depth_lambda.py` | `kmer_depth_lambda.py`: coverage measured from BUSCO read depth and the host depth partition |
| `test_megahit_checkpoint_key.py` | `megahit_checkpoint_key.sh` and the MEGAHIT checkpoint decisions |
| `test_draft_genome_stats_backfill.py` | `draft_genome_stats.py` / `backfill_draft_genome_stats.py`: parsing and the validate-then-apply backfill |

There are no pipeline-level nf-test tests. The only local module nf-test is
`modules/local/bbmap/filterbyname/tests/main.nf.test`; the vendored nf-core modules carry their own.
