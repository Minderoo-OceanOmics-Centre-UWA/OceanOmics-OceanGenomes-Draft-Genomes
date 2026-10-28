# OceanOmics-OceanGenomes-Draft-Genomes: Documentation

- [Usage](usage.md): how to run the pipeline, the samplesheet, BUSCO lineage selection, reference
  data, QC thresholds, skipping stages, resuming MEGAHIT checkpoints and backfilling the database.
- [Output](output.md): every published file and how to interpret the genome size and coverage verdict.
- [Internal SOP](internal_sop.md): the internal OceanGenomes run procedure using
  `nextflow_run_template.sh`, including the mitogenome hand-off and backups.

Elsewhere in the repository:

- [`scripts/backup_scripts/README.md`](../scripts/backup_scripts/README.md): raw and final backups and the output file check.
- [`compute-audit/README.md`](../compute-audit/README.md): Slurm SU accounting per run and per sample.
- [`scripts/taxonomy/README.md`](../scripts/taxonomy/README.md): loading NCBI taxa into the `species` table.
- [`tests/README.md`](../tests/README.md): running the unit tests.
- [`scripts/rerun_NOVA_260909_LA/README.md`](../scripts/rerun_NOVA_260909_LA/README.md): the NOVA_260909_LA remediation re-run.
