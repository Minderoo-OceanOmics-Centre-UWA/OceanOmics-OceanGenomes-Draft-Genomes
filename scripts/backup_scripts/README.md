# Backup scripts

Backup and output-check scripts for a draft genome run. `nextflow_run_template.sh` copies this
directory into the run's outdir and sets `RUN` in the copy of `configfile.txt`, so always run the
copies in `/scratch/pawsey1348/$USER/$RUN/backup_scripts/`, not these originals. The full procedure
is in the [internal SOP](../../docs/internal_sop.md).

| Step | Directory | Run | What it does |
| --- | --- | --- | --- |
| 1. Raw backup (after download) | `raw_data_backup_aws/` | `sbatch 01_aws_raw_backup.sh` | Copies `basespace/$RUN` to `s3:ocom.oceangenomes/illumina-raw/$RUN`, then `02_data_audit.sh` writes `_draftRAWcheck-join.$RUN.TSV` comparing local and S3 counts and sizes. |
| 2. Output check (after the run) | `check_files/` | `bash 01_checkfiles.sh` | Checks each sample for every published output listed in the script (reads, kraken2, coverage verdict, k-mers, assemblies, decontamination, QC, per-sample MultiQC, all meryl shards) and writes `_FileCheck.tsv`: one row per expected file, one column per sample. |
| 3. Final backup | `finall_all_data_backup_acacia/` | `sbatch 01_final_genome_backup.sh` | Copies trimmed reads, fastp JSON, the `.fna` and the per-sample MultiQC report to `s3:ocom.oceangenomes/analysed-data/draft-genomes`, archives `.meryl` directories, then **moves** `draftgenomes/` (except `fastp/`) to `pawsey0964:oceanomics-draftgenomes/genomes.v2/`. `02_data_audit_local.sh` and `03_data_audit_acacia.sh` write `_draftcheck-workflow-{acacia,s3}-join.$RUN.tsv`. |

`configfile.txt` derives the Pawsey project and the run directory (`/scratch/<project>/$USER/$RUN`)
from where the copy sits, so it works for any project as long as the outdir follows that layout.

Both `sbatch` scripts charge `pawsey0964`, the same account the pipeline's own jobs use
([`pawsey_profile.config`](../../pawsey_profile.config)).

To check for a new output, add its pattern to the `expected` list in `01_checkfiles.sh`; the row
labels are generated from the same list.
