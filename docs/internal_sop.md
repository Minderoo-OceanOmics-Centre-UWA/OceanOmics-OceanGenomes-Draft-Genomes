# OceanGenomes Draft Genome Pipeline SOP

This SOP describes the standard internal workflow for running OceanGenomes short-read draft genome assemblies on Pawsey using `nextflow_run_template.sh`.

## 1. Prepare the repository

Clone or update the draft genome pipeline repository on scratch:

```bash
cd /scratch/pawsey1348/$USER/agent
git clone https://github.com/Minderoo-OceanOmics-Centre-UWA/OceanOmics-OceanGenomes-Draft-Genomes.git
cd OceanOmics-OceanGenomes-Draft-Genomes
```

For an existing checkout:

```bash
cd /scratch/pawsey1348/$USER/_NFCORE/OceanOmics-OceanGenomes-Draft-Genomes
git pull
```

The run template assumes the mitogenome pipeline is available at:

```bash
/scratch/pawsey0964/$USER/_NFCORE/Oceanomics-OceanGenomes-Mitogenomes
```

If your clone is somewhere else, update `--mitogenome_nfcore_dir` in the copied run script before launching.
If you havent cloned the mitogenome pipeline, clone it so that the run script will automatically be generated from the draft genome pipeline to run on the samples from this run.

## 2. Create a run-specific launcher

Set the BaseSpace run ID and copy the template to a run-specific script:

```bash
RUN=NEXT_250724_ET
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
```

Review `nextflow_run_${RUN}.sh` before starting. The template is set up for standard OceanGenomes runs and will:

- load Nextflow and Singularity modules
- create `/scratch/pawsey0964/$USER/$RUN`
- copy `scripts/backup_scripts` into the run directory
- update the copied `backup_scripts/configfile.txt` with the selected `RUN`
- run Nextflow from the run output directory (this is what the run script does, you need to launch the run script from the cloned repo directory)
- keep work files under `/scratch/pawsey0964/$USER/$RUN/work/$RUN`

## 3. Run the pipeline

Run the launcher inside `tmux` so the workflow continues if your shell disconnects:

```bash
tmux new -s "$RUN"
bash "nextflow_run_${RUN}.sh"
```

To detach from `tmux`, press `Ctrl-b` then `d`. To reattach:

```bash
tmux attach -t "$RUN"
```

The run output directory is:

```bash
/scratch/pawsey0964/$USER/$RUN
```

The automatically generated samplesheet is written under:

```bash
/scratch/pawsey0964/$USER/$RUN/samplesheet/${RUN}_samplesheet.csv
```

## 4. Rerun using the generated samplesheet

After the first BaseSpace download run, you can rerun from the generated samplesheet instead of downloading reads again.

In `nextflow_run_${RUN}.sh`:

- set `--skip_skip_bs_download true`
- set `--skip_download_reads true`
- uncomment the `--input "$OUT/samplesheet/${RUN}_samplesheet.csv"` line at the end of the command

Keep `--run "$RUN"` so output naming and metadata remain tied to the same sequencing run.

## 5. Raw data backup

Raw data backup is run after BaseSpace download has completed and the raw data are present under:

```bash
/scratch/pawsey0964/$USER/$RUN/basespace/$RUN
```

The template copies backup scripts into the run directory, so run backups from the copied scripts:

```bash
cd /scratch/pawsey0964/$USER/$RUN/backup_scripts/raw_data_backup_aws
sbatch 01_aws_raw_backup.sh
```
Alternatively the backup can be ran in a tmux session.


`01_aws_raw_backup.sh` copies the raw BaseSpace download to:

```bash
s3:oceanomics/OceanGenomes/illumina-raw/$RUN
```

It then runs `02_data_audit.sh`, which compares local and S3 file counts/sizes and writes audit TSVs in the backup working directory. Check the joined audit file:

```bash
_draftRAWcheck-join.$RUN.TSV
```

There is no separate fastq/FastQC backup step. The old FastQ backup scripts are no longer used.

## 6. Check expected output files

When the pipeline has completed, check that expected files are present before final backup:

```bash
cd /scratch/pawsey0964/$USER/$RUN/backup_scripts/check_files
bash 01_checkfiles.sh
```

This writes:

```bash
_FileCheck.tsv
```

Review `_FileCheck.tsv` for missing files. If important outputs are missing, rerun the affected samples or steps before final backup.

## 7. Final genome backup

Run the final backup after the pipeline is complete and output checks have passed:

```bash
cd /scratch/pawsey0964/$USER/$RUN/backup_scripts/finall_all_data_backup_acacia
sbatch 01_final_genome_backup.sh
```
This can also be ran in a tmux session as an alternative.


The final backup script:

- copies selected final per-sample outputs to S3 under `s3:oceanomics/OceanGenomes/analysed-data/draft-genomes`
- archives each `.meryl` directory to `.meryl.tar.gz` and removes the unarchived `.meryl` directory after verifying the archive
- runs `02_data_audit_local.sh` to record the local subsets that should match Acacia and S3
- moves the final `draftgenomes` directory to Acacia at `pawsey0964:oceanomics-draftgenomes/genomes.v2/`
- runs `03_data_audit_acacia.sh` to compare the matching remote Acacia and S3 subsets

Important: the Acacia transfer uses `rclone move`, not `rclone copy`. After the final backup completes successfully, the moved files are no longer retained in the local `draftgenomes` directory.

Check the joined final audit files:

```bash
_draftcheck-workflow-acacia-join.$RUN.tsv
_draftcheck-workflow-s3-join.$RUN.tsv
```

The Acacia audit compares the local final workflow subset excluding `fastp/`, because `fastp/` is not moved to Acacia. The S3 audit compares only the selected files copied to S3: trimmed FASTQs, `fastp` JSON files, final `.fna` assemblies, and MultiQC HTML reports.

Both joined audit files include `ObjectCountMatch` and `ByteCountMatch` columns. These should be `YES` for every sample before the backup is treated as complete.

## 8. Run records

For each run, make sure to upload the complete multiqc report for the run onto labarchives. The link to the location should be included in the asana card for this task.
