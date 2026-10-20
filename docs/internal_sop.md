# OceanGenomes Draft Genome Pipeline SOP

This SOP describes the standard internal workflow for running OceanGenomes short-read draft genome
assemblies on Pawsey (Setonix) using `nextflow_run_template.sh`. Parameter details are in
[usage](usage.md) and output files in [output](output.md).

## 1. Prepare the repositories

The pipeline and the mitogenome pipeline are kept side by side under `/software`:

```bash
cd /software/projects/pawsey1348/$USER/repos
git clone https://github.com/Minderoo-OceanOmics-Centre-UWA/OceanOmics-OceanGenomes-Draft-Genomes.git
```

For an existing checkout:

```bash
cd /software/projects/pawsey1348/$USER/repos/OceanOmics-OceanGenomes-Draft-Genomes
git pull
```

The run template expects the mitogenome pipeline at:

```bash
/software/projects/pawsey1348/$USER/repos/Oceanomics-OceanGenomes-Mitogenomes
```

If your clone is elsewhere, change `MITO_PIPELINE_DIR` in the copied run script. The directory must
contain `nextflow_run.sh`, or the run fails at start. The draft pipeline uses it to write a ready-made
mitogenome launcher for this run (section 4).

## 2. Create a run-specific launcher

Set the BaseSpace run ID and copy the template to a run-specific script:

```bash
RUN=NOVA_260925_LA
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
```

Review `nextflow_run_${RUN}.sh` before starting. It will:

- load the Nextflow and Singularity modules
- create the run directory `/scratch/pawsey1348/$USER/$RUN`
- copy `scripts/backup_scripts` into it and set `RUN` in the copied `backup_scripts/configfile.txt`
- change into the run directory and launch from there, with work files under
  `/scratch/pawsey1348/$USER/$RUN/work/$RUN`
- after the run, write per-sample compute cost to `pipeline_info/cost_per_sample.csv`

The script must be started from the repository directory (it uses relative paths to copy the backup
scripts and find `main.nf`).

## 3. Run the pipeline

Run the launcher inside `tmux` so the workflow continues if your shell disconnects:

```bash
tmux new -s "$RUN"
bash "nextflow_run_${RUN}.sh"
```

To detach, press `Ctrl-b` then `d`. To reattach: `tmux attach -t "$RUN"`.

Run directory and generated samplesheet:

```bash
/scratch/pawsey1348/$USER/$RUN
/scratch/pawsey1348/$USER/$RUN/samplesheet/${RUN}_samplesheet.csv
```

Optional: to check the samplesheet before spending SUs on assembly, add
`--stop_after_samplesheet true` to the first launch. Review (and if needed trim) the samplesheet,
then remove the flag and launch again; `-resume` reuses the download and the edited sheet.

If the run stops because a sample's taxonomy could not be resolved, see
[usage: when a sample's class is unknown](usage.md#when-a-samples-class-is-unknown).

## 4. Mitogenomes

Once fastp has finished, the pipeline writes
`run_${RUN}_mitogenomes_<timestamp>.sh` into the mitogenome pipeline directory. It is not started
automatically. Launch it from that directory (in its own `tmux` session) when ready:

```bash
cd /software/projects/pawsey1348/$USER/repos/Oceanomics-OceanGenomes-Mitogenomes
bash run_${RUN}_mitogenomes_<timestamp>.sh
```

## 5. Rerun using the generated samplesheet

After the first BaseSpace download run, you can rerun from the generated samplesheet instead of
downloading reads again. In `nextflow_run_${RUN}.sh`:

- set `--skip_bs_download true`
- set `--skip_download_reads true`
- uncomment the `--input "$OUT/samplesheet/${RUN}_samplesheet.csv"` line at the end of the command

Keep `--run "$RUN"` so output naming and metadata remain tied to the same sequencing run.

## 6. Raw data backup

Run the raw backup once the BaseSpace download is complete and the raw data are under
`/scratch/pawsey1348/$USER/$RUN/basespace/$RUN`. Use the copies in the run directory:

```bash
cd /scratch/pawsey1348/$USER/$RUN/backup_scripts/raw_data_backup_aws
sbatch 01_aws_raw_backup.sh
```

This can also be run in a `tmux` session.

`01_aws_raw_backup.sh` copies the raw BaseSpace download to `s3:ocom.oceangenomes/illumina-raw/$RUN`,
then runs `02_data_audit.sh`, which compares local and S3 file counts and sizes. Check the joined
audit file `_draftRAWcheck-join.$RUN.TSV`.

There is no separate FASTQ/FastQC backup step.

## 7. Check expected output files

When the pipeline has completed, check that expected files are present before the final backup:

```bash
cd /scratch/pawsey1348/$USER/$RUN/backup_scripts/check_files
bash 01_checkfiles.sh
```

`_FileCheck.tsv` has one row per expected file (`<OG>` standing for the sample ID) and one column per
sample, each cell `Found <path>` or `Missing`. The script also prints how many rows have a missing
file. Review the missing rows; if important outputs are missing, rerun the affected samples or steps
before the final backup. `kraken2/` rows are expected to be missing for runs launched with
`--skip_kraken2_decontamination`.

## 8. Final genome backup

Run the final backup after the pipeline is complete and output checks have passed:

```bash
cd /scratch/pawsey1348/$USER/$RUN/backup_scripts/finall_all_data_backup_acacia
sbatch 01_final_genome_backup.sh
```

(The directory name really is spelled `finall_...`.) This can also be run in a `tmux` session.

The final backup script:

- copies selected per-sample outputs to `s3:ocom.oceangenomes/analysed-data/draft-genomes/<OG>/`:
  trimmed FASTQs and fastp JSON from `fastp/`, the final `.fna` assembly, and the per-sample MultiQC
  HTML report
- archives each `.meryl` directory to `.meryl.tar.gz` and removes the directory once the archive
  verifies (if any archive fails, it stops before the move)
- runs `02_data_audit_local.sh` to record the local subsets that should match Acacia and S3
- **moves** the `draftgenomes` directory, except every `fastp/` folder, to Acacia at
  `pawsey0964:oceanomics-draftgenomes/genomes.v2/`
- runs `03_data_audit_acacia.sh` to compare the remote Acacia and S3 copies with the local record

Because the Acacia transfer uses `rclone move`, everything except `fastp/` is gone from the local
`draftgenomes` directory once it succeeds. The `fastp/` folders (trimmed reads, already on S3) stay on
scratch and can be deleted once the audits pass and the mitogenome run no longer needs them.
Check the joined audit files:

```bash
_draftcheck-workflow-acacia-join.$RUN.tsv
_draftcheck-workflow-s3-join.$RUN.tsv
```

The Acacia audit compares the local workflow subset excluding `fastp/`. The S3 audit compares only the
files copied to S3. Both files have `ObjectCountMatch` and `ByteCountMatch` columns, which should be
`YES` for every sample before the backup is treated as complete.

## 9. Run records

Upload the run's complete MultiQC report (`multiqc/multiqc_report.html`) to LabArchives and put the
link in the Asana card for the run.
