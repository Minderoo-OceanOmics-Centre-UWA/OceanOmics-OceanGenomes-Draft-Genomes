# compute-audit

Slurm service-unit (SU) accounting for a finished Nextflow run, per job and per sample. Vendored
here so the pipeline can report its own cost without depending on any other repository.

## Scripts

- `nf_workflow_cost.sh`: costs a whole Nextflow run. It reads the job list from Nextflow's execution
  database (`nextflow log <session-id>`) rather than the rotating `.nextflow.log` files, so tasks that
  ran once early and were cached on every later resume are still counted. It calls `job_usage.sh` for
  each Slurm job.
- `job_usage.sh`: SU usage for one Slurm job ID, as text, CSV or JSON. Vendored from
  [SarahBeecroft/slurm_post_script](https://github.com/SarahBeecroft/slurm_post_script) at commit
  `4559fe9`; re-sync by hand if the upstream SU calculation or CSV schema changes.

## Usage

```bash
bash compute-audit/nf_workflow_cost.sh [--tag-filter <file>] <nextflow-run-dir> [output_csv]
```

- `<nextflow-run-dir>`: the directory Nextflow was launched from. It must still contain `.nextflow/`
  with its history and `cache/<session-id>/`.
- `output_csv`: per-job usage, default `metrics_summary.csv`.
- `--tag-filter <file>`: only count jobs whose task tag (usually the sample ID) matches a line in
  the file, to cost a subset of samples from a shared run.

Outputs, written next to `output_csv`:

| File | Contents |
| --- | --- |
| `<output_csv>` | one row per Slurm job, with a `tag` column |
| `<output_csv basename>.report.txt` | readable summary |
| `cost_per_sample.csv` | SU per sample: direct SU plus a share of run-level overhead |

Nextflow must be loadable (`module load nextflow/25.04.6`) and `sacct` available.

## How the pipeline uses it

The tail of [`nextflow_run_template.sh`](../nextflow_run_template.sh) runs

```bash
bash compute-audit/nf_workflow_cost.sh "$OUT" "$OUT/pipeline_info/compute_usage.csv"
```

after every launch, producing `pipeline_info/compute_usage.csv`, `compute_usage.report.txt` and
`cost_per_sample.csv`. A failure here is reported but never changes the launcher's exit status. Set
`OCEANOMICS_SKIP_COST=1` to skip it.
