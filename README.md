<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-oceangenomesdraftgenomes_logo_dark.png">
    <img alt="OceanOmics-OceanGenomes-Draft-Genomes" src="docs/images/nf-core-oceangenomesdraftgenomes_logo_light.png">
  </picture>
</h1>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**OceanOmics-OceanGenomes-Draft-Genomes** builds short-read draft genome assemblies for Ocean Omics /
Ocean Genomes samples, fish and invertebrates alike. It downloads an Illumina run from BaseSpace (or
takes existing FASTQs), assembles each sample, removes contamination, checks the assembly and the
genome size estimate, and loads the QC metrics into the OceanOmics database.

1. Download and re-pair reads from Illumina BaseSpace, and build the samplesheet from the database,
   with an NCBI taxdump fallback for taxa not yet curated (optional if `--input` is used)
2. Read QC and trimming with fastp and FastQC
3. Read-level decontamination with kraken2 + KrakenTools (bacteria, archaea and viruses removed
   before assembly)
4. K-mer profiling with meryl and GenomeScope2
5. Assembly with MEGAHIT, checkpointed and fingerprinted so interrupted assemblies resume safely
6. Assembly decontamination with NCBI FCS-GX, FCS-adaptor and Tiara
7. Genome QC with BUSCO (lineage chosen by taxonomic class), BWA-MEM2, Merqury, gfastats and SeqKit
8. Genome size and coverage verdict: GenomeScope fit checks, an independent coverage measurement
   from read depth over BUSCO genes, a host/symbiont depth partition of the assembly, an optional
   reseeded GenomeScope refit, and a cross-check against the assembly
9. Run-level and per-sample MultiQC reports
10. Upload of QC metrics to the OceanOmics PostgreSQL database
11. A ready-to-run launcher for the [mitogenome pipeline](https://github.com/Minderoo-OceanOmics-Centre-UWA/Oceanomics-OceanGenomes-Mitogenomes)
    on the same reads, and per-sample compute (SU) cost

## Usage

For standard OceanGenomes runs on Setonix, copy [`nextflow_run_template.sh`](nextflow_run_template.sh)
to a run-specific launcher, set `RUN`, and launch it from the repository root:

```bash
RUN=NOVA_260925_LA
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
bash "nextflow_run_${RUN}.sh"
```

The launcher writes everything for the run under `/scratch/pawsey1348/$USER/$RUN`, launches Nextflow
from that directory (so several runs can share one checkout), copies the backup scripts into it, and
records per-sample compute cost when the run ends. To start from existing FASTQs instead of
BaseSpace, use `--input <samplesheet.csv> --skip_download_reads true`.

> [!WARNING]
> Provide pipeline parameters on the command line or with `-params-file`. Config files given with
> `-c` are for resources and executor settings only.

## Documentation

- [Usage](docs/usage.md): running the pipeline, the samplesheet, BUSCO lineage selection, reference
  data, QC thresholds, skipping stages, resuming MEGAHIT checkpoints, backfilling the database.
- [Output](docs/output.md): every published file and how to read the coverage verdict.
- [Internal SOP](docs/internal_sop.md): the step-by-step OceanGenomes run and backup procedure.
- [`nextflow_schema.json`](nextflow_schema.json): every parameter (`nextflow run main.nf --help`).
- [`CHANGELOG.md`](CHANGELOG.md): what changed and why.

Supporting tools in this repository:

| Path | What it is |
| --- | --- |
| [`scripts/backup_scripts/`](scripts/backup_scripts/README.md) | Raw and final backups to S3 and Acacia, and the output file check |
| [`compute-audit/`](compute-audit/README.md) | Per-run and per-sample Slurm SU accounting |
| [`scripts/taxonomy/`](scripts/taxonomy/README.md) | Bulk-load NCBI taxa into the `species` table |
| [`bin/backfill_draft_genome_stats.py`](docs/usage.md#backfilling-database-statistics-from-acacia) | Backfill database statistics from archived runs |
| [`scripts/rerun_NOVA_260909_LA/`](scripts/rerun_NOVA_260909_LA/README.md) | Remediation re-run of NOVA_260909_LA and the investigations behind the genome size checks |
| [`tests/`](tests/README.md) | Unit tests for the Python helpers in `bin/` |

## Testing

The Python helpers in `bin/` have unit tests that run with the system Python, no extra packages
needed:

```bash
python3 -m unittest discover -s tests -p "test_*.py"
```

See [`tests/README.md`](tests/README.md) for what each file covers.

## Credits

OceanOmics-OceanGenomes-Draft-Genomes was originally written by Tyler Peirce.

We thank the following people for their extensive assistance in the development of this pipeline:

- Ocean Omics / Ocean Genomes project team and Pawsey collaborators

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
