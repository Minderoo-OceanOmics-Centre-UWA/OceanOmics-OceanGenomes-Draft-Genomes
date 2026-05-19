<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-oceangenomesdraftgenomes_logo_dark.png">
    <img alt="OceanOmics-OceanGenomes-Draft-Genomes" src="docs/images/nf-core-oceangenomesdraftgenomes_logo_light.png">
  </picture>
</h1>

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**OceanOmics-OceanGenomes-Draft-Genomes** builds short-read draft assemblies for Ocean Omics / Ocean Genomes samples. It can download Illumina runs directly from BaseSpace or consume existing FASTQs, performs adapter/quality trimming, estimates genome size and coverage, assembles with MEGAHIT, screens contaminants (NCBI FCS-GX, FCS adaptor, Tiara + BBMap), runs genome QC (BUSCO, BWA-MEM2, Merqury, Gfastats), and collates results with MultiQC. Optional steps push QC summaries into the project SQL database.

Key steps:

- Download and re-pair reads from Illumina BaseSpace (optional if `--skip_download_reads` or `--input` used)
- Read QC and trimming with FastQC / fastp
- K-mer profiling (meryl, GenomeScope2) and coverage estimation
- Assembly with MEGAHIT
- Contamination filtering with NCBI FCS-GX + adaptor screening, Tiara classification, and BBMap cleanup
- Genome QC with BUSCO, BWA-MEM2 alignments, Merqury, and Gfastats
- MultiQC summary and optional SQL upload of QC metrics

## Usage

> [!NOTE]
> If you are new to Nextflow, install it and make sure to test your setup with `-profile test` before running the workflow on actual data.

For standard OceanGenomes runs on Pawsey, use the provided [`nextflow_run_template.sh`](nextflow_run_template.sh). The template is preconfigured with the project database, BaseSpace, BUSCO, contamination-screening, mitogenome, per-run working-directory, and backup-script settings. In most cases you only need to copy the template to a run-specific launcher, update the `RUN` variable, and launch it from the repository root:

```bash
RUN=NEXT_250724_ET
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
bash "nextflow_run_${RUN}.sh"
```

The run-specific launcher writes each run under `/scratch/pawsey0964/$USER/$RUN`, runs Nextflow from that output directory to avoid clashes between runs, and copies the backup helper scripts into the run directory. The plain [`nextflow_run.sh`](nextflow_run.sh) script is kept as an alternative for running from the repository directory, with the work directory under `./work/$RUN`. Use `--input <samplesheet.csv> --skip_download_reads true` if you already have FASTQs instead of pulling from BaseSpace.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_.

For more details and further functionality, please refer to the local [usage documentation](docs/usage.md), [internal SOP](docs/internal_sop.md), and [`nextflow_schema.json`](nextflow_schema.json).

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

OceanOmics-OceanGenomes-Draft-Genomes was originally written by Tyler Peirce.

We thank the following people for their extensive assistance in the development of this pipeline:

- Ocean Omics / Ocean Genomes project team and Pawsey collaborators

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
