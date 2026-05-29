# OceanOmics-OceanGenomes-Draft-Genomes: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

Some folders are only present if the corresponding step is run (for example, `basespace/` when reads are downloaded, or `busco/` when genome QC is enabled).

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Download and repair reads](#download-and-repair-reads-optional) (BaseSpace + BBMap repair)
- [Read QC and trimming](#read-qc-and-trimming-fastqc--fastp)
- [K-mer profiling and coverage](#k-mer-profiling-and-coverage)
- [Assembly](#assembly-megahit)
- [Decontamination](#decontamination-fcs-gx--tiara--bbmap)
- [Genome QC](#genome-qc-busco-bwa-mem2-merqury-gfastats)
- [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

### Download and repair reads (optional)

<details markdown="1">
<summary>Output files</summary>

- `basespace/<RUN>/`
  - Raw FASTQs downloaded per dataset from Illumina BaseSpace plus run JSON metadata.
- `pooled/<RUN>/`
  - Repaired and paired FASTQs from `BBMAP_REPAIR` for each sample (`*.R1.fq.gz`, `*.R2.fq.gz`).

</details>

Reads are pulled via the BaseSpace CLI when `--skip_download_reads` is false, then re-paired to ensure matching R1/R2 files before downstream QC.

### Read QC and trimming (FastQC / fastp)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/fastp/`
  - Adapter/quality-trimmed FASTQs.
  - `*.json`, `*.html`: fastp reports.
- `draftgenomes/<sample>/fastp/fastqc/`
  - `*_fastqc.html`, `*_fastqc.zip`: FastQC reports for the raw/trimmed reads.

</details>

fastp performs adapter/quality trimming and filtering; FastQC provides per-sample QC summaries. Both feed into MultiQC.

### K-mer profiling and coverage

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/coverage/`
  - Per-sample coverage JSON files from `CALCULATE_SEQUENCING_COVERAGE`.
- `coverage_summary/genome_coverage_summary.csv`
  - Combined coverage summary compiled across samples.
- `genomescope2/`
  - GenomeScope2 model fit outputs (`*_summary.txt`, plots) derived from k-mer histograms.
- `meryl/` (or similar, depending on profile)
  - meryl count/unionsum/histogram databases and plots used for genome size estimation.

</details>

K-mer histograms are generated with meryl and modelled by GenomeScope2 to estimate genome size, heterozygosity, and duplication; coverage metrics are derived by combining fastp and GenomeScope outputs.

### Assembly (MEGAHIT)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/assemblies/genome/`
  - MEGAHIT contigs (`*.contigs.fa`) and reformatted FASTA files.

</details>

MEGAHIT assembles the trimmed reads into draft contigs that are passed to decontamination and QC.

### Decontamination (FCS-GX / Tiara / BBMap)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/assemblies/genome/NCBI/`
  - NCBI FCS-GX contamination reports and cleaned/contaminant FASTAs.
- `draftgenomes/<sample>/assemblies/genome/NCBI/adaptor/`
  - Adaptor screening reports from `fcs-adaptor`.
- `draftgenomes/<sample>/assemblies/genome/tiara/`
  - Tiara classification summaries (`*.tiara_filter_summary.txt`).
- `draftgenomes/<sample>/assemblies/genome/`
  - Final cleaned contigs after Tiara/BBMap filtering and counts of short contigs (`*.contig_count_500bp.txt`).

</details>

Assemblies are screened with NCBI FCS-GX for contamination, adapters are removed, organelle sequences flagged with Tiara, and remaining unwanted contigs filtered with BBMap.

### Genome QC (BUSCO, BWA-MEM2, Merqury, Gfastats)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/assemblies/genome/busco/`
  - BUSCO short summaries and full results; extracted BUSCO sequences under `busco_sequences/`.
- `draftgenomes/<sample>/assemblies/genome/bwa/`
  - BWA-MEM2 alignments of reads back to the assembly (BAM + logs).
- `draftgenomes/<sample>/kmers/`
  - Merqury completeness/QV statistics (`*.completeness.stats`, `*.qv`).
- `draftgenomes/<sample>/assemblies/genome/gfastats/`
  - Gfastats assembly summary metrics.

</details>

Quality metrics cover gene content (BUSCO), read mapping (BWA-MEM2), k-mer based QV/completeness (Merqury), and assembly structure (Gfastats). BUSCO lineage selection is based on the samplesheet `class`: `Actinopteri` uses `--busco_acti_db`, `Anthozoa`/`Cnidaria` uses `--busco_metazoa_db`, and other values fall back to `--busco_vert_db`.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
