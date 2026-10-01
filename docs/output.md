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

The per-sample `*_kmer_depth.json` also carries a genome size that does not come from the k-mer histogram at all. The k-mer route collapses below about 10x host coverage, so `SUMMARISE_KMER_COVERAGE` partitions the assembly by per-contig read depth around the depth of the host's own single-copy BUSCO genes: `host_assembly_size` is the mass sitting at host depth, `symbiont_assembly_size` the mass well above it, and `host_assembly_fraction` how much of the assembly is the animal. `partition_status` says whether the partition was possible (`ok`, `no_host_depth`, `no_contig_depth`). Both sizes are **lower bounds**: a symbiont that happens to sit at the host's depth counts as host, and a high-copy host repeat counts as symbiont. Where the partition succeeded, the genome size cross-check divides by the host mass rather than the whole assembly, and records which it used in `assembly_ratio_basis`.

### Assembly (MEGAHIT)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<sample>/assemblies/genome/`
  - MEGAHIT contigs (`*.contigs.fa`) and reformatted FASTA files.
  - `*.megahit_checkpoint.txt` - which checkpoint path the assembly took, and the input fingerprint it used.
- `megahit_checkpoints/<prefix>_megahit_out/` and `<prefix>_megahit_out.key`
  - The resumable checkpoint itself, and the fingerprint of the inputs it was built from.

</details>

MEGAHIT assembles the trimmed reads into draft contigs that are passed to decontamination and QC.

Assembly is checkpointed outside the work directory so that a task killed part way
through (an OOM, a walltime) resumes rather than starting over. That puts the checkpoint
beyond Nextflow's caching, so it carries a fingerprint of its inputs: the reads (by name
and size), the assembly arguments and the megahit version. A checkpoint whose fingerprint
no longer matches is rebuilt rather than reused, which is what stops a change upstream of
assembly, enabling kraken2 for instance, from silently republishing the previous
assembly. Memory and thread counts are excluded from the fingerprint, so a retry at
higher memory still resumes.

`--megahit_checkpoint_unkeyed` decides what happens to checkpoints written before
fingerprinting existed (`invalidate`, the default, or `adopt`), `--megahit_stale_checkpoint`
what happens when a fingerprint no longer matches (`rerun`, `archive` or `fail`), and
`--megahit_checkpoint_cleanup` whether the intermediate assembly graph is pruned once a
checkpoint completes. `*.megahit_checkpoint.txt` records which path each task took, so
whether an assembly was actually rebuilt is a file to read rather than an inference.

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

Quality metrics cover gene content (BUSCO), read mapping (BWA-MEM2), k-mer based QV/completeness (Merqury), and assembly structure (Gfastats). BUSCO lineage selection is based on the samplesheet `class`: ray-finned fish classes use `--busco_acti_db`, other vertebrate classes use `--busco_vert_db`, and everything else (all invertebrates) falls back to `--busco_metazoa_db`.

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
