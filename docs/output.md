# OceanOmics-OceanGenomes-Draft-Genomes: Output

## Introduction

This document describes what the pipeline publishes and how to read it. All paths are relative to
`--outdir`. For a standard launch from [`nextflow_run_template.sh`](../nextflow_run_template.sh) that
is `/scratch/pawsey1348/$USER/<RUN>`.

Publish locations are set in [`conf/modules.config`](../conf/modules.config). A folder only appears
if its step ran: `basespace/` and `pooled/` need the BaseSpace download, `kraken2/` needs
`--kraken2_db`, `push/` needs the SQL upload, and so on.

Throughout, `<RUN>` is the sequencing run (e.g. `NOVA_260925_LA`), `<OG>` the sample ID (e.g.
`OG2622`) and `<prefix>` the samplesheet prefix (e.g. `OG2622.ilmn.260925`). The assembly prefix
is `<prefix>.v129mh` (MEGAHIT v1.2.9).

Most directories also contain two kinds of housekeeping file, which are not described below:

- `versions.yml`: tool versions for the step that published into that directory.
- `NN_<tool>.tool_params_mqcrow.html`: the arguments that step ran with, rendered as a row of the
  MultiQC "tool parameters" table. The number orders the rows.

## Layout at a glance

```text
<outdir>/
├── samplesheet/                 run samplesheet and taxonomy provenance
├── download/taxonkit_dbs/       NCBI taxdump cache (only when --taxonkit_db_dir points here)
├── basespace/<RUN>/             raw BaseSpace download
├── pooled/<RUN>/                lane-merged, re-paired reads
├── draftgenomes/<OG>/           everything per sample (see below)
├── megahit_checkpoints/         resumable MEGAHIT state, outside the work dir
├── coverage_summary/            run-level coverage table
├── push/                        SQL upload receipts
├── multiqc/                     run-level MultiQC report
└── pipeline_info/               Nextflow reports, versions, compute cost
```

```text
draftgenomes/<OG>/
├── fastp/                       trimmed reads + fastp reports
│   └── fastqc/
├── kraken2/                     read-level decontamination reports
├── kmers/                       meryl, Merqury
│   └── <prefix>.genomescope/    GenomeScope2 fit and reseed decision
├── coverage/                    genome size / coverage verdict
├── assemblies/genome/           raw and final assemblies
│   ├── NCBI/                    FCS-GX
│   │   └── adaptor/             FCS-adaptor
│   ├── tiara/
│   ├── busco/
│   │   └── busco_sequences/
│   ├── bwa/
│   ├── gfastats/
│   └── seqkit/
└── multiqc/                     per-sample MultiQC report
```

The deliverable is **`draftgenomes/<OG>/assemblies/genome/<prefix>.v129mh.fna`**, the
decontaminated assembly. Everything else describes it or the reads it was built from.

## Pipeline overview

- [Samplesheet](#samplesheet)
- [Download and repair reads](#download-and-repair-reads-optional)
- [Read QC and trimming](#read-qc-and-trimming-fastp--fastqc)
- [Read-level decontamination](#read-level-decontamination-kraken2)
- [K-mer profiling](#k-mer-profiling-meryl--genomescope2)
- [Assembly](#assembly-megahit)
- [Decontamination](#assembly-decontamination-fcs-gx--fcs-adaptor--tiara)
- [Genome QC](#genome-qc-busco-bwa-mem2-merqury-gfastats-seqkit)
- [Genome size and coverage verdict](#genome-size-and-coverage-verdict)
- [MultiQC](#multiqc)
- [SQL upload receipts](#sql-upload-receipts)
- [Mitogenome launcher](#mitogenome-launcher)
- [Pipeline information and compute cost](#pipeline-information-and-compute-cost)

### Samplesheet

<details markdown="1">
<summary>Output files</summary>

- `samplesheet/`
  - `<RUN>_samplesheet.csv`: the samplesheet the run used, built from the downloaded reads and the
    `sample`/`species` tables.
  - `taxonomy_resolution.tsv`: where each sample's `taxon_id` and `class` came from (`db`,
    `db+taxdump`, `taxdump`, `unresolved`).
- `download/taxonkit_dbs/`: the NCBI taxdump, cached with `storeDir`. It lands here only when
  `--taxonkit_db_dir` is the outdir; the template points it at `/scratch/pawsey1348/$USER`, so it
  is normally shared across runs and not under the outdir at all.

</details>

The published samplesheet is reused automatically on the next launch, so hand edits to it stick.
See [usage: samplesheet input](usage.md#samplesheet-input).

### Download and repair reads (optional)

<details markdown="1">
<summary>Output files</summary>

- `basespace/<RUN>/<dataset>/`: raw FASTQs and the dataset JSON metadata, as downloaded.
- `pooled/<RUN>/`: one pair per sample, `<prefix>.R1.fq.gz` / `<prefix>.R2.fq.gz`, lane-merged and
  re-paired with `BBMAP_REPAIR`.

</details>

Present only when reads are downloaded (`--skip_download_reads false`). The raw backup
(`backup_scripts/raw_data_backup_aws`) copies `basespace/<RUN>/` to S3.

### Read QC and trimming (fastp / FastQC)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/fastp/`
  - `<prefix>.R1.fastq.gz`, `<prefix>.R2.fastq.gz`: trimmed reads. These are the reads every
    later step starts from, and what the mitogenome pipeline is pointed at.
  - `<OG>.ilmn.<RUN>.fastp.json`, `<prefix>.fastp.html`, `<prefix>.fastp.log`: fastp reports.
- `draftgenomes/<OG>/fastp/fastqc/`
  - `*_fastqc.html`, `*_fastqc.zip`: FastQC reports.

</details>

### Read-level decontamination (kraken2)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/kraken2/`
  - `<prefix>.kraken2.report.txt`: the standard kraken2 taxonomic report.
  - `<prefix>.kraken2_retention.txt`: how many read pairs were removed per excluded taxon and the
    fraction kept.

</details>

kraken2 runs on the trimmed reads before meryl and MEGAHIT, and KrakenTools removes pairs assigned
to `--kraken2_exclude_taxids` (Bacteria, Archaea, Viruses and their descendants by default).
Unclassified reads are kept. The filtered reads and the per-read assignment file are intermediates
and are deliberately **not** published; the published `fastp/` reads are the pre-kraken2 reads.
`--kraken2_db` is required unless `--skip_kraken2_decontamination` is set, in which case MEGAHIT and
meryl use the fastp reads directly and this folder is absent.

### K-mer profiling (meryl / GenomeScope2)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/kmers/`
  - `<prefix>.meryl/`: the merged meryl k-mer database (`MERYL_UNIONSUM`). The final backup
    archives it to `<prefix>.meryl.tar.gz`.
  - `<prefix>.meryl.hist`: the k-mer histogram GenomeScope2 fits.
- `draftgenomes/<OG>/kmers/<prefix>.genomescope/`
  - `<prefix>_genomescope2_summary.txt`, `<prefix>_genomescope2_model.txt`: the GenomeScope2 fit.
  - `<prefix>_genomescope2_{linear,log,transformed_linear,transformed_log}_plot.png`.
  - `<prefix>_reseed_decision.json`: whether a reseeded refit was attempted and which fit won
    (`seed`, `lambda_depth`, `winner`, `reason`, and the original and reseeded fit, kmercov and
    genome size).
  - `<prefix>_genomescope2_reseeded_summary.txt`, `<prefix>_genomescope2_reseeded_model.txt`:
    only when a reseed actually ran.

</details>

When `--genomescope2_l` is set the GenomeScope2 files carry an `.l<value>` suffix. When the first
fit is flagged unreliable, `GENOMESCOPE_RESEED` refits with `-l` set to the coverage measured from
read depth, and keeps the refit only if it improves `model_fit_full` by at least
`--genomescope_min_reseed_fit_gain` **and** moves kmercov toward the measurement. Downstream steps
use whichever fit won; `_reseed_decision.json` says which.

### Assembly (MEGAHIT)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/assemblies/genome/`
  - `<prefix>.v129mh.fasta`: raw MEGAHIT contigs, before any decontamination.
  - `<prefix>.megahit_checkpoint.txt`: which checkpoint path the task took (e.g. `fresh`, `skip`,
    `continue`, `invalidated-stale-fresh`, `adopted-unkeyed-skip`, `resume-failed-rebuilt`) and the
    input fingerprint it used.
- `megahit_checkpoints/`
  - `<prefix>_megahit_out/`: the resumable MEGAHIT checkpoint.
  - `<prefix>_megahit_out.key`: the fingerprint of the inputs it was built from.

</details>

Assembly is checkpointed outside the work directory so that a task killed part way through (an
OOM, a walltime) resumes rather than starting over. That puts the checkpoint beyond Nextflow's
caching, so it carries a fingerprint of its inputs: the reads (by name and size), the assembly
arguments and the MEGAHIT version. A checkpoint whose fingerprint no longer matches is rebuilt rather
than reused, which is what stops a change upstream of assembly (enabling kraken2, for instance) from
silently republishing the previous assembly. Memory and thread counts are excluded from the
fingerprint, so a retry at higher memory still resumes.

`--megahit_checkpoint_unkeyed` decides what happens to checkpoints written before fingerprinting
existed (`invalidate`, the default, or `adopt`), `--megahit_stale_checkpoint` what happens when a
fingerprint no longer matches (`rerun`, `archive` or `fail`), and `--megahit_checkpoint_cleanup`
whether the intermediate assembly graph is pruned once a checkpoint completes.

### Assembly decontamination (FCS-GX / FCS-adaptor / Tiara)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/assemblies/genome/NCBI/`
  - `<prefix>.v129mh.fcs_gx_report.txt`, `<prefix>.v129mh.taxonomy.rpt`,
    `<prefix>.v129mh.summary.txt`: FCS-GX reports.
  - `<prefix>.v129mh.rc.fasta`, `<prefix>.v129mh.contam.fasta`: `gx clean-genome` output (EXCLUDE
    and TRIM applied) and the sequence it removed.
  - `<prefix>.review_scaffolds.txt`: the contig names removed in the next step, chosen by
    `--fcs_review_action`.
  - `<prefix>.filter_report.txt`: contigs and bases flagged EXCLUDE, TRIM and REVIEW.
  - `<prefix>.contig_count_500bp.txt`: number of contigs under 500 bp dropped.
  - `<prefix>.v129mh.fa`: assembly after REVIEW removal and the 500 bp floor, passed to adaptor
    screening.
- `draftgenomes/<OG>/assemblies/genome/NCBI/adaptor/`: FCS-adaptor reports.
- `draftgenomes/<OG>/assemblies/genome/`
  - `<prefix>.v129mh.rmadapt.fasta`, `<prefix>.v129mh.adaptor-contam.fasta`: assembly after
    adaptor removal, and what was removed.
  - **`<prefix>.v129mh.fna`: the final decontaminated assembly.**
- `draftgenomes/<OG>/assemblies/genome/tiara/`
  - `<prefix>.v129mh.tiara.txt`, `log_<prefix>.v129mh.tiara.txt`: Tiara classifications.
  - `<prefix>.v129mh.tiara.contig_removal.txt`: contigs removed on Tiara's call.
  - `<prefix>.v129mh.tiara_filter_summary.txt`: counts and bases removed per class.

</details>

Order of operations: FCS-GX screen and clean (`.rc.fasta`), REVIEW removal and 500 bp floor
(`NCBI/*.fa`), FCS-adaptor screen and clean (`.rmadapt.fasta`), Tiara classification of contigs of at
least `--tiara_min_len`, then removal of contigs Tiara calls archaea, bacteria, prokarya, organelle,
mitochondrion or plastid to produce the `.fna`.

### Genome QC (BUSCO, BWA-MEM2, Merqury, gfastats, SeqKit)

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/assemblies/genome/busco/`
  - `<prefix>.v129mh.busco.<lineage>.short_summary.{txt,json}`, `.full_table.tsv`,
    `.missing_busco_list.tsv`, `.batch_summary.txt`, `.logs`.
  - `.busco_sequences.tar.gz` (+ `.md5`): BUSCO's per-gene sequences, archived.
  - `<lineage>` is the first four letters of the lineage used (`acti`, `vert`, `meta`, `crus`,
    `moll`, `anth`, `arth`).
- `draftgenomes/<OG>/assemblies/genome/busco/busco_sequences/`
  - `*.busco_sequences.fasta`, `*.busco_coordinates.bed`: Complete single-copy BUSCOs extracted
    from the assembly, and their coordinates.
- `draftgenomes/<OG>/assemblies/genome/bwa/`
  - `<OG>.sorted.bam`: trimmed reads mapped back to the `.fna`.
  - `<OG>-sn_results.tsv`: `samtools stats` summary numbers.
- `draftgenomes/<OG>/kmers/`
  - `<prefix>.v129mh.merqury.qv`, `<prefix>.v129mh.merqury.completeness.stats`: Merqury QV and
    k-mer completeness, plus `spectra-cn` and `spectra-asm` plots.
- `draftgenomes/<OG>/assemblies/genome/gfastats/<prefix>.v129mh.assembly_summary`: contiguity stats
  (N50, contig count, sizes).
- `draftgenomes/<OG>/assemblies/genome/seqkit/<prefix>.v129mh.seqkit_stats.tsv`: `seqkit stats --all`
  for the final assembly.

</details>

BUSCO lineage is chosen from the samplesheet `class`; see the table in
[usage](usage.md#busco-lineage-selection).

### Genome size and coverage verdict

<details markdown="1">
<summary>Output files</summary>

- `draftgenomes/<OG>/coverage/`
  - `<prefix>_coverage_summary.json`: the per-sample verdict (fields below).
  - `<prefix>_sequencing_coverage.txt`: the same verdict as a readable report.
  - `<prefix>_kmer_depth.json`: coverage measured from read depth over Complete BUSCO genes, and the
    host/symbiont partition of the assembly.
- `coverage_summary/genome_coverage_summary.csv` (or `coverage_summary/<report_subdir>/`): one row
  per sample, compiled across the run.

</details>

The verdict is written twice. `CALCULATE_SEQUENCING_COVERAGE` runs before assembly and writes a
provisional summary with `"assembly_cross_check": "pending"`. `RECHECK_GENOME_SIZE` runs during QC,
adds the checks that need the assembly, and publishes over the same filename. If the recheck was
skipped or failed, the provisional file stays and `pending` says so.

`coverage_summary/genome_coverage_summary.csv` is compiled from the final summaries once QC has
run, so it matches the per-sample JSONs and the database. Its `assembly_cross_check` column says
which rows are final (`done`). Only when QC is skipped and no published QC results are available
is it built from the provisional summaries, and then that column reads `pending`. (Runs before
this change compiled it from the provisional summaries every time.)

Key fields in `<prefix>_coverage_summary.json`:

| Field | Meaning |
| --- | --- |
| `estimated_genome_size`, `genome_size_min/max`, `model_fit_*`, `heterozygosity_*`, `kmercov` | From the GenomeScope2 fit that won (original or reseeded). |
| `theoretical_coverage`, `coverage_before_filtering`, `data_retention_percent` | Read bases over estimated genome size, before and after fastp. |
| `coverage_status` | `EXCELLENT`, `GOOD`, `ADEQUATE`, `LOW` or `INSUFFICIENT`, or `UNRELIABLE_GENOME_SIZE_ESTIMATE` when any check below fails. |
| `coverage_recommendation` | Plain-language next step. |
| `genome_size_reliable`, `genome_size_flags` | Whether the size estimate passed, and the reasons if not (e.g. `no_kmer_peak`, `genome_size_disagrees_with_assembly`, `fitted_coverage_disagrees_with_read_depth`, `assembly_lacks_host_gene_space`). |
| `assembly_cross_check`, `assembly_size`, `assembly_ratio_basis`, `assembly_ratio_denominator` | `pending` or `done`; which size the estimate was compared with (`host` or `total`) and the number used. |
| `lambda_depth`, `kmercov_over_lambda_depth` | Haploid coverage measured from BUSCO read depth, and GenomeScope's fit relative to it. |
| `busco_complete` | BUSCO C% used by the gene-space check. |
| Host partition fields | See below. |

The per-sample `*_kmer_depth.json` also carries a genome size that does not come from the k-mer
histogram at all. The k-mer route collapses below about 10x host coverage, so
`SUMMARISE_KMER_COVERAGE` partitions the assembly by per-contig read depth around the depth of the
host's own single-copy BUSCO genes: `host_assembly_size` is the mass sitting at host depth,
`symbiont_assembly_size` the mass well above it, `low_depth_assembly_size` the mass below it, and
`host_assembly_fraction` how much of the assembly is the animal. `partition_status` says whether the
partition was possible (`ok`, `no_host_depth`, `no_contig_depth`). Both sizes are **lower bounds**: a
symbiont that happens to sit at the host's depth counts as host, and a high-copy host repeat counts
as symbiont. Where the partition succeeded, the genome size cross-check divides by the host mass
rather than the whole assembly, and records which it used in `assembly_ratio_basis`.

Two refinements to that partition:

- **Haplotype-split loci are counted once.** Within the host band, contigs below 0.75x host depth are
  reported as `host_haplotig_size` and the rest as `host_collapsed_size`. Where heterozygosity is
  high, the assembler keeps both haplotypes of a locus as separate contigs at about half depth, so
  `host_assembly_size` counts that locus twice. `host_haploid_size` (collapsed + haplotig / 2) counts
  it once, and is the size the cross-check divides by.
- **The partition says when it cannot be trusted.** `depth_mode` is the length-weighted modal contig
  depth and `anchor_over_mode` is the BUSCO-gene depth divided by it. `partition_confidence` is `low`,
  with a reason in `partition_note`, when there are fewer than `--host_depth_min_regions` Complete
  BUSCO regions (`few_busco_regions`), no depth band holds 40% of the assembly (`no_dominant_mode`),
  or the BUSCO anchor sits outside the host window of the mode (`anchor_off_mode`). A low-confidence
  partition is not used by the cross-check, which falls back to the total assembly. Where the anchor
  is off a dominant mode, `host_assembly_size_mode_anchored` reports the same partition re-centred on
  that mode, as a second figure rather than a replacement.

The final verdict also checks gene space: an assembly whose BUSCO completeness is below
`--genomescope_min_busco_complete` (10%) is flagged `assembly_lacks_host_gene_space`. Without it, a
failed library whose size estimate matches its own tiny assembly can pass every other check.

The thresholds behind each flag are listed in [usage](usage.md#genome-size-and-coverage-checks).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/` (or `multiqc/<report_subdir>/`)
  - `multiqc_report.html`: run-level report across all samples in the launch.
  - `multiqc_data/`: parsed statistics behind the report.
- `draftgenomes/<OG>/multiqc/`
  - `<OG>_multiqc_report.html`, `<OG>_multiqc_report_data/`: the same report for one sample. This is
    the copy the final backup puts on S3.

</details>

Both reports include the tool-parameters table built from the `*.tool_params_mqcrow.html` rows, the
kraken2 retention, coverage verdict and genome size recheck sections, and software versions.

### SQL upload receipts

<details markdown="1">
<summary>Output files</summary>

- `push/<OG>.{fastp,assembly,decontamination,busco,merqury,gfastats}.upload.txt`: one receipt per
  metric family per sample, recording what was written to the OceanOmics database.

</details>

Produced unless `--skip_upload_results`. A missing receipt means that family was not uploaded for
that sample.

### Mitogenome launcher

<details markdown="1">
<summary>Output files</summary>

- `<mitogenome_nfcore_dir>/run_<RUN>_mitogenomes_<YYMMDD_HHMMSS>.sh`: written into the mitogenome
  pipeline checkout, **not** the outdir.

</details>

`TRIGGER_MITOGENOME` writes, but does not run, a launcher for the mitogenome pipeline. It sets the
`OG_*` run values (run ID, `/scratch/pawsey1348/$USER/<RUN>_mitogenomes` as the outdir, the
`draftgenomes/*/fastp/*.fastq.gz` reads as input, BaseSpace and SQL configs) and then appends a
verbatim copy of that repo's `nextflow_run.sh`. Run it from the mitogenome checkout once the
draft run's fastp step has finished.

### Pipeline information and compute cost

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/` (or `pipeline_info/<report_subdir>/` for `software_versions.yml`)
  - `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`,
    `execution_trace_<timestamp>.txt`, `pipeline_dag_<timestamp>.html`: Nextflow reports, one set per
    launch (a resumed run adds a new set rather than overwriting).
  - `params_<timestamp>.json`: parameters for each launch.
  - `software_versions.yml`: collated tool versions.
  - `compute_usage.csv`, `compute_usage.report.txt`, `cost_per_sample.csv`: Slurm SU usage per job,
    a readable summary, and SU per sample with shared overhead amortised. Written by the template's
    tail via [`compute-audit/nf_workflow_cost.sh`](../compute-audit/README.md); absent for runs not
    launched from the template, or when `OCEANOMICS_SKIP_COST` is set.
  - `pipeline_report.{html,txt}`: only when `--email` / `--email_on_fail` is used.

</details>

## What gets backed up

The final backup ([internal SOP section 7](internal_sop.md#7-final-genome-backup)) copies a subset of
each `draftgenomes/<OG>/` to S3 (`fastp/*.fastq.gz`, `fastp/*fastp.json`, `assemblies/genome/*.fna`,
`multiqc/*multiqc_report.html`), archives each `.meryl` directory, then **moves** the `draftgenomes/`
tree, except the `fastp/` folders, to Acacia. After it completes only `draftgenomes/<OG>/fastp/` is
left on scratch.
