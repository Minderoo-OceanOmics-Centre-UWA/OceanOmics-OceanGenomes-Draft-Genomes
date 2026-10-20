# OceanOmics-OceanGenomes-Draft-Genomes: Usage

> Every parameter, with its type and default, is in [`nextflow_schema.json`](../nextflow_schema.json)
> and is printed by `nextflow run main.nf --help` (add `--show_hidden` for all of them). The
> reasoning behind the QC defaults is in the comments in [`nextflow.config`](../nextflow.config).
> This page covers how to run the pipeline and what the important options do.

## Introduction

The pipeline has two entry points:

- **Download from Illumina BaseSpace**: provide `--run` and a BaseSpace CLI config (`--bs_config`).
  The workflow pulls every dataset for that run, merges lanes, re-pairs reads, and builds a
  samplesheet from the `sample` and `species` tables in the database (`--sql_config`).
- **Use existing FASTQs**: provide `--input` with a samplesheet and set `--skip_download_reads true`.

## Standard OceanGenomes run on Pawsey

Copy [`nextflow_run_template.sh`](../nextflow_run_template.sh) to a run-specific launcher, set `RUN`,
and run it **from the repository root**:

```bash
RUN=NOVA_260925_LA
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
bash "nextflow_run_${RUN}.sh"
```

The launcher:

- loads `nextflow/25.04.6` and `singularity/4.1.0-nompi`;
- creates the outdir `/scratch/pawsey1348/$USER/$RUN`, copies `scripts/backup_scripts` into it and
  stamps `RUN` into the copied `backup_scripts/configfile.txt`;
- changes into the outdir and launches from there, so several runs can share one checkout. Work files
  go to `$OUT/work/$RUN` and the log to `$OUT/.nextflow_$RUN.log`;
- passes the project database, BaseSpace, FCS-GX, kraken2, BUSCO, taxdump and mitogenome settings
  (see [Reference data and configs](#reference-data-and-configs));
- after Nextflow exits, fails the launcher if the log shows `Disk quota exceeded` or an uncreatable
  work directory (Nextflow can exit 0 in that case with late steps never run);
- writes per-sample compute cost into `pipeline_info/` (see [Compute cost](#compute-cost)).

Run it inside `tmux` so it survives a disconnect. The internal step-by-step procedure, including
backups, is in the [internal SOP](internal_sop.md).

[`nextflow_run.sh`](../nextflow_run.sh) is an older alternative that launches from the repository
directory with a fixed `OUT` (`/scratch/pawsey1348/$USER/_NFCORE/_OUT_DIR`) and `./work/$RUN`. Edit
`RUN`, `OUT` and `--mitogenome_nfcore_dir` before using it. It does not copy backup scripts or record
compute cost.

## Using existing FASTQs

```bash
nextflow run main.nf \
  -profile singularity \
  -c pawsey_profile.config \
  --input samplesheet.csv \
  --skip_download_reads true \
  --run NOVA_260925_LA \
  --outdir /scratch/pawsey1348/$USER/NOVA_260925_LA \
  <reference data options as in the template>
```

To rerun a run that was downloaded earlier, keep the template and change the tail of the command:
set `--skip_bs_download true` and `--skip_download_reads true`, and uncomment the
`--input "$OUT/samplesheet/${RUN}_samplesheet.csv"` line. Keep `--run` so names stay tied to the run.

## Samplesheet input

With `--run` and no `--input`, the pipeline builds the samplesheet and publishes it as
`<outdir>/samplesheet/<RUN>_samplesheet.csv`. On a later launch it **reuses** that file if it exists,
so hand edits stick. Delete it to have it regenerated. With `--input`, that file is used instead.

Columns are validated against [`assets/schema_input.json`](../assets/schema_input.json). Header, in
this order:

```csv title="samplesheet.csv"
sample,run,date,prefix,nom_species_id,taxon_id,class,fastq_1,fastq_2
OG1408,NOVA_251031_LMAD,251031,OG1408.ilmn.251031,Balaenoptera acutorostrata,9767,Mammalia,/path/OG1408.ilmn.NOVA_251031_LMAD.R1.fq.gz,/path/OG1408.ilmn.NOVA_251031_LMAD.R2.fq.gz
OG2622,NOVA_260925_LA,260925,OG2622.ilmn.260925,Amphinomidae (family),104701,Polychaeta,/path/OG2622.ilmn.NOVA_260925_LA.R1.fq.gz,/path/OG2622.ilmn.NOVA_260925_LA.R2.fq.gz
```

| Column | Description |
| --- | --- |
| `sample` | OG identifier (`OG` + digits, optional letter / `-n` / `_HIC` suffix). Rows sharing a `sample` are merged. |
| `run` | Sequencing run ID. |
| `date` | Six-digit date code from the run (`YYMMDD`). |
| `prefix` | File-name prefix for everything downstream, usually `OG###.ilmn.<date>`. |
| `nom_species_id` | Nominal species name from the `sample` table (a scientific name, not an ID). |
| `taxon_id` | NCBI taxid, or `unknown`/`None`. Passed to FCS-GX. |
| `class` | Taxonomic class; selects the BUSCO lineage. `unknown` or empty aborts the run. |
| `fastq_1`, `fastq_2` | Full paths; must end `.R1.fastq.gz`/`.R1.fq.gz` and `.R2.fastq.gz`/`.R2.fq.gz`. |

An [example samplesheet](../assets/samplesheet.csv) is provided (replace the FASTQ paths).

### BUSCO lineage selection

`class` picks the BUSCO database (`buscoDbForClass` in
[`subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline/main.nf`](../subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline/main.nf)).
Rows are tested top to bottom:

| `class` | Database parameter |
| --- | --- |
| `Actinopteri`, `Actinopterygii`, `Teleostei` | `--busco_acti_db` |
| `Chondrichthyes`, `Mammalia`, `Aves`, `Reptilia`, `Amphibia`, `Myxini`, `Hyperoartia`, `Coelacanthimorpha`, `Dipneusti`, `Lepidosauria`, `Testudines` | `--busco_vert_db` |
| `Malacostraca`, `Thecostraca`, `Branchiopoda`, `Copepoda`, `Ostracoda`, `Maxillopoda` | `--busco_crustacea_db` |
| `Gastropoda`, `Bivalvia`, `Polyplacophora`, `Cephalopoda`, `Scaphopoda`, `Monoplacophora` | `--busco_mollusca_db` |
| `Anthozoa` | `--busco_anthozoa_db` |
| `Pycnogonida`, `Merostomata`, `Chilopoda`, `Diplopoda`, `Symphyla`, `Pauropoda` | `--busco_arthropoda_db` |
| anything else (all other invertebrates) | `--busco_metazoa_db` |
| `unknown` or empty | run aborts |

The four invertebrate lineages are optional: if one is not set, its classes fall back to metazoa.
There is no echinoderm or sponge lineage in odb12, so Asteroidea, Ophiuroidea, Demospongiae,
Hexactinellida and Porifera are scored against metazoa by necessity. To route a new class, add it to
the matching list in `buscoDbForClass`.

### When a sample's class is unknown

A sample whose `class` is `unknown` (or empty) aborts the run rather than being guessed at. The
samplesheet and its `taxonomy_resolution.tsv` are still published to `<outdir>/samplesheet/`, and the
run stops before any assembly task is submitted, so no SUs are spent. (The stop happens in the
pipeline rather than in `bin/create_samplesheet.py` because an `unknown` reaching FCS-GX as
`--tax-id unknown` fails only after MEGAHIT has already run.)

To fix, either correct the `nominal_species_id` in the `sample` table and delete the published
samplesheet so it is regenerated, or edit `taxon_id` and `class` in the published samplesheet and
re-run.

### Where taxonomy comes from

`taxon_id` and `class` are resolved in two steps, and the per-sample outcome is recorded in
`taxonomy_resolution.tsv`, next to the samplesheet:

| `source` | Meaning |
| --- | --- |
| `db` | the curated `species` table had the answer |
| `db+taxdump` | the table had part of it, the taxdump filled the rest |
| `taxdump` | the table had no row; resolved entirely from NCBI |
| `unresolved` | neither source knew the name; the run aborts on these |

The curated table is authoritative and always tried first. It only holds taxa someone has loaded, so
pass `--taxonkit_db_dir <dir>` to enable the NCBI taxdump fallback. The dump is downloaded once and
cached there with `storeDir`; point it at the same directory the mitogenome pipeline uses and the two
share one copy (the template uses `/scratch/pawsey1348/$USER`). Without the flag, any sample the
`species` table does not carry aborts the run.

The fallback matches `nominal_species_id` against NCBI scientific names at any rank within Metazoa,
so genus, family, order and class names all resolve, as do the `Asteroidea (Class)` style values the
collection records. Where NCBI has no class rank for a lineage (e.g. Porifera), the phylum is used as
the class, matching what `scripts/taxonomy/load_taxonomy.py` writes into `species.class`.

It cannot rescue a name that is not a taxon at all: a misspelling (`Actinaria` for `Actiniaria`) or a
placeholder (`Larval fish`). Those are listed by the failing step and need the `nominal_species_id`
corrected in the `sample` table.

### Loading taxonomy for new clades

Use this when the name is correct and you want it curated permanently, or when running without
`--taxonkit_db_dir`. [`scripts/taxonomy/load_taxonomy.py`](../scripts/taxonomy/load_taxonomy.py)
bulk-loads NCBI species-rank taxa into the `species` table from the NCBI `new_taxdump`:

```bash
# dry-run: writes a CSV preview, touches nothing
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata,Arthropoda,Porifera

# after inspecting the preview
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata --apply ~/postgresql_details/oceanomics.cfg
```

Use `--phylum` rather than `--class` for invertebrates: NCBI leaves the class column empty for many
invertebrate lineages, and a class-only load silently misses them. Where a matched row has no NCBI
class, the phylum name is written into `species.class`, which routes to the metazoa BUSCO database.
Existing rows are never modified (`ON CONFLICT (species) DO NOTHING`).

### Reviewing the samplesheet before assembly

`--stop_after_samplesheet true` downloads the reads, builds and publishes the samplesheet (checking
its taxonomy) and then stops without assembling anything. Edit
`<outdir>/samplesheet/<RUN>_samplesheet.csv` (for example to drop samples not worth the SUs), then
launch again without the flag and with `-resume`; the edited sheet and the download are reused.

## Reference data and configs

Values used by the template on Setonix:

| Parameter | Purpose | Template value |
| --- | --- | --- |
| `--bs_config` | BaseSpace CLI config | `~/.basespace/default.cfg` |
| `--sql_config` | OceanOmics PostgreSQL config (samplesheet lookup, result upload) | `~/postgresql_details/oceanomics.cfg` |
| `--gxdb` | FCS-GX database | `/scratch/references/Foreign_Contamination_Screening` |
| `--ramdisk_path` | where FCS-GX stages its database in memory | `/tmp/gxdb/` |
| `--kraken2_db` | kraken2 database (`hash.k2d`, `opts.k2d`, `taxo.k2d`) | `/scratch/references/kraken2/pluspfp_20230605` |
| `--busco_acti_db`, `--busco_vert_db` | fish / other vertebrate lineages | `/scratch/references/busco_db/{actinopterygii,vertebrata}_odb10` |
| `--busco_metazoa_db` | default invertebrate lineage | `/software/projects/pawsey0964/busco_db/metazoa_odb12` |
| `--busco_{crustacea,mollusca,anthozoa,arthropoda}_db` | optional invertebrate lineages | `/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/*_odb12` |
| `--taxonkit_db_dir` | NCBI taxdump cache | `/scratch/pawsey1348/$USER` |
| `--mitogenome_nfcore_dir` | mitogenome pipeline checkout | `/software/projects/pawsey1348/$USER/repos/Oceanomics-OceanGenomes-Mitogenomes` |
| `--kvalue` | k for meryl/GenomeScope | `21` |
| `--tempdir` | accepted for compatibility; not currently used by any step | `$BASE/tmp` |
| `--binddir` | extra Singularity bind path (`-B`), if a database lives outside the default binds | not set |

Slurm job settings (account, queue size, Singularity cache) are in
[`pawsey_profile.config`](../pawsey_profile.config).

## Read-level decontamination (kraken2)

kraken2 classifies the trimmed reads and KrakenTools drops pairs assigned to the excluded taxa
before meryl and MEGAHIT, so both the k-mer profile and the assembly are built from host reads.
Screening contigs after assembly cannot undo a host/symbiont co-assembly.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--kraken2_db` | none | Required unless kraken2 is skipped. |
| `--kraken2_confidence` | `0.05` | kraken2 `--confidence`. Its own default of 0 assigns on one k-mer and deletes host reads from poorly covered taxa. |
| `--kraken2_exclude_taxids` | `2 2157 10239` | Taxa removed with their descendants (Bacteria, Archaea, Viruses). Unclassified reads are kept. |
| `--skip_kraken2_decontamination` | `false` | Assemble the fastp reads directly. |

Changing any of these changes the reads MEGAHIT sees, so existing MEGAHIT checkpoints are rebuilt
(see below).

## Assembly decontamination

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--fcs_review_action` | `exclude` | FCS-GX REVIEW contigs: `exclude` removes all, `short-only` removes only those of 1000 bp or less (pre-Sept 2026 behaviour), `keep` retains them. |
| `--tiara_min_len` | `1000` | Tiara `--min_len`; shorter contigs are not classified. Tiara's own 3000 leaves most of a fragmented assembly unscreened. |

Contigs under 500 bp are always dropped after FCS-GX.

## Genome size and coverage checks

The GenomeScope2 genome size estimate is gated before coverage is graded. When any check fails the
verdict is `UNRELIABLE_GENOME_SIZE_ESTIMATE` and the reasons are listed in `genome_size_flags` (see
[output](output.md#genome-size-and-coverage-verdict)).

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--genomescope2_m` | `1000` | GenomeScope2 `-m` (max k-mer coverage in the fit). Raising it inflates repeat length without recovering host signal. |
| `--genomescope2_l` | `false` | Optional GenomeScope2 `-l` seed; adds `.l<value>` to file names. |
| `--genomescope_min_model_fit` | `50` | Minimum model fit (upper value, %). |
| `--genomescope_max_model_bias` | `5` | Maximum bias coefficient. |
| `--genomescope_min_peak_kmercov_ratio`, `--genomescope_max_peak_kmercov_ratio` | `0.5`, `3` | Allowed ratio of the observed histogram mode to fitted kmercov. |
| `--genomescope_peak_search_max` | `1000` | Coverage ceiling for k-mer peak detection. |
| `--genomescope_min_assembly_ratio`, `--genomescope_max_assembly_ratio` | `0.3`, `2.0` | Allowed estimate / assembly size (host mass where the partition is trusted). |
| `--genomescope_min_busco_complete` | `10` | BUSCO C% below which the assembly lacks host gene space. |
| `--kmer_depth_min_mapq` | `20` | MAPQ floor for read depth over Complete BUSCOs. |
| `--kmer_depth_min_regions` | `20` | Minimum usable BUSCO regions to report a measured coverage. |
| `--kmer_depth_max_ratio` | `1.5` | Allowed fitted kmercov / measured coverage (either direction). |
| `--kmer_depth_min_seedable` | `5.0` | Below this measured coverage no reseed is attempted. |
| `--host_depth_window` | `2.0` | Host depth band half-width (multiplicative). |
| `--host_depth_min_contig` | `500` | Contigs shorter than this are left out of the partition. |
| `--host_depth_min_regions` | `50` | Fewer BUSCO regions and the partition is low confidence. |
| `--genomescope_reseed_when_flagged` | `true` | Refit GenomeScope seeded with the measured coverage when flagged. |
| `--genomescope_min_reseed_fit_gain` | `2` | Model-fit points a refit must gain (and it must move kmercov toward the measurement) to be kept. |

## Resuming a run that already has MEGAHIT checkpoints

MEGAHIT keeps a resumable checkpoint per sample under `<outdir>/megahit_checkpoints` so a task
killed by an OOM or a walltime does not start over. Each one carries a fingerprint of the reads,
assembly arguments and MEGAHIT version it was built from, and a checkpoint whose fingerprint no
longer matches is rebuilt rather than reused. That stops a change upstream of assembly (enabling
kraken2, for example) from silently republishing the previous assembly.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--megahit_checkpoint_unkeyed` | `invalidate` | Checkpoints written before fingerprinting: `invalidate` reassembles, `adopt` claims them for the current inputs. |
| `--megahit_stale_checkpoint` | `rerun` | Fingerprint mismatch: `rerun`, `archive` (move aside to `<dir>.stale.<timestamp>` first) or `fail`. |
| `--megahit_checkpoint_cleanup` | `true` | Prune the intermediate graph once a checkpoint completes (tens of GB per sample). |

**On the first resume of an older run whose reads have not changed, pass
`--megahit_checkpoint_unkeyed adopt`** to avoid a full reassembly. That is an assertion, so only make
it when nothing upstream of assembly has changed. After that pass every checkpoint is fingerprinted
and the flag stops mattering. `<prefix>.megahit_checkpoint.txt` in each sample's assembly directory
records which path the task took.

## Skipping stages

Each stage can be skipped; the pipeline then reads that stage's results from the outdir using the
`--precomputed_*` globs, whose defaults point at the standard publish locations. Normally you only
set the `--skip_*` flags and leave the globs alone.

| Flag | Reads instead (default glob under `--outdir`) |
| --- | --- |
| `--skip_bs_download` | `--bs_fastq_glob`: `basespace/<RUN>/*/*.fastq.gz` |
| `--skip_download_reads` | `--precomputed_download_reads_results`: `pooled/<RUN>/*.{R1,R2}.fq.gz` (or the samplesheet given by `--input`) |
| `--skip_fastp_fastqc` | `--precomputed_fastp_fastqc_results`, `--precomputed_fastp_json`: `draftgenomes/*/fastp/` |
| `--skip_genome_assembly` | `--precomputed_genome_assembly_results` (`*mh.fasta`), `--precomputed_genomescope_summary`, `--precomputed_genomescope_model`, `--precomputed_meryl_hist`, `--precomputed_meryl_results`, `--precomputed_coverage_summary` |
| `--skip_genome_decontamination` | `--precomputed_genome_decontamination_results` (`*.fna`), `--precomputed_filter_report`, `--precomputed_contigs_under_500bp`, `--precomputed_tiara_filter_summary` |
| `--skip_genome_qc` | `--precomputed_genome_qc_results`, `--precomputed_busoco_short_summary_results`, `--precomputed_merqury_results_results`, `--precomputed_gfastats_results_results` |
| `--skip_upload_results` | nothing; the database is not written |
| `--skip_kraken2_decontamination` | nothing; assembles the fastp reads |

An empty glob logs a warning rather than failing, so check the log when skipping stages.

## Running several launches against one outdir

`--report_subdir <name>` publishes the run-level aggregates (`multiqc/`, `coverage_summary/` and
`pipeline_info/software_versions.yml`) under a subdirectory. Use it when launches that each cover
only some samples share an outdir, as the tiered re-runs in
[`scripts/rerun_NOVA_260909_LA`](../scripts/rerun_NOVA_260909_LA/README.md) do; otherwise each launch
overwrites the last one's partial report.

## Mitogenome hand-off

Once fastp has finished for all samples, `TRIGGER_MITOGENOME` writes
`run_<RUN>_mitogenomes_<timestamp>.sh` into `--mitogenome_nfcore_dir`. It is not run automatically.
The script sets the run-specific `OG_*` values (outdir `/scratch/pawsey1348/$USER/<RUN>_mitogenomes`,
input `<outdir>/draftgenomes/*/fastp/*.fastq.gz`, BaseSpace and SQL configs) and then appends a copy
of the mitogenome repo's own `nextflow_run.sh`, so launcher changes are made once, in that repo.
`--mitogenome_nfcore_dir` must contain `nextflow_run.sh` or the run fails at start.

## Compute cost

After Nextflow exits, the template runs
[`compute-audit/nf_workflow_cost.sh`](../compute-audit/README.md) against the outdir and writes
`pipeline_info/compute_usage.csv`, `compute_usage.report.txt` and `cost_per_sample.csv` (SU per
sample). It is best-effort and never changes the launcher's exit status. Set `OCEANOMICS_SKIP_COST=1`
to skip it (the automation poller does, and records cost itself).

## Backfilling database statistics from Acacia

`bin/backfill_draft_genome_stats.py` reads only small report files from mounted object-storage
archives. The `genomes.v2` Acacia backup intentionally excludes `fastp/`; mount
`s3:ocom.oceangenomes/analysed-data/draft-genomes` as a second read-only tree and pass it with
`--fastp-root`. The command does not stage FASTQs, assemblies, or Meryl databases. Run it inside a
Slurm compute allocation where both mounts and the PostgreSQL service are reachable.

Create a tab-separated manifest with one unambiguous database key per row. A template is available
at [`assets/backfill_manifest.tsv`](../assets/backfill_manifest.tsv).

```tsv
og_id	seq_date
OG123	250101
OG456	250205
```

First validate the complete set. Validation is the default and does not connect to or modify
PostgreSQL:

```bash
singularity exec \
  -B /path/to/acacia-mount,/path/to/s3-mount,/path/to/repository,/path/to/db-config,/path/to/reports \
  docker://tylerpeirce/psycopg2:0.1 \
  python3 /path/to/repository/bin/backfill_draft_genome_stats.py \
    --archive-root /path/to/acacia-mount/genomes.v2 \
    --fastp-root /path/to/s3-mount/draft-genomes \
    --manifest /path/to/manifest.tsv \
    --db-config /path/to/db-config/oceanomics.cfg \
    --report-dir /path/to/reports/validation
```

Do not add `--apply` until `inventory.tsv` contains only `PASS` rows. Make a one-row canary
manifest, apply it, and inspect `upload.tsv` and `verification.tsv` before applying the full
manifest:

```bash
singularity exec \
  -B /path/to/acacia-mount,/path/to/s3-mount,/path/to/repository,/path/to/db-config,/path/to/reports \
  docker://tylerpeirce/psycopg2:0.1 \
  python3 /path/to/repository/bin/backfill_draft_genome_stats.py \
    --archive-root /path/to/acacia-mount/genomes.v2 \
    --fastp-root /path/to/s3-mount/draft-genomes \
    --manifest /path/to/canary.tsv \
    --db-config /path/to/db-config/oceanomics.cfg \
    --report-dir /path/to/reports/canary \
    --apply
```

All six metric families for one `(og_id, seq_date)` are written and verified in one transaction. A
parsing, SQL, or source-to-database comparison failure rolls back that sample. Successful and failed
samples can safely be rerun because the writes use `ON CONFLICT (og_id, seq_date) DO UPDATE`.

## Parameter files

Parameters can be given in a YAML or JSON file with `-params-file`:

```bash
nextflow run main.nf -profile singularity -c pawsey_profile.config -params-file params.yaml
```

```yaml title="params.yaml"
run: 'NOVA_260925_LA'
outdir: '/scratch/pawsey1348/<user>/NOVA_260925_LA'
skip_download_reads: true
input: '/scratch/pawsey1348/<user>/NOVA_260925_LA/samplesheet/NOVA_260925_LA_samplesheet.csv'
```

> [!WARNING]
> Do not use `-c <file>` to set pipeline parameters. Config files given with `-c` are for resources,
> executor settings and module `ext.args` only.

## Updating the pipeline and reproducibility

The pipeline runs from a local checkout; update it with `git pull`. The version printed at start-up
and recorded in `pipeline_info/` comes from `git describe --tags --always`, so record it (or
`git rev-parse --short HEAD`) for each production run. Each launch also writes
`pipeline_info/params_<timestamp>.json`.

## Core Nextflow arguments

> [!NOTE]
> These are Nextflow options and take a _single_ hyphen.

- `-profile`: use `singularity` on Setonix, together with `-c pawsey_profile.config`. `test` and
  `test_full` profiles exist ([`conf/test.config`](../conf/test.config),
  [`conf/test_full.config`](../conf/test_full.config)); the generic nf-core container profiles
  (`docker`, `apptainer`, `conda`, ...) are defined in `nextflow.config` but are not used in
  production.
- `-resume`: reuse cached tasks whose inputs are unchanged. The template always passes it. Use
  `nextflow log` to list previous run names.
- `-c`: extra config for resources or executor settings.

## Custom configuration

- **Resources**: defaults and retry behaviour are in [`conf/base.config`](../conf/base.config). Most
  per-sample steps retry on resource kills and are then ignored so one failed genome does not stop
  the rest; run-level steps (`CREATE_SAMPLESHEET`, `DOWNLOAD_TAXONKIT_DB`, `BASESPACE`,
  `COMPILE_JSON_TO_CSV`, `MULTIQC`) stop the run. Override with a `-c` config using `withName`
  selectors. On Setonix memory, not CPUs, drives the SU charge.
- **Tool arguments**: module `ext.args` and publish locations are in
  [`conf/modules.config`](../conf/modules.config). Anything added to MEGAHIT's `ext.args` changes the
  checkpoint fingerprint.

## Nextflow memory

If the Nextflow JVM itself uses too much memory, limit it in your environment:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
