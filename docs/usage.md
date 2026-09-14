# OceanOmics-OceanGenomes-Draft-Genomes: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

This pipeline supports two entry points:

- **Download from Illumina BaseSpace**: provide `--run` and a BaseSpace CLI config (`--bs_config`). The workflow pulls datasets for that run ID, re-pairs lanes, and builds a samplesheet automatically using metadata looked up via `--sql_config`.
- **Use existing FASTQs**: provide `--input` with a validated samplesheet and set `--skip_download_reads true`. This is useful when FASTQs are already staged or come from outside BaseSpace.

For standard OceanGenomes runs on Pawsey, the simplest route is to copy [`nextflow_run_template.sh`](../nextflow_run_template.sh) to a run-specific launcher, update the `RUN` variable in that copied script, and run it from the repository root:

```bash
RUN=NEXT_250724_ET
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
bash "nextflow_run_${RUN}.sh"
```

The template is set up for the OceanGenomes project defaults: it loads Nextflow and Singularity modules, creates `/scratch/pawsey0964/$USER/$RUN`, copies the backup helper scripts into the run directory, stamps the backup config with the selected run ID, changes into the run output directory, and launches the workflow with the project BaseSpace, SQL, BUSCO, contamination-screening, mitogenome, temp-directory, and Pawsey profile settings.

When using pre-existing FASTQs:

```bash
nextflow run main.nf \
  -profile singularity \
  -c pawsey_profile.config \
  --input assets/samplesheet.csv \
  --skip_download_reads true \
  --outdir /scratch/pawsey0964/$USER/oceangenomesdraftgenomes
```

If your environment uses different contamination/BUSCO databases or a different temporary directory, set `--gxdb`, `--busco_acti_db`, `--busco_vert_db`, and `--tempdir` accordingly (see `nextflow_run*.sh` for a full example).

## Samplesheet input

If you supply `--run`, the pipeline will create a samplesheet for you and place it under `samplesheet/<RUN>_samplesheet.csv` in your `--outdir` using metadata from `--sql_config`. To supply your own, point `--input` at a CSV that matches `assets/schema_input.json`.

Required header (order fixed):

```
sample,run,date,prefix,nom_species_id,taxon_id,class,fastq_1,fastq_2
```

```csv title="samplesheet.csv"
sample,run,date,prefix,nom_species_id,taxon_id,class,fastq_1,fastq_2
OG747,NOVA_250131_AD,250131,OG747.ilmn.250131,70868,70868,Actinopteri,/data/OG747.ilmn.250131.L002_R1.fastq.gz,/data/OG747.ilmn.250131.L002_R2.fastq.gz
OG747,NOVA_250131_AD,250131,OG747.ilmn.250131,70868,70868,Actinopteri,/data/OG747.ilmn.250131.L003_R1.fastq.gz,/data/OG747.ilmn.250131.L003_R2.fastq.gz
OG846,NOVA_250131_AD,250131,OG846.ilmn.250131,13397,13397,Chondrichthyes,/data/OG846.ilmn.250131.L002_R1.fastq.gz,/data/OG846.ilmn.250131.L002_R2.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | OG identifier; use the same value for multiple lanes/runs of the same specimen.                                                                                                        |
| `run` / `date` / `prefix` | Run metadata (prefix usually follows `OG###.ilmn.<date>` and is used to name downstream files).                                                                        |
| `nom_species_id`, `taxon_id`, `class` | Taxonomic metadata used for BUSCO lineage selection and reporting; populate from your LIMS/SQL source or set to `unknown` if not available.                  |
| `fastq_1` | Full path to R1 FASTQ (`.R1.fastq.gz`/`.R1.fq.gz`).                                                                                                                                    |
| `fastq_2` | Full path to R2 FASTQ (`.R2.fastq.gz`/`.R2.fq.gz`).                                                                                                                                    |

The pipeline concatenates multiple rows with the same `sample` before processing.

The `class` value also controls BUSCO lineage selection during genome QC:

| `class` value | BUSCO database parameter |
| ------------- | ------------------------ |
| `Actinopteri`, `Actinopterygii`, `Teleostei` | `--busco_acti_db` |
| Other vertebrate classes (`Chondrichthyes`, `Mammalia`, `Aves`, `Reptilia`, `Amphibia`, ...) | `--busco_vert_db` |
| Anything else, i.e. all invertebrates | `--busco_metazoa_db` |
| `unknown` or empty | run aborts |

The default is metazoa, not vertebrata, so a new invertebrate class needs no code change to
be scored against a sensible lineage. The two vertebrate lists live at the top of the BUSCO
selection block in `subworkflows/local/genome_qc/main.nf`; adding a class is a one-line edit
there.

A sample whose `class` is `unknown` (or empty) aborts the run rather than being guessed at.

The samplesheet is still written and published in that case. `bin/create_samplesheet.py`
marks the offending rows `unknown` and exits successfully, so the sheet and its
`taxonomy_resolution.tsv` land in `<outdir>/samplesheet/` where you can see them; the run is
then stopped by the pipeline before any assembly task is submitted. Nothing is assembled and
no SUs are spent until the rows are fixed. (The stop has to happen at that point rather than
in the script: an `unknown` reaching FCS-GX as `--tax-id unknown` fails only after MEGAHIT
has already run.)

To fix, either correct the `nominal_species_id` in the `sample` table and delete the
published samplesheet so it is regenerated, or edit the `taxon_id` and `class` columns in the
published samplesheet directly and re-run. A published samplesheet is reused automatically on
the next run, so an edited one is picked up without needing `--input`.

### Resuming a run that already has MEGAHIT checkpoints

MEGAHIT keeps a resumable checkpoint per sample under `<outdir>/megahit_checkpoints` so a
task killed by an OOM or a walltime does not start over. Each one now carries a
fingerprint of the reads, assembly arguments and megahit version it was built from, and a
checkpoint whose fingerprint no longer matches is rebuilt rather than reused. That is what
stops a change upstream of assembly, enabling kraken2 for example, from silently
republishing the previous assembly with every downstream metric describing it.

Checkpoints written before fingerprinting existed carry no fingerprint, and their inputs
cannot be known. The default `--megahit_checkpoint_unkeyed invalidate` discards them and
reassembles, which is the safe reading but costs a full reassembly. **On the first resume
of a run whose reads have not changed, pass `--megahit_checkpoint_unkeyed adopt`** to claim
those checkpoints instead. That is an assertion, so only make it when nothing upstream of
assembly has changed. After that first pass every checkpoint is fingerprinted and the flag
stops mattering.

`<prefix>.megahit_checkpoint.txt` in each sample's assembly directory records which path
the task took, so whether an assembly was rebuilt is something to read rather than infer.

### Where taxonomy comes from

`taxon_id` and `class` are resolved in two steps, and the per-sample outcome is recorded in
`taxonomy_resolution.tsv`, published next to the samplesheet:

| `source` | meaning |
| --- | --- |
| `db` | the curated `species` table had the answer |
| `db+taxdump` | the table had part of it, the taxdump filled the rest |
| `taxdump` | the table had no row; resolved entirely from NCBI |
| `unresolved` | neither source knew the name -- the run aborts on these |

The curated table is authoritative and always tried first. It only holds taxa someone has
loaded, though, and the invertebrate runs draw from most of Metazoa, so pass
`--taxonkit_db_dir <dir>` to enable the NCBI taxdump fallback. The dump is downloaded once
and cached there with `storeDir`; point it at the same directory the mitogenome pipeline
uses and the two share one copy. Without the flag, any sample the `species` table does not
carry aborts the run, which is the pre-fallback behaviour.

The fallback matches the `nominal_species_id` against NCBI scientific names at any rank
within Metazoa, so genus, family, order and class names all resolve, as do the
`Asteroidea (Class)` style values the collection records. Where NCBI has no class rank for a
lineage (e.g. Porifera), the phylum is used as the class, matching what
`scripts/taxonomy/load_taxonomy.py` writes into `species.class`.

What it cannot rescue is a name that is not a taxon at all -- a misspelling (`Actinaria` for
`Actiniaria`) or a placeholder (`Larval fish`). Those are listed by the failing step and need
the `nominal_species_id` corrected in the `sample` table.

### Loading taxonomy for new clades

Still the right tool when the name is correct and you want it curated permanently, and the
only option if you are running without `--taxonkit_db_dir`.

`scripts/taxonomy/load_taxonomy.py` bulk-loads NCBI species-rank taxa into the `species`
table from the NCBI `new_taxdump`. Target one or more classes and/or phyla:

```bash
# dry-run: writes a CSV preview, touches nothing
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata,Arthropoda,Porifera

# after inspecting the preview
python3 scripts/taxonomy/load_taxonomy.py --phylum Mollusca,Echinodermata --apply ~/postgresql_details/oceanomics.cfg
```

Use `--phylum` rather than `--class` for invertebrates: NCBI leaves the class column empty
for many invertebrate lineages, and a class-only load silently misses them. Where a matched
row has no NCBI class, the phylum name is written into `species.class`, which routes to the
metazoa BUSCO database as intended. Existing rows are never modified
(`ON CONFLICT (species) DO NOTHING`).

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline (fill in FASTQ paths before use).

## Backfilling database statistics from Acacia

`bin/backfill_draft_genome_stats.py` reads only small report files from mounted
object-storage archives. The `genomes.v2` Acacia backup intentionally excludes
`fastp/`; mount `s3:oceanomics/OceanGenomes/analysed-data/draft-genomes` as a
second read-only tree and pass it with `--fastp-root`. The command does not stage
FASTQs, assemblies, or Meryl databases. Run it inside a Slurm compute allocation
where both mounts and the PostgreSQL service are reachable.

Create a tab-separated manifest with one unambiguous database key per row. A
template is available at [`assets/backfill_manifest.tsv`](../assets/backfill_manifest.tsv).

```tsv
og_id	seq_date
OG123	250101
OG456	250205
```

First validate the complete set. Validation is the default and does not connect
to or modify PostgreSQL:

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

Do not add `--apply` until `inventory.tsv` contains only `PASS` rows. Make a
one-row canary manifest, apply it, and inspect `upload.tsv` and
`verification.tsv` before applying the full manifest:

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

All six metric families for one `(og_id, seq_date)` are written and verified in
one transaction. A parsing, SQL, or source-to-database comparison failure rolls
back that sample. Successful and failed samples can safely be rerun because the
writes use `ON CONFLICT (og_id, seq_date) DO UPDATE`.

## Running the pipeline

The recommended command for OceanGenomes BaseSpace runs is:

```bash
RUN=NEXT_250724_ET
cp nextflow_run_template.sh "nextflow_run_${RUN}.sh"
sed -i "s/^RUN=.*/RUN=${RUN}/" "nextflow_run_${RUN}.sh"
bash "nextflow_run_${RUN}.sh"
```

Use the plain [`nextflow_run.sh`](../nextflow_run.sh) script as an alternative when you want to launch from the repository directory and keep the Nextflow work directory under `./work/$RUN`. In that mode, update `RUN` and `OUT` in `nextflow_run.sh` before launching.

To reuse existing FASTQs instead of downloading:

```bash
nextflow run main.nf \
  -profile singularity \
  -c pawsey_profile.config \
  --input ./samplesheet.csv \
  --skip_download_reads true \
  --outdir ./results
```

See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for process resource settings, other infrastructural tweaks, or module arguments (args). Use `-params-file` for pipeline parameters.

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run main.nf -profile singularity -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
skip_download_reads: true
<...>
```

### Updating the pipeline

This pipeline is intended to be run from a local clone or working copy. To update it, update the repository checkout itself before launching Nextflow:

```bash
git pull
```

### Reproducibility

It is a good idea to record the pipeline commit or tag used for each production run. This ensures that a specific version of the pipeline code and software can be recovered later if results need to be reproduced.

You can record the current commit before launching a run:

```bash
git rev-parse --short HEAD
```

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share a parameter file, make sure it does not include private credentials, cluster-specific paths, or other local-only settings.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes test settings so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). Use this for infrastructure and executor settings, not pipeline parameters. See the [Nextflow config documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For many pipeline steps, failed jobs are automatically retried with higher resource requests according to the retry settings in [`conf/base.config`](../conf/base.config). If a step still fails after the configured retries, the pipeline execution is stopped.

To change resource requests, provide a custom Nextflow config with `-c` and override the relevant `process` selectors or labels.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline step for a particular tool. Many modules use containers and software from the [BioContainers](https://biocontainers.pro/) or [Bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline-specified version may be out of date.

To use a different container from the default container or conda environment specified in the pipeline, override the relevant process settings in a custom Nextflow config.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in the workflow. Where modules expose `ext.args` or similar settings, you can provide additional arguments via a custom Nextflow config.

See the local [`conf/modules.config`](../conf/modules.config) file for module-specific argument hooks already used by this pipeline.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
