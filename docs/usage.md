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
| `Actinopteri` | `--busco_acti_db` |
| `Anthozoa` or `Cnidaria` | `--busco_metazoa_db` |
| Any other value, including `unknown` | `--busco_vert_db` |

If you manually edit or provide a samplesheet, check that coral/cnidarian samples use `Anthozoa` or `Cnidaria` so they are not assigned the default vertebrate BUSCO database.

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline (fill in FASTQ paths before use).

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
