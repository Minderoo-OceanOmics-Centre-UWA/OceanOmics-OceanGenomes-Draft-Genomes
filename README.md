<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-oceangenomesdraftgenomes_logo_dark.png">
    <img alt="nf-core/oceangenomesdraftgenomes" src="docs/images/nf-core-oceangenomesdraftgenomes_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/oceangenomesdraftgenomes/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/oceangenomesdraftgenomes/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/oceangenomesdraftgenomes/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/oceangenomesdraftgenomes/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/oceangenomesdraftgenomes/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/oceangenomesdraftgenomes)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23oceangenomesdraftgenomes-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/oceangenomesdraftgenomes)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/oceangenomesdraftgenomes** builds short-read draft assemblies for Ocean Omics / Ocean Genomes samples. It can download Illumina runs directly from BaseSpace or consume existing FASTQs, performs adapter/quality trimming, estimates genome size and coverage, assembles with MEGAHIT, screens contaminants (NCBI FCS-GX, FCS adaptor, Tiara + BBMap), runs genome QC (BUSCO, BWA-MEM2, Merqury, Gfastats), and collates results with MultiQC. Optional steps push QC summaries into the project SQL database.

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
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
  -profile singularity \
  -c pawsey_profile.config \
  --run NEXT_250724_ET \
  --bs_config ~/.basespace/default.cfg \
  --sql_config ~/postgresql_details/oceanomics.cfg \
  --outdir /scratch/pawsey0964/$USER/oceangenomesdraftgenomes
```

Use `--input <samplesheet.csv> --skip_download_reads true` if you already have FASTQs instead of pulling from BaseSpace. The pipeline needs paths for the contamination DB (`--gxdb`), BUSCO DBs (`--busco_acti_db` / `--busco_vert_db`), and a temp directory suited to your system; see the provided `nextflow_run*.sh` templates for a full set of Pawsey-friendly flags.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/oceangenomesdraftgenomes/usage) and the [parameter documentation](https://nf-co.re/oceangenomesdraftgenomes/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/oceangenomesdraftgenomes/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/oceangenomesdraftgenomes/output).

## Credits

nf-core/oceangenomesdraftgenomes was originally written by Tyler Peirce.

We thank the following people for their extensive assistance in the development of this pipeline:

- Ocean Omics / Ocean Genomes project team and Pawsey collaborators

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#oceangenomesdraftgenomes` channel](https://nfcore.slack.com/channels/oceangenomesdraftgenomes) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

The Zenodo DOI will be added at the first tagged release (update the badge above accordingly).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
