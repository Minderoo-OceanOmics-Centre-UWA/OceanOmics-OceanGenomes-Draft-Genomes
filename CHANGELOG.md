# OceanOmics-OceanGenomes-Draft-Genomes: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

Generalised from "fish + corals" to all invertebrates.

### `Added`

- `busco_metazoa_db` registered in `nextflow_schema.json`, `conf/test.config`,
  `conf/test_full.config` and `nextflow_run.sh`. It was used by the code and passed by
  `nextflow_run_template.sh` but had never been declared anywhere else.
- `--phylum` targeting in the taxonomy loader, and support for multiple comma-separated or
  repeated targets in one invocation. NCBI leaves the class column empty for many
  invertebrate lineages, so a class-only load silently missed them; where a matched row has
  no class, the phylum name is written into `species.class` so the column is never null.

### `Changed`

- BUSCO lineage selection now defaults to metazoa instead of vertebrata. Ray-finned fish
  classes get `--busco_acti_db`, an explicit list of other vertebrate classes gets
  `--busco_vert_db`, and everything else gets `--busco_metazoa_db`. Previously any class
  other than `Actinopteri`/`Anthozoa`/`Cnidaria` was silently scored against the vertebrate
  database, which reads as a bad assembly rather than a bad lineage.
- `scripts/taxonomy/load_anthozoa_taxonomy.py` renamed to `scripts/taxonomy/load_taxonomy.py`.

### `Fixed`

- `bin/create_samplesheet.py` no longer writes `unknown` placeholders for samples whose
  taxonomy cannot be resolved. It now lists the offending samples and exits non-zero without
  writing a samplesheet, because `unknown` was passed through to FCS-GX as
  `--tax-id unknown` and failed only after MEGAHIT had already burned the SUs.
- Genome QC aborts on a sample with an `unknown` or empty class rather than guessing a
  BUSCO lineage for it.

## v1.1.0 - 2026-08-21

Coral support, per-sample reporting, and backup/cost tooling.

### `Added`

- Coral/cnidarian support in genome QC: samplesheet `class` values `Anthozoa` and `Cnidaria` now select `--busco_metazoa_db`, `Actinopteri` selects `--busco_acti_db`, and anything else falls back to `--busco_vert_db`.
- `--nt_blast_db` so invertebrate runs can be screened against an appropriate BLAST database, and MITOS reference database parameters in the auto-generated mitogenome run script.
- `SEQKIT_STATS` module for assembly statistics.
- Per-sample MultiQC report, uploaded to S3 alongside the assembly so each genome carries a record of what was run, including module parameters and tool versions.
- `bin/draft_genome_stats.py` and `bin/backfill_draft_genome_stats.py` to backfill draft-genome statistics into the SQL database from mounted Acacia/S3 archives, with a manifest template (`assets/backfill_manifest.tsv`), validate-before-apply workflow, and tests (`tests/test_draft_genome_stats_backfill.py`).
- `scripts/taxonomy/load_anthozoa_taxonomy.py` to load Anthozoa taxonomy.
- Self-contained per-genome compute cost accounting under `compute-audit/`.
- `--skip_bs_download` and `--bs_fastq_glob` to reuse locally available BaseSpace FASTQs.
- `docs/internal_sop.md` describing the internal operating procedure.

### `Changed`

- Pipeline version reported in the manifest is now derived from `git describe --tags` rather than being hardcoded.
- The samplesheet now supplies all required metadata; workflows read the meta map instead of re-deriving it.
- `nextflow_run_template.sh` copies the backup scripts into the run directory and stamps the backup config with the run ID, so the backup scripts can be run directly from there.
- The mitogenome trigger module copies the backup scripts into the mitogenome run directory.
- Contamination BLAST steps now use a MITOS-only database instead of the full `nt` database, which was too slow.
- The `db_used` filename variable takes the first four letters of the BUSCO database name, so metazoa runs are labelled correctly alongside acti/vert.
- Tiara filtering tightened to remove more contaminants; MEGAHIT initial memory raised to 80 GB.
- Backup scripts reorganised and renumbered, with extra checks around tar creation and Meryl directory deletion, and now back up to S3.
- Paths updated for the new S3 locations.
- BUSCO database schema entries corrected to `directory-path`.
- Merqury and Gfastats precomputed-result globs updated to match the current output layout.
- Documentation and README de-nf-core-ised and pointed at the local `docs/`.

### `Fixed`

- `BBMAP_FILTERBYNAME` was emitting its input, so contigs below 500 bp were not actually removed from the final assembly.
- Removed the redundant `bbmap/reformat` module, which duplicated filtering already done (correctly) in `bbmap/filterbyname`.
- MEGAHIT no longer errors when re-run over a previously completed output directory.
- BUSCO `--lineage_dataset` variable fixed, and the checkpoints directory removed now that BUSCO v6 no longer needs it.
- `UPLOAD_RESULTS` works when all other steps are skipped.
- Check files are removed automatically instead of needing manual deletion.
- Coverage modules clean up intermediates that previously created thousands of files.
- Run IDs handled correctly when multiple runs are given in any order.
- Gfastats `--nstar-report` added; `--output-format` dropped so the genome is not regenerated when only stats are needed.

## v1.0.0 - 2025-12-11

Initial release.
