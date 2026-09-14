# OceanOmics-OceanGenomes-Draft-Genomes: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

Generalised from "fish + corals" to all invertebrates.

### `Added`

- MEGAHIT checkpoints are fingerprinted on their inputs. The checkpoint directory under
  `${params.outdir}/megahit_checkpoints` exists so a killed assembly can resume, which
  puts it outside Nextflow's caching; it was keyed on the sample prefix alone, so a
  checkpoint was reused even when the reads that produced it had changed. With kraken2
  now sitting upstream of assembly that was a live hazard: all 80 NOVA_260909_LA
  checkpoints are present and marked `done`, so a kraken2 re-run would have republished
  the old contaminated assembly while every downstream metric claimed to describe a clean
  one. `bin/megahit_checkpoint_key.sh` hashes the reads (basename and byte size), the
  assembly arguments and the megahit version; `-m`/`-t` are excluded so a retry at higher
  memory still resumes, and path and mtime are excluded so a fresh work directory or
  restaged reads do not force a needless reassembly. New `--megahit_checkpoint_unkeyed`
  (`invalidate`/`adopt`, for checkpoints predating the fingerprint) and
  `--megahit_stale_checkpoint` (`rerun`/`archive`/`fail`).
- `MEGAHIT` publishes `*.megahit_checkpoint.txt` and records the decision in its MultiQC
  tool_params row. The skip path was otherwise invisible: the trace showed MEGAHIT
  running for seconds and emitting contigs, with nothing to say whether they were built.
- `--megahit_checkpoint_cleanup` (default true) prunes the intermediate assembly graph
  from a checkpoint once it completes. A completed checkpoint can never be continued, so
  its intermediates were tens of GB per sample retained indefinitely under the outdir.

- kraken2 read-level decontamination before assembly, as `KRAKEN2_KRAKEN2` +
  `KRAKENTOOLS_EXTRACTREADS` in `GENOME_ASSEMBLY`. Read pairs assigned to Bacteria,
  Archaea or Viruses (`--kraken2_exclude_taxids`, descendants included) are dropped
  before meryl and megahit; unclassified reads are kept, which is where an
  under-referenced invertebrate host's genome mostly sits. Screening contigs after
  assembly cannot undo a host/symbiont co-assembly: in NOVA_260909_LA, 25 sponges had
  FCS-GX infer a *prokaryotic* primary division, and their assemblies ran to 150k-660k
  contigs at N50 ~1000. New `--kraken2_db`, `--kraken2_confidence` (default 0.05, since
  kraken2's own default of 0 assigns on a single k-mer hit and would delete host reads),
  and `--skip_kraken2_decontamination`. Per-sample retained-read fractions go to MultiQC.
- Per-class BUSCO lineages: `--busco_crustacea_db`, `--busco_mollusca_db`,
  `--busco_anthozoa_db`, `--busco_arthropoda_db`, selected off `meta.class` ahead of the
  metazoa fallback. Note there is no echinoderm or sponge lineage in odb12, so
  Asteroidea, Ophiuroidea, Demospongiae, Hexactinellida and Porifera stay on metazoa by
  necessity. Unset lineages fall back to metazoa, so this is backwards compatible.
- A reliability gate on the GenomeScope size estimate. `CALCULATE_SEQUENCING_COVERAGE`
  now takes the GenomeScope model file and the meryl histogram, and flags an estimate as
  unreliable on inverted model-fit bounds, a fit below `--genomescope_min_model_fit`,
  implausible heterozygosity, a degenerate single-solution fit, no genomic k-mer peak, or
  a peak that disagrees with the fitted `kmercov` by more than
  `--genomescope_max_peak_kmercov_ratio`. When flagged it reports
  `coverage_status: UNRELIABLE_GENOME_SIZE_ESTIMATE` with reasons in `genome_size_flags`
  instead of a confident grade. Validated against all 80 samples of NOVA_260909_LA: 18
  pass, 62 are flagged with a specific reason, and the five best fishes all pass.
- `scripts/rerun_NOVA_260909_LA/`: tiered re-run launchers, per-tier samplesheets, a
  baseline snapshot script and a before/after comparison for remediating that run.

- NCBI taxdump fallback for samplesheet taxonomy. `bin/create_samplesheet.py` now resolves a
  sample's `taxon_id` and `class` from the NCBI taxdump when the curated `species` table has
  no row for its `nominal_species_id`, via a new `bin/taxdump_lineage.py` (kept
  byte-identical to the mitogenome pipeline's copy). Enabled with the new
  `--taxonkit_db_dir` parameter, whose `storeDir` cache is shared with the mitogenome
  pipeline; without it the step behaves exactly as before. The curated table stays
  authoritative and is only topped up where it is blank.
- `taxonomy_resolution.tsv`, published alongside the samplesheet, recording per-sample
  taxonomy provenance (`db`, `db+taxdump`, `taxdump`, `unresolved`). Without it a
  taxdump-derived class is indistinguishable from a curated one. Written on the failure path
  too, since that is when it is most needed.
- `busco_metazoa_db` registered in `nextflow_schema.json`, `conf/test.config`,
  `conf/test_full.config` and `nextflow_run.sh`. It was used by the code and passed by
  `nextflow_run_template.sh` but had never been declared anywhere else.
- `--phylum` targeting in the taxonomy loader, and support for multiple comma-separated or
  repeated targets in one invocation. NCBI leaves the class column empty for many
  invertebrate lineages, so a class-only load silently missed them; where a matched row has
  no class, the phylum name is written into `species.class` so the column is never null.

### `Changed`

- `megahit` moved from `modules/nf-core/` to `modules/local/` and dropped from
  `modules.json`. It had diverged far enough (renamed output, tool_params row, the whole
  checkpoint mechanism) that `nf-core modules update megahit` would have silently
  reverted all of it, and nothing recorded it as patched. Its nf-core test went with it:
  it asserted on `k_contigs`, `addi_contigs`, `local_contigs`, `kfinal_contigs` and
  `log`, none of which this module has emitted for some time.

- BUSCO lineage selection now defaults to metazoa instead of vertebrata. Ray-finned fish
  classes get `--busco_acti_db`, an explicit list of other vertebrate classes gets
  `--busco_vert_db`, and everything else gets `--busco_metazoa_db`. Previously any class
  other than `Actinopteri`/`Anthozoa`/`Cnidaria` was silently scored against the vertebrate
  database, which reads as a bad assembly rather than a bad lineage.
- FCS-GX REVIEW contigs are now removed by default. `BBMAP_FILTERBYNAME` previously
  removed only REVIEW contigs of 1000 bp or less, and tiara (the intended safety net) ran
  at its default `--min_len 3000`, so REVIEW contigs between 1 and 3 kb escaped both
  filters: 818 Mb across NOVA_260909_LA, up to 12% of individual assemblies. Those
  contigs are not ambiguous host sequence -- OG2630's REVIEW set is taxonomically
  indistinguishable from its EXCLUDE set and near-entirely prokaryotic. FCS-GX downgrades
  EXCLUDE to REVIEW when its confidence is low, which for taxa the GX database covers
  poorly (`agg-cvg` 0.08-0.19 here, against 0.33-0.39 for the fishes) is most of what it
  finds. `--fcs_review_action` takes `exclude` (default), `short-only` (the old rule) or
  `keep`. Output renamed `*.review_scaffolds_1kb.txt` -> `*.review_scaffolds.txt`.
- tiara now runs at `--min_len 1000` (`--tiara_min_len`) instead of its 3000 default. On
  these megahit assemblies the old floor screened as little as 1.3% of an assembly's
  bases (OG3037), versus 84% for a well-covered fish.
- GenomeScope2's k-mer coverage ceiling is now `--genomescope2_m` (default 10000), up
  from a hard-coded `-m 1000`. 48 of 78 NOVA_260909_LA samples had more than 30% of their
  total k-mer mass above 1000x and therefore excluded from the fit, which broke the model
  and pushed Genome Haploid Length down by up to 10x.
- meryl now counts the kraken2-filtered reads rather than the raw fastp reads, so the
  k-mer profile, the GenomeScope estimate and the merqury QV all describe the same
  sequence as the assembly.
- `scripts/taxonomy/load_anthozoa_taxonomy.py` renamed to `scripts/taxonomy/load_taxonomy.py`.

### `Fixed`

- `MEGAHIT` treated a checkpoint as complete on a `done` marker *or* a non-empty contigs
  file. A checkpoint with `done` but a missing or truncated contigs file surfaced several
  lines later as an unexplained `cp` failure; it now requires both and otherwise treats
  the checkpoint as partial.

- GenomeScope2 was hard-coded to `-k 21` while meryl used `--kvalue`. Changing `--kvalue`
  would have had GenomeScope silently model a different k than the histogram was built
  with. It now reads `--kvalue` too.
- `contig_count_500bp.txt` counted wrapped FASTA sequence *lines* shorter than 500
  characters, not contigs. OG2630 reported 4,775,697 for an assembly of 167,330 contigs
  of which 2 are under 500 bp. The figure was published to MultiQC and pushed to the
  database.
- GenomeScope heterozygosity and model fit were parsed with `([0-9.]+)%`, which drops a
  minus sign. OG2653's `-100%` heterozygosity, a total model failure, was read as `100%`.
- `modules.json` was not valid JSON (two trailing commas), which broke
  `nf-core modules` commands against this repo.
- An unresolvable sample no longer costs you the whole samplesheet. `bin/create_samplesheet.py`
  used to exit non-zero without writing anything, and because `publishDir` does not run on a
  failed task, `<outdir>/samplesheet/` stayed empty: there was nothing to inspect and nothing
  to correct. The sheet is now always written, with the offending rows marked `unknown`, and
  the run is stopped instead by a new `validateTaxonomy()` in
  `subworkflows/local/prepare_samplesheet` after the sheet is published and before any
  assembly task is submitted. The guarantee that `--tax-id unknown` never reaches FCS-GX is
  unchanged; it is now enforced in the workflow rather than by refusing to write the file.
- The samplesheet-reuse path never matched. `main.nf` passed `"<run>_samplesheet"` as the
  prefix and the subworkflow appended `_samplesheet.csv`, so it looked for
  `<run>_samplesheet_samplesheet.csv` while the module publishes `<run>_samplesheet.csv`. A
  re-run therefore regenerated from the database and republished over a hand-corrected sheet,
  losing the edits.
- The samplesheet is read by column name instead of by position. `prepare_samplesheet` rebuilt
  its meta map from `sample_record[4]`/`[5]`/`[6]`, so reordering a column in
  `schema_input.json` would have silently reassigned every taxonomy field. All the meta columns
  are now annotated with `meta` in `assets/schema_input.json`, as the mitogenome pipeline does.
  Safe to do here because every column is `required`, so nf-schema can never substitute the
  empty-list placeholder that makes this change hazardous when columns are optional; the
  resulting meta map is key-for-key and type-for-type identical to the hand-built one, leaving
  every `join(..., by: 0)` and task hash untouched.
- The `taxon_id` guard before FCS-GX rejected only null, not the string `unknown`. A non-empty
  string is truthy, so `if (!meta.taxon_id)` waved a placeholder straight through.
- Invertebrate runs no longer abort on uncurated taxa. `NOVA_260909_LA` lost 67 of its
  samples to `missing taxon_id, class` in one go; with the taxdump fallback the same run
  resolves 65 of them and stops on the 2 whose `nominal_species_id` is not a real taxon
  (`Actinaria`, a misspelling of `Actiniaria`, and `Larval fish`). Bulk-loading every
  invertebrate the runs touch was never going to scale, since they draw from most of
  Metazoa.
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
