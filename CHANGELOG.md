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
- A reliability gate on the GenomeScope size estimate, in `bin/genomescope_reliability.py`
  so it can be unit tested and replayed offline over published results.
  `CALCULATE_SEQUENCING_COVERAGE` takes the GenomeScope model file and the meryl
  histogram, and flags an estimate as unreliable on inverted model-fit bounds, a fit below
  `--genomescope_min_model_fit`, implausible heterozygosity, a degenerate fit, a bias
  coefficient above `--genomescope_max_model_bias`, no genomic k-mer peak, a peak outside
  `--genomescope_min/max_peak_kmercov_ratio` of the fitted `kmercov`, or a size outside
  `--genomescope_min/max_assembly_ratio` of the decontaminated assembly. When flagged it
  reports `coverage_status: UNRELIABLE_GENOME_SIZE_ESTIMATE` with reasons in
  `genome_size_flags` instead of a confident grade. Across all 80 NOVA_260909_LA samples:
  20 pass, 60 are flagged with a specific reason, seven of the nine fishes pass.
- `RECHECK_GENOME_SIZE` (`modules/local/coverage/recheck`), which runs during QC and adds
  the one check the k-mer histogram cannot make: whether the size estimate agrees with the
  assembly built from the same reads. `CALCULATE_SEQUENCING_COVERAGE` runs before MEGAHIT,
  so its summary is provisional (`"assembly_cross_check": "pending"`) and this republishes
  it as the final verdict. It is what catches a model that converged cleanly onto
  something that is not the genome: OG3037 fitted at 86% and reported 209 Mb against a
  1.26 Gb assembly, and no fit statistic objected.
- `MEASURE_KMER_COVERAGE` (`modules/local/coverage/kmer_depth`), which measures haploid
  k-mer coverage from read depth over Complete BUSCO genes using the alignments the
  pipeline already makes. It is an estimate of the quantity GenomeScope fits, arrived at
  independently of the k-mer histogram and of the assembly's total size, so the two can
  finally be compared. On NOVA_260909_LA it agrees with GenomeScope within 8% on all eight
  samples whose fits were independently trustworthy (ratio 0.92-1.00). New flags
  `fitted_coverage_disagrees_with_read_depth`, `assembly_not_dominated_by_host` and
  `host_coverage_too_low`, plus `--kmer_depth_*` parameters.
- The module also reports host single-copy depth against assembly-wide depth. That ratio
  separates "GenomeScope fitted badly" from "the host is a minority of this assembly":
  OG3037's BUSCO genes sit at 5.8x against an assembly-wide 31.6x while 95-99% of reads
  map, so its fitted kmercov of 39 was not modelling the host at all. No fit statistic can
  see that.
- `GENOMESCOPE_RESEED` (`modules/local/genomescope_reseed`), which refits GenomeScope with
  the measured lambda as `-l` when the first fit was flagged, and keeps the result only if
  `model_fit_full` improves AND `kmercov` moves toward that measurement. Both conditions,
  because a reseed that lands on the measured coverage while the fit degrades has found a
  different optimum rather than a better one.

  An earlier version of this entry credited OG2617 with "96.7% at kmercov 8.3 ... against
  a measured 9.2", which overstated what reseeding achieves. Per the l-sweep, 96.7% is the
  `lhalf` seed (8.298); the `ldepth` seed at the measured lambda gave **77.8%**, worse than
  the original 83.5%. The rule this module implements seeds at the measured lambda, so it
  would have DECLINED OG2617's reseed even once the reseed was attempted. What the module
  buys is the attempt and the recorded reason, not a rescued fit.
- **Host genome size from read depth**, in `SUMMARISE_KMER_COVERAGE`. The k-mer route to
  a genome size collapses below about 10x host coverage: of the 32 NOVA_260909_LA samples
  with measured host lambda under 10, not one produced an estimate within 2x of its own
  assembly, and for many the assembly is largely not the animal (OG2634's is 46% symbiont
  by depth). Partitioning the assembly by per-contig depth around the host's own depth
  answers both questions, needs no k-mer peak and runs on data the pipeline already has:
  `MEASURE_KMER_COVERAGE` adds one `samtools coverage` pass and `bin/kmer_depth_lambda.py`
  does the arithmetic. New keys `host_assembly_size`, `host_assembly_fraction`,
  `symbiont_assembly_size`, `low_depth_assembly_size`, `host_contig_count` and
  `partition_status`, new `--host_depth_window` (2.0) and `--host_depth_min_contig` (500),
  a MultiQC row, and two new `draft_genomes` columns.

  The clean control reproduces GenomeScope to 4%: OG2906 gives 344.6 Mb at 91.8% host
  against a fitted 358.9 Mb. OG2634, where the fit has nothing to work with, gives
  207.0 Mb at 63.2% host with 13.6% of the assembly above host depth.

  Those OG2634 figures supersede the 168.2 Mb / 51.3% of the design note, which was
  computed with `idxstats` (reads x read length / contig length) rather than
  `samtools coverage`. That approximation credits a contig the FULL length of every read
  recorded against it, but only 55.3% of a mapped read's bases actually align on this
  sample -- the rest is soft-clipped -- so it overstates per-contig depth by 1.81x on
  average and 2.2x on the sub-700 bp contigs this 460k-contig assembly is mostly made of.
  Worse, it overstated only the CONTIG depths: the host depth defining the band came from
  `bedcov`, which counts aligned bases, so the two sides of the comparison were in
  different units and 100 Mb of host mass was pushed above the band. Reproducing the
  approximation exactly reproduces the design note's numbers (170.5 Mb / 52.0% host,
  147.9 Mb / 45.1% above), which is how this was confirmed. MAPQ filtering accounts for
  only 7 Mb of the difference, and secondary alignments and coverage gaps for none of it
  -- the clean control barely moves (345.7 -> 344.6 Mb) because a well-assembled sample
  has little soft-clipping to over-count.

  It is a **lower bound** in both directions and is described as one everywhere it
  surfaces: a symbiont sitting at the host's own depth counts as host, and a high-copy
  host repeat counts as symbiont.
- `scripts/rerun_NOVA_260909_LA/95_depth_lambda.sh` and `94_genomescope_l_sweep.sh`, the
  offline investigations those two modules came out of.
- `scripts/rerun_NOVA_260909_LA/96_gate_dryrun.sh` replays the gate over every published
  sample in seconds, so a threshold change can be evaluated against a whole run without
  re-running it. `97_genomescope_control.sh` refits GenomeScope from existing histograms
  and is the GenomeScope half of the tier-3 regression control, which
  `tier3_fish_regression.sh` cannot provide because it skips assembly.
  `98_genomescope_m_sweep.sh` sweeps `--max_kmercov` over existing histograms.
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

- **`--kmer_depth_max_ratio` lowered from 2.0 to 1.5.** 2.0 is the worst possible place for
  that boundary: fitting the homozygous peak as heterozygous or the reverse is
  GenomeScope's characteristic failure, so a threshold at exactly 2.0 sits on top of the
  most common failure mode, and four NOVA_260909_LA samples cluster just underneath it
  (OG3000 1.99, OG2617 1.95, OG2648 1.95, OG3065 1.90). The eight independently trusted
  fits agree with the measurement to within 8.3% (0.917-0.996), so 1.5 keeps a wide margin
  over everything known to be right. Blast radius: eight samples newly flag, six of them
  already `UNRELIABLE_GENOME_SIZE_ESTIMATE` for other reasons; only OG2617 and OG3065
  change verdict, and 21 reliable becomes 19.
- **The assembly cross-check compares against the host mass, not the whole assembly.**
  GenomeScope estimates the host genome, so `genome_size_disagrees_with_assembly` was
  dividing a host-only estimate by a total assembly that can be half symbiont -- the wrong
  denominator in exactly the samples the check exists for. Where the depth partition
  succeeded it now uses `host_assembly_size`, and the summary records which was used in
  `assembly_ratio_basis` (`host`/`total`) with the number itself in
  `assembly_ratio_denominator`, so a published ratio can always be reproduced.
- **A reseed is attempted on the lambda ratio as well as the verdict.**
  `genomescope_reseed_decide.py plan` gated only on `genome_size_reliable`, taken from the
  PROVISIONAL summary -- computed before the assembly exists, so blind to both the
  assembly cross-check and the lambda ratio. OG2617 recorded "no reseed attempted" while
  its fitted/measured coverage was 1.954. It now also fires when
  `kmercov_over_lambda_depth` (in either direction) exceeds `--kmer_depth_max_ratio`.
- `gfastats` moved from `modules/nf-core/` to `modules/local/` and dropped from
  `modules.json`, for the same reason megahit was: it has diverged (genome size parsed
  from a staged GenomeScope summary, `meta.assembly_prefix`, a tool_params row) and
  `nf-core modules update` would have reverted all of it silently. Its nf-core test went
  with it -- it drove an eight-input interface this module has never had.
- `fitted_kmer_coverage` renamed to `kmercov` in the coverage summary, one name for one
  number. The old spelling is written alongside it for one release and still read by
  `parse_coverage_summary`, which is the single parser: summaries published before the
  rename carry only the old key, and the `--skip_genome_assembly` tiers read exactly those.

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
- GenomeScope2's k-mer coverage ceiling is now the `--genomescope2_m` parameter, still
  defaulting to the previous hard-coded `-m 1000`. It was briefly raised to 10000 on the
  theory that the low ceiling was discarding host k-mer mass; a sweep over 27 histograms
  (`scripts/rerun_NOVA_260909_LA/98_genomescope_m_sweep.sh`) showed that it is not. The
  estimate climbs monotonically with the ceiling and never plateaus (median 1.13x from
  1000 to 10000, up to 1.81x), the whole gain lands in Genome Repeat Length while Genome
  Unique Length does not move, the same change moved repeat length +54.8% on raw reads and
  -1.6% on kraken2-cleaned reads, and one sample (OG3043) lost its fit entirely above 1000
  with `kmercov` collapsing from 40.7 to 0.6. Decontamination is the lever; the ceiling is
  not. Raise it only for a specific sample with evidence its genomic signal is truncated.
- meryl now counts the kraken2-filtered reads rather than the raw fastp reads, so the
  k-mer profile, the GenomeScope estimate and the merqury QV all describe the same
  sequence as the assembly.
- `scripts/taxonomy/load_anthozoa_taxonomy.py` renamed to `scripts/taxonomy/load_taxonomy.py`.

### `Fixed`

- **GFASTATS died on a collapsed GenomeScope fit.** OG2653's model returned `-1` for every
  property at Model Fit 0%; the module parsed that `-1` out of the summary and passed it as
  gfastats' expected-size positional, and gfastats aborted with
  `basic_string: construction from null` (exit 134), three attempts, then ignored. Any
  sample with a collapsed fit did the same. When the parsed value is not a positive
  integer the positional is now omitted entirely: gfastats reports N50 and skips the NG
  statistics, which is honest, where substituting the assembly size would have published an
  NG50 that is really an N50. The MultiQC tool_params row says why the NG columns are
  empty.
- **`kmercov` was null for 16 `draft_genomes` rows.** All 16 come from
  `--skip_genome_assembly` tiers -- tier3 (9), tier6 (4), tier2 (2), plus OG3009 -- where
  `CALCULATE_SEQUENCING_COVERAGE` never runs, so the summary is a precomputed published
  file that predates the field, while `kmer_depth_lambda.py` parsed the same model file
  successfully every time (OG2941 carried `genomescope_kmercov: 36.0` beside
  `kmercov: null`). `recheck_against_assembly` read that value into a local and never wrote
  it back; it now merges it in wherever the summary has none of its own.

- **`RECHECK_GENOME_SIZE` collected no output and silently corrupted its own input.** Its
  output glob `*_coverage_summary.json` also matched the provisional summary staged as its
  input, and Nextflow excludes input files from output matching, so every task failed with
  `MissingFileException` under `errorStrategy ignore` -- invisibly, at exit 0 with empty
  stderr. All 68 re-run samples of NOVA_260909_LA were affected: the final verdict never
  published, and `PUSH_ASSEMBLY_RESULTS` never ran at all because the channel feeding it
  stayed empty. Worse, `open("<prefix>_coverage_summary.json", "w")` resolved through the
  input symlink and overwrote `CALCULATE_SEQUENCING_COVERAGE`'s own work-directory copy,
  which is the only reason four samples published a correct summary despite failing. The
  provisional summary is now staged under `provisional/`, so the output name is free and the
  write cannot reach upstream. `scripts/rerun_NOVA_260909_LA/93_clear_stale_provisional_summaries.sh`
  clears the 64 cached summaries the old behaviour overwrote.
- **No assembly row reached the database once `RECHECK_GENOME_SIZE` was in the graph.**
  `UPLOAD_RESULTS` joined the GenomeScope summary to the coverage summary with `by: 0`,
  which compares whole meta maps. The two channels sit on opposite sides of
  decontamination, and `GENOME_DECONTAMINATION` adds `assembly_prefix` to the meta, so the
  maps were never equal and the join matched nothing -- no error, just an empty channel and
  a process that started and did nothing. Joined on `meta.id` now. It had worked only
  because the coverage summary used to come from `CALCULATE_SEQUENCING_COVERAGE`, which also
  runs before decontamination.
- **Each k-mer-depth flag was stored twice.** `recheck_against_assembly` seeds its flag list
  from the summary it is given and then appended `fitted_coverage_disagrees_with_read_depth`,
  `assembly_not_dominated_by_host` and `host_coverage_too_low` without checking whether they
  were already there. Any run whose provisional input had already been rechecked doubled
  them, which is what all 23 tier-4 samples pushed. The function is idempotent now, with a
  test that rechecks its own output and asserts nothing moves.
- **Only the right-hand column of every GenomeScope property reached the database.** For
  the length and heterozygosity rows the two columns are genuine bounds, so half of each
  interval was dropped -- OG2949 spans 331-1113 Mb and only 1113 was stored. New
  `genomesize_min`, `repeatsize_min`, `uniquesize_min`, `heterozygosity_min` and
  `homozygosity_min` columns.
- **GenomeScope's two "Model Fit" values are not a minimum and a maximum**, though they are
  printed under a `min   max` header. Its R source emits `allscore` ("Percent Kmers Modeled
  (All Kmers)") then `fullscore` ("Percent Kmers Modeled (Full Model)") -- unrelated
  statistics, which is why they can invert (OG2647 82/25, OG2624 70/15). Only `fullscore`
  was stored, so a sample could be recorded as an 84% fit while the model accounted for
  17% of its k-mers, which is OG2644 exactly; the gap exceeds 30 points in 48 of the 80
  baseline fits. The new column is `modelfit_allkmers`, named for what it is rather than as
  a `modelfit_min` it is not, and the same rename runs through
  `bin/genomescope_reliability.py` and the coverage modules.
- **A diverged fit lost its entire database row.** GenomeScope writes `Inf bp` when the
  model diverges; `_float` returned `None`, the completeness check raised, and OG3043 and
  OG2647 got no GenomeScope record at all. The finite bound is now stored with the other
  left `NULL`.
- `draft_genomes` also gains `genome_size_reliable`, `genome_size_flags`, `kmercov` and
  `lambda_depth`, so a stored genome size can be audited rather than taken on trust.
- The GenomeScope percentage columns were `NUMERIC(5,2)`, which stored 0.381807% and
  0.421044% both as ~0.4 and could not distinguish the low-heterozygosity fishes. Widened
  to `NUMERIC(8,4)`.
- **The precomputed BUSCO glob could attach the wrong lineage's scores.** A sample re-run
  against a different lineage keeps both result sets, so `*short_summary.json` matches more
  than one file for it -- 30 of 80 samples on NOVA_260909_LA. Those went straight into
  `join()`, which consumes one item per key and silently drops the rest, with arrival order
  deciding which. Selection is now by the lineage tag BUSCO writes into the filename,
  derived from the sample's class through a `buscoDbForClass` helper shared with
  `GENOME_QC` so the two cannot drift, and a missing lineage fails loudly instead of
  substituting a neighbour's scores -- BUSCO scores are not comparable across lineages.
- The reliability gate thresholded the **lower** of GenomeScope's two Model Fit values.
  That value is a residual over whatever `--max_kmercov` admitted, so it falls whenever
  the ceiling rises: the gate was measuring its own `-m` setting. It flagged 57 of 64
  re-run NOVA_260909_LA samples, including OG2906 whose real fit was 98.6498% at every
  cutoff tested. It now uses the upper value, which was identical to two decimal places at
  every `-m` for 20 of 27 histograms and moved only where a fit genuinely broke.
- A non-finite genome size bound was silently read as the opposite bound. GenomeScope
  writes `Inf bp` when the model diverges; the parser matched `([0-9,]+) bp`, so taking
  the last match returned the *lower* bound. OG3043 was published with a 2,499,788 bp
  genome size, and every coverage figure divided by it. Non-finite bounds now raise
  `non_finite_genome_size`.
- The k-mer peak detector invented a peak on histograms that do not have one, then flagged
  the disagreement it had just created -- on 22 of 23 tier-4 samples, reporting ratios of
  11x to 267x against a threshold of 3. Its trough walk re-seated the trough every time
  the count fell, so on a monotonically decreasing histogram (a coverage-starved library)
  it ran out into the sparse tail and the guard meant to return "no peak" never tripped.
  Detection now requires a genuine mode above a genuine trough, within the coverage range
  that still holds k-mer mass, and returns `no_kmer_peak` when there isn't one -- which
  for those libraries is the actionable answer.
- Peak detection searched `1..params.genomescope2_m`, so changing the model cutoff changed
  which mode was found and therefore the flag. It now uses its own
  `--genomescope_peak_search_max`.
- Degeneracy was detected as `model_fit_min == model_fit_max`, which stopped firing for
  OG2983 the moment `-m` changed, letting a 10.9 Mb "genome" at `kmercov` 569 and `bias`
  37 be graded EXCELLENT. It is now detected from the model parameters themselves.
- `scripts/rerun_NOVA_260909_LA/99_compare_to_baseline.sh` reported the lower Model Fit
  bound (the one that moves with `-m`) and silently printed the lower size bound for a
  diverged fit. It now reports the upper fit, shows `NA` for a non-finite bound, and adds
  Genome Unique Length -- without which a size estimate inflating purely in its repeat
  component looks like a size estimate improving.

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
