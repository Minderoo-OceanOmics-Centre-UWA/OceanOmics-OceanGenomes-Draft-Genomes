# NOVA_260909_LA remediation

This directory holds the re-run plan for NOVA_260909_LA, the first mostly-invertebrate
run through this pipeline. You run these; nothing here submits anything on its own.

## Why

The run completed but the outputs are largely unusable: median BUSCO completeness (metazoa_odb12)
around 15%, and only 5 of 80 samples clear C:45% - all 5 of those are the fishes. Three
separate causes, all now fixed in the pipeline:

1. **FCS-GX REVIEW contigs were retained.** The filter only removed REVIEW contigs of
   1000 bp or less, and tiara (the intended safety net) ran at its default
   `--min_len 3000`. REVIEW contigs between 1 and 3 kb escaped both: **818 Mb across the
   run**, up to 12% of individual assemblies. They were not ambiguous - OG2630's REVIEW
   set is taxonomically indistinguishable from its EXCLUDE set, near-entirely prokaryotic.
2. **The genome size estimates were not trustworthy, and the gate on them was worse.**
   OG2617 was published as `"coverage_status": "EXCELLENT"` on what the pipeline reported
   as a 48% model fit.

   The first attempt at this raised GenomeScope's `-m` from 1000 to 10000, on the theory
   that the low ceiling was discarding host k-mer mass. **That was wrong and has been
   reverted** -- see "The -m 10000 detour" below. The real defects were in the reliability
   gate, and they are fixed: it now thresholds the fit GenomeScope actually reports rather
   than a residual that moves with `-m`, detects a k-mer peak instead of inventing one,
   and cross-checks the size against the assembly.
3. **The sponges are microbially dominated.** For 25 of them FCS-GX inferred a
   *prokaryotic* primary division, i.e. it judged bacteria to be the dominant organism in
   the assembly. Screening contigs after assembly cannot undo a host/symbiont
   co-assembly, which is why those assemblies have 150k-660k contigs at N50 ~1000.

## Order

```
./00_snapshot_baseline.sh          # required first: the tiers publish over the originals
./tier3_fish_regression.sh        # the control. STOP HERE if it regresses.
./99_compare_to_baseline.sh tier3
./97_genomescope_control.sh tier3  # the other half of the control: tier 3 skips GenomeScope
./tier2_resume.sh                 # the two that never finished
./6_kraken2_sanity_check.sh       # required before tier 4
./tier4_sponges_full.sh           # the expensive tier, where the gains are
./99_compare_to_baseline.sh tier4
./tier5_refit_and_decon.sh        # cheap: refit sizes, re-filter, rescore BUSCO
./99_compare_to_baseline.sh tier5
./tier6_failed_libraries_qc.sh    # the four failed libraries: verdicts only, no re-filter
```

Tier 1 is not a script: see `tier1_no_rerun.txt`.

## Tiers

| tier | n | what it does | cost |
|---|---|---|---|
| 1 | 4 | **No re-assembly and no re-filtering.** Failed libraries, request resequencing. See `tier1_no_rerun.txt`. Tier 6 now runs QC over them so the resequencing case is recorded as numbers rather than prose. | none |
| 1b | 1 | **Nothing yet.** OG3009 needs a sample-identity check first. | none |
| 2 | 2 | Resume OG2625 and OG3036, which died mid-decontamination. | low |
| 3 | 9 | Fishes, decontamination + QC only. The regression control. | low |
| 4 | 23 | Sponges + OG2917: full re-run from fastp reads **with kraken2**. | high |
| 5 | 41 | Refit GenomeScope, re-filter, rescore BUSCO. No re-assembly. | medium |
| 6 | 4 | The failed libraries from tier 1. QC and gate only, no re-filter, no re-assembly: gives them the reliability verdict and measured coverage every other sample now has. OG3009 stays out until its identity is settled. | low |

## Running tiers concurrently

Every tier launches from its own directory, `${OUT}/.nf_<tier>`, seeded once from the
run's original `.nextflow`. That gives each tier its own history, cache database and
session, which is what makes running them all at once safe -- Nextflow opens the cache for
writing, so two tiers sharing one was a hazard independent of anything else.

It also fixes resume. `_run_tier.sh` used to pass a bare `-resume`, which resolves to the
last run in the launch directory's history; with every tier launching from `$OUT` that was
usually a *different tier*. Tier 3 relaunched against tier 4's index and recomputed all
nine samples from scratch, while tier 5 landed on a cache that happened to fit and reused
860 tasks. Each tier now resumes its own last successful run by name, and prints which one
in the banner before starting:

```
 resuming:     NOVA_260909_LA_tier3_20260914_142712
```

Override with `RESUME_RUN=<run name>` to resume something specific, or `RESUME_RUN=none`
to force a clean run.

Three outputs are run-level rather than per-sample and would otherwise be written by every
tier at once: `multiqc/`, `coverage_summary/` and `pipeline_info/software_versions.yml`.
Each tier writes them under its own `--report_subdir <tier>`, so nothing races. They each
describe only that tier's samples, so regenerate a single run-wide set over all 80 once the
tiers are done.

## Two things that will bite you

**MEGAHIT checkpoint reuse.** MEGAHIT resumes from
`$OUT/megahit_checkpoints/<prefix>_megahit_out` and skips entirely when it finds a
finished assembly there. All 80 checkpoints from the original run are present and marked
`done`. That is what makes tier 5 cheap, and it used to be a trap: any tier that changed
the *reads* would get the old contaminated assembly republished, with every downstream
metric describing it.

Checkpoints are now fingerprinted on their inputs (`bin/megahit_checkpoint_key.sh`), so a
changed read set invalidates one automatically and tier 4 needs nothing cleared by hand.
The one thing the fingerprint cannot decide is a checkpoint written *before*
fingerprinting existed, which is all 80 of the ones on disk: their inputs are unknowable.
Each tier therefore declares what it wants:

| tier | `MEGAHIT_UNKEYED` | why |
|---|---|---|
| 2, 3 | `invalidate` | assembly is skipped outright, so it never comes up |
| 4 | `invalidate` | reads change (kraken2), the old assemblies must not be reused |
| 5 | `adopt` | same fastp reads as the originals, so reuse is correct and saves 41 reassemblies |

`adopt` is an assertion that the reads have not changed. If you change anything upstream
of assembly in tier 5, put it back to `invalidate`. After that first pass every checkpoint
carries a fingerprint and this choice stops mattering.

Each MEGAHIT task now publishes `<prefix>.megahit_checkpoint.txt` recording which path it
took (`fresh`, `continue`, `skip`, `invalidated-*`) and the key it used, so "was this
assembly actually rebuilt" is a file to read rather than an inference from the trace.

**The tiers publish into the original outdir.** They have to: the `precomputed_*` globs
resolve against `$OUT`. So the numbers this review is based on are overwritten as soon as
tier 3 starts. `00_snapshot_baseline.sh` copies the small summary files to `$OUT/rerun_baseline` first (a few MB,
no FASTA/FASTQ) and `_run_tier.sh` refuses to run without that snapshot.

## What "success" looks like

- **tier 3**: BUSCO C within ~1% of baseline (OG2941 73.8%, OG3078 72.3%) and assembly
  size down no more than ~1%. A material drop means `--fcs_review_action exclude` is too
  aggressive for high-`agg-cvg` samples and should be gated by class or `agg-cvg`.
- **tier 4**: contig count falls several-fold and N50 rises. If N50 stays near 1000 bp
  after clean reads go in, that sample is coverage-limited as well and belongs on the
  resequencing list rather than in another re-run.
- **tier 5**: samples whose fit is bad report `coverage_status:
  UNRELIABLE_GENOME_SIZE_ESTIMATE` with a specific reason in `genome_size_flags` instead
  of a confident grade. Do **not** expect the size estimates themselves to move much:
  nothing in this remediation improves a GenomeScope fit, it only stops the pipeline
  presenting a bad one as fact.

Validated against all 80 samples with `./96_gate_dryrun.sh`, the gate passes 20 and flags
60, each with a specific reason: `no_kmer_peak` (57 -- these libraries have no genomic
mode at all, which means resequence), `genome_size_disagrees_with_assembly` (27),
`implausible_model_bias` (10), `degenerate_model_fit` (9), `model_fit_below_50pc` (5),
`non_finite_genome_size` (4), `inverted_model_fit_bounds` (4). Seven of the nine fishes
pass; the two that do not are OG2647 (a diverged fit reporting 36 Mb against a 904 Mb
assembly) and OG3074 (bias 19.2 at kmercov 71.2, no k-mer peak).

## The -m 10000 detour

Worth reading before anyone proposes it again. `98_genomescope_m_sweep.sh` refit 27
histograms at 1000, 3000, 10000 and adaptive multiples of each sample's fitted coverage:

- **The size estimate never plateaus.** It climbs monotonically with the ceiling (median
  1.13x from 1000 to 10000, up to 1.81x for OG2940). There is no coverage above which the
  extra mass stops arriving, so there is no principled place to stop. That is a continuous
  contamination and error tail, not repeat structure.
- **All of it lands in repeat length.** Genome Unique Length, the component the model
  constrains, does not move: across tier 5 it went -2.4% while repeat length went +54.8%.
- **Cleaning the reads is the lever, not the ceiling.** The identical `-m` change moved
  repeat length +54.8% on raw reads and -1.6% on kraken2-cleaned reads.
- **It broke one sample outright.** OG3043's kmercov collapsed from 40.7 to 0.6 and its
  model fit from 88% to 34% above `-m 1000`.
- **It made the gate meaningless.** The gate thresholded the *lower* of GenomeScope's two
  Model Fit values, which is a residual over whatever `-m` admitted and so falls whenever
  the ceiling rises. 57 of 64 re-run samples were flagged, including OG2906 at a real fit
  of 98.6%. The upper value -- what the gate uses now -- was identical to two decimal
  places at every `-m` for 20 of 27 histograms, and moved only where a fit genuinely broke.

`genomescope2_m` is back to 1000 and remains a parameter. Raise it only for a specific
sample with evidence its genomic signal is truncated, not as a default.

## Open decisions for you

- **BUSCO generation mismatch.** The four new lineages
  (`crustacea/mollusca/anthozoa/arthropoda_odb12`) are OrthoDB **12.1** (built 2026-05-22).
  The installed `metazoa_odb12` is OrthoDB **12.0** (2025-07-01). Sponges and echinoderms
  have no odb12 lineage of their own so they stay on metazoa, which means this run will
  mix generations. Scores are not comparable across lineages anyway, but if you want them
  internally consistent, download `metazoa_odb12` at 2026-05-22 and repoint
  `--busco_metazoa_db`. That changes the sponge/echinoderm baseline, so it is a deliberate
  choice, not a cleanup.
- **BUSCO lineage location.** `/software/projects/pawsey0964/busco_db/` is owned by another
  user and not writable, so the new lineages went to
  `/software/projects/pawsey1348/$USER/busco_db/busco_downloads/lineages/`. Move them
  somewhere shared if other people need them.
- **Resequencing.** Tier 1's four libraries, and whichever tier 4/5 samples stay
  coverage-limited.
