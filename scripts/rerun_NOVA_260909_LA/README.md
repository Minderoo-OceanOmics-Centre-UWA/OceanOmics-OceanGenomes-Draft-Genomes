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
2. **GenomeScope's `-m 1000` truncated the histogram.** 48 of 78 samples had more than
   30% of their k-mer mass above that cutoff and thrown away, which broke the fit and
   pushed the size estimate down by up to 10x. The bad size then propagated into
   `theoretical_coverage`, and OG2617 was published as `"coverage_status": "EXCELLENT"` on
   a 48% model fit.
3. **The sponges are microbially dominated.** For 25 of them FCS-GX inferred a
   *prokaryotic* primary division, i.e. it judged bacteria to be the dominant organism in
   the assembly. Screening contigs after assembly cannot undo a host/symbiont
   co-assembly, which is why those assemblies have 150k-660k contigs at N50 ~1000.

## Order

```
./00_snapshot_baseline.sh          # required first: the tiers publish over the originals
./tier3_fish_regression.sh        # the control. STOP HERE if it regresses.
./99_compare_to_baseline.sh tier3
./tier2_resume.sh                 # the two that never finished
./6_kraken2_sanity_check.sh       # required before tier 4
./tier4_sponges_full.sh           # the expensive tier, where the gains are
./99_compare_to_baseline.sh tier4
./tier5_refit_and_decon.sh        # cheap: refit sizes, re-filter, rescore BUSCO
./99_compare_to_baseline.sh tier5
```

Tier 1 is not a script: see `tier1_no_rerun.txt`.

## Tiers

| tier | n | what it does | cost |
|---|---|---|---|
| 1 | 4 | **Nothing.** Failed libraries, request resequencing. See `tier1_no_rerun.txt`. | none |
| 1b | 1 | **Nothing yet.** OG3009 needs a sample-identity check first. | none |
| 2 | 2 | Resume OG2625 and OG3036, which died mid-decontamination. | low |
| 3 | 9 | Fishes, decontamination + QC only. The regression control. | low |
| 4 | 23 | Sponges + OG2917: full re-run from fastp reads **with kraken2**. | high |
| 5 | 41 | Refit GenomeScope, re-filter, rescore BUSCO. No re-assembly. | medium |

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
- **tier 5**: GenomeScope sizes rise substantially (OG3037's 123 Mb estimate against a
  1.27 Gb assembly is the clearest test), and samples whose fit is still bad now report
  `coverage_status: UNRELIABLE_GENOME_SIZE_ESTIMATE` with a reason in
  `genome_size_flags` instead of a confident grade.

Validated against all 80 samples of the original run, the reliability gate passes 18 and
flags 62, including every sample identified as broken, with a specific reason each:
`inverted_model_fit_bounds` (OG2624, OG2647), `degenerate_model_fit` (OG2639, OG2959,
OG2983), `implausible_heterozygosity` (OG2653), `model_fit_below_50pc`, and
`kmer_peak_disagrees_with_fitted_coverage`. All five of the best fishes pass.

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
