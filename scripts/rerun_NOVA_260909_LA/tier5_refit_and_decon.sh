#!/usr/bin/env bash
# TIER 5 -- the remaining 41 samples: crustaceans, molluscs, echinoderms, anthozoans.
#
# These have almost no detectable contamination (FCS EXCLUDE+REVIEW under 2% for most;
# OG3037 is 0.3%) yet BUSCO C of 3-20% at N50 665-950 bp, on assemblies of 0.9-1.8 Gb
# against GenomeScope estimates of 100-130 Mb. The size estimate is the thing that was
# wrong, not the assembly: most of their k-mer mass sits above the old -m 1000 ceiling,
# and real per-haplotype coverage is roughly 10-20x.
#
# So this refits GenomeScope with -m 10000, re-runs decontamination with the new filters
# and rescores BUSCO against the per-class lineages. It does NOT re-assemble: kraken2 is
# off and the reads are the same fastp reads these assemblies were built from, so
# MEGAHIT_UNKEYED=adopt tells MEGAHIT to claim the existing checkpoints rather than
# discard them, and only meryl, GenomeScope and everything downstream re-run.
#
# adopt is an assertion that the reads have not changed. It holds here. If you change
# anything upstream of assembly for this tier, drop it back to invalidate or you will
# publish metrics describing an assembly built from different reads.
#
# After this completes the open question is resequencing depth, not another pipeline
# pass. Use the corrected genome sizes against the delivered bases to make that case.
set -euo pipefail
TIER=tier5
SKIP_FASTP=true
SKIP_ASSEMBLY=false
SKIP_KRAKEN=true
SKIP_DECON=false
SKIP_QC=false
MEGAHIT_UNKEYED=adopt
source "$(dirname "${BASH_SOURCE[0]}")/_run_tier.sh"
