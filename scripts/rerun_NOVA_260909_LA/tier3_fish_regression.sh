#!/usr/bin/env bash
# TIER 3 -- run this FIRST.
#
# The nine fishes are the only samples in NOVA_260909_LA that worked. This re-runs
# decontamination and QC on their existing assemblies with the new filters (all FCS-GX
# REVIEW contigs removed, tiara --min_len 1000) and the fixed GenomeScope gate, and
# changes nothing else: no kraken2, no re-assembly.
#
# It exists to answer one question before any compute goes into the sponges: do the new
# filters damage a well-covered vertebrate assembly? Check with
# 99_compare_to_baseline.sh tier3 -- BUSCO C should stay within ~1% and assembly size
# should not drop more than ~1%. If either moves materially, --fcs_review_action is too
# aggressive for high-agg-cvg samples and should be gated by class or agg-cvg before
# tiers 4 and 5 run.
#
# IT DOES NOT CONTROL GENOMESCOPE. SKIP_ASSEMBLY=true below means meryl, GENOMESCOPE2 and
# CALCULATE_SEQUENCING_COVERAGE never execute, so every GenomeScope column in the
# comparison is identical before and after by construction -- which is why this tier
# reported a clean bill of health while the -m 1000 -> 10000 change was inflating genome
# sizes by up to 1.8x on the other tiers. For that half of the control, run:
#
#   ./97_genomescope_control.sh tier3
#
# It refits the nine fishes from their existing histograms at the current config and runs
# the reliability gate, in seconds, without touching the pipeline.
set -euo pipefail
TIER=tier3
SKIP_FASTP=true
SKIP_ASSEMBLY=true
SKIP_KRAKEN=true
SKIP_DECON=false
SKIP_QC=false
MEGAHIT_UNKEYED=invalidate
source "$(dirname "${BASH_SOURCE[0]}")/_run_tier.sh"
