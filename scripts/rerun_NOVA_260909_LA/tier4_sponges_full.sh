#!/usr/bin/env bash
# TIER 4 -- the 23 microbially dominated samples. This is the expensive tier and the one
# expected to actually improve assemblies.
#
# 22 sponges (Demospongiae, Hexactinellida, Porifera) plus OG2917, the one crustacean
# with the same signature. In these libraries FCS-GX EXCLUDE+REVIEW ran to 35-240% of
# the final assembly size, tiara removed 21-527 Mb per sample, and for most of them
# FCS-GX inferred a PROKARYOTIC primary division: it judged the dominant organism in the
# assembly to be bacteria. OG2630 inferred 13 separate prokaryotic divisions.
#
# Full re-run from the fastp reads with kraken2 read-level decontamination, so the
# bacterial reads are gone BEFORE megahit and meryl rather than screened out of the
# contigs afterwards. The reads change, so the existing assemblies must not be reused;
# MEGAHIT invalidates its own checkpoints here, both because these were written before
# checkpoint fingerprinting existed and because the kraken2-filtered reads no longer
# match the fingerprint. Nothing has to be deleted by hand.
#
# Demosponges and hexactinellids are genuinely high-microbial-abundance hosts, so part
# of this is biology rather than a pipeline defect. If a sample's N50 stays around
# 1000 bp after clean reads go in, it is coverage-limited too and belongs on the
# resequencing list, not in another re-run.
#
# Run 6_kraken2_sanity_check.sh before this. Do not skip that step: it is what catches a
# confidence threshold that would delete host reads.
set -euo pipefail
TIER=tier4
SKIP_FASTP=true
SKIP_ASSEMBLY=false
SKIP_KRAKEN=false
SKIP_DECON=false
SKIP_QC=false
MEGAHIT_UNKEYED=invalidate
source "$(dirname "${BASH_SOURCE[0]}")/_run_tier.sh"
