#!/usr/bin/env bash
# TIER 2 -- the two samples that never finished.
#
# OG2625 (Ophiuroidea) has a GenomeScope summary but no FCS-GX report and no assembly
# QC: decontamination never ran. OG3036 (Polyplacophora) has an FCS-GX report
# (135,603 EXCLUDE / 179,044 REVIEW) but no tiara, seqkit, BUSCO, merqury or .fna: it
# died mid-decontamination. Neither failed loudly enough to be noticed.
#
# This resumes them from their existing assemblies with the new filters applied. Once
# they complete they belong in tier 5 (both are Mode B, coverage-limited), so re-run
# them there if you want their GenomeScope estimates refitted as well.
set -euo pipefail
TIER=tier2
SKIP_FASTP=true
SKIP_ASSEMBLY=true
SKIP_KRAKEN=true
SKIP_DECON=false
SKIP_QC=false
MEGAHIT_UNKEYED=invalidate
source "$(dirname "${BASH_SOURCE[0]}")/_run_tier.sh"
