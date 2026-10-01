#!/usr/bin/env bash
# Shared launcher for the NOVA_260909_LA remediation tiers.
#
# Not run directly. Each tierN_*.sh sets the variables below and sources this file.
#
# Required from the caller:
#   TIER                     tier label, e.g. "tier3"
#   SKIP_FASTP               true/false  --skip_fastp_fastqc
#   SKIP_ASSEMBLY            true/false  --skip_genome_assembly
#   SKIP_KRAKEN              true/false  --skip_kraken2_decontamination
#   SKIP_DECON               true/false  --skip_genome_decontamination
#   SKIP_QC                  true/false  --skip_genome_qc
#   MEGAHIT_UNKEYED          invalidate/adopt  what to do with the pre-existing MEGAHIT
#                            checkpoints, which carry no input fingerprint (see below)
set -euo pipefail

module load nextflow/25.04.6
module load singularity/4.1.0-nompi

RUN=NOVA_260909_LA
BASE="/scratch/pawsey1348/$USER"
OUT="${BASE}/${RUN}"
MITO_PIPELINE_DIR="/software/projects/pawsey1348/$USER/repos/Oceanomics-OceanGenomes-Mitogenomes-inverts"
BUSCO_LINEAGES="/software/projects/pawsey1348/$USER/busco_db/busco_downloads/lineages"

# Repo root, two levels up from scripts/rerun_NOVA_260909_LA/
RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
TIER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
SAMPLESHEET="${TIER_DIR}/samplesheets/${TIER}_samplesheet.csv"

if [ ! -f "$SAMPLESHEET" ]; then
    echo "ERROR: no samplesheet for ${TIER} at ${SAMPLESHEET}" >&2
    exit 1
fi
if [ ! -d "$OUT" ]; then
    echo "ERROR: run directory $OUT does not exist." >&2
    exit 1
fi
if [ ! -d "${OUT}/rerun_baseline" ]; then
    echo "ERROR: no baseline snapshot found at ${OUT}/rerun_baseline." >&2
    echo "ERROR: run 00_snapshot_baseline.sh first -- these tiers publish over the" >&2
    echo "ERROR: existing results, and without a baseline there is nothing to compare to." >&2
    exit 1
fi

N=$(($(wc -l < "$SAMPLESHEET") - 1))
echo "=============================================================="
echo " ${RUN} remediation: ${TIER}"
echo " samples:      ${N}  (${SAMPLESHEET})"
echo " outdir:       ${OUT}"
echo " skip fastp:   ${SKIP_FASTP}      skip assembly: ${SKIP_ASSEMBLY}"
echo " skip kraken2: ${SKIP_KRAKEN}      skip decon:    ${SKIP_DECON}    skip qc: ${SKIP_QC}"

# MEGAHIT keeps a resumable checkpoint per sample under ${OUT}/megahit_checkpoints, and
# skips assembly outright when it finds a completed one. That is what makes tier 5 cheap
# and what used to be a trap for any tier that changes the READS going in: with kraken2
# enabled and a stale checkpoint present, MEGAHIT would silently republish the old
# contaminated assembly and every downstream metric would describe it.
#
# The module now fingerprints each checkpoint against its inputs, so a changed read set
# invalidates the checkpoint by itself and tier 4 needs no manual clearing. The one case
# it cannot decide is a checkpoint written BEFORE fingerprinting existed, which is all 80
# of the ones on disk now: their inputs are unknowable, so the default is to discard and
# reassemble. Tiers that deliberately reuse those assemblies (tier 5) set
# MEGAHIT_UNKEYED=adopt to assert that the reads have not changed.

# Each tier gets its own launch directory, so each gets its own .nextflow: its own
# history, its own cache database, its own session id, its own execution report.
#
# This is what makes running several tiers at once safe. They used to all launch from
# $OUT, sharing one history and one cache -- and a bare `-resume` resolves to the LAST run
# in that shared history, not to the tier's own. Tier 3 relaunched against tier 4's index
# and recomputed all 9 samples from scratch while tier 5, launched a minute later, landed
# on a cache that happened to fit and reused 860 tasks. Nextflow also opens the cache for
# writing, so two concurrent runs sharing one LevelDB is a hazard in its own right.
#
# Seed the new directory from the original shared .nextflow the first time, or every tier
# starts with an empty cache and re-runs everything -- the opposite of the point. It is
# about 14 MB.
# Same path tier5_standalone.sh established, so tier 5's existing session lineage and its
# cached tasks carry over rather than starting again from nothing.
LAUNCH="$OUT/.nf_${TIER}"
mkdir -p "$LAUNCH"
if [ ! -d "$LAUNCH/.nextflow" ] && [ -d "$OUT/.nextflow" ]; then
    echo "seeding ${TIER} launch dir with the existing Nextflow cache"
    cp -a "$OUT/.nextflow" "$LAUNCH/.nextflow"
fi
cd "$LAUNCH"

# Resume this tier's own last SUCCESSFUL run, by name.
#
# A bare `-resume` picks whatever ran last in this directory's history, which after the
# copy above still includes every other tier -- that is the bug being fixed. Naming the run
# is the fix; picking the successful one is what makes the name the right one.
#
# Status matters because a failed run's index holds only what that run executed. Tier 3's
# most recent run is an ERR with a 391-byte index, because it resumed tier 4 by mistake and
# recomputed everything; the run worth resuming is the OK one before it, at 3247 bytes.
# Falling back to the most recent of any status covers a run that was interrupted rather
# than failed, where the partial cache is still worth having.
#
# Override with RESUME_RUN=<run name> to resume something specific, or RESUME_RUN=none to
# force a clean run.
RESUME_ARGS=()
HISTORY="$LAUNCH/.nextflow/history"
if [ "${RESUME_RUN:-}" = "none" ]; then
    RESUME_TARGET=""
    echo " resuming:     nothing (RESUME_RUN=none)"
else
    RESUME_TARGET="${RESUME_RUN:-}"
    if [ -z "$RESUME_TARGET" ] && [ -f "$HISTORY" ]; then
        RESUME_TARGET="$(awk -F'\t' -v prefix="^${RUN}_${TIER}_" \
            '$3 ~ prefix && $4 == "OK" { last = $3 } END { if (last != "") print last }' \
            "$HISTORY")"
        [ -z "$RESUME_TARGET" ] && RESUME_TARGET="$(awk -F'\t' -v prefix="^${RUN}_${TIER}_" \
            '$3 ~ prefix { last = $3 } END { if (last != "") print last }' "$HISTORY")"
    fi
fi

if [ -n "$RESUME_TARGET" ]; then
    if compgen -G "$LAUNCH/.nextflow/cache/*/index.$RESUME_TARGET" > /dev/null; then
        RESUME_ARGS=(-resume "$RESUME_TARGET")
        echo " resuming:     ${RESUME_TARGET}"
    else
        echo " resuming:     nothing (no cache index for ${RESUME_TARGET})" >&2
    fi
elif [ "${RESUME_RUN:-}" != "none" ]; then
    echo " resuming:     nothing (no previous ${TIER} run found)"
fi
echo "=============================================================="

nextflow -log ".nextflow_${RUN}_${TIER}.log" \
    run "$RUN_DIR/main.nf" \
    "${RESUME_ARGS[@]}" \
    -name "${RUN}_${TIER}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir "$OUT/work/${RUN}_${TIER}" \
    -c "$RUN_DIR/pawsey_profile.config" \
    -profile singularity \
    -with-report \
    --run "$RUN" \
    --outdir "$OUT" \
    --mitogenome_nfcore_dir "$MITO_PIPELINE_DIR" \
    --kvalue "21" \
    --genomescope2_l false \
    --genomescope2_m 1000 \
    --tiara_min_len 1000 \
    --fcs_review_action exclude \
    --megahit_checkpoint_unkeyed "$MEGAHIT_UNKEYED" \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --kraken2_db "/scratch/references/kraken2/pluspfp_20230605" \
    --kraken2_confidence 0.05 \
    --kraken2_exclude_taxids "2 2157 10239" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --busco_metazoa_db "/software/projects/pawsey0964/busco_db/metazoa_odb12" \
    --busco_crustacea_db "${BUSCO_LINEAGES}/crustacea_odb12" \
    --busco_mollusca_db "${BUSCO_LINEAGES}/mollusca_odb12" \
    --busco_anthozoa_db "${BUSCO_LINEAGES}/anthozoa_odb12" \
    --busco_arthropoda_db "${BUSCO_LINEAGES}/arthropoda_odb12" \
    --tempdir "$BASE/tmp" \
    --refresh-modules \
    --skip_bs_download true \
    --skip_download_reads true \
    --skip_fastp_fastqc "$SKIP_FASTP" \
    --skip_kraken2_decontamination "$SKIP_KRAKEN" \
    --skip_genome_assembly "$SKIP_ASSEMBLY" \
    --skip_genome_decontamination "$SKIP_DECON" \
    --skip_genome_qc "$SKIP_QC" \
    --skip_upload_results false \
    --report_subdir "$TIER" \
    --taxonkit_db_dir "$BASE" \
    --input "$SAMPLESHEET"
NF_STATUS=$?

# Nextflow exits 0 even when a task was never dispatched because its work directory
# could not be created (/scratch quota exhausted), leaving a stale report published as
# this run's output. Same guard as the main launcher.
NF_LOG="$LAUNCH/.nextflow_${RUN}_${TIER}.log"
if [ -f "$NF_LOG" ] && grep -qE 'Disk quota exceeded|Unable to create directory=' "$NF_LOG"; then
    echo "ERROR: filesystem errors in $NF_LOG (disk quota or unwritable work dir)." >&2
    echo "ERROR: the run is incomplete -- free space on /scratch and resume." >&2
    NF_STATUS=1
fi

echo
echo "${TIER} finished with status ${NF_STATUS}."
echo "Compare against the baseline with: ${TIER_DIR}/99_compare_to_baseline.sh ${TIER}"
exit $NF_STATUS
