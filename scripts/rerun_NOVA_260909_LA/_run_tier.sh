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
MITO_PIPELINE_DIR="/software/projects/pawsey1348/$USER/Oceanomics-OceanGenomes-Mitogenomes"
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
echo "=============================================================="

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

cd "$OUT"

nextflow -log ".nextflow_${RUN}_${TIER}.log" \
    run "$RUN_DIR/main.nf" \
    -resume \
    -name "${RUN}_${TIER}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir "./work/${RUN}_${TIER}" \
    -c "$RUN_DIR/pawsey_profile.config" \
    -profile singularity \
    -with-report \
    --run "$RUN" \
    --outdir "$OUT" \
    --mitogenome_nfcore_dir "$MITO_PIPELINE_DIR" \
    --kvalue "21" \
    --genomescope2_l false \
    --genomescope2_m 10000 \
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
    --taxonkit_db_dir "$BASE" \
    --input "$SAMPLESHEET"
NF_STATUS=$?

# Nextflow exits 0 even when a task was never dispatched because its work directory
# could not be created (/scratch quota exhausted), leaving a stale report published as
# this run's output. Same guard as the main launcher.
NF_LOG="$OUT/.nextflow_${RUN}_${TIER}.log"
if [ -f "$NF_LOG" ] && grep -qE 'Disk quota exceeded|Unable to create directory=' "$NF_LOG"; then
    echo "ERROR: filesystem errors in $NF_LOG (disk quota or unwritable work dir)." >&2
    echo "ERROR: the run is incomplete -- free space on /scratch and resume." >&2
    NF_STATUS=1
fi

echo
echo "${TIER} finished with status ${NF_STATUS}."
echo "Compare against the baseline with: ${TIER_DIR}/99_compare_to_baseline.sh ${TIER}"
exit $NF_STATUS
