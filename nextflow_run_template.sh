module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Per-run OceanGenomes launcher.
# Copy this template to nextflow_run_<RUN>.sh, set RUN below, then run the copied script.
RUN=XXXX_0000_XX
BASE="/scratch/pawsey1348/$USER"
MITO_PIPELINE_DIR="/software/projects/pawsey1348/$USER/Oceanomics-OceanGenomes-Mitogenomes"

# Outdir is made inside the base directory
OUT="${BASE}/${RUN}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd -P)"
cp -r scripts/backup_scripts "$OUT"

# Stamp RUN into copied backup config
sed -i "s/^RUN=.*/RUN=$RUN/" "$OUT/backup_scripts/configfile.txt"

# Get the absolute path to the current directory
RUN_DIR="$(pwd -P)"

# Create output directory and then change into this directory to run the nextflow pipeline.
# This allows you to run multiple runs from the same repository without conflicts.

cd "$OUT"

# To reuse an existing samplesheet, uncomment the --input line at the end
# and set --skip_bs_download and --skip_download_reads true.
nextflow -log ".nextflow_${RUN}.log" \
    run "$RUN_DIR/main.nf" \
    -resume \
    -name "${RUN}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir "./work/$RUN" \
    -c "$RUN_DIR/pawsey_profile.config" \
    -profile singularity \
    -with-report \
    --run "$RUN" \
    --outdir "$OUT" \
    --mitogenome_nfcore_dir "$MITO_PIPELINE_DIR" \
    --kvalue "21" \
    --genomescope2_l false \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --kraken2_db "/scratch/references/kraken2/pluspfp_20230605" \
    --kraken2_confidence 0.05 \
    --kraken2_exclude_taxids "2 2157 10239" \
    --taxonkit_db_dir "$BASE" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --busco_metazoa_db "/software/projects/pawsey0964/busco_db/metazoa_odb12" \
    --busco_crustacea_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/crustacea_odb12" \
    --busco_mollusca_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/mollusca_odb12" \
    --busco_anthozoa_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/anthozoa_odb12" \
    --busco_arthropoda_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/arthropoda_odb12" \
    --tempdir "$BASE/tmp" \
    --refresh-modules \
    --skip_bs_download false \
    --skip_download_reads false \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_kraken2_decontamination false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_upload_results false \
    # --input "$OUT/samplesheet/${RUN}_samplesheet.csv"
NF_STATUS=$?

# --- Filesystem guard ------------------------------------------------------
# Nextflow exits 0 even when a task was never dispatched because its work
# directory could not be created (/scratch quota exhausted). The run then looks
# successful while late processes -- MULTIQC above all -- never ran, leaving a
# stale report from an earlier attempt published as this run's output. Catch
# that here so the launcher's exit status reflects it.
NF_LOG="$OUT/.nextflow_${RUN}.log"
if [ -f "$NF_LOG" ] && grep -qE 'Disk quota exceeded|Unable to create directory=' "$NF_LOG"; then
    echo "ERROR: filesystem errors in $NF_LOG (disk quota or unwritable work dir)." >&2
    echo "ERROR: the run is incomplete -- free space on /scratch and resume." >&2
    NF_STATUS=1
fi

# --- Per-genome compute cost (best-effort, self-contained) -----------------
# Repo-local script: writes pipeline_info/cost_per_sample.csv (SU per OG sample)
# for this run. Skipped when the poller runs this (it records cost, incl. the
# central ledger, via post_draft.py). The `|| echo` keeps a cost failure from
# affecting the run's exit status. Inherits this script's NXF_HOME + modules.
if [ -z "${OCEANOMICS_SKIP_COST:-}" ]; then
    COST_SCRIPT="$RUN_DIR/compute-audit/nf_workflow_cost.sh"
    if [ -f "$COST_SCRIPT" ]; then
        mkdir -p "$OUT/pipeline_info"
        bash "$COST_SCRIPT" "$OUT" "$OUT/pipeline_info/compute_usage.csv" \
            || echo "compute cost: accounting failed (non-fatal)"
    fi
fi

exit $NF_STATUS
