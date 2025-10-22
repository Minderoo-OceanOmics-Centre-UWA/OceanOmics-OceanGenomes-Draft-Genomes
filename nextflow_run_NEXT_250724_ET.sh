module load nextflow/24.10.0
module load singularity/4.1.0-nompi

# Define the run name and base directory
RUN=NEXT_250724_ET
BASE="/scratch/pawsey0964/$USER"
# Outdir is made inside the base directory
OUT="${BASE}/${RUN}"
# Get the absolute path to the current directory
RUN_DIR="$(realpath .)"
# Create output directory and then change into this directory to run the nextflow pipeline.
# This allows yout to run multiple runs from the same nf-core directory without conflicts.
mkdir -p $OUT
cd $OUT

nextflow -log .nextflow_$RUN.log \
    run $RUN_DIR/main.nf \
    -resume \
    -name "${RUN}_$(date +%Y%m%d_%H%M%S)" \
    -work-dir ./work/$RUN \
    -c $RUN_DIR/pawsey_profile.config \
    -profile singularity \
    -with-report \
    --run $RUN \
    --outdir $OUT \
    --mitogenome_nfcore_dir $BASE/_NFCORE/Oceanomics-OceanGenomes-Mitogenomes \
    --kvalue "21" \
    --genomescope2_l false \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --tempdir $BASE/tmp \
    --refresh-modules \
    --skip_download_reads false \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_upload_results false \

    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
