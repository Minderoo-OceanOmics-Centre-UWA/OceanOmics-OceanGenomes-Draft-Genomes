module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Define the run name and base directory
RUN=XXXX_0000_XX
BASE="/scratch/pawsey0964/$USER"
# Outdir is made inside the base directory
OUT="${BASE}/${RUN}"
mkdir -p $OUT
OUT="$(cd "$OUT" && pwd -P)"
cp -r scripts/backup_scripts $OUT
# Stamp RUN into copied backup config
sed -i "s/^RUN=.*/RUN=$RUN/" "$OUT/backup_scripts/configfile.txt"
# Get the absolute path to the current directory
RUN_DIR="$(pwd -P)"
# Create output directory and then change into this directory to run the nextflow pipeline.
# This allows yout to run multiple runs from the same nf-core directory without conflicts.

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
    --busco_metazoa_db "/software/projects/pawsey0964/busco_db/metazoa_odb12" \
    --tempdir $BASE/tmp \
    --refresh-modules \
    --skip_bs_download false \
    --skip_download_reads false \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_upload_results false \
    # --input $OUT/samplesheet/${RUN}_samplesheet.csv
    
    
