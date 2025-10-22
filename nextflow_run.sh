module load nextflow/24.10.0
module load singularity/4.1.0-nompi

RUN=NOVA_250606_ADJP

nextflow run main.nf \
    -work-dir ./work/$RUN \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --run "$RUN" \
    --outdir /scratch/pawsey0964/$USER/_NFCORE/_OUT_DIR \
    --mitogenome_nfcore_dir /scratch/pawsey0964/$USER/_NFCORE/Oceanomics-OceanGenomes-Mitogenomes \
    --kvalue "21" \
    --genomescope2_l false \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --tempdir /scratch/pawsey0964/$USER/tmp \
    --refresh-modules \
    --skip_download_reads true \
    --skip_fastp_fastqc true \
    --skip_genome_assembly true \
    --skip_genome_decontamination true \
    --skip_genome_qc false \
    --skip_upload_results false \


    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
