module load nextflow/24.10.0
module load singularity/4.1.0-nompi

nextflow run main.nf \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --run 'NOVA_250606_ADJP' \
    --outdir /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR_testparams \
    --mitogenome_nfcore_dir /scratch/pawsey0964/tpeirce/_NFCORE/Oceanomics-OceanGenomes-Mitogenomes \
    --kvalue "21" \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    --refresh-modules \
    --skip_download_reads true \
    --skip_fastp_fastqc true \
    --skip_genome_assembly false \
    --skip_genome_decontamination true \
    --skip_genome_qc true \


# include the corresponding precomputed results if skipping an analysis    
    # --precomputed_download_reads_results '' \
    # --precomputed_fastp_fastqc_results '' \
    # --precomputed_genome_assembly_results '' \        
    # --precomputed_genomescope_summary '' \ 
    # --precomputed_meryl_results '' \
    # --precomputed_genome_decontamination_results '' \
    # --precomputed_genome_qc_results '' \
    # --input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
