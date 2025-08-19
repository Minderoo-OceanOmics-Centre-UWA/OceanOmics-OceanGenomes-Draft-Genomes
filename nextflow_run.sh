module load nextflow/24.10.0
module load singularity/4.1.0-nompi

nextflow run main.nf \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --outdir /scratch/pawsey0964/tpeirce/_NFCORE/_OUT_DIR \
    --curated_blast_db /scratch/pawsey0964/pbayer/OceanGenomes.CuratedNT.NBDLTranche1and2.CuratedBOLD.fasta \
    --organelle_type "animal_mt" \
    --kvalue "21" \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --binddir /scratch \
    --tempdir /scratch/pawsey0964/tpeirce/tmp \
    --refresh-modules \
    --skip_download_reads true \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_mitogenome_assembly false \
    --skip_mitogenome_annotation false \
    --skip_upload_results false \
    
   
    #--input assets/samplesheet.csv \  # include a samplesheet if you are not downloading sample.
    
