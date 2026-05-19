module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Alternative repo-local launcher.
# This runs from the repository and keeps Nextflow work/log files under this directory.
# For standard OceanGenomes production runs, prefer copying nextflow_run_template.sh
# to nextflow_run_<RUN>.sh so each run launches from its own output directory.
RUN=XXXX_0000_XX
OUT="/scratch/pawsey1348/$USER/_NFCORE/_OUT_DIR"

# To reuse an existing samplesheet, uncomment the --input line at the end
# and set --skip_download_reads true.
nextflow -log ".nextflow_${RUN}.log" \
    run main.nf \
    -work-dir "./work/$RUN" \
    -c pawsey_profile.config \
    -resume \
    -profile singularity \
    -with-report \
    --run "$RUN" \
    --outdir "$OUT" \
    --mitogenome_nfcore_dir "/scratch/pawsey1348/$USER/Oceanomics-OceanGenomes-Mitogenomes" \
    --kvalue "21" \
    --genomescope2_l false \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg \
    --gxdb "/scratch/references/Foreign_Contamination_Screening" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --tempdir "/scratch/pawsey1348/$USER/tmp" \
    --refresh-modules \
    --skip_bs_download false \
    --skip_download_reads false \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_upload_results false \
    # --input assets/samplesheet.csv
    
