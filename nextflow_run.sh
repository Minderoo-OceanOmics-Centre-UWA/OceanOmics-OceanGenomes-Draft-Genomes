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
    --kraken2_db "/scratch/references/kraken2/pluspfp_20230605" \
    --kraken2_confidence 0.05 \
    --kraken2_exclude_taxids "2 2157 10239" \
    --taxonkit_db_dir "/scratch/pawsey1348/$USER" \
    --ramdisk_path "/tmp/gxdb/" \
    --busco_acti_db "/scratch/references/busco_db/actinopterygii_odb10" \
    --busco_vert_db "/scratch/references/busco_db/vertebrata_odb10" \
    --busco_metazoa_db "/software/projects/pawsey0964/busco_db/metazoa_odb12" \
    --busco_crustacea_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/crustacea_odb12" \
    --busco_mollusca_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/mollusca_odb12" \
    --busco_anthozoa_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/anthozoa_odb12" \
    --busco_arthropoda_db "/software/projects/pawsey1348/tpeirce/busco_db/busco_downloads/lineages/arthropoda_odb12" \
    --tempdir "/scratch/pawsey1348/$USER/tmp" \
    --refresh-modules \
    --skip_bs_download false \
    --skip_download_reads false \
    --skip_fastp_fastqc false \
    --skip_genome_assembly false \
    --skip_kraken2_decontamination false \
    --skip_genome_decontamination false \
    --skip_genome_qc false \
    --skip_upload_results false \
    # --input assets/samplesheet.csv
    
NF_STATUS=$?

# --- Filesystem guard ------------------------------------------------------
# Nextflow exits 0 even when a task was never dispatched because its work
# directory could not be created (/scratch quota exhausted). The run then looks
# successful while late processes -- MULTIQC above all -- never ran, leaving a
# stale report from an earlier attempt published as this run's output. Catch
# that here so the launcher's exit status reflects it.
NF_LOG=".nextflow_${RUN}.log"
if [ -f "$NF_LOG" ] && grep -qE 'Disk quota exceeded|Unable to create directory=' "$NF_LOG"; then
    echo "ERROR: filesystem errors in $NF_LOG (disk quota or unwritable work dir)." >&2
    echo "ERROR: the run is incomplete -- free space on /scratch and resume." >&2
    NF_STATUS=1
fi

exit $NF_STATUS
