process TRIGGER_MITOGENOME {
    tag "trigger_mitogenome"
    
    input:
    val all_samples // Used to trigger this process after all samples are collected
    
    output:
    path "run_*.sh", emit: mito_run_script
    
    
script:
def input_pattern = "${params.outdir}/draftgenomes/*/fastp/*.fastq.gz"

"""
# Generate timestamp once
TIMESTAMP=\$(date +%y%m%d_%H%M%S)
SCRIPT_NAME="run_${params.run}_mitogenomes_\${TIMESTAMP}.sh"

cat > \$SCRIPT_NAME << "EOF"
module load nextflow/25.04.6
module load singularity/4.1.0-nompi

# Get the absolute path to the current directory
RUN_DIR="\$(pwd -P)"
BASE="/scratch/pawsey1348/\$USER"
OUT_DIR="\$BASE/${params.run}_mitogenomes"

mkdir -p "\$OUT_DIR"

# Get the absolute path to the output directory. This ensures that the script can be run from any location and still correctly reference the output directory.
OUT_DIR="\$(cd "\$OUT_DIR" && pwd -P)"

# Stage the backup scripts alongside the results and bake this run's outdir into
# them, so the backup is a no-argument `sbatch $OUT_DIR/backup_scripts/backup.sh`
# once the pipeline finishes.
mkdir -p "\$OUT_DIR/backup_scripts"
command cp -r "\$RUN_DIR/backup_scripts/." "\$OUT_DIR/backup_scripts/"
sed -i "s|^RUNDIR_DEFAULT=.*|RUNDIR_DEFAULT=\\\"\$OUT_DIR\\\"|" "\$OUT_DIR/backup_scripts/backup.sh"

if ! grep -q "^RUNDIR_DEFAULT=\\\"\$OUT_DIR\\\"\$" "\$OUT_DIR/backup_scripts/backup.sh"; then
    echo "WARNING: could not bake RUNDIR into \$OUT_DIR/backup_scripts/backup.sh;" \\
         "run it with -r \$OUT_DIR"
fi

# Change to output directory to run Nextflow there
cd \$OUT_DIR

nextflow -log \$OUT_DIR/.nextflow_${params.run}.log \\
    run \$RUN_DIR/main.nf \\
    -work-dir ./work/${params.run}_mitogenomes \\
    -c \$RUN_DIR/pawsey_profile.config \\
    -resume \\
    -profile singularity \\
    -with-report \\
    --input_dir \"${input_pattern}\" \\
    --outdir "\$OUT_DIR" \\
    --blast_db_dir \$(realpath ../blast_dbs) \\
    --taxonkit_db_dir \$(realpath ../) \\
    --curated_blast_db /software/projects/pawsey0964/curated_db/OceanGenomes.CuratedNT.NBDLTranche1and2and3.CuratedBOLD.NoDuplicate.fasta \\
    --nt_blast_db /scratch/references/blastdb_update/blast-2026-02-01/db/mito \\
    --mitos_refdb /software/projects/pawsey0964/mitos_refdb \\
    --mitos_refseq_ver refseq89m  \\
    --organelle_type \"animal_mt\" \\
    --kvalue \"21\" \\
    --bs_config "${params.bs_config}" \\
    --sql_config "${params.sql_config}" \\
    --enable_oatk_fallback true \
    --oatk_mito_db /software/projects/pawsey0964/oatk_db/actinopterygii_mito.fam \
    --binddir /scratch \\
    --tempdir /scratch/pawsey0964/\$USER/tmp \\
    --refresh-modules \\
    --skip_mitogenome_assembly_getorg false \\
    --skip_mitogenome_assembly_hifi false \\
    --skip_mitogenome_annotation false \\
    --skip_upload_results false \\
    --samplesheet_prefix \"samplesheet\" \\
    --template_sbt \"/home/\$USER/template.sbt\" \\
    --force_db_overwrite false \\
    --translation_table \"2\" \\
    --ena_webin_validate true \\
    --ena_study "PRJEB110568" \\
    --ena_validation_attempt "initial"
EOF

chmod +x \$SCRIPT_NAME
"""
}