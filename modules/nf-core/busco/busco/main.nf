process BUSCO_BUSCO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c607f319867d96a38c8502f751458aa78bbd18fe4c7c4fa6b9d8350e6ba11ebe/data'
        : 'community.wave.seqera.io/library/busco_sepp:f2dbc18a2f7a5b64'}"

    input:
    tuple val(meta), path(fasta, stageAs:'tmp_input/*'), path(busco_db) // 
    val mode                              // Required:    One of genome, proteins, or transcriptome                      
   

    output:
    tuple val(meta), path("*busco*.batch_summary.txt")                                        , emit: batch_summary
    tuple val(meta), path("*.short_summary.txt")                                              , emit: short_summaries_txt , optional: true
    tuple val(meta), path("*.short_summary.json")                                             , emit: short_summaries_json, optional: true
    tuple val(meta), path("*busco*.logs")                                                      , emit: log                 , optional: true
    tuple val(meta), path("*busco*.full_table.tsv")                                   , emit: full_table          , optional: true
    tuple val(meta), path("*busco*.missing_busco_list.tsv")                           , emit: missing_busco_list  , optional: true
    // tuple val(meta), path("*busco*/*/run_*/single_copy_proteins.faa")                         , emit: single_copy_proteins, optional: true
    tuple val(meta), path("*busco*.busco_sequences.tar.gz*")                                  , emit: seq_dir             , optional: true
    // tuple val(meta), path("*busco*/*/translated_proteins")                                    , emit: translated_dir      , optional: true
    // tuple val(meta), path("*busco*")                                                          , emit: busco_dir
    // tuple val(meta), path("busco_downloads/lineages/*")                                       , emit: downloaded_lineages , optional: true
    // tuple val(meta), path("*busco*/*/run_*/busco_sequences/single_copy_busco_sequences/*.faa"), emit: single_copy_faa     , optional: true
    // tuple val(meta), path("*busco*/*/run_*/busco_sequences/single_copy_busco_sequences/*.fna"), emit: single_copy_fna     , optional: true

    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (mode !in ['genome', 'proteins', 'transcriptome']) {
        error("Mode must be one of 'genome', 'proteins', or 'transcriptome'.")
    }
    def args = task.ext.args ?: ''
    def db_used = meta.class == 'Actinopteri' ? "acti" : "vert"
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}.busco.${db_used}"
    busco_lineage = "--lineage_dataset ${busco_db}" // Function in workflow determines busco databased used in OceanOmics pipeline.
    // Use persistent output directory for checkpoint/resume capability
    def checkpoint_base = "${params.outdir}/busco_checkpoints"
    def output_dir = "${checkpoint_base}/${prefix}"    

    """
    # Fix Augustus for Apptainer
    ENV_AUGUSTUS=/opt/conda/etc/conda/activate.d/augustus.sh
    set +u
    if [ -z "\${AUGUSTUS_CONFIG_PATH}" ] && [ -f "\${ENV_AUGUSTUS}" ]; then
        source "\${ENV_AUGUSTUS}"
    fi
    set -u

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    mkdir -p $output_dir

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out_path ${output_dir} \\
        --mode ${mode} \\
        --restart \\
        ${busco_lineage} \\
        ${args}

    # clean up
    rm -rf "\$INPUT_SEQS"

    # find and remove broken symlinks from the cleanup
    find . -xtype l -delete

    
    # Move files to avoid staging/publishing issues
    cp ${output_dir}/*/batch_summary.txt ${prefix}.batch_summary.txt
    cp ${output_dir}/*/*/*/short_summary.txt ${prefix}.short_summary.txt
    cp ${output_dir}/*/*/*/short_summary.json ${prefix}.short_summary.json
    cp ${output_dir}/*/*/*/full_table.tsv ${prefix}.full_table.tsv
    cp -r ${output_dir}/*/*/*/busco_sequences ${prefix}.busco_sequences
    cp ${output_dir}/*/*/*/missing_busco_list.tsv ${prefix}.missing_busco_list.tsv
    cp -r ${output_dir}/*/*/logs ${prefix}.logs
    cp ${output_dir}/*/*/*/logs/* ${prefix}.logs

    tar -czvf ${prefix}.busco_sequences.tar.gz ${prefix}.busco_sequences
    md5sum ${prefix}.busco_sequences.tar.gz > ${prefix}.busco_sequences.tar.gz.md5 &&  rm -rf ${prefix}.busco_sequences
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def fasta_name = files(fasta).first().name - '.gz'
    """
    touch ${prefix}-busco.batch_summary.txt
    mkdir -p ${prefix}-busco/${fasta_name}/run_${lineage}/busco_sequences

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
