process BUSCO_BUSCO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/41/4137d65ab5b90d2ae4fa9d3e0e8294ddccc287e53ca653bb3c63b8fdb03e882f/data'
        : 'community.wave.seqera.io/library/busco:6.0.0--a9a1426105f81165'}"

    input:
    tuple val(meta), path(fasta, stageAs:'tmp_input/*'), val(busco_db)
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
    tuple val(meta), path("30_busco.tool_params_mqcrow.html")                                  , emit: tool_params
    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (mode !in ['genome', 'proteins', 'transcriptome']) {
        error("Mode must be one of 'genome', 'proteins', or 'transcriptome'.")
    }
    def args = task.ext.args ?: ''
    // Determine database prefix for output naming using the first 4 letters of the database name.
    def db_used = busco_db.tokenize('/').last().take(4)
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}.busco.${db_used}"
    def effective_args = ["--cpu ${task.cpus}", "--mode ${mode}", "--lineage_dataset ${busco_db}", args].findAll { it?.trim() }.join(' ')
    def busco_db_name = busco_db.tokenize('/').last()
    def note = "Output prefix ${prefix}; lineage dataset ${busco_db_name}."


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

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix} \\
        --mode ${mode} \\
        --lineage_dataset ${busco_db} \\
        ${args}

    # find and remove broken symlinks from the cleanup
    find . -xtype l -delete

    # Move files to avoid staging/publishing issues
    mv ${prefix}/batch_summary.txt ${prefix}.batch_summary.txt
    mv ${prefix}/*/*/short_summary.txt ${prefix}.short_summary.txt
    mv ${prefix}/*/*/short_summary.json ${prefix}.short_summary.json
    mv ${prefix}/*/*/full_table.tsv ${prefix}.full_table.tsv
    mv ${prefix}/*/*/busco_sequences ${prefix}.busco_sequences
    mv ${prefix}/*/*/missing_busco_list.tsv ${prefix}.missing_busco_list.tsv
    mv ${prefix}/logs ${prefix}.logs
    mv ${prefix}/*/logs/* ${prefix}.logs

    tar -czvf ${prefix}.busco_sequences.tar.gz ${prefix}.busco_sequences
    md5sum ${prefix}.busco_sequences.tar.gz > ${prefix}.busco_sequences.tar.gz.md5 &&  rm -rf ${prefix}.busco_sequences

    if grep 'Run failed; check logs' ${prefix}.batch_summary.txt > /dev/null
    then
        echo "Busco run failed"
        exit 1
    fi

    # clean up
    rm -rf "\$INPUT_SEQS"
    rm -rf ${prefix}/*
    rm -rf tmp.*

    cat <<-END_TOOL_PARAMS > 30_busco.tool_params_mqcrow.html
    <tr><td>BUSCO</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}.busco.${db_used}"
    def fasta_name = files(fasta).first().name - '.gz'
    def effective_args = ["--cpu ${task.cpus}", "--mode ${mode}", "--lineage_dataset ${busco_db}", args].findAll { it?.trim() }.join(' ')
    def busco_db_name = busco_db.tokenize('/').last()
    def note = "Output prefix ${prefix}; lineage dataset ${busco_db_name}."
    """
    touch ${prefix}-busco.batch_summary.txt
    mkdir -p ${prefix}-busco/${fasta_name}/run_${lineage}/busco_sequences
    cat <<-END_TOOL_PARAMS > 30_busco.tool_params_mqcrow.html
    <tr><td>BUSCO</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
