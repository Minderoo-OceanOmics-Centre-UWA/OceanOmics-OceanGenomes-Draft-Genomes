process MERYL_COUNT {
    tag "$meta.id"
    label 'process_extra_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_1':
        'biocontainers/meryl:1.4.1--h4ac6f70_1' }"

    input:
    tuple val(meta), path(reads)
    val kvalue

    output:
    tuple val(meta), path("*.meryl")    , emit: meryl_dbs
    tuple val(meta), path("17_meryl_count.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def reduced_mem = task.memory.multiply(0.9).toGiga()
    def read_list = reads.collect { it.getName() }.join(' and ')
    def effective_args = ["k=${kvalue}", "threads=${task.cpus}", "memory=${reduced_mem}", args].findAll { it?.trim() }.join(' ')
    def note = "Counts ${kvalue}-mers separately for ${read_list}."
    """
    for READ in ${reads}; do
        meryl count \\
            k=${kvalue} \\
            threads=${task.cpus} \\
            memory=${reduced_mem} \\
            ${args} \\
            \$READ \\
            output \${READ%.f*}.meryl
    done
    cat <<-END_TOOL_PARAMS > 17_meryl_count.tool_params_mqcrow.html
    <tr><td>Meryl Count</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& awk '{ print \$2 }' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reduced_mem = task.memory.multiply(0.9).toGiga()
    def read_list = reads.collect { it.getName() }.join(' and ')
    def effective_args = ["k=${kvalue}", "threads=${task.cpus}", "memory=${reduced_mem}", args].findAll { it?.trim() }.join(' ')
    def note = "Counts ${kvalue}-mers separately for ${read_list}."
    """
    for READ in ${reads}; do
        touch ${prefix}.\${READ%.f*}.meryl
    done
    cat <<-END_TOOL_PARAMS > 17_meryl_count.tool_params_mqcrow.html
    <tr><td>Meryl Count</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& awk '{ print \$2 }' )
    END_VERSIONS
    """
}
