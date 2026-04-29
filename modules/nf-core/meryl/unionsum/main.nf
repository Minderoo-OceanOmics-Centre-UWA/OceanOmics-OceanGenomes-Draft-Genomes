process MERYL_UNIONSUM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_0':
        'biocontainers/meryl:1.4.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(meryl_dbs)
    val kvalue

    output:
    tuple val(meta), path("*.meryl"), emit: meryl_db
    tuple val(meta), path("19_meryl_unionsum.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def input_list = meryl_dbs.collect { it.getName() }.join(', ')
    def effective_args = ["threads=${task.cpus}", "memory=${task.memory.toGiga()}", args].findAll { it?.trim() }.join(' ')
    def note = "Merges ${input_list} into ${prefix}.meryl."
    """
    meryl union-sum \\
        threads=$task.cpus \\
        memory=${task.memory.toGiga()} \\
        $args \\
        output ${prefix}.meryl \\
        $meryl_dbs
    cat <<-END_TOOL_PARAMS > 19_meryl_unionsum.tool_params_mqcrow.html
    <tr><td>Meryl Union-sum</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = meryl_dbs.collect { it.getName() }.join(', ')
    def effective_args = ["threads=${task.cpus}", "memory=${task.memory.toGiga()}", args].findAll { it?.trim() }.join(' ')
    def note = "Merges ${input_list} into ${prefix}.meryl."
    """
    touch ${prefix}.unionsum.meryl
    cat <<-END_TOOL_PARAMS > 19_meryl_unionsum.tool_params_mqcrow.html
    <tr><td>Meryl Union-sum</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
