process MERYL_HISTOGRAM {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.4.1--h4ac6f70_0':
        'biocontainers/meryl:1.4.1--h4ac6f70_0' }"

    input:
    tuple val(meta), path(meryl_db)
    val kvalue

    output:
    tuple val(meta), path("*.hist"), emit: hist
    tuple val(meta), path("18_meryl_histogram.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def effective_args = ["threads=${task.cpus}", "memory=${task.memory.toGiga()}", args].findAll { it?.trim() }.join(' ')
    def note = "Builds ${prefix}.meryl.hist from ${meryl_db.getName()}."
    """
    meryl histogram \\
        threads=$task.cpus \\
        memory=${task.memory.toGiga()} \\
        $args \\
        $meryl_db > ${prefix}.meryl.hist
    cat <<-END_TOOL_PARAMS > 18_meryl_histogram.tool_params_mqcrow.html
    <tr><td>Meryl Histogram</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def effective_args = ["threads=${task.cpus}", "memory=${task.memory.toGiga()}", args].findAll { it?.trim() }.join(' ')
    def note = "Builds ${prefix}.meryl.hist from ${meryl_db.getName()}."
    """
    touch ${prefix}.hist
    cat <<-END_TOOL_PARAMS > 18_meryl_histogram.tool_params_mqcrow.html
    <tr><td>Meryl Histogram</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
