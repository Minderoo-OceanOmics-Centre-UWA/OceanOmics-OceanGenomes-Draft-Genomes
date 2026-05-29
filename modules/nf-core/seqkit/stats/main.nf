process SEQKIT_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0':
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.seqkit_stats.tsv")                  , emit: stats
    tuple val(meta), path("80_seqkit_stats.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    def effective_args = [args, "--tabular", "--threads ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def note = "Per-sequence summary statistics computed with seqkit stats; writes ${prefix}.seqkit_stats.tsv."
    """
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${args} \\
        ${reads} > ${prefix}.seqkit_stats.tsv

    cat <<-END_TOOL_PARAMS > 80_seqkit_stats.tool_params_mqcrow.html
    <tr><td>Seqkit stats</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/^.*v//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: '--all'
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    def effective_args = [args, "--tabular", "--threads ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def note = "Per-sequence summary statistics computed with seqkit stats; writes ${prefix}.seqkit_stats.tsv."
    """
    touch ${prefix}.seqkit_stats.tsv

    cat <<-END_TOOL_PARAMS > 80_seqkit_stats.tool_params_mqcrow.html
    <tr><td>Seqkit stats</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/^.*v//' )
    END_VERSIONS
    """
}
