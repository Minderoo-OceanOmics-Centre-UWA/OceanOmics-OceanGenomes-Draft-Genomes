process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py311r42hdfd78af_6':
        'biocontainers/genomescope2:2.0--py311r42hdfd78af_6' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("${prefix}_linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("${prefix}_transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("${prefix}_log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("${prefix}_transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("${prefix}_model.txt")                  , emit: model
    tuple val(meta), path("${prefix}_summary.txt")                , emit: summary
    // tuple val(meta), path("${prefix}_lookup_table.txt")           , emit: lookup_table, optional: true
    // tuple val(meta), path("${prefix}_fitted_hist.png")            , emit: fitted_histogram_png, optional: true
    tuple val(meta), path("16_genomescope2.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-k 21 -m 1000'
    prefix = task.ext.prefix ?: "${meta.prefix}_genomescope2"
    def effective_args = ["--input ${histogram}", args, '--output .', "--name_prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    def note = 'K-mer spectrum model fit from the sample meryl histogram.'
    """
    genomescope2 \\
        --input $histogram \\
        $args \\
        --output . \\
        --name_prefix ${prefix}

    test -f "fitted_hist.png" && mv fitted_hist.png ${prefix}_fitted_hist.png
    test -f "lookup_table.txt" && mv lookup_table.txt ${prefix}_lookup_table.txt
    cat <<-END_TOOL_PARAMS > 16_genomescope2.tool_params_mqcrow.html
    <tr><td>Genomescope2</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def effective_args = ["--input ${histogram}", args, '--output .', "--name_prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    def note = 'K-mer spectrum model fit from the sample meryl histogram.'
    """
    touch ${prefix}_linear_plot.png
    touch ${prefix}_transformed_linear_plot.png
    touch ${prefix}_log_plot.png
    touch ${prefix}_transformed_log_plot.png
    touch ${prefix}_model.txt
    touch ${prefix}_summary.txt
    touch ${prefix}_fitted_hist.png
    touch ${prefix}_lookup_table.txt
    cat <<-END_TOOL_PARAMS > 16_genomescope2.tool_params_mqcrow.html
    <tr><td>Genomescope2</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """
}
