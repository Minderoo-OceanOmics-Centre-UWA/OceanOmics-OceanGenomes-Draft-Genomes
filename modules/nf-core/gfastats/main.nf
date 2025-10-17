process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.6--hdcf5f25_3':
        'biocontainers/gfastats:1.3.6--hdcf5f25_3' }"

    input:
    tuple val(meta), path(assembly), path(summary)
    val out_fmt
 
    output:
    tuple val(meta), path("*.assembly_summary"), emit: assembly_summary
    tuple val(meta), path("*.${out_fmt}")   , emit: assembly        , optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    def output_sequences = out_fmt ? "--out-format ${prefix}.${out_fmt}" : ""
    """
    # Get genomesize from $summary file
    genome_size=\$(cat $summary | grep 'Genome Haploid Length' | grep -o 'bp.*' | sed 's/bp//g' | sed 's/ //g' | sed 's/,//g')

    gfastats \\
        $args \\
        --threads $task.cpus \\
        $output_sequences \\
        $assembly \\
        \$genome_size \\
        > ${prefix}.assembly_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.${out_fmt}.gz
    touch ${prefix}.assembly_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """
}
