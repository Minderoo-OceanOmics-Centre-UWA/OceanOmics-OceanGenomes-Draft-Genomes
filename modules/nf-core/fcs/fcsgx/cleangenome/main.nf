process FCSGX_CLEANGENOME {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

    input:
    tuple val(meta), path(fasta), path(action_report)
    val(cleaned_suffix)
    val(contam_suffix)

    output:
    tuple val(meta), path("*.${cleaned_suffix}.fasta")     , emit: cleaned
    tuple val(meta), path("*.${contam_suffix}.fasta"), emit: contaminants
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    """
    gx \\
        clean-genome \\
        --input ${fasta} \\
        --action-report ${action_report} \\
        --output ${prefix}.${cleaned_suffix}.fasta \\
        --contam-fasta-out ${prefix}.${contam_suffix}.fasta \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cleaned.fasta
    touch ${prefix}.contaminants.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """
}
