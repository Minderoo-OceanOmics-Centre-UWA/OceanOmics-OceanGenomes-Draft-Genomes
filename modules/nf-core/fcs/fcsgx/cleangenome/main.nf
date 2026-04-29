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
    tuple val(meta), path("*.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    def row_filename = cleaned_suffix == 'rmadapt' ? '27_fcsgx_cleangenome_adaptor.tool_params_mqcrow.html' : '24_fcsgx_cleangenome.tool_params_mqcrow.html'
    def tool_label = cleaned_suffix == 'rmadapt' ? 'FCSGX Clean Genome Adaptor' : 'FCSGX Clean Genome'
    def note = cleaned_suffix == 'rmadapt' ? 'Applies the adaptor-screen action report to generate the adaptor-cleaned assembly.' : 'Applies the FCS-GX action report to generate cleaned and contaminant assemblies.'
    def effective_args = ['gx clean-genome', "--input ${fasta}", "--action-report ${action_report}", "--output ${prefix}.${cleaned_suffix}.fasta", "--contam-fasta-out ${prefix}.${contam_suffix}.fasta", args].findAll { it?.trim() }.join(' ')
    """
    gx \\
        clean-genome \\
        --input ${fasta} \\
        --action-report ${action_report} \\
        --output ${prefix}.${cleaned_suffix}.fasta \\
        --contam-fasta-out ${prefix}.${contam_suffix}.fasta \\
        ${args}
    cat <<-END_TOOL_PARAMS > ${row_filename}
    <tr><td>${tool_label}</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def row_filename = cleaned_suffix == 'rmadapt' ? '27_fcsgx_cleangenome_adaptor.tool_params_mqcrow.html' : '24_fcsgx_cleangenome.tool_params_mqcrow.html'
    def tool_label = cleaned_suffix == 'rmadapt' ? 'FCSGX Clean Genome Adaptor' : 'FCSGX Clean Genome'
    def note = cleaned_suffix == 'rmadapt' ? 'Applies the adaptor-screen action report to generate the adaptor-cleaned assembly.' : 'Applies the FCS-GX action report to generate cleaned and contaminant assemblies.'
    def effective_args = ['gx clean-genome', "--input ${fasta}", "--action-report ${action_report}", "--output ${prefix}.${cleaned_suffix}.fasta", "--contam-fasta-out ${prefix}.${contam_suffix}.fasta", args].findAll { it?.trim() }.join(' ')
    """
    touch ${prefix}.${cleaned_suffix}.fasta
    touch ${prefix}.${contam_suffix}.fasta
    cat <<-END_TOOL_PARAMS > ${row_filename}
    <tr><td>${tool_label}</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """
}
