process FASTQC {
    tag "${meta.id}"
    //label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    tuple val(meta), path("10_fastqc.tool_params_mqcrow.html"), emit: tool_params
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when ?: true

    script:
    def args          = task.ext.args ?: ''
    //def prefix        = task.ext.prefix ?: "${meta.id}.ilmn.${params.run}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    //def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    //def rename_to     = old_new_pairs*.join(' ').join(' ')
    //def renamed_files = old_new_pairs.collect{ _old_name, new_name -> new_name }.join(' ')

    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpu allows to stick to requested amount of RAM in the label
     def memory_in_mb = task.memory ? (task.memory.toUnit('MB').toFloat() / task.cpus).toInteger() : null
    // // FastQC memory value allowed range (100 - 10000)
     def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    def effective_args = [args, "--threads ${task.cpus}", "--memory ${fastqc_memory}"].findAll { it?.trim() }.join(' ')
    def input_count = reads instanceof List ? reads.size() : 1

    """

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        --memory ${fastqc_memory} \\
        ${reads.join(' ')}

    cat <<-END_TOOL_PARAMS > 10_fastqc.tool_params_mqcrow.html
    <tr><td>FASTQC</td><td><samp>${effective_args}</samp></td><td>${input_count} input FASTQ file(s).</td></tr>
    END_TOOL_PARAMS


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory_in_mb = task.memory ? (task.memory.toUnit('MB').toFloat() / task.cpus).toInteger() : null
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    def effective_args = [args, "--threads ${task.cpus}", "--memory ${fastqc_memory}"].findAll { it?.trim() }.join(' ')
    def input_count = reads instanceof List ? reads.size() : 1
    """
    touch ${prefix}.html
    touch ${prefix}.zip
    cat <<-END_TOOL_PARAMS > 10_fastqc.tool_params_mqcrow.html
    <tr><td>FASTQC</td><td><samp>${effective_args}</samp></td><td>${input_count} input FASTQ file(s).</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
