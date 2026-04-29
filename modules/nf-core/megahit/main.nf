process MEGAHIT {
    tag "${meta.id}"
    label 'process_extra_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2cb827988dca7067ff8096c37cb20bc841c878013da52ad47a50865d54efe83/data' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.v129mh.fasta")                            , emit: contigs
    // tuple val(meta), path("intermediate_contigs/k*.contigs.fa")      , emit: k_contigs
    // tuple val(meta), path("intermediate_contigs/k*.addi.fa")         , emit: addi_contigs
    // tuple val(meta), path("intermediate_contigs/k*.local.fa")        , emit: local_contigs
    // tuple val(meta), path("intermediate_contigs/k*.final.contigs.fa"), emit: kfinal_contigs
    // tuple val(meta), path('*.log')                                      , emit: log
    tuple val(meta), path("22_megahit.tool_params_mqcrow.html")        , emit: tool_params
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def memory = task.memory.toBytes()
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def reads_command = meta.single_end || !reads[1] ? "-r ${reads[0].join(',')}" : "-1 ${reads[0].join(',')} -2 ${reads[1].join(',')}"
    def checkpoint_base = "${params.outdir}/megahit_checkpoints"
    def output_dir = "${checkpoint_base}/${prefix}_megahit_out"
    def effective_args = [args, "-m ${memory}", "-t ${task.cpus}", reads_command, "--out-prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    def note = "Assembles the trimmed reads to ${prefix}.v129mh.fasta and resumes from ${output_dir} when checkpoints are present."
    """
    # Always use the same output directory for this sample
    # This way MEGAHIT can resume if it exists
    if [ -d "${output_dir}" ] && [ -f "${output_dir}/checkpoints.txt" ]; then
        echo "Found existing MEGAHIT checkpoint, resuming..."
        megahit \\
            ${args} \\
            -m ${memory} \\
            -t ${task.cpus} \\
            --continue \\
            -o ${output_dir}
    else
        echo "Starting fresh MEGAHIT assembly..."
        mkdir -p ${checkpoint_base}
        megahit \\
            ${args} \\
            -m ${memory} \\
            -t ${task.cpus} \\
            ${reads_command} \\
            --out-prefix ${prefix} \\
            -o ${output_dir}
    fi

    mv ${output_dir}/${prefix}.contigs.fa ${prefix}.v129mh.fasta

    cat <<-END_TOOL_PARAMS > 22_megahit.tool_params_mqcrow.html
    <tr><td>MEGAHIT</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end || !reads[1] ? "-r ${reads[0].join(',')}" : "-1 ${reads[0].join(',')} -2 ${reads[1].join(',')}"
    def memory = task.memory.toBytes()
    def output_dir = "${params.outdir}/megahit_checkpoints/${prefix}_megahit_out"
    def effective_args = [args, "-m ${memory}", "-t ${task.cpus}", reads_command, "--out-prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    def note = "Assembles the trimmed reads to ${prefix}.v129mh.fasta and resumes from ${output_dir} when checkpoints are present."
    """
    touch ${prefix}.v129mh.fasta
    cat <<-END_TOOL_PARAMS > 22_megahit.tool_params_mqcrow.html
    <tr><td>MEGAHIT</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
