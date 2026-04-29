process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/88/889a182b8066804f4799f3808a5813ad601381a8a0e3baa4ab8d73e739b97001/data' :
        'community.wave.seqera.io/library/fastp:0.24.0--62c97b06e8447690' }"

    input:
    tuple val(meta), path(reads)
    path  adapter_fasta
    val   discard_trimmed_pass
    val   save_trimmed_fail
    val   save_merged

    output:
    tuple val(meta), path('*.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    tuple val(meta), path("20_fastp.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def base_prefix = task.ext.prefix ?: "${meta.id}.ilmn"
    def prefix = task.ext.prefix ?: "${base_prefix}.${meta.date}"
    def prefix_json = task.ext.prefix ?: "${base_prefix}.${meta.run}"
    def effective_args = ['--dedup', '--cut_tail', '--verbose', '--max_len1 300', '--max_len2 300', '--length_required 100', "--thread ${task.cpus}"].join(' ')
    def note = args ? "Paired-end trimming with JSON/HTML reports. Configured ext.args <samp>${args}</samp> is defined but not interpolated by this module." : 'Paired-end trimming with JSON/HTML reports.'
   
    """
    fastp \\
        --dedup \\
        --cut_tail \\
        --in1 ${reads[0]} \\
        --out1 ${prefix}.R1.fastq.gz \\
        --in2 ${reads[1]} \\
        --out2 ${prefix}.R2.fastq.gz \\
        --verbose \\
        --max_len1 300 \\
        --max_len2 300 \\
        --length_required 100 \\
        --json '${prefix_json}.fastp.json' \\
        --html '${prefix}.fastp.html' \\
        --report_title="${prefix} fastp" \\
        --thread $task.cpus \\
        2>&1 | tee ${prefix}.fastp.log

    cat <<-END_TOOL_PARAMS > 20_fastp.tool_params_mqcrow.html
    <tr><td>FASTP</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """


    stub:
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def effective_args      = ['--dedup', '--cut_tail', '--verbose', '--max_len1 300', '--max_len2 300', '--length_required 100', "--thread ${task.cpus}"].join(' ')
    def note                = args ? "Paired-end trimming with JSON/HTML reports. Configured ext.args <samp>${args}</samp> is defined but not interpolated by this module." : 'Paired-end trimming with JSON/HTML reports.'
    def touch_reads         = (discard_trimmed_pass) ? "" : (is_single_output) ? "echo '' | gzip > ${prefix}.fastp.fastq.gz" : "echo '' | gzip > ${prefix}_1.fastp.fastq.gz ; echo '' | gzip > ${prefix}_2.fastp.fastq.gz"
    def touch_merged        = (!is_single_output && save_merged) ? "echo '' | gzip >  ${prefix}.merged.fastq.gz" : ""
    def touch_fail_fastq    = (!save_trimmed_fail) ? "" : meta.single_end ? "echo '' | gzip > ${prefix}.fail.fastq.gz" : "echo '' | gzip > ${prefix}.paired.fail.fastq.gz ; echo '' | gzip > ${prefix}_1.fail.fastq.gz ; echo '' | gzip > ${prefix}_2.fail.fastq.gz"
    """
    $touch_reads
    $touch_fail_fastq
    $touch_merged
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"
    cat <<-END_TOOL_PARAMS > 20_fastp.tool_params_mqcrow.html
    <tr><td>FASTP</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
