// Classify reads against a kraken2 database so non-target reads can be removed before
// assembly. Adapted from nf-core/modules kraken2/kraken2 to follow this pipeline's
// conventions: a versions.yml file rather than a versions topic channel, and a
// tool_params row for the per-sample MultiQC report.
process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    tuple val(meta), path(reads)
    path  db

    output:
    tuple val(meta), path("*.kraken2.report.txt")            , emit: report
    tuple val(meta), path("*.kraken2.classifiedreads.txt.gz"), emit: classified_reads_assignment
    tuple val(meta), path("14_kraken2.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def paired = meta.single_end ? '' : '--paired'
    def read_list = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    def effective_args = ["kraken2 --db ${db}", "--threads ${task.cpus}", '--gzip-compressed', paired, args].findAll { it?.trim() }.join(' ')
    def note = "Classifies the trimmed reads against ${db} so reads from non-target clades can be removed before k-mer counting and assembly."
    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --output ${prefix}.kraken2.classifiedreads.txt \\
        --gzip-compressed \\
        $paired \\
        $args \\
        $read_list

    # The per-read assignment file carries one line per read pair, so it is large for
    # a whole library. Compress it for publishing; the extraction step reads it back
    # with zcat.
    pigz -p ${task.cpus} ${prefix}.kraken2.classifiedreads.txt

    cat <<-END_TOOL_PARAMS > 14_kraken2.tool_params_mqcrow.html
    <tr><td>Kraken2</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version 2>&1 | head -1 | sed 's/^.*Kraken version //; s/ .*//')
        pigz: \$(pigz --version 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def paired = meta.single_end ? '' : '--paired'
    def effective_args = ["kraken2 --db ${db}", "--threads ${task.cpus}", '--gzip-compressed', paired, args].findAll { it?.trim() }.join(' ')
    def note = "Classifies the trimmed reads against ${db} so reads from non-target clades can be removed before k-mer counting and assembly."
    """
    touch ${prefix}.kraken2.report.txt
    echo | gzip > ${prefix}.kraken2.classifiedreads.txt.gz

    cat <<-END_TOOL_PARAMS > 14_kraken2.tool_params_mqcrow.html
    <tr><td>Kraken2</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: 2.1.3
        pigz: 2.8
    END_VERSIONS
    """
}
