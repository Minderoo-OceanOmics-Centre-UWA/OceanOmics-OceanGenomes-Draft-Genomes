process BWAMEM2_INDEX {
    tag "$fasta"
    // NOTE Requires 28N GB memory where N is the size of the reference sequence
    // source: https://github.com/bwa-mem2/bwa-mem2/issues/9
    memory { 28.B * fasta.size() }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9ac054213e67b3c9308e409b459080bbe438f8fd6c646c351bc42887f35a42e7/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    tuple val(meta), path("40_bwamem2_index.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    def effective_args = ['index', args, "-p bwamem2/${prefix}"].findAll { it?.trim() }.join(' ')
    def note = "Reference input ${fasta.name}."
    """
    mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwamem2/${prefix}

    cat <<-END_TOOL_PARAMS > 40_bwamem2_index.tool_params_mqcrow.html
    <tr><td>BWA-MEM2 Index</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta}"
    def args = task.ext.args ?: ''
    def effective_args = ['index', args, "-p bwamem2/${prefix}"].findAll { it?.trim() }.join(' ')
    def note = "Reference input ${fasta.name}."

    """
    mkdir bwamem2
    touch bwamem2/${prefix}.0123
    touch bwamem2/${prefix}.ann
    touch bwamem2/${prefix}.pac
    touch bwamem2/${prefix}.amb
    touch bwamem2/${prefix}.bwt.2bit.64
    cat <<-END_TOOL_PARAMS > 40_bwamem2_index.tool_params_mqcrow.html
    <tr><td>BWA-MEM2 Index</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
