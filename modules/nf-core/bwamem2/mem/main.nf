process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9a/9ac054213e67b3c9308e409b459080bbe438f8fd6c646c351bc42887f35a42e7/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:e1f420694f8e42bd' }"

    input:
    tuple val(meta), path(reads), path(index), path(fasta)


    output:
    tuple val(meta), path("*.sam")  , emit: sam , optional:true
    tuple val(meta), path("*.bam")  , emit: bam , optional:true
    tuple val(meta), path("*.cram") , emit: cram, optional:true
    tuple val(meta), path("*.crai") , emit: crai, optional:true
    tuple val(meta), path("*.csi")  , emit: csi , optional:true
    tuple val(meta), path("50_bwamem2_mem.tool_params_mqcrow.html"), emit: tool_params
    path "${meta.id}-sn_results.tsv"  , emit: results
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    def mem_args = [args, "-t ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def sort_args = ["-@ ${task.cpus}", reference].findAll { it?.trim() }.join(' ')
    def effective_args = "bwa-mem2 mem ${mem_args} | samtools view -b | samtools sort ${sort_args}".trim()
    def note = args2 ? "Produces sorted BAM and samtools stats. Configured ext.args2 <samp>${args2}</samp> is defined but not interpolated by this module." : 'Produces sorted BAM and samtools stats.'

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools view -b \\
        | samtools sort -@ $task.cpus ${reference} -o ${meta.id}.sorted.bam -
    
    samtools index ${meta.id}.sorted.bam

    samtools stats ${meta.id}.sorted.bam | grep "^SN" | cut -f 2- > ${meta.id}-sn_results.tsv
    cat <<-END_TOOL_PARAMS > 50_bwamem2_mem.tool_params_mqcrow.html
    <tr><td>BWA-MEM2 Mem</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    def reference = fasta && extension=="cram" ? "--reference ${fasta}" : ""
    def mem_args = [args, "-t ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def sort_args = ["-@ ${task.cpus}", reference].findAll { it?.trim() }.join(' ')
    def effective_args = "bwa-mem2 mem ${mem_args} | samtools view -b | samtools sort ${sort_args}".trim()
    def note = args2 ? "Produces sorted BAM and samtools stats. Configured ext.args2 <samp>${args2}</samp> is defined but not interpolated by this module." : 'Produces sorted BAM and samtools stats.'

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }

    """
    touch ${prefix}.${extension}
    ${create_index}
    cat <<-END_TOOL_PARAMS > 50_bwamem2_mem.tool_params_mqcrow.html
    <tr><td>BWA-MEM2 Mem</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
