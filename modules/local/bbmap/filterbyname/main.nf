process BBMAP_FILTERBYNAME {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(action_report), path(reads)
    val(output_format)
 
    output:
    tuple val(meta), path("$fully_filtered_reads")  , emit: fully_filtered_reads
    tuple val(meta), path ("$filter_report")                         , emit: filter_report
    tuple val(meta), path("$names_to_filter")                        , emit: names_to_filter 
    tuple val(meta), path("${meta.prefix}.contig_count_500bp.txt")   , emit: contigs_under_500bp
    tuple val(meta), path("25_bbmap_filterbyname.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    input  = "in=${reads}"
    first_filtered_reads = "${prefix}.rf.fa"
    fully_filtered_reads = "${prefix}.${output_format}"
    filter_report = "${meta.prefix}.filter_report.txt"
    names_to_filter = "${meta.prefix}.review_scaffolds_1kb.txt"
    contigs_under_500bp = "${meta.prefix}.contig_count_500bp.txt"
    
    def avail_mem = 3
    if (!task.memory) {
        log.info '[filterbyname] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def effective_args = [
        "filterbyname.sh -Xmx${avail_mem}g",
        input,
        "out=${first_filtered_reads}",
        "names=${names_to_filter} exclude",
        args
    ].findAll { it?.trim() }.join(' ')
    def effective_reformat = "reformat.sh in=${first_filtered_reads} out=${fully_filtered_reads} minlength=500"
    def note = 'Removes FCS-GX flagged scaffolds, records review and trim counts, and drops contigs shorter than 500 bp.'

    """
    # count the number of contigs and the number of base pairs being removed across EXCLUDE and TRIM 

    exclude_lines=\$(grep -w EXCLUDE "${action_report}" || true)
    if [[ -n "\$exclude_lines" ]]; then
        count=\$(echo "\$exclude_lines" | cut -f 1 | sort -u | wc -l)
        bp=\$(echo "\$exclude_lines" | awk '{sum+=\$3-\$2+1}END{print sum}')
    else
        count=0
        bp=0
    fi
    echo "EXCLUDE \$count \$bp" | tee -a $filter_report


    trim_lines=\$(grep -w TRIM "${action_report}" || true)
    if [[ -n "\$trim_lines" ]]; then
        count=\$(echo "\$trim_lines" | cut -f 1 | sort -u | wc -l)
        bp=\$(echo "\$trim_lines" | awk '{sum+=\$3-\$2+1}END{print sum}')
    else
        count=0
        bp=0
    fi
    echo "TRIM \$count \$bp" | tee -a $filter_report

    review_lines=\$(grep -w REVIEW "${action_report}" || true)
    if [[ -n "\$review_lines" ]]; then
        count=\$(echo "\$review_lines" | cut -f 1 | sort -u | wc -l)
        bp=\$(echo "\$review_lines" | awk '{sum+=\$3-\$2+1}END{print sum}')
    else
        count=0
        bp=0
    fi
    echo "REVIEW \$count \$bp" | tee -a $filter_report
  
    #generate a txt file with the name of the contigs that are in review that are less that 1000bp.
    if [[ -n "\$review_lines" ]]; then
        echo "\$review_lines" | awk '\$4 <= 1000 {print \$1}' > $names_to_filter
    else
        : > $names_to_filter
    fi

    # First pass of filtering - remove all EXCLUDE and TRIM contigs
    # REVIEW contigs are retained for manual review
    filterbyname.sh \\
        -Xmx${avail_mem}g \\
        $input \\
        out=$first_filtered_reads \\
        names=$names_to_filter exclude \\
        $args
     
    # Wait for the first bbmap script to complete before moving on
    wait


    grep -v '^>' "$first_filtered_reads" \\
        | awk 'length(\$0) < 500 {count++} END {print "Number of contigs less than 500bp:", count}' \\
        > $contigs_under_500bp


    #remove the contigs that are less than 500bp from the assembly 
    reformat.sh \\
        in="$first_filtered_reads" \\
        out="$fully_filtered_reads" \\
        minlength=500
    cat <<-END_TOOL_PARAMS > 25_bbmap_filterbyname.tool_params_mqcrow.html
    <tr><td>BBMap FilterByName</td><td><samp>${effective_args}; ${effective_reformat}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    input  = "in=${reads}"
    first_filtered_reads = "${prefix}.rf.fa"
    fully_filtered_reads = "${prefix}.${output_format}"
    filter_report = "${meta.prefix}.filter_report.txt"
    names_to_filter = "${meta.prefix}.review_scaffolds_1kb.txt"
    contigs_under_500bp = "${meta.prefix}.contig_count_500bp.txt"
    def avail_mem = task.memory ? task.memory.giga : 3
    def effective_args = [
        "filterbyname.sh -Xmx${avail_mem}g",
        input,
        "out=${first_filtered_reads}",
        "names=${names_to_filter} exclude",
        args
    ].findAll { it?.trim() }.join(' ')
    def effective_reformat = "reformat.sh in=${first_filtered_reads} out=${fully_filtered_reads} minlength=500"
    def note = 'Removes FCS-GX flagged scaffolds, records review and trim counts, and drops contigs shorter than 500 bp.'

    """
    touch $first_filtered_reads
    touch $filter_report
    touch $names_to_filter
    touch $contigs_under_500bp
    cat <<-END_TOOL_PARAMS > 25_bbmap_filterbyname.tool_params_mqcrow.html
    <tr><td>BBMap FilterByName</td><td><samp>${effective_args}; ${effective_reformat}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

}
