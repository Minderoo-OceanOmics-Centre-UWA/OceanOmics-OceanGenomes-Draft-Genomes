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
    // How to treat FCS-GX REVIEW contigs: exclude (default) drops all of them,
    // short-only keeps the historical "<= 1000 bp" rule, keep retains them all.
    def review_action = params.fcs_review_action ?: 'exclude'
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    input  = "in=${reads}"
    first_filtered_reads = "${prefix}.rf.fa"
    fully_filtered_reads = "${prefix}.${output_format}"
    filter_report = "${meta.prefix}.filter_report.txt"
    names_to_filter = "${meta.prefix}.review_scaffolds.txt"
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
    def note = "Removes FCS-GX flagged scaffolds (REVIEW action: ${review_action}), records review and trim counts, and drops contigs shorter than 500 bp."

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
  
    # Build the list of REVIEW contigs to remove.
    #
    # FCS-GX assigns REVIEW rather than EXCLUDE when its confidence is low. For taxa
    # the GX database covers poorly -- sponges, deep-sea crustaceans, echinoderms --
    # that is most of the contamination it finds, and the REVIEW sets come back
    # taxonomically indistinguishable from the EXCLUDE sets (near-entirely
    # prokaryotic). Nothing downstream actually reviews them, so retaining them just
    # publishes bacterial contigs. Default to dropping them all.
    case "${review_action}" in
        exclude)
            if [[ -n "\$review_lines" ]]; then
                echo "\$review_lines" | cut -f 1 | sort -u > $names_to_filter
            else
                : > $names_to_filter
            fi
            ;;
        short-only)
            # Historical behaviour, kept so earlier runs can be reproduced.
            if [[ -n "\$review_lines" ]]; then
                echo "\$review_lines" | awk '\$4 <= 1000 {print \$1}' | sort -u > $names_to_filter
            else
                : > $names_to_filter
            fi
            ;;
        keep)
            : > $names_to_filter
            ;;
        *)
            echo "ERROR: unknown fcs_review_action '${review_action}' (expected exclude, short-only or keep)" >&2
            exit 1
            ;;
    esac

    # EXCLUDE and TRIM were already applied upstream by FCSGX_CLEANGENOME, so this
    # pass removes only the REVIEW contigs selected above.
    filterbyname.sh \\
        -Xmx${avail_mem}g \\
        $input \\
        out=$first_filtered_reads \\
        names=$names_to_filter exclude \\
        $args
     
    # Wait for the first bbmap script to complete before moving on
    wait


    # Count whole contigs below 500 bp, not wrapped sequence lines: every line of a
    # wrapped FASTA is shorter than 500 characters, which made this a count of lines.
    awk '/^>/ {if (len > 0 && len < 500) count++; len = 0; next}
         {len += length(\$0)}
         END {if (len > 0 && len < 500) count++; print "Number of contigs less than 500bp:", count+0}' \\
        "$first_filtered_reads" \\
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
    names_to_filter = "${meta.prefix}.review_scaffolds.txt"
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
