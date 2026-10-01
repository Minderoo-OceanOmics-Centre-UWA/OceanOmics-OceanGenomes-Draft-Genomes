// Local, and deliberately not an nf-core module any more.
//
// It has diverged: the genome size comes from a GenomeScope summary staged as a second
// input, the prefix is meta.assembly_prefix, and it writes a tool_params row for MultiQC.
// `nf-core modules update` would revert all of that silently, the same trap the CHANGELOG
// records for megahit, so it lives here and is out of modules.json.
process GFASTATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.6--hdcf5f25_3':
        'biocontainers/gfastats:1.3.6--hdcf5f25_3' }"

    input:
    tuple val(meta), path(assembly), path(summary)
    val out_fmt
 
    output:
    tuple val(meta), path("*.assembly_summary"), emit: assembly_summary
    tuple val(meta), path("*.${out_fmt}")   , emit: assembly        , optional: true
    tuple val(meta), path("70_gfastats.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    def effective_args = [args, "--threads ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def note = "Genome size is parsed from ${summary.name}; writes ${prefix}.assembly_summary."
    // def output_sequences = out_fmt ? "--out-format ${prefix}.${out_fmt}" : ""
    """
    # Get genomesize from $summary file
    genome_size=\$(cat $summary | grep 'Genome Haploid Length' | grep -o 'bp.*' | sed 's/bp//g' | sed 's/ //g' | sed 's/,//g')

    # A collapsed GenomeScope fit writes -1 for every property (OG2653: Model Fit 0%), and
    # gfastats takes the expected size as a positional -- it aborts on a non-positive one
    # with "basic_string: construction from null", exit 134, three attempts, then ignored.
    #
    # When the size is not a positive integer, omit the positional entirely. gfastats then
    # reports N50 and skips the NG statistics, which is honest. Substituting the assembly
    # size would publish an NG50 that is really an N50.
    NOTE_SUFFIX=""
    case "\$genome_size" in
        ''|*[!0-9]*|0)
            echo "GFASTATS: genome size '\$genome_size' is not a positive integer; omitting the expected-size argument, so NG statistics will be absent" >&2
            genome_size=""
            NOTE_SUFFIX=" GenomeScope reported no usable genome size, so the expected-size argument was omitted and the NG statistics are absent."
            ;;
    esac

    gfastats \\
        $args \\
        --threads $task.cpus \\
        $assembly \\
        \$genome_size \\
        > ${prefix}.assembly_summary

    cat <<-END_TOOL_PARAMS > 70_gfastats.tool_params_mqcrow.html
    <tr><td>Gfastats</td><td><samp>${effective_args}</samp></td><td>${note}\${NOTE_SUFFIX}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def effective_args = [args, "--threads ${task.cpus}"].findAll { it?.trim() }.join(' ')
    def note = "Genome size is parsed from ${summary.name}; writes ${prefix}.assembly_summary."
    """
    echo | gzip > ${prefix}.${out_fmt}.gz
    touch ${prefix}.assembly_summary
    cat <<-END_TOOL_PARAMS > 70_gfastats.tool_params_mqcrow.html
    <tr><td>Gfastats</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
    END_VERSIONS
    """
}
