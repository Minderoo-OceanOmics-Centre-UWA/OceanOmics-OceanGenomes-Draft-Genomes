// Turn BUSCO-region read depth into a haploid k-mer coverage, and compare it to the one
// GenomeScope fitted.
//
// Split from MEASURE_KMER_COVERAGE because that half needs samtools and this half needs
// python, and no single container has both. All the arithmetic and file parsing lives
// here, in bin/kmer_depth_lambda.py, where it is unit tested -- the shell substitutes
// that a one-container version would have forced are exactly where the bugs were.
//
// Emits a row for every sample, including those with no measurement. A sample whose
// assembly is too fragmented to yield 20 Complete BUSCOs gets lambda_depth null and a
// status saying so, rather than a failed task that would remove it from the downstream
// joins and from the database.
process SUMMARISE_KMER_COVERAGE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(bedcov), path(status), path(contig_depth), path(fastp_json), path(seqkit_stats), path(genomescope_model)

    output:
    tuple val(meta), path("*_kmer_depth.json"), emit: kmer_depth
    tuple val(meta), path("33_summarise_kmer_coverage.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    STATUS=\$(cat ${status})
    BEDCOV_ARG=""
    if [ "\$STATUS" = ok ]; then
        BEDCOV_ARG="--bedcov ${bedcov}"
    fi

    # An empty contig-depth file is the too_few_busco_regions branch, not a failure: pass
    # it only when it has content, so the partition reports no_contig_depth rather than
    # parsing nothing and calling it a zero-sized host.
    DEPTH_ARG=""
    if [ -s "${contig_depth}" ]; then
        DEPTH_ARG="--contig-depth ${contig_depth}"
    fi

    kmer_depth_lambda.py \\
        \$BEDCOV_ARG \\
        \$DEPTH_ARG \\
        --host-depth-window ${params.host_depth_window} \\
        --host-depth-min-contig ${params.host_depth_min_contig} \\
        --status "\$STATUS" \\
        --fastp-json ${fastp_json} \\
        --seqkit-stats ${seqkit_stats} \\
        --model ${genomescope_model} \\
        --kmer ${params.kvalue ?: 21} \\
        --sample ${meta.id} > ${prefix}_kmer_depth.json

    cat <<-END_TOOL_PARAMS > 33_summarise_kmer_coverage.tool_params_mqcrow.html
    <tr><td>Summarise K-mer Coverage</td><td><samp>kmer_depth_lambda.py --kmer ${params.kvalue ?: 21} --host-depth-window ${params.host_depth_window}</samp></td><td>Haploid k-mer coverage from read depth, the ratio of host single-copy depth to assembly-wide depth, and a lower bound on host genome size from per-contig depth partitioning.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    """
    echo '{"sample_id": "${meta.id}", "lambda_depth": 30.0, "genomescope_kmercov": 29.7, "host_assembly_size": 345700000, "host_assembly_fraction": 0.92, "partition_status": "ok"}' > ${prefix}_kmer_depth.json
    touch 33_summarise_kmer_coverage.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.9
    END_VERSIONS
    """
}
