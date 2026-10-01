// Refit GenomeScope seeded with the measured k-mer coverage, where the first fit failed.
//
// -l is an initial kmercov estimate for the optimiser, not a constraint, so it only helps
// where the default start converged somewhere wrong. A sweep over the NOVA_260909_LA
// histograms set the rules this module follows:
//
//   * Seeding a CORRECT lambda is free. All eight samples with independently trustworthy
//     fits returned bit-identical results when seeded with their measured lambda.
//   * Seeding the raw histogram mode is harmful. The mode sits at 2 x lambda, and seeding
//     it cost those same eight samples 7-20 points of fit. The seed here is therefore the
//     measured lambda, never a peak.
//   * Below ~5x measured lambda there is nothing to fit. Seeds of 3-8 on samples with no
//     k-mer peak produced fits of 5-8% with non-finite genome lengths, and half-mode seeds
//     produced NEGATIVE fit values (-344.9, -431.4). Those samples are skipped and flagged
//     host_coverage_too_low instead.
//   * -l takes an integer only ("invalid int value: '28.50'"), so the seed is rounded.
//
// Where it works it is worth a lot: OG2617 went from a 83.5% fit at kmercov 17.9 claiming
// a 363 Mb genome, to a 96.7% fit at kmercov 8.3 -- against a measured 9.2 -- claiming
// 945 Mb, on an 829 Mb assembly.
process GENOMESCOPE_RESEED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py311r42hdfd78af_6':
        'biocontainers/genomescope2:2.0--py311r42hdfd78af_6' }"

    input:
    tuple val(meta), path(histogram), path(coverage_json), path(kmer_depth_json), path(summary), path(model)

    output:
    tuple val(meta), path("${prefix}_reseed_decision.json"), emit: decision
    tuple val(meta), path("${prefix}_genomescope2_reseeded_summary.txt"), emit: summary, optional: true
    tuple val(meta), path("${prefix}_genomescope2_reseeded_model.txt"), emit: model, optional: true
    tuple val(meta), path("33_genomescope_reseed.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.prefix}"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Decide whether to attempt a reseed at all, and with what.
    SEED=\$(genomescope_reseed_decide.py plan \\
        --coverage-json ${coverage_json} \\
        --kmer-depth-json ${kmer_depth_json} \\
        --min-seedable-lambda ${params.kmer_depth_min_seedable} \\
        --max-lambda-ratio ${params.kmer_depth_max_ratio} \\
        --enabled ${params.genomescope_reseed_when_flagged})

    if [ "\$SEED" != "skip" ]; then
        genomescope2 --input ${histogram} \\
            -k ${params.kvalue ?: 21} \\
            -m ${params.genomescope2_m} \\
            -l "\$SEED" \\
            --output reseed \\
            --name_prefix reseed || true
    fi

    genomescope_reseed_decide.py choose \\
        --coverage-json ${coverage_json} \\
        --kmer-depth-json ${kmer_depth_json} \\
        --original-summary ${summary} \\
        --original-model ${model} \\
        --reseeded-summary reseed/reseed_summary.txt \\
        --reseeded-model reseed/reseed_model.txt \\
        --seed "\$SEED" \\
        --min-reseed-fit-gain ${params.genomescope_min_reseed_fit_gain} \\
        --out ${prefix}_reseed_decision.json

    # Publish the reseeded fit only when it won, so there is never a better-looking file
    # on disk than the one the pipeline actually used.
    if [ "\$(python3 -c "import json;print(json.load(open('${prefix}_reseed_decision.json'))['winner'])")" = "reseeded" ]; then
        cp reseed/reseed_summary.txt ${prefix}_genomescope2_reseeded_summary.txt
        cp reseed/reseed_model.txt ${prefix}_genomescope2_reseeded_model.txt
    fi

    cat <<-END_TOOL_PARAMS > 33_genomescope_reseed.tool_params_mqcrow.html
    <tr><td>GenomeScope Reseed</td><td><samp>genomescope2 -l &lt;measured lambda&gt;, when the first fit was flagged and lambda >= ${params.kmer_depth_min_seedable}</samp></td><td>Refits with the k-mer coverage measured from read depth as the optimiser's starting point, and keeps the result only if the fit improves and kmercov moves toward that measurement.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.prefix}"
    """
    echo '{"winner": "original", "reason": "stub"}' > ${prefix}_reseed_decision.json
    touch 33_genomescope_reseed.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomescope2: 2.0
    END_VERSIONS
    """
}
