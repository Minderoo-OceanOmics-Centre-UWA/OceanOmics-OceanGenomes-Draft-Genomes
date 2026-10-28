// Upload read depth and the host partition from *_kmer_depth.json, and nothing else.
//
// Used by the depth-only entry (workflows/depth_only.nf). push_kmer_depth_to_sqldb.py writes
// DEPTH_COLUMNS only, so the GenomeScope numbers and every other family already in
// draft_genomes stay as they are, and it refuses a draft that has no row.
process PUSH_KMER_DEPTH_RESULTS {
    tag "$meta.id"
    label 'process_upload'

    conda "${moduleDir}/../assembly/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(kmer_depth)
    path config

    output:
    path "${meta.id}.kmer_depth.upload.txt", emit: upload
    path "versions.yml"                    , emit: versions

    script:
    """
    push_kmer_depth_to_sqldb.py \\
        -c $config \\
        -k ${kmer_depth} \\
        > ${meta.id}.kmer_depth.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upload: "push_kmer_depth_to_sqldb.py"
    END_VERSIONS
    """

    stub:
    """
    : > ${meta.id}.kmer_depth.upload.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upload: "stub"
    END_VERSIONS
    """
}
