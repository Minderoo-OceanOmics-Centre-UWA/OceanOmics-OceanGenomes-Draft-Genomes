process SPECIES_VALIDATION {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(blast_files), path(lca_files)
    path config

    output:
    path "lca_results.${meta.id}.tsv", emit: summary
    tuple val(meta), path("lca_combined.${meta.id}.tsv"), path("blast_combined.${meta.id}.tsv"), emit: full
    path "versions.yml"                   , emit: versions

    script:
    """
    species_validation.py \\
        $config \\
        ${meta.id} \\
        "${lca_files.join(',')}" \\
        "${blast_files.join(',')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control this species validation script.
    END_VERSIONS
    """
    }