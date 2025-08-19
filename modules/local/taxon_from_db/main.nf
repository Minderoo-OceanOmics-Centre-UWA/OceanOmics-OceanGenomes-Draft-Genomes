process TAXON {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(repaired_files)
    path config

    output:
    tuple val(meta), path(repaired_files), path("${meta.id}_taxon_meta.csv")

    script:
    """
    taxon_id.py \\
        $config \\
        ${meta.id} > ${meta.id}_taxon_meta.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control this species validation script.
    END_VERSIONS
    """
    }