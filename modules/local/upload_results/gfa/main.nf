process PUSH_GFA_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(gfastats_results)
    path config

    output:
    path "${meta.id}.gfastats.upload.txt", emit: upload
    path "versions.yml"                   , emit: versions

    script:
    """
    # Push the results to SQL database
    push_gfa_results_to_sqldb.py \\
        -c $config \\
        -f ${gfastats_results} \\
        > ${meta.id}.gfastats.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    
    stub:
    """
    : > ${meta.id}.lca_blast.upload.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upload: "stub"
    END_VERSIONS
    """
    }


    
