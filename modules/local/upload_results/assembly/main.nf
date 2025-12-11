process PUSH_ASSEMBLY_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(genomescope_summary)
    path config

    output:
    path "${meta.id}.assembly.upload.txt", emit: upload
    path "versions.yml"                   , emit: versions

    script:
    """
    # Push the results to SQL database
    push_assembly_results_to_sqldb.py \\
        -c $config \\
        -f ${genomescope_summary} \\
        > ${meta.id}.assembly.upload.txt 

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


    
