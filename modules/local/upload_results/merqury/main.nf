process PUSH_MERQURY_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(completeness_stats), path(qv_tsv)
    path config

    output:
    path "${meta.id}.lca_blast.upload.txt", emit: upload
    path "versions.yml"                   , emit: versions

    script:
    """
    # Push the results to SQL database
    # QV-only
    push_merqury_results_to_sqldb.py -c $config --qv ${qv_tsv}

    # Completeness-only
    push_merqury_results_to_sqldb.py -c $config --comp ${completeness_stats}

    # Both
    push_merqury_results_to_sqldb.py -c $config \\
        --qv ${qv_tsv} \\
        --comp ${completeness_stats}
        ${meta.id}.merqury_results.txt \\
        > ${meta.id}.merqury.upload.txt

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


    
