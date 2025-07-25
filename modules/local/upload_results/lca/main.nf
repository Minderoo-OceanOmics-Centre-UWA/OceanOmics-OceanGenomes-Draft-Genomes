process PUSH_LCA_BLAST_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple path (lca_results), path (blast_results)
    path config

    output:
    path "${meta.id}.lca_blast.upload.txt", emit: upload

    script:
    """
    # Push the results to SQL database
    push_emma_annotation_results.py \\
        $config \\
        ${meta.id} \\
        ${lca_results} \\
        ${blast_results} \\
        > ${meta.id}.lca_blast.upload.txt

    """
    }


    