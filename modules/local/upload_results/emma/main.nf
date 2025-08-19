process PUSH_MTDNA_ANNOTATION_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(annotations) 
    path config

    output:
    path "${meta.id}.annotation.upload.txt", emit: upload
    path "${meta.id}.annotation_stats.csv", emit: stats
    path "versions.yml"                   , emit: versions

    script:
    """
    # Compile the statistics 
    annotation_stats.py \\
        *.gff \\
        proteins
    
    wait

    # Push the results to SQL database
    push_emma_annotation_results.py \\
        $config \\
        ${meta.id} \\
        ${meta.id}.annotation_stats.csv \\
        > ${meta.id}.annotation.upload.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    }


    