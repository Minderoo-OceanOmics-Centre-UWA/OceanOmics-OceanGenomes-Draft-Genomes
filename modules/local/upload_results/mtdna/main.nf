process PUSH_MTDNA_ASSM_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(log), val(assembly_prefix) 
    path config

    output:
    path "${meta.id}.upload.txt", emit: upload

    script:
    """
    push_mtdna_assm_results.py \\
        $config \\
        $assembly_prefix \\
        $log \\
        $fasta\\
        > ${meta.id}.mtdna.upload.txt

    """
    }