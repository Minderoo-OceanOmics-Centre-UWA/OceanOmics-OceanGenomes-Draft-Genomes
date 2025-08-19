process PUSH_MTDNA_ASSM_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(fasta), path(out_log)
    path config

    output:
    path "${meta.id}.mtdna.upload.txt", emit: upload
    path "versions.yml"         , emit: versions

    script:
    """
    push_mtdna_assm_results.py \\
        $config \\
        ${meta.mt_assembly_prefix} \\
        $out_log \\
        $fasta \\
        > ${meta.id}.mtdna.upload.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    }