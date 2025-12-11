process CREATE_SAMPLESHEET {
    tag "$meta"
    label 'process_low'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    val reads_list
    path config

    output:
    path("${params.run}_samplesheet.csv"), emit: samplesheet

    when:
    !params.input
    
    script:
    // Turn the Groovy list into a valid Python literal string
    // e.g. ['OG1323', ['/path/R1', '/path/R2'], 'OG1336', ...]
    def reads_literal = reads_list.inspect()

    """
    create_samplesheet.py \\
        $config \\
        $params.run \\
        $params.outdir/pooled/$params.run

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control this species validation script.
    END_VERSIONS
    """
    }