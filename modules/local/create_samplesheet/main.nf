process CREATE_SAMPLESHEET {
    tag "$meta"
    label 'process_low'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"
    conda "${moduleDir}/environment.yml"

    input:
    val reads_list
    path config
    path taxdump

    output:
    path("${params.run}_samplesheet.csv"), emit: samplesheet
    path "taxonomy_resolution.tsv"        , emit: taxonomy_resolution
    path "versions.yml"                   , emit: versions

    when:
    !params.input
    
    script:
    // Turn the Groovy list into a valid Python literal string
    // e.g. ['OG1323', ['/path/R1', '/path/R2'], 'OG1336', ...]
    def reads_literal = reads_list.inspect()
    // The curated species table only holds taxa someone has loaded, and the
    // invertebrate runs draw from most of Metazoa. Without the taxdump fallback
    // every uncurated sample aborts the run at this step.
    def taxdump_arg = taxdump ? "--taxdump-dir ${taxdump}" : ''

    """
    create_samplesheet.py \\
        $config \\
        $params.run \\
        $params.outdir/pooled/$params.run \\
        ${taxdump_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control this species validation script.
    END_VERSIONS
    """
    }
