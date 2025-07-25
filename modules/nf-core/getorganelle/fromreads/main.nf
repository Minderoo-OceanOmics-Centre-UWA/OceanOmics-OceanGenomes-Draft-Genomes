process GETORGANELLE_FROMREADS {
    tag "$meta.id"
    label 'process_high'

    publishDir = [
        path: { "${params.outdir}/mitogenomes/${meta.id}/${meta.id}.ilmn.${meta.date}.getorg${version.replaceAll('\\.', '')}" },
        mode: params.publish_dir_mode
    ]

    conda "bioconda::getorganelle=1.7.7.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/getorganelle:1.7.7.0--pyh7cba7a3_0':
        'biocontainers/getorganelle:1.7.7.0--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fastq), val(organelle_type), path(db), val(version)  // getOrganelle has a database and config file

    output:
    tuple val(meta), path("mtdna/${meta.id}.ilmn.${meta.date}.getorg${version.replaceAll('\\.', '')}.fasta"), emit: fasta
    tuple val(meta), path("mtdna/*get_org.log.txt")     , emit: log
    path("mtdna/*")                                     , emit: etc // the rest of the result files
    path "versions.yml"                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.ilmn.${meta.date}.getorg${version.replaceAll('\\.', '')}"

    """
    get_organelle_from_reads.py \\
        $args \\
        --prefix ${prefix}. \\
        -F $organelle_type \\
        --config-dir $db \\
        -t $task.cpus \\
        -1 ${fastq[0]} \\
        -2 ${fastq[1]} \\
        -o mtdna

    wait
    
    # Move and rename the output file to include the version in the filename
    mv mtdna/${prefix}.*1.1.*.fasta mtdna/${prefix}.fasta
    sed -i "/^>/s/.*/>${prefix}/g" mtdna/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_from_reads.py --version | sed 's/^GetOrganelle v//g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def base_prefix = task.ext.prefix ?: "${meta.id}.ilmn.${meta.date}.getorg${version.replaceAll('\\.', '')}"
    """
    mkdir -p mtdna
       
    touch "mtdna/${prefix}.fasta"
    touch "mtdna/stub_output.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getorganelle: \$(get_organelle_config.py --version | sed 's/^GetOrganelle v//g')
    END_VERSIONS
    """
}