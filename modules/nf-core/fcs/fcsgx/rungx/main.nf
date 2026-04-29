process FCSGX_RUNGX {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.5--h9948957_0':
        'biocontainers/ncbi-fcs-gx:0.5.5--h9948957_0' }"

    input:
    tuple val(meta), path(fasta)
    path gxdb
    val ramdisk_path

    output:
    tuple val(meta), path("*.fcs_gx_report.txt"), emit: fcsgx_report
    tuple val(meta), path("*.taxonomy.rpt")     , emit: taxonomy_report
    tuple val(meta), path("*.summary.txt")      , emit: log
    // tuple val(meta), path("*.hits.tsv.gz")      , emit: hits, optional: true
    tuple val(meta), path("23_fcsgx_rungx.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--debug '
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    // def mv_database_to_ram = ramdisk_path ? "rclone copy $gxdb $ramdisk_path" : '' // have just added in custom pawsey code that works
    def database = ramdisk_path ? "$ramdisk_path" : gxdb // Use task.index to make memory location unique
    def effective_args = ["GX_NUM_CORES=${task.cpus}; run_gx.py", "--fasta ${fasta}", "--gx-db ${database}", "--tax-id ${meta.taxon_id}", "--generate-logfile true", "--out-basename ${prefix}", "--out-dir .", args].findAll { it?.trim() }.join(' ')
    def note = 'Copies the GX database to the ramdisk when configured, then screens the assembly for contamination.'
    """
    # Copy DB to RAM-disk when supplied. Otherwise, the tool is very slow.
    mkdir $ramdisk_path
    cp -v $gxdb/gxdb/all.gxi $ramdisk_path
                cp -v $gxdb/gxdb/all.gxs $ramdisk_path
                cp -v $gxdb/gxdb/all.meta.jsonl $ramdisk_path
                cp -v $gxdb/gxdb/all.blast_div.tsv.gz $ramdisk_path 
                cp -v $gxdb/gxdb/all.taxa.tsv $ramdisk_path 
            echo ‘done copying’ 
            ls -l $ramdisk_path

    export GX_NUM_CORES=${task.cpus}
    run_gx.py \\
        --fasta ${fasta} \\
        --gx-db ${database} \\
        --tax-id ${meta.taxon_id} \\
        --generate-logfile true \\
        --out-basename ${prefix} \\
        --out-dir . \\
        ${args}
    cat <<-END_TOOL_PARAMS > 23_fcsgx_rungx.tool_params_mqcrow.html
    <tr><td>FCSGX Rungx</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """

    stub:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: '--debug '
    def database = ramdisk_path ? "$ramdisk_path" : gxdb
    def effective_args = ["GX_NUM_CORES=${task.cpus}; run_gx.py", "--fasta ${fasta}", "--gx-db ${database}", "--tax-id ${meta.taxon_id}", "--generate-logfile true", "--out-basename ${prefix}", "--out-dir .", args].findAll { it?.trim() }.join(' ')
    def note = 'Copies the GX database to the ramdisk when configured, then screens the assembly for contamination.'
    """
    touch ${prefix}.fcs_gx_report.txt
    touch ${prefix}.taxonomy.rpt
    touch ${prefix}.summary.txt
    echo "" | gzip > ${prefix}.hits.tsv.gz
    cat <<-END_TOOL_PARAMS > 23_fcsgx_rungx.tool_params_mqcrow.html
    <tr><td>FCSGX Rungx</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcsgx: \$( gx --help | sed '/build/!d; s/.*:v//; s/-.*//' )
    END_VERSIONS
    """
}
