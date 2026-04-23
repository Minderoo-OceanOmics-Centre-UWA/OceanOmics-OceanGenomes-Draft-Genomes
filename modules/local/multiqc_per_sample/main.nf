process MULTIQC_PER_SAMPLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.27--pyhdfd78af_0' :
        'biocontainers/multiqc:1.27--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*")
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    tuple val(meta), path("*multiqc_report.html"), emit: report
    tuple val(meta), path("*_data")              , emit: data
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    while IFS= read -r -d '' versions_file; do
        cp -L "\$versions_file" "\${versions_file%.yml}_mqc_versions.yml"
    done < <(find . \\( -type f -o -type l \\) -name 'software_versions.yml' -print0)

    while IFS= read -r -d '' versions_file; do
        cp -L "\$versions_file" "\${versions_file%.yaml}_mqc_versions.yaml"
    done < <(find . \\( -type f -o -type l \\) -name 'software_versions.yaml' -print0)

    multiqc \\
        --force \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
