process BBMAP_REPAIR {
    tag "$ogid"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(ogid), path(reads)
    val(interleave)

    output:
    tuple val(ogid), path("*.{R1,R2}.fq.gz")         , emit: repaired
    path  "versions.yml"                                 , emit: versions
    path  "*paicheck.log"                                        , emit: log

    when:
    task.ext.when ?: true

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${ogid}.ilmn.${params.run}"
    in_reads  = ( interleave )  ?: "in=cat/${prefix}.cat.R1.fq.gz in2=cat/${prefix}.cat.R2.fq.gz"
    out_reads = ( interleave )  ?: "out=${prefix}.R1.fq.gz out2=${prefix}.R2.fq.gz"
      
 
    """
    mkdir -p cat
    cat ${ogid}*R1*.gz > cat/${prefix}.cat.R1.fq.gz && echo "cat R1 completed"
    cat ${ogid}*R2*.gz > cat/${prefix}.cat.R2.fq.gz && echo "cat R2 completed" 

    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    repair.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        threads=${task.cpus} \\
        ${args} 
    
    cp .command.log ${prefix}.paicheck.log
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}