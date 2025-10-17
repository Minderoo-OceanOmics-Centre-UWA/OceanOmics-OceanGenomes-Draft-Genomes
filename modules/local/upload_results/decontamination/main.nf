process PUSH_DECONTAMINATION_RESULTS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://tylerpeirce/psycopg2:0.1' :
        'tylerpeirce/psycopg2:0.1' }"

    input:
    tuple val(meta), path(filter_report), path(contigs_under_500bp), path(tiara_filter_summary)
    path config

    output:
    path "${meta.id}.lca_blast.upload.txt", emit: upload
    path "versions.yml"                   , emit: versions

    script:
    """
    # Push the results to SQL database
    # Tiara-only upsert for one sample
    # python push_decontamination_results_to_sqldb.py -c ../configfile.txt --tiara ${tiara_sample_tsv}

    # NCBI-only upsert for one sample
    # python push_decontamination_results_to_sqldb.py -c ../configfile.txt --ncbi ${ncbi_sample_tsv}

    # Do both in one call
    python push_decontamination_results_to_sqldb.py -c $config \\
        --tiara ${tiara_filter_summary} \\
        --filter ${filter_report} \\
        --contigs ${contigs_under_500bp} \\
        ${meta.id}.decontamination_results.txt \
        > ${meta.id}.decontamination.upload.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version 1 - need to version control these upload scripts.
    END_VERSIONS
    """
    
    stub:
    """
    : > ${meta.id}.lca_blast.upload.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        upload: "stub"
    END_VERSIONS
    """
    }


    
