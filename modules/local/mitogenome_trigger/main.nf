process TRIGGER_MITOGENOME {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(assembly_files)
    
    output:
    tuple val(sample_id), path("trigger_${sample_id}.done"), emit: triggered
    
    when:
    params.enable_mitogenome_trigger
    
    script:
    """
    # Create sample-specific directories
    mkdir -p ${params.mitogenome_outdir}/inputs/${sample_id}
    mkdir -p ${params.mitogenome_outdir}/logs
    mkdir -p ${params.mitogenome_outdir}/${sample_id}
    
    # Copy/link relevant files for mitogenome assembly
    ln -s \$(readlink -f ${assembly_files}) ${params.mitogenome_outdir}/inputs/${sample_id}/
    
    # Create sample sheet for this specific sample
    echo "sample_id,assembly_path" > ${params.mitogenome_outdir}/inputs/${sample_id}/samplesheet.csv
    echo "${sample_id},${params.mitogenome_outdir}/inputs/${sample_id}/" >> ${params.mitogenome_outdir}/inputs/${sample_id}/samplesheet.csv
    
    # Create tmux session name
    SESSION_NAME="mitogenome_${sample_id}"
    
    # Kill existing session if it exists
    if tmux has-session -t "\$SESSION_NAME" 2>/dev/null; then
        tmux kill-session -t "\$SESSION_NAME"
    fi
    
    # Launch mitogenome pipeline in background tmux session (fire-and-forget)
    tmux new-session -d -s "\$SESSION_NAME" -c "\$(pwd)" \\
        "nextflow -log ${params.mitogenome_outdir}/logs/${sample_id}_nextflow.log \\
            run ./mitogenome_assembly.nf \\
            --input ${params.mitogenome_outdir}/inputs/${sample_id}/samplesheet.csv \\
            --outdir ${params.mitogenome_outdir}/${sample_id} \\
            --sample_id ${sample_id} \\
            -work-dir ${params.mitogenome_outdir}/${sample_id}/work \\
            -c ${params.mitogenome_config} \\
            -profile docker \\
            -resume 2>&1 | tee ${params.mitogenome_outdir}/logs/${sample_id}_mitogenome.log; \\
         echo 'Pipeline finished. Session will close in 10 seconds...'; \\
         sleep 10"
    
    # Wait a moment to ensure tmux session is created, then immediately continue
    sleep 2
    
    # Create completion marker (indicates trigger was successful, not pipeline completion)
    touch trigger_${sample_id}.done
    echo "Successfully triggered mitogenome assembly for ${sample_id}" > trigger_${sample_id}.done
    echo "tmux session: \$SESSION_NAME"
    echo "Nextflow log: ${params.mitogenome_outdir}/logs/${sample_id}_nextflow.log"
    echo "Pipeline log: ${params.mitogenome_outdir}/logs/${sample_id}_mitogenome.log"
    
    echo "Main workflow will continue without waiting for mitogenome completion"
    """
}