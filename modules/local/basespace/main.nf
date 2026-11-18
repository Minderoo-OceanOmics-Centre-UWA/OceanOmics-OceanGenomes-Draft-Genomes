process BASESPACE {
    tag "$run_id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/csgenetics/basespace-cli:0.1' :
        'quay.io/csgenetics/basespace-cli@sha256:514422008a5e32ba3795c2d61992dd3c9aaf589156440cd378e5a3157d82eb7d' }"

    input:
    val run_id
    path config

    output:
    path "*/*fastq.gz", emit: fastqs
    path "*/*json"    , emit: jsons
    path "versions.yml", emit: versions

    script:
    // Lanes download into a staging dir outside the work dir so a failed attempt,
    // a walltime kill or a later -resume carries on from what is already there:
    // `bs download` skips files that exist at the same size and re-fetches partial
    // ones. Staging sits on /scratch like the work dir, so the final mv is a rename.
    def stage = params.bs_staging_dir ?: "${params.outdir}/basespace_staging/${run_id}"
    """
    cp /bin/bs .
    STAGE="${stage}"
    mkdir -p "\$STAGE"

    # BaseSpace API calls fail transiently (e.g. "TLS handshake timeout", exit 3).
    # --retry gives bs's own 5 retries; on top of that retry the whole call with a
    # growing wait so a short outage doesn't fail a multi-hour download.
    bs_retry() {
        local attempt=1 wait=60
        until ./bs --retry "\$@"; do
            if [ \$attempt -ge 5 ]; then
                echo "bs \$* failed after \$attempt attempts" >&2
                return 1
            fi
            echo "bs \$* failed (attempt \$attempt), retrying in \$wait s" >&2
            sleep \$wait
            attempt=\$((attempt + 1)); wait=\$((wait * 2))
        done
    }

    RUNID=\$(bs_retry list run | { grep "$run_id " || true; } | awk '{print \$4}')
    if [ -z "\$RUNID" ]; then
        echo "Run $run_id not found in BaseSpace" >&2
        exit 1
    fi

    # List the run's datasets once (Name,Id) and look IDs up from this file, instead
    # of one API call per lane.
    bs_retry list dataset --input-run "\$RUNID" -f csv -F Name -F Id \\
        | tail -n +2 | grep '^OG' > ${run_id}.datasets.csv
    cut -d, -f1 ${run_id}.datasets.csv > ${run_id}.prefix.txt
    if [ ! -s ${run_id}.prefix.txt ]; then
        echo "No OG datasets found for run $run_id (\$RUNID)" >&2
        exit 1
    fi

    while IFS=, read -r PREFIX ID; do
        echo "\$PREFIX \$ID >> $run_id"
        bs_retry download dataset -i "\$ID" -o "\$STAGE/\$PREFIX" --no-progress-bars
    done < ${run_id}.datasets.csv

    # Only hand the lanes over once every one has R1, R2 and its json; otherwise keep
    # staging so the next attempt resumes.
    missing=0
    while read -r PREFIX; do
        n_r1=\$(ls "\$STAGE/\$PREFIX"/*_R1_*fastq.gz 2>/dev/null | wc -l)
        n_r2=\$(ls "\$STAGE/\$PREFIX"/*_R2_*fastq.gz 2>/dev/null | wc -l)
        n_js=\$(ls "\$STAGE/\$PREFIX"/*json 2>/dev/null | wc -l)
        if [ \$n_r1 -lt 1 ] || [ \$n_r2 -lt 1 ] || [ \$n_js -lt 1 ]; then
            echo "Incomplete lane \$PREFIX (R1=\$n_r1 R2=\$n_r2 json=\$n_js)" >&2
            missing=\$((missing + 1))
        fi
    done < ${run_id}.prefix.txt
    if [ \$missing -gt 0 ]; then
        echo "\$missing lane(s) incomplete; staging kept at \$STAGE for the retry" >&2
        exit 1
    fi

    while read -r PREFIX; do
        mv "\$STAGE/\$PREFIX" .
    done < ${run_id}.prefix.txt
    rmdir "\$STAGE" 2>/dev/null || echo "Staging dir \$STAGE not empty, left in place" >&2

    bs_version=\$(./bs --version 2>/dev/null | head -n 1 | awk '{print \$2}')
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bs: "\${bs_version:-unknown}"
    END_VERSIONS
    """
}
