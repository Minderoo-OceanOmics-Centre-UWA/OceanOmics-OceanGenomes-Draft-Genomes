// Adapted from nf-core/megahit. It has diverged enough that it no longer tracks
// upstream: it publishes a single renamed contigs file, emits a tool_params MultiQC
// row, and keeps a resumable checkpoint directory outside the work directory so a
// failed assembly does not have to start from zero.
//
// That checkpoint directory is the reason this module needs care. It lives under
// params.outdir, so Nextflow's caching cannot see it, and a checkpoint keyed on the
// sample prefix alone gets reused even when the inputs that produced it have changed.
// It is therefore fingerprinted: see bin/megahit_checkpoint_key.sh.
process MEGAHIT {
    tag "${meta.id}"
    label 'process_extra_high'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f2cb827988dca7067ff8096c37cb20bc841c878013da52ad47a50865d54efe83/data' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.v129mh.fasta")                     , emit: contigs
    tuple val(meta), path("*.megahit_checkpoint.txt")           , emit: checkpoint_status
    tuple val(meta), path("22_megahit.tool_params_mqcrow.html") , emit: tool_params
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def memory = task.memory.toBytes()
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def reads_command = meta.single_end || !reads[1] ? "-r ${reads[0].join(',')}" : "-1 ${reads[0].join(',')} -2 ${reads[1].join(',')}"
    def read_files = reads.flatten().join(' ')
    def checkpoint_base = "${params.outdir}/megahit_checkpoints"
    def output_dir = "${checkpoint_base}/${prefix}_megahit_out"
    def unkeyed_policy = params.megahit_checkpoint_unkeyed ?: 'invalidate'
    def stale_policy = params.megahit_stale_checkpoint ?: 'rerun'
    def cleanup = params.megahit_checkpoint_cleanup == null ? true : params.megahit_checkpoint_cleanup
    def effective_args = [args, "-m ${memory}", "-t ${task.cpus}", reads_command, "--out-prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    """
    set -o pipefail

    # Validate the policies before doing any work. The schema flags a bad value as a
    # warning rather than an error, and a policy that is only checked when a stale
    # checkpoint happens to turn up is a typo that surfaces days later.
    case "${unkeyed_policy}" in
        invalidate|adopt) ;;
        *) echo "ERROR: unknown megahit_checkpoint_unkeyed '${unkeyed_policy}' (expected invalidate or adopt)" >&2; exit 1 ;;
    esac
    case "${stale_policy}" in
        rerun|archive|fail) ;;
        *) echo "ERROR: unknown megahit_stale_checkpoint '${stale_policy}' (expected rerun, archive or fail)" >&2; exit 1 ;;
    esac

    megahit_version=\$(megahit -v 2>&1 | sed 's/MEGAHIT v//')

    # Fingerprint this assembly's inputs. A checkpoint is reusable only if it was
    # started from the same reads, the same assembly args and the same megahit.
    #
    # The key lives beside the checkpoint directory rather than inside it, because
    # megahit refuses to write into an -o that already exists. Keeping it outside lets
    # the key be written BEFORE the assembly starts, which is what makes a half-built
    # checkpoint identifiable as a checkpoint of these inputs and therefore resumable.
    key_file="${output_dir}.key"
    key=\$(megahit_checkpoint_key.sh \\
        --version "\$megahit_version" \\
        --args '${args}' \\
        --prefix '${prefix}' \\
        ${read_files})
    stored_key=""
    if [ -f "\$key_file" ]; then
        stored_key=\$(cat "\$key_file")
    fi

    # Decide before acting. The key comparison has to come BEFORE the completion
    # check: testing for "done" first is what let a stale assembly be republished.
    if [ ! -d "${output_dir}" ]; then
        decision="fresh"
    elif [ -z "\$stored_key" ]; then
        # A checkpoint predating fingerprinting. Its inputs are unknowable.
        case "${unkeyed_policy}" in
            invalidate) decision="invalidated-unkeyed" ;;
            adopt)      decision="adopt-unkeyed" ;;
        esac
    elif [ "\$stored_key" != "\$key" ]; then
        case "${stale_policy}" in
            rerun|archive) decision="invalidated-stale" ;;
            fail)
                echo "ERROR: the MEGAHIT checkpoint at ${output_dir} was built from different inputs." >&2
                echo "  stored key:  \$stored_key" >&2
                echo "  current key: \$key" >&2
                echo "  megahit_stale_checkpoint = 'fail', so refusing to either reuse or discard it." >&2
                echo "  Inspect the inputs with: megahit_checkpoint_key.sh --manifest ..." >&2
                exit 1
                ;;
        esac
    elif [ -f "${output_dir}/done" ] && [ -s "${output_dir}/${prefix}.contigs.fa" ]; then
        # Both, not either: a bare "done" with a missing or truncated contigs file used
        # to surface several lines later as an unexplained cp failure.
        decision="skip"
    elif [ -f "${output_dir}/checkpoints.txt" ]; then
        decision="continue"
    else
        decision="fresh"
    fi

    # An operator asserting that the inputs have not changed takes ownership of the
    # checkpoint, after which it is an ordinary keyed checkpoint.
    if [ "\$decision" = "adopt-unkeyed" ]; then
        echo "Adopting unkeyed MEGAHIT checkpoint at ${output_dir} (megahit_checkpoint_unkeyed = adopt)."
        echo "\$key" > "\$key_file"
        if [ -f "${output_dir}/done" ] && [ -s "${output_dir}/${prefix}.contigs.fa" ]; then
            decision="adopted-unkeyed-skip"
        elif [ -f "${output_dir}/checkpoints.txt" ]; then
            decision="adopted-unkeyed-continue"
        else
            decision="invalidated-unkeyed"
        fi
    fi

    # Clear or set aside anything unusable before assembling.
    case "\$decision" in
        invalidated-*)
            if [ "${stale_policy}" = "archive" ]; then
                archive_dir="${output_dir}.stale.\$(date +%Y%m%dT%H%M%S)"
                echo "Setting aside unusable MEGAHIT checkpoint as \$archive_dir"
                mv "${output_dir}" "\$archive_dir"
                if [ -n "\$stored_key" ]; then
                    mv "\$key_file" "\$archive_dir.key"
                fi
            else
                echo "Discarding unusable MEGAHIT checkpoint at ${output_dir}"
                rm -rf "${output_dir}" "\$key_file"
            fi
            decision="\${decision}-fresh"
            ;;
    esac

    case "\$decision" in
        skip|adopted-unkeyed-skip)
            echo "Found a completed MEGAHIT assembly for the current inputs, skipping assembly."
            ;;
        continue|adopted-unkeyed-continue)
            echo "Found a MEGAHIT checkpoint for the current inputs, resuming."
            if ! megahit \\
                    ${args} \\
                    -m ${memory} \\
                    -t ${task.cpus} \\
                    --continue \\
                    -o ${output_dir}; then
                echo "Resume failed, wiping the checkpoint and restarting fresh."
                rm -rf "${output_dir}" "\$key_file"
                decision="resume-failed-rebuilt"
                mkdir -p ${checkpoint_base}
                echo "\$key" > "\$key_file"
                megahit \\
                    ${args} \\
                    -m ${memory} \\
                    -t ${task.cpus} \\
                    ${reads_command} \\
                    --out-prefix ${prefix} \\
                    -o ${output_dir}
            fi
            ;;
        *-fresh|fresh)
            echo "Starting a fresh MEGAHIT assembly."
            mkdir -p ${checkpoint_base}
            echo "\$key" > "\$key_file"
            megahit \\
                ${args} \\
                -m ${memory} \\
                -t ${task.cpus} \\
                ${reads_command} \\
                --out-prefix ${prefix} \\
                -o ${output_dir}
            ;;
        *)
            echo "ERROR: internal error, unhandled checkpoint decision '\$decision'" >&2
            exit 1
            ;;
    esac

    cp ${output_dir}/${prefix}.contigs.fa ${prefix}.v129mh.fasta

    # A checkpoint that is done can never be continued, so its intermediate graph is
    # dead weight: tens of GB per sample retained indefinitely under params.outdir.
    # Pruning it keeps the skip path working at a fraction of the disk.
    pruned="no"
    if [ "${cleanup}" = "true" ] && [ -f "${output_dir}/done" ]; then
        rm -rf ${output_dir}/intermediate_contigs ${output_dir}/tmp
        pruned="yes"
    fi

    # The skip path is otherwise invisible: the trace shows MEGAHIT running for thirty
    # seconds and emitting contigs, with nothing to say whether they were rebuilt.
    {
        printf 'sample\\t%s\\n' "${prefix}"
        printf 'checkpoint_decision\\t%s\\n' "\$decision"
        printf 'checkpoint_key\\t%s\\n' "\$key"
        printf 'checkpoint_key_stored\\t%s\\n' "\${stored_key:-none}"
        printf 'checkpoint_dir\\t%s\\n' "${output_dir}"
        printf 'checkpoint_pruned\\t%s\\n' "\$pruned"
        printf 'megahit_version\\t%s\\n' "\$megahit_version"
    } > ${prefix}.megahit_checkpoint.txt

    note="Assembles the trimmed reads to ${prefix}.v129mh.fasta. Checkpoint \$decision, key \${key:0:12}, at ${output_dir}."
    printf '<tr><td>MEGAHIT</td><td><samp>%s</samp></td><td>%s</td></tr>\\n' \\
        '${effective_args}' "\$note" > 22_megahit.tool_params_mqcrow.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$megahit_version
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end || !reads[1] ? "-r ${reads[0].join(',')}" : "-1 ${reads[0].join(',')} -2 ${reads[1].join(',')}"
    def memory = task.memory.toBytes()
    def output_dir = "${params.outdir}/megahit_checkpoints/${prefix}_megahit_out"
    def effective_args = [args, "-m ${memory}", "-t ${task.cpus}", reads_command, "--out-prefix ${prefix}"].findAll { it?.trim() }.join(' ')
    def note = "Assembles the trimmed reads to ${prefix}.v129mh.fasta and resumes from ${output_dir} when the checkpoint fingerprint matches."
    """
    touch ${prefix}.v129mh.fasta
    printf 'sample\\t${prefix}\\ncheckpoint_decision\\tstub\\n' > ${prefix}.megahit_checkpoint.txt
    cat <<-END_TOOL_PARAMS > 22_megahit.tool_params_mqcrow.html
    <tr><td>MEGAHIT</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
