#!/usr/bin/env bash
# Fingerprint the inputs of a MEGAHIT assembly.
#
# MEGAHIT checkpoints live under ${params.outdir}/megahit_checkpoints so they can
# survive a failed task, which means Nextflow's own caching cannot see them. Without
# a fingerprint, a checkpoint is keyed on the sample prefix alone: change the reads
# (insert kraken2 read-level decontamination, say) or change the assembly args, and
# the module short-circuits on the old "done" marker and republishes an assembly
# built from inputs that no longer exist. Every downstream metric then describes an
# assembly nobody asked for.
#
# This prints a sha256 over the things that change the assembly. Two exclusions are
# deliberate:
#
#   -m / -t (memory and threads)
#       These change on every retry via conf/base.config, and they do not change the
#       contigs. Including them would wipe the checkpoint on exactly the retry it
#       exists to serve.
#
#   realpath and mtime
#       A re-run with a fresh work directory, or reads restaged from the
#       precomputed_* published globs, gives byte-identical reads at a new path with
#       a new mtime. Keying on either would buy a full reassembly of every sample in
#       exchange for nothing.
#
# What is keyed instead is basename plus byte size, which is what actually separates
# the cases that matter: kraken2-filtered reads arrive under a different basename and
# a different size. The residual blind spot is a content change that preserves the
# exact byte length of a gzip stream, which does not happen in practice.

set -euo pipefail

usage() {
    cat >&2 <<'USAGE_EOF'
Usage: megahit_checkpoint_key.sh --version VERSION --args ARGS --prefix PREFIX READS...

  --version VERSION  megahit version string
  --args ARGS        assembly-affecting megahit arguments (-m/-t are stripped)
  --prefix PREFIX    --out-prefix passed to megahit
  READS...           the read files handed to this assembly

Prints a hex sha256 of the manifest. With --manifest, prints the manifest instead,
which is what you want when debugging why a key changed.
USAGE_EOF
}

version=""
args=""
prefix=""
show_manifest=false
reads=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --version) version="${2-}"; shift 2 ;;
        --args)    args="${2-}";    shift 2 ;;
        --prefix)  prefix="${2-}";  shift 2 ;;
        --manifest) show_manifest=true; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; reads+=("$@"); break ;;
        -*) echo "ERROR: unknown option '$1'" >&2; usage; exit 2 ;;
        *)  reads+=("$1"); shift ;;
    esac
done

if [[ -z "$prefix" ]]; then
    echo "ERROR: --prefix is required" >&2
    exit 2
fi

if [[ ${#reads[@]} -eq 0 ]]; then
    echo "ERROR: at least one read file is required" >&2
    exit 2
fi

# Strip the resource flags and collapse whitespace so that a reformatted but
# equivalent ext.args does not invalidate a checkpoint.
normalise_args() {
    local -a kept=()
    local -a tokens
    read -r -a tokens <<<"${1-}"
    local i=0
    while [[ $i -lt ${#tokens[@]} ]]; do
        case "${tokens[$i]}" in
            -m|-t|--memory|--num-cpu-threads)
                i=$((i + 2))
                continue
                ;;
            -m=*|-t=*|--memory=*|--num-cpu-threads=*)
                i=$((i + 1))
                continue
                ;;
            *)
                kept+=("${tokens[$i]}")
                i=$((i + 1))
                ;;
        esac
    done
    if [[ ${#kept[@]} -gt 0 ]]; then
        printf '%s\n' "${kept[*]}"
    else
        printf '\n'
    fi
}

# Collect the read entries before any output is produced. Doing this inside the
# manifest pipeline below would let a missing file emit a partial manifest to stdout
# before it failed, and a partial manifest hashes to a perfectly plausible key.
read_entries=()
for read_file in "${reads[@]}"; do
    if [[ ! -e "$read_file" ]]; then
        echo "ERROR: read file does not exist: $read_file" >&2
        exit 1
    fi
    read_entries+=("$(printf 'read=%s\t%s' "$(basename "$read_file")" "$(stat -Lc '%s' "$read_file")")")
done

manifest() {
    printf 'megahit_version=%s\n' "$version"
    printf 'args=%s\n' "$(normalise_args "$args")"
    printf 'out_prefix=%s\n' "$prefix"
    # Sorted so that channel ordering cannot move the key.
    printf '%s\n' "${read_entries[@]}" | sort
}

if [[ "$show_manifest" == true ]]; then
    manifest
else
    manifest | sha256sum | cut -d' ' -f1
fi
