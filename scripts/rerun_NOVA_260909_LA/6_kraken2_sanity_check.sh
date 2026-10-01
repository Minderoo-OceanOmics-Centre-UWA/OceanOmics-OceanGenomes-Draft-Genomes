#!/usr/bin/env bash
# Run BEFORE tier 4. Classifies one sponge and one fish and prints what kraken2 would
# remove from each.
#
# What to look for:
#   OG2630 (Porifera)     -- expect a large Bacteria fraction. Its assembly already
#                            inferred 13 prokaryotic primary divisions and tiara removed
#                            223 Mb from it, so if kraken2 reports almost nothing the
#                            database or the taxid list is wrong.
#   OG2941 (Actinopteri)  -- expect nearly everything unclassified. This is the control.
#                            If a meaningful fraction of the fish reads is being assigned
#                            to Bacteria, --kraken2_confidence is too low: STOP and raise
#                            it before running tier 4, because that same threshold would
#                            be deleting real host sequence from every sponge.
set -euo pipefail

module load singularity/4.1.0-nompi

RUN=NOVA_260909_LA
OUT="/scratch/pawsey1348/$USER/${RUN}"
DB="/scratch/references/kraken2/pluspfp_20230605"
CONF="${KRAKEN2_CONFIDENCE:-0.05}"
WORK="${OUT}/kraken2_sanity"
# The same image the KRAKEN2_KRAKEN2 module uses, taken from the pipeline's shared
# nextflow singularity cache. Do not put the https:// URL here: SINGULARITY_CACHEDIR
# points at pawsey0964, which is over its group quota, so the download fails only at
# close() and leaves a 0-byte file that then reports "image format not recognized".
CONTAINER="${KRAKEN2_SIF:-/software/projects/pawsey1348/singularity/nextflow_cache/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-0f-0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82-data.img}"

for f in hash.k2d opts.k2d taxo.k2d; do
    if [ ! -f "$DB/$f" ]; then
        echo "ERROR: $DB is not a loadable kraken2 database ($f missing)." >&2
        echo "ERROR: /scratch/references/kraken_dec2025 fails this check -- it has a" >&2
        echo "ERROR: kraken2 hash but no opts.k2d/taxo.k2d and is a Centrifuger build." >&2
        exit 1
    fi
done

mkdir -p "$WORK"
cd "$WORK"

for S in OG2630 OG2941; do
    R1=$(ls "$OUT"/draftgenomes/$S/fastp/*.R1.fastq.gz 2>/dev/null | head -1)
    R2=$(ls "$OUT"/draftgenomes/$S/fastp/*.R2.fastq.gz 2>/dev/null | head -1)
    if [ -z "$R1" ] || [ -z "$R2" ]; then
        echo "SKIP $S: fastp reads not found"; continue
    fi
    echo "=============================================================="
    echo " kraken2 sanity check: $S  (--confidence $CONF)"
    echo "=============================================================="
    # highmem, not work. hash.k2d is 147.5 GiB and kraken2 reads it into anonymous
    # memory, so on a 245000 MB work node the 230G cgroup has no headroom left for the
    # Lustre page cache of the same 158 GB file and the job gets OOM-killed. This is the
    # tier conf/base.config already retries KRAKEN2_KRAKEN2 on.
    srun --account=pawsey0964 --partition=highmem --cpus-per-task=16 --mem=460G --time=04:00:00 \
      singularity exec "$CONTAINER" \
        kraken2 --db "$DB" --threads 16 --gzip-compressed --paired \
                --confidence "$CONF" \
                --report "${S}.sanity.report.txt" \
                --output /dev/null \
                "$R1" "$R2"

    python3 - "$S" <<'PY'
import sys
s = sys.argv[1]
total = removed = 0
rows = []
excluded = {"2", "2157", "10239"}
for line in open(f"{s}.sanity.report.txt"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 6:
        continue
    n, rank, taxid, name = int(p[1]), p[3].strip(), p[4].strip(), p[5].strip()
    if rank == "U" or taxid == "1":
        total += n
    if taxid in excluded:
        removed += n
        rows.append((name, taxid, n))
unclass = next((int(l.split("\t")[1]) for l in open(f"{s}.sanity.report.txt")
                if len(l.split("\t")) > 3 and l.split("\t")[3].strip() == "U"), 0)
print(f"  total read pairs      : {total:,}")
print(f"  unclassified          : {unclass:,} ({unclass/total*100:.2f}%)" if total else "  unclassified: n/a")
for name, taxid, n in rows:
    print(f"  would remove {name} ({taxid}): {n:,} ({n/total*100:.2f}%)" if total else "")
print(f"  WOULD REMOVE TOTAL    : {removed:,} ({removed/total*100:.2f}%)" if total else "")
print(f"  WOULD RETAIN          : {total-removed:,} ({(total-removed)/total*100:.2f}%)" if total else "")
PY
    echo
done

echo "Reports left in $WORK for reference."
echo
echo "Decision: if OG2941 (the fish control) loses more than a few percent to Bacteria,"
echo "raise --kraken2_confidence and re-run this check before tier 4."
