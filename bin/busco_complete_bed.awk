#!/usr/bin/awk -f
#
# BUSCO "Complete" genes -> BED, validated against the contigs an alignment actually has.
#
#   busco_complete_bed.awk contigs.tsv full_table.tsv > out.bed
#
# contigs.tsv is "name<TAB>length" per line, from `samtools view -H | awk '/^@SQ/...'`.
#
# This is awk rather than python because it runs in the samtools container, which has
# /usr/bin/awk and no python at all -- a task that needs both tools needs the one container
# that has them.
#
# Coordinate handling, all of it load-bearing:
#   * full_table is 1-based inclusive; BED is 0-based half-open, so the start shifts by one.
#   * Start > End on the minus strand, so the pair is normalised before use.
#   * Regions naming a contig the alignment does not have are DROPPED, not passed through:
#     samtools bedcov fails the entire sample on the first unresolvable line, and a sample
#     re-run against a different BUSCO lineage keeps a stale table whose contigs belong to
#     an assembly that no longer exists.
#   * Regions running past a contig's end are clipped to it.
#
# Writes the usable region count to stderr, and to the file named by -v count_out= if given.
BEGIN { FS = "\t"; OFS = "\t"; kept = 0; skipped = 0; total = 0 }

# First file: the contig lengths.
NR == FNR { len[$1] = $2; next }

# Second file: the BUSCO full table.
/^#/ { next }
$2 != "Complete" { next }
$3 == "" { next }
{
    total++
    if (!($3 in len)) { skipped++; next }

    start = $4 + 0
    end = $5 + 0
    if (start > end) { tmp = start; start = end; end = tmp }
    if (start < 1) start = 1
    if (end > len[$3]) end = len[$3]
    if (end < start) { skipped++; next }

    print $3, start - 1, end, $1
    kept++
}

END {
    printf "%d usable of %d Complete BUSCOs (%d unusable against this BAM)\n",
           kept, total, skipped > "/dev/stderr"
    if (count_out != "") printf "%d\n", kept > count_out
}
