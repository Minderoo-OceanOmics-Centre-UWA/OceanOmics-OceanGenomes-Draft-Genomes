// Drop read pairs kraken2 assigned to non-target clades (by default Bacteria, Archaea
// and Viruses) and keep everything else, including unclassified reads.
//
// Order matters: this runs before meryl and megahit. Screening contigs after assembly
// cannot undo a co-assembly of host and symbiont, which is what produced 150k-660k
// contig assemblies at N50 ~1000 for the sponges in NOVA_260909_LA.
//
// Reads are kept unless their assignment falls inside an excluded clade, so anything
// kraken2 could not classify -- which for an under-referenced invertebrate host is most
// of the target genome -- survives.
process KRAKENTOOLS_EXTRACTREADS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krakentools:1.2.1--pyh7e72e81_0':
        'quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(reads), path(classified_reads_assignment), path(report)
    val taxids

    output:
    tuple val(meta), path("*.kraken2filt.R{1,2}.fastq.gz")     , emit: reads
    tuple val(meta), path("*.kraken2_retention.txt")           , emit: retention
    tuple val(meta), path("*_kraken2_retention_mqc.yaml")      , emit: multiqc
    tuple val(meta), path("15_krakentools_extractreads.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def effective_args = ["extract_kraken_reads.py -t ${taxids}", '--exclude', '--include-children', '--fastq-output', args].findAll { it?.trim() }.join(' ')
    def note = "Removes read pairs assigned to taxids ${taxids} (and their descendants), keeping unclassified reads and every other assignment."
    """
    # extract_kraken_reads.py reads the per-read assignment file with plain open(), so
    # it has to be decompressed first. It is deleted again below to keep the work
    # directory from holding two copies of a whole library's assignments.
    zcat ${classified_reads_assignment} > kraken2_assignments.txt

    extract_kraken_reads.py \\
        -t ${taxids} \\
        --exclude \\
        --include-children \\
        --fastq-output \\
        -k kraken2_assignments.txt \\
        -r ${report} \\
        -s1 ${reads[0]} \\
        -s2 ${reads[1]} \\
        -o ${prefix}.kraken2filt.R1.fastq \\
        -o2 ${prefix}.kraken2filt.R2.fastq \\
        ${args}

    rm -f kraken2_assignments.txt

    gzip -f ${prefix}.kraken2filt.R1.fastq
    gzip -f ${prefix}.kraken2filt.R2.fastq

    # How much of the library was removed. Read this from the kraken2 report rather
    # than by counting FASTQ records: the report already carries cumulative read counts
    # per clade, so this costs nothing and stays exact. If a sample loses most of its
    # reads here, that is the headline result for that sample and needs to be visible
    # rather than buried in a log.
    python3 - <<'PY' > ${prefix}.kraken2_retention.txt
excluded = set("${taxids}".split())
total = 0
removed = 0
rows = []
with open("${report}") as fh:
    for line in fh:
        parts = line.rstrip("\\n").split("\\t")
        if len(parts) < 6:
            continue
        clade_reads = int(parts[1])
        rank = parts[3].strip()
        taxid = parts[4].strip()
        name = parts[5].strip()
        # "U" (unclassified) and "R" (root) together account for every read.
        if rank == "U" or taxid == "1":
            total += clade_reads
        if taxid in excluded:
            removed += clade_reads
            rows.append((name, taxid, clade_reads))

kept = total - removed
pct_removed = (removed / total * 100) if total else 0.0
pct_kept = (kept / total * 100) if total else 0.0

print("KRAKEN2 READ-LEVEL DECONTAMINATION")
print("=" * 50)
print()
print(f"Sample ID: ${meta.id}")
print(f"Excluded taxids: ${taxids} (including descendants)")
print()
print(f"Total read pairs (kraken2 fragments): {total:,}")
for name, taxid, n in rows:
    share = (n / total * 100) if total else 0.0
    print(f"  removed {name} (taxid {taxid}): {n:,} pairs ({share:.2f}%)")
print(f"Read pairs removed: {removed:,} ({pct_removed:.2f}%)")
print(f"Read pairs retained: {kept:,} ({pct_kept:.2f}%)")
PY

    python3 - <<'PY' > ${prefix}_kraken2_retention_mqc.yaml
import html, re

text = open("${prefix}.kraken2_retention.txt").read()
def grab(pattern):
    m = re.search(pattern, text)
    return m.group(1) if m else "NA"

rows = [
    ("Sample ID", "${meta.id}"),
    ("Excluded taxids", "${taxids}"),
    ("Total read pairs", grab(r"Total read pairs \\(kraken2 fragments\\): ([\\d,]+)")),
    ("Read pairs removed", grab(r"Read pairs removed: ([\\d,]+ \\([\\d.]+%\\))")),
    ("Read pairs retained", grab(r"Read pairs retained: ([\\d,]+ \\([\\d.]+%\\))")),
]

lines = [
    "id: 'nf-core-oceangenomesdraftgenomes-kraken2-retention'",
    "description: 'Read pairs removed before assembly because kraken2 assigned them to an excluded clade.'",
    "section_name: 'Kraken2 Read Decontamination'",
    "plot_type: 'html'",
    "data: |",
    "    <table class=\\"table table-condensed\\">",
    "    <tbody>",
]
for label, value in rows:
    lines.append(f"    <tr><th>{html.escape(str(label))}</th><td>{html.escape(str(value))}</td></tr>")
lines += ["    </tbody>", "    </table>"]
print("\\n".join(lines))
PY

    cat <<-END_TOOL_PARAMS > 15_krakentools_extractreads.tool_params_mqcrow.html
    <tr><td>KrakenTools ExtractReads</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def effective_args = ["extract_kraken_reads.py -t ${taxids}", '--exclude', '--include-children', '--fastq-output', args].findAll { it?.trim() }.join(' ')
    def note = "Removes read pairs assigned to taxids ${taxids} (and their descendants), keeping unclassified reads and every other assignment."
    """
    echo | gzip > ${prefix}.kraken2filt.R1.fastq.gz
    echo | gzip > ${prefix}.kraken2filt.R2.fastq.gz
    touch ${prefix}.kraken2_retention.txt
    touch ${prefix}_kraken2_retention_mqc.yaml

    cat <<-END_TOOL_PARAMS > 15_krakentools_extractreads.tool_params_mqcrow.html
    <tr><td>KrakenTools ExtractReads</td><td><samp>${effective_args}</samp></td><td>${note}</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2.1
    END_VERSIONS
    """
}
