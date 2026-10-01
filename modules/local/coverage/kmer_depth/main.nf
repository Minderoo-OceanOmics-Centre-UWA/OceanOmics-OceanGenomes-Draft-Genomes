// Read depth over single-copy BUSCO genes, as an independent estimate of the k-mer
// coverage GenomeScope fits.
//
// GenomeScope's fitted kmercov is the denominator of every genome size this pipeline
// publishes, and until now nothing checked it. Depth over Complete BUSCO genes estimates
// the same quantity from the alignments, independently of the k-mer histogram and of the
// assembly's total size, so the two can be compared.
//
// This half is samtools only. The container has awk and no python whatsoever, which is
// why the BED comes from bin/busco_complete_bed.awk and every file that needs real
// parsing is handed to SUMMARISE_KMER_COVERAGE instead. An earlier version called a
// python script from here and every task died with exit 127.
process MEASURE_KMER_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1':
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam), path(full_table)

    output:
    tuple val(meta), path("*_bedcov.tsv"), path("status.txt"), emit: bedcov
    tuple val(meta), path("*_contig_depth.tsv"), emit: contig_depth
    tuple val(meta), path("*_busco_complete.bed"), emit: regions
    tuple val(meta), path("32_measure_kmer_coverage.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    def min_mq = params.kmer_depth_min_mapq
    """
    # Contigs this BAM actually knows about, with their lengths. Regions naming anything
    # else are dropped by the BED step: bedcov fails the whole sample on the first
    # unresolvable line, and a stale BUSCO table names contigs from an assembly that no
    # longer exists.
    samtools view -H ${bam} \\
        | awk '/^@SQ/{split(\$2,a,":"); split(\$3,b,":"); print a[2]"\\t"b[2]}' > contigs.tsv

    awk -v count_out=region_count.txt -f \$(which busco_complete_bed.awk) \\
        contigs.tsv ${full_table} > ${prefix}_busco_complete.bed

    # Too few usable regions is a finding, not a failure. 12 of the 80 samples on
    # NOVA_260909_LA have fewer than 20 Complete BUSCOs because their assemblies are too
    # fragmented to recover them. Failing here would drop every one of them out of the
    # downstream joins and so out of the database -- strictly worse than an unmeasured
    # row that says why.
    REGIONS=\$(cat region_count.txt)
    if [ "\$REGIONS" -lt ${params.kmer_depth_min_regions} ]; then
        echo "too_few_busco_regions" > status.txt
        : > ${prefix}_bedcov.tsv
        # Empty, not absent. A sample with no measurement still has to reach the database
        # with a row that says so -- the same reasoning as the bedcov branch above.
        : > ${prefix}_contig_depth.tsv
    else
        echo "ok" > status.txt
        # Index into the task directory, never beside the published BAM.
        samtools index -@ ${task.cpus} -o ${prefix}.bam.bai ${bam}

        # bedcov's default filter already drops unmapped, secondary, QC-fail and duplicate
        # reads. MAPQ filtering is safe here: Q0 and Q20 agreed within 3% on every sample
        # tested, fragmented assemblies included.
        samtools bedcov -X -Q ${min_mq} ${prefix}_busco_complete.bed ${bam} ${prefix}.bam.bai \\
            > ${prefix}_bedcov.tsv

        # Exact per-contig mean depth, for the host/symbiont partition. `samtools coverage`
        # computes meandepth from the pileup; the reads x read-length approximation from
        # idxstats that the hand validation used over-counts soft-clipped and
        # partially-mapped reads, so it is not used here.
        samtools coverage -q ${min_mq} ${bam} > ${prefix}_contig_depth.tsv
    fi

    cat <<-END_TOOL_PARAMS > 32_measure_kmer_coverage.tool_params_mqcrow.html
    <tr><td>Measure K-mer Coverage</td><td><samp>samtools bedcov -Q ${min_mq} over Complete BUSCO genes; samtools coverage -q ${min_mq}</samp></td><td>Read depth over conserved single-copy genes, the input to an independent estimate of GenomeScope's fitted k-mer coverage, plus per-contig mean depth for the host/symbiont partition.</td></tr>
    END_TOOL_PARAMS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.prefix}"
    """
    : > ${prefix}_bedcov.tsv
    : > ${prefix}_contig_depth.tsv
    echo ok > status.txt
    touch ${prefix}_busco_complete.bed
    touch 32_measure_kmer_coverage.tool_params_mqcrow.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.20
    END_VERSIONS
    """
}
