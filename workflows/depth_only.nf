/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEPTH_ONLY: measure read depth for drafts assembled before MEASURE_KMER_COVERAGE existed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Run with: nextflow run main.nf -entry DEPTH_ONLY --depth_samplesheet <csv> --outdir <dir>
              --sql_config <cfg>

    Drafts up to 260909 have no busco_read_depth, so ENA-draft-genomes has nothing to send as
    COVERAGE (genomes_db_review Finding 17.2; plan docs/r6_gate_plan.md, A.4). Everything the
    measurement needs is already archived beside each draft: the assembly, the BAM of its reads
    mapped back to it, the BUSCO full table and the GenomeScope model. So this runs only the
    two halves of the depth module, and SEQKIT_STATS where an archived run has no stats.

    Nothing else is recomputed: no re-mapping, BUSCO, reseed or RECHECK_GENOME_SIZE. The
    genome size, the reliability verdict and every other QC value in draft_genomes stay as they
    are; the upload writes the depth and host partition columns only. A full GENOME_QC re-run
    of these drafts is a separate decision.

    Samplesheet columns (CSV, header row):
      og_id, seq_date            identify the draft row in draft_genomes
      assembly                   the published .fna
      bam                        assemblies/genome/bwa/*.sorted.bam, mapped to that assembly
      busco_full_table           the BUSCO *full_table.tsv for that assembly
      genomescope_model          the GenomeScope *_model.txt
      fastp_json                 the fastp report (read length)
      seqkit_stats               optional; computed from the assembly when blank

    If the BAM was mapped to a different version of the assembly, its contig names will not
    match the BUSCO table. MEASURE_KMER_COVERAGE drops regions naming contigs the BAM does not
    have, and reports too_few_busco_regions when that leaves under kmer_depth_min_regions.
    Check the first samples' status before running a whole batch.
*/

include { SEQKIT_STATS              } from '../modules/nf-core/seqkit/stats/main'
include { MEASURE_KMER_COVERAGE     } from '../modules/local/coverage/kmer_depth'
include { SUMMARISE_KMER_COVERAGE   } from '../modules/local/coverage/kmer_depth_summarise'
include { PUSH_KMER_DEPTH_RESULTS   } from '../modules/local/upload_results/kmer_depth'

workflow DEPTH_ONLY {

    take:
    samplesheet // path: the CSV described above
    sql_config  // path: database config for the upload

    main:
    ch_rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            for (column in ['og_id', 'seq_date', 'assembly', 'bam', 'busco_full_table',
                            'genomescope_model', 'fastp_json']) {
                if (!row[column]) {
                    error "depth samplesheet: ${row.og_id ?: '?'}/${row.seq_date ?: '?'} has no ${column}"
                }
            }
            def meta = [id: row.og_id, seq_date: row.seq_date, prefix: "${row.og_id}.ilmn.${row.seq_date}"]
            [meta, row]
        }

    // Archived runs before 260514 have no seqkit stats. SUMMARISE_KMER_COVERAGE needs the
    // assembly size from them, and seqkit on one assembly is cheap.
    ch_split = ch_rows.branch { meta, row ->
        has_stats: row.seqkit_stats
        needs_stats: true
    }
    SEQKIT_STATS (
        ch_split.needs_stats.map { meta, row -> [meta, file(row.assembly, checkIfExists: true)] }
    )
    ch_seqkit = ch_split.has_stats
        .map { meta, row -> [meta.prefix, file(row.seqkit_stats, checkIfExists: true)] }
        .mix(SEQKIT_STATS.out.stats.map { meta, stats -> [meta.prefix, stats] })

    MEASURE_KMER_COVERAGE (
        ch_rows.map { meta, row ->
            [meta, file(row.bam, checkIfExists: true), file(row.busco_full_table, checkIfExists: true)]
        }
    )

    // Every input is keyed on meta.prefix, so a sample that fails one step drops out rather
    // than being paired with another sample's files.
    ch_summarise_input = MEASURE_KMER_COVERAGE.out.bedcov
        .map { meta, bedcov, status -> [meta.prefix, meta, bedcov, status] }
        .join(MEASURE_KMER_COVERAGE.out.contig_depth.map { meta, depth -> [meta.prefix, depth] })
        .join(ch_rows.map { meta, row ->
            [meta.prefix, file(row.fastp_json, checkIfExists: true), file(row.genomescope_model, checkIfExists: true)]
        })
        .join(ch_seqkit)
        .map { prefix, meta, bedcov, status, depth, fastp_json, model, stats ->
            [meta, bedcov, status, depth, fastp_json, stats, model]
        }

    SUMMARISE_KMER_COVERAGE ( ch_summarise_input )

    if (!params.skip_upload_results) {
        PUSH_KMER_DEPTH_RESULTS ( SUMMARISE_KMER_COVERAGE.out.kmer_depth, sql_config )
    }

    emit:
    kmer_depth = SUMMARISE_KMER_COVERAGE.out.kmer_depth // channel: tuple val(meta), path("*_kmer_depth.json")
}
