/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME_SIZE_ONLY: re-check the GenomeScope genome size of drafts that have a measured depth
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Run with: nextflow run main.nf -entry GENOME_SIZE_ONLY --genome_size_samplesheet <csv>
              --outdir <dir> --sql_config <cfg>

    Drafts assembled before 260909 never went through GENOMESCOPE_RESEED or
    RECHECK_GENOME_SIZE, so their stored genomesize has no reliability verdict. The depth-only
    re-run (workflows/depth_only.nf) measured their read depth, and for 183 of 542 the fitted
    GenomeScope coverage disagrees with it by more than 25% -- OG1076 stored 364 Mb against a
    2.68 Gb host assembly (genomes_db_review docs/r6_gate_a3_a5_record.md, A.4 and A.6).

    This runs the genome-size half of GENOME_QC, unchanged, on archived files:
      CALCULATE_SEQUENCING_COVERAGE  the provisional verdict, as assembly writes it
      GENOMESCOPE_RESEED             the refit seeded with the measured depth, where flagged
      RECHECK_GENOME_SIZE            the final verdict against the assembly
    and uploads the chosen fit with its verdict through PUSH_ASSEMBLY_RESULTS.

    No BAM or assembly is needed: the depth is already in *_kmer_depth.json and the assembly
    size in the seqkit stats. Nothing else is recomputed. BUSCO, Merqury, gfastats and the
    decontamination numbers do not depend on the genome size and stay as they are.

    Samplesheet columns (CSV, header row), all required:
      og_id, seq_date            identify the draft row in draft_genomes
      fastp_json                 the fastp report
      genomescope_summary        the GenomeScope *_summary.txt
      genomescope_model          the GenomeScope *_model.txt
      meryl_hist                 the *.meryl.hist GenomeScope was fitted to
      kmer_depth_json            the depth-only *_kmer_depth.json, status ok
      seqkit_stats               the assembly's *.seqkit_stats.tsv
      busco_short_summary        BUSCO *short_summary.json, the lineage whose full table the
                                 depth was measured over
*/

include { CALCULATE_SEQUENCING_COVERAGE } from '../modules/local/coverage/calculations'
include { GENOMESCOPE_RESEED            } from '../modules/local/genomescope_reseed'
include { RECHECK_GENOME_SIZE           } from '../modules/local/coverage/recheck'
include { PUSH_ASSEMBLY_RESULTS         } from '../modules/local/upload_results/assembly'

workflow GENOME_SIZE_ONLY {

    take:
    samplesheet // path: the CSV described above
    sql_config  // path: database config for the upload

    main:
    def columns = ['og_id', 'seq_date', 'fastp_json', 'genomescope_summary', 'genomescope_model',
                   'meryl_hist', 'kmer_depth_json', 'seqkit_stats', 'busco_short_summary']
    ch_rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            for (column in columns) {
                if (!row[column]) {
                    error "genome size samplesheet: ${row.og_id ?: '?'}/${row.seq_date ?: '?'} has no ${column}"
                }
            }
            def meta = [id: row.og_id, seq_date: row.seq_date, prefix: "${row.og_id}.ilmn.${row.seq_date}"]
            def files = columns.drop(2).collectEntries { [(it): file(row[it], checkIfExists: true)] }
            [meta, files]
        }

    CALCULATE_SEQUENCING_COVERAGE (
        ch_rows.map { meta, f -> [meta, f.fastp_json, f.genomescope_summary, f.genomescope_model, f.meryl_hist] }
    )

    // Every input is keyed on meta.prefix, so a sample that fails one step drops out rather
    // than being paired with another sample's files.
    ch_by_prefix = ch_rows.map { meta, f -> [meta.prefix, meta, f] }

    ch_reseed_input = CALCULATE_SEQUENCING_COVERAGE.out.coverage_json_meta
        .map { meta, coverage -> [meta.prefix, coverage] }
        .join(ch_by_prefix)
        .map { prefix, coverage, meta, f ->
            [meta, f.meryl_hist, coverage, f.kmer_depth_json, f.genomescope_summary, f.genomescope_model]
        }

    GENOMESCOPE_RESEED ( ch_reseed_input )

    ch_recheck_input = CALCULATE_SEQUENCING_COVERAGE.out.coverage_json_meta
        .map { meta, coverage -> [meta.prefix, coverage] }
        .join(GENOMESCOPE_RESEED.out.decision.map { meta, decision -> [meta.prefix, decision] })
        .join(GENOMESCOPE_RESEED.out.chosen.map { meta, summary, model -> [meta.prefix, summary, model] })
        .join(ch_by_prefix)
        .map { prefix, coverage, decision, summary, model, meta, f ->
            [meta, coverage, f.seqkit_stats, f.kmer_depth_json, f.busco_short_summary,
             decision, summary, model, f.meryl_hist]
        }

    RECHECK_GENOME_SIZE ( ch_recheck_input )

    // The fit that was chosen goes up with the verdict about it, as PUSH_ASSEMBLY_RESULTS
    // does at the end of a full run.
    ch_push_input = GENOMESCOPE_RESEED.out.chosen
        .map { meta, summary, model -> [meta.prefix, meta, summary] }
        .join(RECHECK_GENOME_SIZE.out.coverage_json.map { meta, coverage -> [meta.prefix, coverage] })
        .map { prefix, meta, summary, coverage -> [meta, summary, coverage] }

    if (!params.skip_upload_results) {
        PUSH_ASSEMBLY_RESULTS ( ch_push_input, sql_config )
    }

    emit:
    coverage_summary = RECHECK_GENOME_SIZE.out.coverage_json // channel: tuple val(meta), path("*_coverage_summary.json")
    reseed_decision = GENOMESCOPE_RESEED.out.decision         // channel: tuple val(meta), path("*_reseed_decision.json")
}
