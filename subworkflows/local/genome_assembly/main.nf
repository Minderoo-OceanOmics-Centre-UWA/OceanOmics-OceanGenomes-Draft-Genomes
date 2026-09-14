/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//assembly
include { MERYL_COUNT               } from '../../../modules/nf-core/meryl/count'
include { MERYL_UNIONSUM            } from '../../../modules/nf-core/meryl/unionsum'
include { MERYL_HISTOGRAM           } from '../../../modules/nf-core/meryl/histogram'
include { GENOMESCOPE2              } from '../../../modules/nf-core/genomescope2'
include { CALCULATE_SEQUENCING_COVERAGE  } from '../../../modules/local/coverage/calculations'
include { COMPILE_JSON_TO_CSV       } from '../../../modules/local/coverage/compile'
include { MEGAHIT                   } from '../../../modules/local/megahit'
include { KRAKEN2_KRAKEN2           } from '../../../modules/local/kraken2/kraken2'
include { KRAKENTOOLS_EXTRACTREADS  } from '../../../modules/local/krakentools/extractreads'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOME_ASSEMBLY {

    take:
    fastp_reads // tuple val(meta), path('*.fastq.gz')
    fastp_json
    //samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_sample_multiqc_inputs = Channel.empty()

    //
    // MODULE: Run kraken2 read-level decontamination
    //
    // This has to happen before meryl and megahit, not after assembly. When half a
    // library is bacterial, megahit co-assembles host and symbiont and no amount of
    // contig screening afterwards recovers the contiguity that was lost, and the k-mer
    // histogram has no separable genomic peak for GenomeScope to fit.
    //
    if (!params.skip_kraken2_decontamination) {
        KRAKEN2_KRAKEN2 (
            fastp_reads, // tuple val(meta), path(reads)
            file(params.kraken2_db, checkIfExists: true) // path db
        )
        ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

        ch_extract_input = fastp_reads
            .join(KRAKEN2_KRAKEN2.out.classified_reads_assignment, by: 0)
            .join(KRAKEN2_KRAKEN2.out.report, by: 0)

        KRAKENTOOLS_EXTRACTREADS (
            ch_extract_input, // tuple val(meta), path(reads), path(assignment), path(report)
            params.kraken2_exclude_taxids // val taxids
        )
        ch_versions = ch_versions.mix(KRAKENTOOLS_EXTRACTREADS.out.versions.first())

        ch_assembly_reads = KRAKENTOOLS_EXTRACTREADS.out.reads

        ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(KRAKEN2_KRAKEN2.out.tool_params)
        ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(KRAKENTOOLS_EXTRACTREADS.out.tool_params)
        ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(KRAKENTOOLS_EXTRACTREADS.out.multiqc)
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.tool_params.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(KRAKENTOOLS_EXTRACTREADS.out.tool_params.collect { it[1] })
    } else {
        ch_assembly_reads = fastp_reads
    }

    //
    // MODULE: Run Meryl count
    //
    // Counted from the decontaminated reads so the k-mer profile, the GenomeScope size
    // estimate and the merqury QV all describe the same sequence as the assembly.
    
    MERYL_COUNT (
        ch_assembly_reads, // tuple val(meta), path(reads)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())
    
    //
    // MODULE: Run Meryl unionsum
    //

    MERYL_UNIONSUM (
        MERYL_COUNT.out.meryl_dbs, // tuple val(meta), path(meryl_dbs)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())

    //
    // MODULE: Run Meryl histogram
    //

    MERYL_HISTOGRAM (
        MERYL_UNIONSUM.out.meryl_db, // tuple val(meta), path(meryl_dbs)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run Genomescope
    //

    GENOMESCOPE2 (
        MERYL_HISTOGRAM.out.hist // tuple val(meta), path(histogram)
    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())

    // The meryl histogram goes in alongside the GenomeScope summary so the coverage
    // step can tell a real genomic k-mer peak from an error shoulder, and refuse to
    // grade coverage when GenomeScope had nothing to fit.
    ch_coverage_calc = fastp_json
        .join(GENOMESCOPE2.out.summary, by:0)
        .join(GENOMESCOPE2.out.model, by:0)
        .join(MERYL_HISTOGRAM.out.hist, by:0)
    
    //
    // MODULE: Calculate Sequencing Coverage
    //
    
    CALCULATE_SEQUENCING_COVERAGE(
        ch_coverage_calc
    )
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(CALCULATE_SEQUENCING_COVERAGE.out.multiqc)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(CALCULATE_SEQUENCING_COVERAGE.out.tool_params)
    
    // Collect all JSON outputs from coverage calculation to compile into a single CSV
    collected_jsons = CALCULATE_SEQUENCING_COVERAGE.out.coverage_json.collect()
    
    //
    // MODULE: Compile samples coverage calculations into one CSV
    //

    COMPILE_JSON_TO_CSV(collected_jsons)
    ch_multiqc_files = ch_multiqc_files.mix(COMPILE_JSON_TO_CSV.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(COMPILE_JSON_TO_CSV.out.tool_params)

    //
    // MODULE: Run Megahit
    //

    MEGAHIT (
        ch_assembly_reads // tuple val(meta), path(reads1), path(reads2)
    )
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

    //
    // Collect files
    //

    // ch_multiqc_files = ch_multiqc_files.mix(BASESPACE.out.json.collect{it[1]})
    // ch_multiqc_files = ch_multiqc_files.mix(BBMAP_REPAIR.out.log.collect{it})
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MERYL_COUNT.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MERYL_UNIONSUM.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MERYL_HISTOGRAM.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(GENOMESCOPE2.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MEGAHIT.out.tool_params)
    ch_multiqc_files = ch_multiqc_files.mix(MERYL_COUNT.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MERYL_UNIONSUM.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MERYL_HISTOGRAM.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GENOMESCOPE2.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_SEQUENCING_COVERAGE.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MEGAHIT.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())
    ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())
    ch_versions = ch_versions.mix(CALCULATE_SEQUENCING_COVERAGE.out.versions.first())
    ch_versions = ch_versions.mix(COMPILE_JSON_TO_CSV.out.versions.first())
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    megahit_assembled_contigs = MEGAHIT.out.contigs // pass to dweconatmination - fcs-gx
    meryl_db = MERYL_UNIONSUM.out.meryl_db  // need to check COUNT output to make sure im passing the right one
    genomescope_summary = GENOMESCOPE2.out.summary // pass into genome QC for size of genome (unique length)
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    multiqc_inputs = ch_sample_multiqc_inputs    // channel: [ tuple(meta), path(multiqc_file) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
