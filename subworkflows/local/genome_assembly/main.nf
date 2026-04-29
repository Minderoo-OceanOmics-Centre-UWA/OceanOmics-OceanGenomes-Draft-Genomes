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
include { MEGAHIT                   } from '../../../modules/nf-core/megahit'


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
    // MODULE: Run Meryl count
    //
    
    MERYL_COUNT (
        fastp_reads, // tuple val(meta), path(reads)
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

    ch_coverage_calc = fastp_json.join(GENOMESCOPE2.out.summary, by:0)
    
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
        fastp_reads // tuple val(meta), path(reads1), path(reads2)
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
