/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../../../modules/nf-core/fastqc'
include { FASTP                  } from '../../../modules/nf-core/fastp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN FASTP AND FASQ WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTP_FASTQC {

    take:
    samplesheet // tuple(meta, reads) meta: id, run, date, prefix
    
    main:

    // Debugging lines to check the input
    // println "samplesheet is: ${samplesheet}"

    // Create empty channels to collect version and multiqc files
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_inputs = Channel.empty()

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        samplesheet
    )

    //
    // MODULE: Run FastQC
    //

    FASTP (
        samplesheet, // tuple val(meta), path(reads)
        [], // path  adapter_fasta
        [], // val   discard_trimmed_pass
        [], // val   save_trimmed_fail
        [], // val   save_merged
    )

    //
    // Collect files
    //

    ch_multiqc_inputs = ch_multiqc_inputs.mix(FASTQC.out.zip)
    ch_multiqc_inputs = ch_multiqc_inputs.mix(FASTP.out.json)
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_inputs.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(FASTP.out.versions.first())


    //
    // Emit outputs
    //

    emit:
    fastp_reads = FASTP.out.reads // channel to the mitogenome pipeline, and genome assembly tuple val(meta), path('*.fastq.gz')
    fastp_json = FASTP.out.json // channel to the mitogenome pipeline, tuple val(meta), path('*.json')
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    multiqc_inputs = ch_multiqc_inputs           // channel: [ tuple(meta), path(multiqc_file) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]
}
