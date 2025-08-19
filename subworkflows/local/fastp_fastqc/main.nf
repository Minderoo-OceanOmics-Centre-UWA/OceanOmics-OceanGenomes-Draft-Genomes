/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../../../modules/nf-core/fastqc'
include { FASTP                  } from '../../../modules/nf-core/fastp'

// Helper functions
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'

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

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(FASTP.out.versions.first())


    //
    // Collate and save software versions
    //

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'oceangenomes_draftgenomes_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    
    //
    // Emit outputs
    //

    emit:
    fastp_reads = FASTP.out.reads // channel to the mitogenome pipeline, and genome assembly tuple val(meta), path('*.fastq.gz')
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
}
