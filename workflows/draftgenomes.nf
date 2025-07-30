/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//download files
include { BASESPACE               } from '../modules/local/basespace/main'

include { FASTP                   } from './subworkflows/local/fastqc'
include { GENOME_ASSEMBLY         } from './subworkflows/local/genome_assembly'
include { GENOME_DECONTAMINATION  } from './subworkflows/local/genome_decontamination'
include { GENOME_QC               } from './subworkflows/local/genome_qc'

include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DRAFTGENOMES {

    take:
    repaired_reads
    //samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    START DRAFT GENOME ASSEMBLY FROM REPAIRED READS
        YOU CAN PROVIDE A SAMPLE SHEET
        OR META CAN BE DEFINED BY THE FILES PROVIDED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Make a new channel from bbmap_repair output or from samplesheet
    samplesheetHybrid(
        repaired_reads
    )
    
    // View samplesheet structure
    samplesheetHybrid.out.samplesheet.view()

    samplesheet_ch = samplesheetHybrid.out.samplesheet
    // MODULE: Run FastQC
    
    FASTP (
        samplesheet_ch
    )

    GENOME_ASSEMBLY (

    )

    GENOME_DECONTAMINATION (

    )

    GENOME_QC (

    )

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

    emit:
    fastp_reads = FASTP.out.reads // channel to the mitogenome pipeline
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
}
