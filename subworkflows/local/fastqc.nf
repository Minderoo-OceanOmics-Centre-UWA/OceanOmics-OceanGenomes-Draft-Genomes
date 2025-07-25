    /*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//download files
include { BASESPACE              } from '../modules/local/basespace/main'
include { samplesheetHybrid      } from '../subworkflows/local/samplesheetHybrid'
include { BBMAP_REPAIR           } from '../modules/nf-core/bbmap/repair/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTP {

    take:
    samplesheet // channel: samplesheet read in from --input og_id, run, R1, R2
    repaired_reads // output from bbmap.
    
    main:

    println "Run ID is: ${run_id}"  // ✅ Works fine here
    println "bs_config is: ${bs_config}"  // ✅ Works fine here
    //samplesheet.view { sheet -> "Samplesheet contents: $sheet" }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


// Make a new channel from bbmap_repair output or from samplesheet
    samplesheetHybrid(
        repaired_ch
    )
    
    // View samplesheet structure
    samplesheetHybrid.out.samplesheet.view()

    samplesheet_ch = samplesheetHybrid.out.samplesheet
    // MODULE: Run FastQC
    
    FASTQC (
        samplesheet_ch
    )

    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
   // sample_ch = params.input ?
    //     Channel.fromPath(params.input).map { parse_csv(it) } :
    //     generate_samplesheet_from_fastqs()//



    FASTP (
        samplesheet_ch, // tuple val(meta), path(reads)
        [], // path  adapter_fasta
        [], // val   discard_trimmed_pass
        [], // val   save_trimmed_fail
        [], // val   save_merged
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

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
