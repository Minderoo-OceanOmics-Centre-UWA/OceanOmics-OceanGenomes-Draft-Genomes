#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/oceangenomes_draftgenomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/oceangenomes_draftgenomes
    Website: https://nf-co.re/oceangenomes_draftgenomes
    Slack  : https://nfcore.slack.com/channels/oceangenomes_draftgenomes
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { OCEANGENOMES_DRAFTGENOMES } from './workflows/draftgenomes'

// Draft genome assembly subworkflows
include { DOWNLOAD_READS            } from './subworkflows/local/download'
include { samplesheetHybrid         } from './subworkflows/local/samplesheetHybrid'

// Pipeline subworkflows
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//

workflow NFCORE_OCEANGENOMES_DRAFTGENOMES {

    take:
        run_id // params.run
        bs_config // params.bs_config
        sql_config // params.sql_config
        
    main:
    
    ch_multiqc_files = Channel.empty()
    samplesheet = params.input ? Channel.fromPath(params.input) : Channel.empty()
    ch_samplesheetHybrid_results = Channel.empty()

    //
    // WORKFLOW: Run pipeline
    //

    if (!params.skip_download_reads) {
        DOWNLOAD_READS(
            run_id, 
            bs_config
        )
        ch_download_reads_results = DOWNLOAD_READS.out.repaired_reads // tuple val(ogid), path("${prefix}.R*.fq.gz") 
    } else if (params.precomputed_download_reads_results) {
        // Use precomputed results if analysis is skipped
        ch_download_reads_results = Channel.fromFilePairs(params.precomputed_download_reads_results, checkIfExists: true)
            .map { sample_id, reads ->
                def meta_id = sample_id.split('\\.')[0] // Extract meta_id which is the first part of sample name seperated by .
                return tuple(meta_id, reads) // or tuple(meta, reads) depending on your downstream processes
            }
    } else {
        ch_download_reads_results = Channel.empty()
    }

    //
    // Run the samplesheetHybrid to process input file or the DOWNLOAD_READS output
    //

    if (!params.skip_fastp_fastqc) {
        samplesheetHybrid(
            ch_download_reads_results,
            samplesheet
        )
        ch_samplesheetHybrid_results = samplesheetHybrid.out.samplesheet // tuple(meta, reads) meta: id, run, date, prefix
    }

    //
    // WORKFLOW: Run main Draft Genomes workflow
    //
    
    OCEANGENOMES_DRAFTGENOMES (
        ch_samplesheetHybrid_results,
        sql_config
    )


    //
    // Emit outputs
    //

    emit:
    multiqc_report = OCEANGENOMES_DRAFTGENOMES.out.multiqc_report // channel: /path/to/multiqc_report.html 

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        // params.input { optional true } // this is now inncliuded in the samplesheetHybrid wf
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_OCEANGENOMES_DRAFTGENOMES (
        params.run,
        params.bs_config,
        params.sql_config,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_OCEANGENOMES_DRAFTGENOMES.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
