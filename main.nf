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

include { DRAFTGENOMES  } from './workflows/draftgenomes'
include { MITOGENOMES  } from './subworkflows/local/oceangenomes_mitogenomes'
include { MITOGENOME_ANNOTATION  } from './subworkflows/local/oceangenomes_mitogenome_annotation_LCA'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'

include { MULTIQC                } from './modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow OCEANGENOMES_DRAFTGENOMES {

    take:
    //samplesheet // channel: samplesheet read in from --input
        run_id // params.run
        bs_config // params.bs_config
        curated_blast_db // params.curated_blast_db
        sql_config // params.sql_config
        
    main:
    
    ch_multiqc_files = Channel.empty()

    //
    // WORKFLOW: Run pipeline
    //
    DRAFTGENOMES (
        //samplesheet
        run_id, 
        bs_config
    )

    println "DRAFTGENOMES emits: ${DRAFTGENOMES.out.fastp_reads}"


    MITOGENOMES (
        fastp_reads    = DRAFTGENOMES.out.fastp_reads,
        organelle_type = "animal_mt"  // << pass it in here
    )

    println "MITOGENOMES emits: ${MITOGENOMES.out.mito_assembly}"

    MITOGENOME_ANNOTATION (
        mito_assembly   = MITOGENOMES.out.mito_assembly,
        curated_blast_db,
        sql_config // params.sql_config
    )

    // println "MITOGENOMES_ANNOTATION emits: ${MITOGENOMES_ANNOTATION.out.  }"

    // MITOGENOME_QC (

    // )
    // Collect all MultiQC files from all subworkflows
    ch_multiqc_files = ch_multiqc_files.mix(DRAFTGENOMES.out.multiqc_files)
    ch_multiqc_files = ch_multiqc_files.mix(MITOGENOMES.out.multiqc_files)
    ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ANNOTATION.out.multiqc_files)
    
    // Collect all versions
    ch_collated_versions = Channel.empty()
    ch_collated_versions = ch_collated_versions.mix(DRAFTGENOMES.out.versions)
    ch_collated_versions = ch_collated_versions.mix(MITOGENOMES.out.versions)
    ch_collated_versions = ch_collated_versions.mix(MITOGENOME_ANNOTATION.out.versions)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html 

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
        params.outdir
        // params.input // this is now inncliuded in the samplesheetHybrid wf
    )

    //
    // WORKFLOW: Run main workflow
    //
    OCEANGENOMES_DRAFTGENOMES (
        params.run,
        params.bs_config,
        params.curated_blast_db,
        params.sql_config
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
        OCEANGENOMES_DRAFTGENOMES.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
