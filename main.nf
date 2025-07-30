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
include { MITOGENOME_ASSEMBLY  } from './subworkflows/local/mitogenome_assembly/getorganelle'
include { MITOGENOME_ANNOTATION  } from './subworkflows/local/mitogenome_annotation_lca'
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
    if (!params.skip_download_reads) {
        DOWNLOAD_READS(
            run_id, 
            bs_config
        )
    }
// if i am providing path to precomputed results i need code to re make meta and tuples
// meta is just: 
//  id = og_id
//  run = params.run - maybe this can be determined from file name if passing in info as may not be running all from the same run
//  date = parames.date - or as above, pull from file name
//  prefix = ${meta.id}.ilmn.${run}

    //
    // FASTP
    //
    if (!params.skip_fastp) {
        FASTP (
            samplesheet_ch
        )
        ch_fastp_results = 
    } else if (params.precomputed_fastp_results) {
        // Use precomputed results if analysis is skipped
        ch_fastp_results = Channel.fromPath(params.precomputed_fastp_results)
    } else {
        ch_fastp_results = Channel.empty()
    }
    
    //
    // GENOME_ASSEMBLY
    //
    if (!params.skip_genome_assembly) {
        GENOME_ASSEMBLY (

        )
        ch_genome_assembly_results = 
    } else if (params.precomputed_genome_assembly_results) {
        // Use precomputed results if analysis is skipped
        ch_genome_assembly_results = Channel.fromPath(params.precomputed_genome_assembly_results)
    } else {
        ch_genome_assembly_results = Channel.empty()
    }

    // GENOME_DECONTAMINATION
    if (!params.skip_genome_decontamination) {
        GENOME_DECONTAMINATION (

        )
        ch_genome_decontamination_results = 
    } else if (params.precomputed_genome_decontamination_results) {
        // Use precomputed results if analysis is skipped
        ch_genome_decontamination_results = Channel.fromPath(params.precomputed_genome_decontamination_results
    } else {
        ch_genome_decontamination_results = Channel.empty()
    }

    //
    // GENOME_QC
    //
    if (!params.skip_genome_qc) {
        GENOME_QC (

        )
        ch_genome_qc_results = // not sure if there is output from this one that will go anywhere, actually mauybe data upload
    } else if (params.precomputed_genome_qc_results) {
        // Use precomputed results if analysis is skipped
        ch_genome_qc_results = Channel.fromPath(params.precomputed_genome_qc_results)
    } else {
        ch_genome_qc_results = Channel.empty()
    }

    //
    // MITOGENOME_ASSEMBLY
    //
    if (!params.skip_upload_results) {
        MITOGENOME_ASSEMBLY (
            ch_fastp_results,
            organelle_type = "animal_mt"  // << pass it in here
        )
        ch_mitogenome_assembly_results = 
    } else if (params.precomputed_mitogenome_assembly_results) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_assembly_results = Channel.fromPath(params.precomputed_mitogenome_assembly_results)
    } else {
        ch_mitogenome_assembly_results = Channel.empty()
    }

    //
    // MITOGENOME_ANNOTATION
    //
    if (!params.skip_upload_results) {
        MITOGENOME_ANNOTATION (
            ch_mitogenome_assembly_results,
            curated_blast_db,
            sql_config // params.sql_config
        )
        ch_mitogenome_annotation_results = MITOGENOME_ANNOTATION.out.annotation_results
        ch_mitogenome_blast_results = MITOGENOME_ANNOTATION.out.blast_filtered_results
        ch_mitogenome_lca_results = MITOGENOME_ANNOTATION.out.lca_results
    } else if (params.precomputed_mitogenome_annotation_results) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_annotation_results = Channel.fromPath(params.precomputed_mitogenome_annotation_results)
        ch_mitogenome_blast_results = Channel.fromPath(params.precomputed_mitogenome_blast_results)
        ch_mitogenome_lca_results = Channel.fromPath(params.precomputed_mitogenome_lca_results)
    } else {
        ch_mitogenome_annotation_results = Channel.empty()
        ch_mitogenome_blast_results = Channel.empty()
        ch_mitogenome_lca_results = Channel.empty()
    }

    //
    // UPLOAD_RESULTS
    //
    // Conditional uploading of results to SQL and species check - only run if not skipped
    // All these processes access the OceanOmics PostgreSQL database.
    if (!params.skip_upload_results) {
        UPLOAD_RESULTS (
            assembly_results        = MITOGENOME_ASSEMBLY.out.assembly_results
            annotation_results      = MITOGENOME_ANNOTATION.out.annotation_results
            blast_filtered_results  = MITOGENOME_ANNOTATION.out.blast_filtered_results
            lca_results             = MITOGENOME_ANNOTATION.out.lca_results
            sql_config // params.sql_config
        )
    }


    // If the LCA validation is correct, then run the QC to prepare for submission to GenBank
    // Need to add this into the pipeline.
    // Now that protein lengths are being added to the database it could provide a list of 
    // non submitted mitogenomes they can be grouped with to submit and then say when there is a 
    // group of similar mitogenomes they can be submitted as a batch.
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
