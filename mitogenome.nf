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

// Mitogenome assembly subworkflows
include { MITOGENOME_ASSEMBLY       } from './subworkflows/local/mitogenome_assembly/getorganelle'
include { MITOGENOME_ANNOTATION     } from './subworkflows/local/mitogenome_annotation_lca'
include { UPLOAD_RESULTS            } from './subworkflows/local/upload_results_mito'
// include { MITOGENOME_QC  } from './subworkflows/local/mitogenome_qc'

// Pipeline subworkflows
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { MULTIQC                   } from './modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from './subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//

workflow OCEANGENOMES_MITOGENOMES {

    take:
        run_id // params.run
        bs_config // params.bs_config
        curated_blast_db // params.curated_blast_db
        sql_config // params.sql_config
        organelle_type // params.organelle_type "animal_mt"
        
    main:
    
    ch_multiqc_files = Channel.empty()
    samplesheet = params.input ? Channel.fromPath(params.input) : Channel.empty()
    /* The if and else statements in this workflow are for when steps are skipped in the nextflow_run script.
        What it is doing is 'if' this processes isnt skipped then run the subworkflow and provide the standard outputs.
        
        Then 'else if' is if the process has been skipped then we will create the required input for the following
        subworkflows using the predefined file paths for the output in the nextflow.config, the "precomputed_*" file
        paths. 
            The relevant information is extracted to create the meta map from the parts of the file name. This is
            assuming that files are named with the sample id, then the type of sequencing, the date and then the other
            information in the file name, all seperated by a '.'
            The mitogenome sections assume that the files are named as above, $sample_id.$sequencing_technology.$date with
            additional information added to this with each proccess. The assembly process will add a .getorganelle${version}
            after the inital prefix and befor the file extension. Then the annotation process will add a .emma${version} to 
            the previous prefix.
        
        The 'else' statement is then just creating and empty channel if neither of the previous steps worked.
    */

    //
    // WORKFLOW: Run pipeline
    //
    
    ch_fastp_fastqc_results = Channel.fromFilePairs(params.precomputed_fastp_fastqc_results, checkIfExists: true)
    .map { sample_id, reads ->
        def basename = reads.baseName
        def meta_id = sample_id.split('\\.')[0]
        def meta = [
            id: meta_id,
            run: params.run,
            date: params.date,
            assembly_prefix: sample_id
        ]
        return tuple(meta, reads)
    }

    // 
    // SUBWORKFLOW: MITOGENOME_ASSEMBLY
    //

    if (!params.skip_mitogenome_assembly) {
        MITOGENOME_ASSEMBLY (
            ch_fastp_fastqc_results,
            organelle_type // params.organelle_type
        )
        ch_mitogenome_assembly_fasta = MITOGENOME_ASSEMBLY.out.assembly_fasta
        ch_mitogenome_assembly_log = MITOGENOME_ASSEMBLY.out.assembly_log
    } else if (params.precomputed_mitogenome_assembly_results) {
        // Use precomputed results if analysis is skipped
        ch_mitogenome_assembly_fasta = Channel.fromPath(params.precomputed_mitogenome_assembly_fasta, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without .fasta extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..2].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
        ch_mitogenome_assembly_log = Channel.fromPath(params.precomputed_mitogenome_assembly_log, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the filename (without extension)
            def filename = file.baseName  // Gets filename without extension
            def parts = filename.split('\\.')
            def meta_id = parts[0]
            def date = parts.length > 2 ? parts[2] : null
            def sample_id = parts.length > 2 ? parts[0..2].join('.') : filename
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                mt_assembly_prefix: sample_id
            ]
            return tuple(meta, file)
        }
    } else {
        ch_mitogenome_assembly_fasta = Channel.empty()
        ch_mitogenome_assembly_log = Channel.empty()
    }

    //
    // SUBWORKFLOW: MITOGENOME_ANNOTATION
    //
    if (!params.skip_mitogenome_annotation) {
        MITOGENOME_ANNOTATION (
            ch_mitogenome_assembly_fasta,
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
    // Combine outputs for data uploads
    //

    ch_mitogenome_assembly_results = ch_mitogenome_assembly_fasta.join(ch_mitogenome_assembly_log, by: 0)

    //
    // SUBWORKFLOW: UPLOAD_RESULTS
    //
    // Conditional uploading of results to SQL and species check - only run if not skipped
    // All these processes access the OceanOmics PostgreSQL database.
    if (!params.skip_upload_results) {
        UPLOAD_RESULTS (
            ch_mitogenome_assembly_results,
            ch_mitogenome_annotation_results,
            ch_mitogenome_blast_results,
            ch_mitogenome_lca_results,
            sql_config // params.sql_config
        )
    }


    // If the LCA validation is correct, then run the QC to prepare for submission to GenBank
    // Need to add this into the pipeline.
    // Now that protein lengths are being added to the database it could provide a list of 
    // non submitted mitogenomes they can be grouped with to submit and then say when there is a 
    // group of similar mitogenomes they can be submitted as a batch.
    // MITOGENOME_QC (

    //
    // Collect all MultiQC files from all subworkflows
    //

    // if (!params.skip_download_reads) {ch_multiqc_files = ch_multiqc_files.mix(DOWNLOAD_READS.out.multiqc_files)}
    // if (!params.skip_fastp_fastqc) {ch_multiqc_files = ch_multiqc_files.mix(FASTP_FASTQC.out.multiqc_files)}
    // if (!params.skip_genome_assembly) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_ASSEMBLY.out.multiqc_files)}
    // if (!params.skip_genome_decontamination) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_DECONTAMINATION.out.multiqc_files)}
    // if (!params.skip_genome_qc) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_QC.out.multiqc_files)}
    if (!params.skip_mitogenome_assembly) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ASSEMBLY.out.multiqc_files)}
    if (!params.skip_mitogenome_annotation) {ch_multiqc_files = ch_multiqc_files.mix(MITOGENOME_ANNOTATION.out.multiqc_files)}

    // 
    // Collect all versions
    //

    ch_collated_versions = Channel.empty()
    // if (!params.skip_download_reads) {ch_collated_versions = ch_collated_versions.mix(DOWNLOAD_READS.out.versions)}
    // if (!params.skip_fastp_fastqc) {ch_collated_versions = ch_collated_versions.mix(FASTP_FASTQC.out.versions)}
    // if (!params.skip_genome_assembly) {ch_collated_versions = ch_collated_versions.mix(GENOME_ASSEMBLY.out.versions)}
    // if (!params.skip_genome_decontamination) {ch_collated_versions = ch_collated_versions.mix(GENOME_DECONTAMINATION.out.versions)}
    // if (!params.skip_genome_qc) {ch_collated_versions = ch_collated_versions.mix(GENOME_QC.out.versions)}
    if (!params.skip_mitogenome_assembly) {ch_collated_versions = ch_collated_versions.mix(MITOGENOME_ASSEMBLY.out.versions)}
    if (!params.skip_mitogenome_annotation) {ch_collated_versions = ch_collated_versions.mix(MITOGENOME_ANNOTATION.out.versions)}

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

    
    //
    // Emit outputs
    //

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
        params.outdir,
        // params.input { optional true } // this is now inncliuded in the samplesheetHybrid wf
    )

    //
    // WORKFLOW: Run main workflow
    //
    OCEANGENOMES_MITOGENOMES (
        params.run,
        params.bs_config,
        params.curated_blast_db,
        params.sql_config,
        params.organelle_type
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
        OCEANGENOMES_MITOGENOMES.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
