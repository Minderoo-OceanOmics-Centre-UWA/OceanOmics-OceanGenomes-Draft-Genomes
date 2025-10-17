/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Draft genome assembly subworkflows
include { FASTP_FASTQC              } from '../subworkflows/local/fastp_fastqc'
include { GENOME_ASSEMBLY           } from '../subworkflows/local/genome_assembly'
include { GENOME_DECONTAMINATION    } from '../subworkflows/local/genome_decontamination'
include { GENOME_QC                 } from '../subworkflows/local/genome_qc'
include { UPLOAD_RESULTS            } from '../subworkflows/local/upload_results'

// Modules
include { TRIGGER_MITOGENOME        } from '../modules/local/trigger_mitogenome'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { TAXON                     } from '../modules/local/taxon_from_db'

// Functions
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_oceangenomes_draftgenomes_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow OCEANGENOMES_DRAFTGENOMES {

    take:
        samplesheet
        sql_config // params.sql_config
        
    main:
    
    ch_multiqc_files = Channel.empty()
    
    /*
    The workflow uses conditional logic to handle steps that may be skipped in the nextflow_run script.
    If a process is not skipped, the corresponding subworkflow is executed and standard outputs are provided.
    If a process is skipped, predefined file paths (from nextflow.config, "precomputed_*") are used to create the required input for subsequent subworkflows.
    Metadata is extracted from file names, which are assumed to follow the format: sample_id.sequencing_technology.date.[additional_info].
    Mitogenome files are named with additional suffixes for each process (e.g., .getorganelle${version}, .emma${version}).
    If neither condition is met, an empty channel is created for the step.
    */

    //
    // SUBWORKFLOW: FASTP_FASTQC
    //

    if (!params.skip_fastp_fastqc) {
        FASTP_FASTQC (
            samplesheet
        )
        ch_fastp_fastqc_results = FASTP_FASTQC.out.fastp_reads // tuple val(meta), path('*.fastq.gz')
        ch_fastp_json = FASTP_FASTQC.out.fastp_json // tuple val(meta), path('*.fastq.gz')
    } else if (params.precomputed_fastp_fastqc_results) {
        // Use precomputed results if analysis is skipped
        ch_fastp_fastqc_results = Channel.fromFilePairs(params.precomputed_fastp_fastqc_results, checkIfExists: true)
            .map { reads_id, reads ->
                def parts = reads_id.split('\\.')
                def meta_id = parts[0]
                def date = parts[2]
                def sample_id = parts[0..2].join('.')  // Parts 0, 1, and 2 joined with dots
                def meta = [
                    id: meta_id,
                    run: params.run,
                    date: date,
                    prefix: sample_id
                ]
                return tuple(meta, reads)
            }

        ch_fastp_json = Channel.fromFilePairs(params.precomputed_fastp_json, checkIfExists: true)
            .map { json_id, json ->
                def parts = json_id.split('\\.')
                def meta_id = parts[0]
                def date = parts[2]
                def sample_id = parts[0..2].join('.')  // Parts 0, 1, and 2 joined with dots
                def meta = [
                    id: meta_id,
                    run: params.run,
                    date: date,
                    prefix: sample_id
                ]
                return tuple(meta, json)
            }
    } else {
        ch_fastp_fastqc_results = Channel.empty()
        ch_fastp_json = Channel.empty()
    }
    

    //
    // MODULE: TRIGGER_MITOGENOME
    //

    TRIGGER_MITOGENOME(
        ch_fastp_fastqc_results.collect()
    )


    //
    // SUBWORKFLOW: GENOME_ASSEMBLY
    //

    if (!params.skip_genome_assembly) {
        GENOME_ASSEMBLY (
            ch_fastp_fastqc_results,
            ch_fastp_json
        )
        ch_genome_assembly_results = GENOME_ASSEMBLY.out.megahit_assembled_contigs
        ch_meryl_db = GENOME_ASSEMBLY.out.meryl_db
        ch_genomescope_summary = GENOME_ASSEMBLY.out.genomescope_summary
    } else if (params.precomputed_genome_assembly_results) {
        // Use precomputed results if analysis is skipped
        ch_genome_assembly_results = Channel.fromPath(params.precomputed_genome_assembly_results, checkIfExists: true) 
        .map { file -> 
            def assembly_id = file.baseName 
            def parts = assembly_id.split('\\.') 
            def meta_id = parts[0] 
            def date = parts[2] 
            def sample_id = parts[0..2].join('.') // Parts 0, 1, and 2 joined with dots 
            def assembly_prefix = parts[0..3].join('.') 
            def meta = [ 
                id: meta_id, 
                run: params.run, // not using the params.run so if running files across multiple runs it doesnt matter what this value is. 
                date: date, 
                prefix: sample_id, 
                assembly_prefix: assembly_prefix 
            ] 
            return tuple(meta, file) 
        }

        ch_genomescope_summary = Channel.fromPath(params.precomputed_genomescope_summary, checkIfExists: true)
        .map { file ->
            // Extract sample_id from the file basename
            def file_name = file.baseName
            def sample_id = file_name.split('\\_')[0] // removes the '_genomescope2_summary' from the name
            def parts = sample_id.split('\\.')
            def meta_id = parts[0]
            def date = parts.size() > 2 ? parts[2] : params.run.tokenize('_')[1]
            
            def meta = [
                id: meta_id,
                run: params.run, // not using the params.run so if running files across multiple runs it doesnt matter what this value is.
                date: date,
                prefix: sample_id
            ]
            
            return tuple(meta, file)
        }

        ch_meryl_db = Channel.fromPath(params.precomputed_meryl_results, type: 'dir', checkIfExists: true)
        .map { meryl_dir ->
            // Extract sample_id from the .meryl directory name (the * before .meryl)
            def sample_id = meryl_dir.name.replace('.meryl', '')
            def parts = sample_id.split('\\.')
            def meta_id = parts[0]
            def date = parts.size() > 2 ? parts[2] : params.run.tokenize('_')[1]
            
            def meta = [
                id: meta_id,
                run: params.run, // not using the params.run so if running files across multiple runs it doesnt matter what this value is.
                date: date,
                prefix: sample_id
            ]
            
            return tuple(meta, meryl_dir)
        }

    } else {
        ch_genome_assembly_results = Channel.empty()
        ch_meryl_db = Channel.empty()
    }


    //
    // MODULE: TAXON - Run once if either decontamination or QC is not skipped
    //

    if (!params.skip_genome_decontamination || !params.skip_genome_qc) {
    // Determine which assembly data to use for taxon lookup
    ch_assembly_for_taxon = !params.skip_genome_decontamination ? 
        ch_genome_assembly_results : 
        Channel.fromPath(params.precomputed_genome_decontamination_results ?: [], checkIfExists: true)
                .map { file ->
                    def assembly_id = file.baseName
                    def parts = assembly_id.split('\\.')
                    def meta_id = parts[0]
                    def date = parts[2]
                    def sample_id = parts[0..2].join('.')
                    def assembly_prefix = parts[0..3].join('.')
                    def meta = [
                        id: meta_id,
                        run: params.run,
                        date: date,
                        prefix: sample_id,
                        assembly_prefix: assembly_prefix
                    ]
                    return tuple(meta, file)
                }

        TAXON (
            ch_assembly_for_taxon,
            params.sql_config
        )
        .map { meta, fasta, taxon_csv_file ->
            def taxon_row = taxon_csv_file
                .splitCsv(header: true)
                .first()

            def updated_meta = meta + [
                nom_species_id: taxon_row.nominal_species_id,
                taxon_id: taxon_row.taxon_id,
                class   : taxon_row.class
            ]

            tuple(updated_meta, fasta)
        }
        .set { ch_assembly_with_taxon }

        // Error handling for missing taxon_id
        ch_assembly_with_taxon.map { meta, assembly_file ->
            if (!meta.taxon_id) error "❗ taxon_id not found for sample ${meta.id}"
            tuple(meta, assembly_file)
        }
    } else {
        ch_assembly_with_taxon = Channel.empty()
    }


    //
    // SUBWORKFLOW: GENOME_DECONTAMINATION
    //

    if (!params.skip_genome_decontamination) {
        GENOME_DECONTAMINATION (
            ch_assembly_with_taxon
        )
        ch_genome_decontamination_assembly = GENOME_DECONTAMINATION.out.decontamined_assembled_reads
        ch_filter_report = GENOME_DECONTAMINATION.out.filter_report
        ch_contigs_under_500bp = GENOME_DECONTAMINATION.out.contigs_under_500bp
        ch_tiara_filter_summary = GENOME_DECONTAMINATION.out.tiara_filter_summary
    } else if (params.precomputed_genome_decontamination_results) {
        // If decontamination is skipped, use the taxon-enriched precomputed data
        // that was already processed in the TAXON section above
        ch_genome_decontamination_assembly = ch_assembly_with_taxon
        ch_filter_report = Channel.fromPath(params.precomputed_filter_report, checkIfExists: true)
        ch_contigs_under_500bp = Channel.fromPath(params.precomputed_contigs_under_500bp, checkIfExists: true)
        ch_tiara_filter_summary = Channel.fromPath(params.precomputed_tiara_filter_summary, checkIfExists: true)
    } else {
        ch_genome_decontamination_assembly = Channel.empty()
        ch_filter_report = Channel.empty()
        ch_contigs_under_500bp = Channel.empty()
        ch_tiara_filter_summary = Channel.empty()
    }


    //
    // SUBWORKFLOW: GENOME_QC
    //

    if (!params.skip_genome_qc) {
        GENOME_QC (
            ch_fastp_fastqc_results,
            ch_genomescope_summary,
            ch_genome_decontamination_assembly,
            ch_meryl_db,
            
        )
        ch_busoco_short_summary = GENOME_QC.out.busco_short_summary
        ch_merqury_results = GENOME_QC.out.merqury_results
        ch_gfastats_results = GENOME_QC.out.gfastats_results
    } else if (params.precomputed_genome_qc_results) {
        // Use precomputed results if analysis is skipped
        /* For the genome QC there is going to be multiple outputs of the results from the different modules run within this subworkflow.
            These results will need to be pushed up to the SQL database. So more modules could be added to the upload subworkflow for these.
            Or have a seperate draft genome results subwokflow which is probably better, to keep the mitogenome pipeline seperate.
            Will need to add in all the outputs and the if skipped paths to files evenbtually.
        */
        ch_busoco_short_summary = Channel.fromPath(params.precomputed_busoco_short_summary_results, checkIfExists: true)
        ch_merqury_results = Channel.fromPath(params.precomputed_merqury_results_results, checkIfExists: true)
        ch_gfastats_results = Channel.fromPath(params.precomputed_gfastats_results_results, checkIfExists: true)
    } else {
        ch_busoco_short_summary = Channel.empty()
        ch_merqury_results = Channel.empty()
        ch_gfastats_results = Channel.empty()
    }


    //
    // SUBWORKFLOW: Upload results to SQL database
    //

    // UPLOAD_RESULTS (
    //     ch_fastp_json,
    //     ch_genomescope_summary,
    //     ch_filter_report,
    //     ch_contigs_under_500bp,
    //     ch_tiara_filter_summary,
    //     ch_busoco_short_summary,
    //     ch_merqury_results,
    //     ch_gfastats_results,
    //     sql_config
    // )

    //
    // Collect all MultiQC files from all subworkflows
    //

    if (!params.skip_fastp_fastqc) {ch_multiqc_files = ch_multiqc_files.mix(FASTP_FASTQC.out.multiqc_files)}
    if (!params.skip_genome_assembly) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_ASSEMBLY.out.multiqc_files)}
    if (!params.skip_genome_decontamination) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_DECONTAMINATION.out.multiqc_files)}
    if (!params.skip_genome_qc) {ch_multiqc_files = ch_multiqc_files.mix(GENOME_QC.out.multiqc_files)}
    // if (!params.skip_upload_results) {ch_multiqc_files = ch_multiqc_files.mix(UPLOAD_RESULTS.out.multiqc_files)}

    // 
    // Collect all versions
    //

    ch_collated_versions = Channel.empty()
    if (!params.skip_fastp_fastqc) {ch_collated_versions = ch_collated_versions.mix(FASTP_FASTQC.out.versions)}
    if (!params.skip_genome_assembly) {ch_collated_versions = ch_collated_versions.mix(GENOME_ASSEMBLY.out.versions)}
    if (!params.skip_genome_decontamination) {ch_collated_versions = ch_collated_versions.mix(GENOME_DECONTAMINATION.out.versions)}
    if (!params.skip_genome_qc) {ch_collated_versions = ch_collated_versions.mix(GENOME_QC.out.versions)}
    // if (!params.skip_upload_results) {ch_collated_versions = ch_collated_versions.mix(UPLOAD_RESULTS.out.versions)}

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
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
}
