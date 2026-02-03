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
    ch_meta_by_prefix = samplesheet
        .map { meta, reads -> tuple(meta.prefix ?: meta.id, meta) }
        .distinct()
    
    def warnIfEmpty = { ch, label ->
        ch.ifEmpty {
            log.warn "No files match pattern `${label}`; continuing with empty channel."
            null
        }.filter { it != null }
    }


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
        ch_precomputed_fastp_fastqc = warnIfEmpty(Channel.fromFilePairs(params.precomputed_fastp_fastqc_results, checkIfExists: false), params.precomputed_fastp_fastqc_results)
            .map { reads_id, reads ->
                def parts = reads_id.split('\\.')
                def sample_id = parts[0..2].join('.')
                def nova_match = reads_id =~ /NOVA_(\d{6})_/
                if (nova_match.find()) {
                    sample_id = "${parts[0]}.${parts[1]}.${nova_match.group(1)}"
                }
                tuple(sample_id, reads)
            }
        ch_fastp_fastqc_results = ch_meta_by_prefix
            .join(ch_precomputed_fastp_fastqc)
            .map { key, meta, reads -> tuple(meta, reads) }

        ch_precomputed_fastp_json = warnIfEmpty(Channel.fromPath(params.precomputed_fastp_json, checkIfExists: false), params.precomputed_fastp_json)
            .map { json ->
                def parts = json.baseName.split('\\.')
                def sample_id = parts[0..2].join('.')
                def nova_match = json.baseName =~ /NOVA_(\d{6})_/
                if (nova_match.find()) {
                    sample_id = "${parts[0]}.${parts[1]}.${nova_match.group(1)}"
                }
                tuple(sample_id, json)
            }
        ch_fastp_json = ch_meta_by_prefix
            .join(ch_precomputed_fastp_json)
            .map { key, meta, json -> tuple(meta, json) }
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
        ch_precomputed_genome_assembly_results = warnIfEmpty(Channel.fromPath(params.precomputed_genome_assembly_results, checkIfExists: false), params.precomputed_genome_assembly_results)
        .map { file ->
            def assembly_id = file.baseName
            def parts = assembly_id.split('\\.')
            def sample_id = parts[0..2].join('.')
            def assembly_prefix = parts[0..3].join('.')
            tuple(sample_id, assembly_prefix, file)
        }
        ch_genome_assembly_results = ch_meta_by_prefix
            .join(ch_precomputed_genome_assembly_results)
            .map { key, meta, assembly_prefix, file ->
                tuple(meta + [assembly_prefix: assembly_prefix], file)
            }

        ch_precomputed_genomescope_summary = warnIfEmpty(Channel.fromPath(params.precomputed_genomescope_summary, checkIfExists: false), params.precomputed_genomescope_summary)
        .map { file ->
            // Extract sample_id from the file basename
            def file_name = file.baseName
            def sample_id = file_name.split('\\_')[0] // removes the '_genomescope2_summary' from the name
            tuple(sample_id, file)
        }
        ch_genomescope_summary = ch_meta_by_prefix
            .join(ch_precomputed_genomescope_summary)
            .map { key, meta, file -> tuple(meta, file) }

        ch_precomputed_meryl_db = warnIfEmpty(Channel.fromPath(params.precomputed_meryl_results, type: 'dir', checkIfExists: false), params.precomputed_meryl_results)
        .map { meryl_dir ->
            // Extract sample_id from the .meryl directory name (the * before .meryl)
            def sample_id = meryl_dir.name.replace('.meryl', '')
            tuple(sample_id, meryl_dir)
        }
        ch_meryl_db = ch_meta_by_prefix
            .join(ch_precomputed_meryl_db)
            .map { key, meta, meryl_dir -> tuple(meta, meryl_dir) }

    } else {
        ch_genome_assembly_results = Channel.empty()
        ch_meryl_db = Channel.empty()
        ch_genomescope_summary = Channel.empty()
    }

    ch_assembly_with_taxon = ch_genome_assembly_results.map { meta, assembly_file ->
        if (!meta.taxon_id) error "❗ taxon_id not found for sample ${meta.id}"
        tuple(meta, assembly_file)
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
        ch_precomputed_genome_decontamination_results = warnIfEmpty(Channel.fromPath(params.precomputed_genome_decontamination_results, checkIfExists: false), params.precomputed_genome_decontamination_results)
            .map { file ->
                def assembly_id = file.baseName
                def parts = assembly_id.split('\\.')
                def sample_id = parts[0..2].join('.')
                def assembly_prefix = parts.size() > 3 ? parts[0..3].join('.') : sample_id
                tuple(sample_id, assembly_prefix, file)
            }
        ch_genome_decontamination_assembly = ch_meta_by_prefix
            .join(ch_precomputed_genome_decontamination_results)
            .map { key, meta, assembly_prefix, file ->
                tuple(meta + [assembly_prefix: assembly_prefix], file)
            }
        ch_precomputed_filter_report = warnIfEmpty(Channel.fromPath(params.precomputed_filter_report, checkIfExists: false), params.precomputed_filter_report)
            .map { file ->
                def file_name = file.baseName
                def sample_id = file_name.split('\\.')[0..2].join('.')  // Parts 0, 1, and 2 joined with dots
                tuple(sample_id, file)
            }
        ch_filter_report = ch_meta_by_prefix
            .join(ch_precomputed_filter_report)
            .map { key, meta, file -> tuple(meta, file) }

        ch_precomputed_contigs_under_500bp = warnIfEmpty(Channel.fromPath(params.precomputed_contigs_under_500bp, checkIfExists: false), params.precomputed_contigs_under_500bp)
            .map { file ->
                def file_name = file.baseName
                def sample_id = file_name.split('\\.')[0..2].join('.')  // Parts 0, 1, and 2 joined with dots
                tuple(sample_id, file)
            }
        ch_contigs_under_500bp = ch_meta_by_prefix
            .join(ch_precomputed_contigs_under_500bp)
            .map { key, meta, file -> tuple(meta, file) }

        ch_precomputed_tiara_filter_summary = warnIfEmpty(Channel.fromPath(params.precomputed_tiara_filter_summary, checkIfExists: false), params.precomputed_tiara_filter_summary)
            .map { file ->
                def file_name = file.baseName
                def sample_id = file_name.split('\\.')[0..2].join('.')  // Parts 0, 1, and 2 joined with dots
                tuple(sample_id, file)
            }
        ch_tiara_filter_summary = ch_meta_by_prefix
            .join(ch_precomputed_tiara_filter_summary)
            .map { key, meta, file -> tuple(meta, file) }
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
        ch_precomputed_busoco_short_summary = warnIfEmpty(Channel.fromPath(params.precomputed_busoco_short_summary_results, checkIfExists: false), params.precomputed_busoco_short_summary_results)
            .map { file ->
                def sample_id = file.baseName.split('\\.')[0..2].join('.')
                tuple(sample_id, file)
            }
        ch_busoco_short_summary = ch_meta_by_prefix
            .join(ch_precomputed_busoco_short_summary)
            .map { key, meta, file -> tuple(meta, file) }

        ch_precomputed_merqury_results = warnIfEmpty(Channel.fromPath(params.precomputed_merqury_results_results, checkIfExists: false), params.precomputed_merqury_results_results)
            .map { file ->
                def sample_id = file.baseName.split('\\.')[0..2].join('.')
                tuple(sample_id, file)
            }
            .groupTuple()
            .map { sample_id, files ->
                def completeness_stats = files.find { it.name.endsWith('completeness.stats') }
                def qv_tsv = files.find { it.name.endsWith('.qv') }
                tuple(sample_id, completeness_stats, qv_tsv)
            }
            .filter { sample_id, completeness_stats, qv_tsv -> completeness_stats && qv_tsv }
        ch_merqury_results = ch_meta_by_prefix
            .join(ch_precomputed_merqury_results)
            .map { key, meta, completeness_stats, qv_tsv -> tuple(meta, completeness_stats, qv_tsv) }

        ch_precomputed_gfastats_results = warnIfEmpty(Channel.fromPath(params.precomputed_gfastats_results_results, checkIfExists: false), params.precomputed_gfastats_results_results)
            .map { file ->
                def sample_id = file.baseName.split('\\.')[0..2].join('.')
                tuple(sample_id, file)
            }
        ch_gfastats_results = ch_meta_by_prefix
            .join(ch_precomputed_gfastats_results)
            .map { key, meta, file -> tuple(meta, file) }
    } else {
        ch_busoco_short_summary = Channel.empty()
        ch_merqury_results = Channel.empty()
        ch_gfastats_results = Channel.empty()
    }


    //
    // SUBWORKFLOW: Upload results to SQL database
    //

    if (!params.skip_upload_results) {
        UPLOAD_RESULTS (
            ch_fastp_json,
            ch_genomescope_summary,
            ch_filter_report,
            ch_contigs_under_500bp,
            ch_tiara_filter_summary,
            ch_busoco_short_summary,
            ch_merqury_results,
            ch_gfastats_results,
            sql_config
        )
    }

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
