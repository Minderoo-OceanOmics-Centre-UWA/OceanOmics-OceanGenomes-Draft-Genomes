/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data upload modules
include { PUSH_FASTP_RESULTS            } from '../../../modules/local/upload_results/fastp'
include { PUSH_ASSEMBLY_RESULTS         } from '../../../modules/local/upload_results/assembly'
include { PUSH_DECONTAMINATION_RESULTS  } from '../../../modules/local/upload_results/decontamination'
include { PUSH_BUSCO_RESULTS            } from '../../../modules/local/upload_results/busco'
include { PUSH_MERQURY_RESULTS          } from '../../../modules/local/upload_results/merqury'
include { PUSH_GFA_RESULTS              } from '../../../modules/local/upload_results/gfa'

// Helper functions
include { softwareVersionsToYAML        } from '../../nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DATE UPLOAD AND SPECIES CHECK WORKFLOW

        - Specific OceanOmics code using the PostgreSQL database.
        - Check the LCA results against nominal species ID and push results to SQL db
        - Determines species that have validated species ID to proceed with QC to prepare
          the sample for submittion to Genbank.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPLOAD_RESULTS {

    take:
    fastp_json // from FASTP_FASTQC.out.json - fastp workflow
    genomescope_summary // from GENOMESCOPE2.out.summary - assembly workflow
    filter_report // from BBMAP_FILTERBYNAME.out.filter_report - decontamination workflow
    contigs_under_500bp // from BBMAP_FILTERBYNAME.out.contigs_under_500bp - decontamination workflow
    tiara_filter_summary // from TIARA_TIARA.out.summary - decontamination workflow
    busco_short_summary // from BUSCO_BUSCO.out.short_summaries_json - genome qc workflow
    merqury_results // from MERQURY_MERQURY.out.stats and MERQURY_MERQURY.out.assembly_qv - genome qc workflow
    gfastats_results // from GFASTATS.out.assembly_summary - genome qc workflow
    sql_config // params.sql_config
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 
    //
    // MODULE: Fastp results upload to SQL database
    //

    PUSH_FASTP_RESULTS (
        fastp_json, // tuple val(meta), path('*.json')
        params.sql_config
    )

    
    //
    // MODULE: Assembly results upload to SQL database
    //

    PUSH_ASSEMBLY_RESULTS (
        genomescope_summary, // tuple val(meta), path(blast_filtered), path(lca_filtered)
        sql_config // params.sql_config
    )

    //
    // MODULE: Decontamination results upload to SQL database
    //
    ch_decontamination = filter_report
        .join(contigs_under_500bp, by: 0)
        .join(tiara_filter_summary, by: 0)

    PUSH_DECONTAMINATION_RESULTS (
        ch_decontamination, // tuple val(meta), path(filter_report), path(contigs_under_500bp), path(tiara_filter_summary)
        sql_config // params.sql_config
    )

    //
    // MODULE: Updating the SQL database with the LCA and filtered BLAST results
    //

    PUSH_BUSCO_RESULTS (
        busco_short_summary, // tuple val(meta), path('*.busco.short_summary.json')
        sql_config // params.sql_config
    )


    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    PUSH_MERQURY_RESULTS (
        merqury_results, // tuple val(meta), path("*.completeness.stats"), path("${prefix}.qv")
        sql_config // params.sql_config
    )


    //
    // MODULE: evaluating the results to determine if to process the sample through QC
    //

    PUSH_GFA_RESULTS (
        gfastats_results, // tuple val(meta), path("*.assembly_summary")
        sql_config // params.sql_config
    )

   
    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    //  - Species validation outputs (per-sample TSVs)
    //  - Annotation statistics CSVs
    //  - QC evaluation flags (for quick visibility in report) and summary table
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_FASTP_RESULTS.out.upload.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_ASSEMBLY_RESULTS.out.upload.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_DECONTAMINATION_RESULTS.out.upload.collect { it[1] })     
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_BUSCO_RESULTS.out.upload.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_MERQURY_RESULTS.out.upload.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(PUSH_GFA_RESULTS.out.upload.collect { it[1] })
    // Collect versions
    ch_versions = ch_versions.mix(PUSH_FASTP_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_ASSEMBLY_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_DECONTAMINATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_BUSCO_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MERQURY_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_GFA_RESULTS.out.versions.first())

    //
    // Emit outputs
    //

    emit:
    multiqc_files = ch_multiqc_files            // channel: [ path(multiqc_files) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
}
