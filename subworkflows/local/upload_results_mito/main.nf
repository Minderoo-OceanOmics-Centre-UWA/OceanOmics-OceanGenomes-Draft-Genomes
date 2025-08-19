/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Data upload modules
include { PUSH_MTDNA_ASSM_RESULTS       } from '../../../modules/local/upload_results/mtdna'
include { SPECIES_VALIDATION            } from '../../../modules/local/species_validation'
include { PUSH_MTDNA_ANNOTATION_RESULTS } from '../../../modules/local/upload_results/emma'
include { PUSH_LCA_BLAST_RESULTS        } from '../../../modules/local/upload_results/lca'

// Helper functions
include { softwareVersionsToYAML        } from '../../nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DATE UPLOAD AND SPECIES CHECK WORKFLOW

        Specific OceanOmics code using the PostgreSQL database.
        Check the LCA results against nominal species ID and push results to SQL db

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow UPLOAD_RESULTS {

    take:
    assembly_results
    annotation_results
    blast_filtered_results
    lca_results
    sql_config // params.sql_config
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 
    //
    // Pulling out the statistics from the assembly files and updating the SQL database
    //

    PUSH_MTDNA_ASSM_RESULTS (
        assembly_results,
        params.sql_config
    )

    //
    // Re grouping the CO1, 12s and 16s BLAST and LCA results per sample
    //

    grouped_lca = lca_results
        .groupTuple(by: 0, size: 3)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }

    grouped_blast = blast_filtered_results
        .groupTuple(by: 0, size: 3)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }
    
    grouped_blast_lca = grouped_blast.join(grouped_lca, by: 0)
    
    //
    // Checking the LCA results against the nominal species ID in the SQL database
    //

    SPECIES_VALIDATION (
        grouped_blast_lca, // tuple val(meta), path(blast_filtered), path(lca_filtered)
        sql_config // params.sql_config
    )

    //
    // Calculating the statistics of the annotations and updating the SQL database
    //

    PUSH_MTDNA_ANNOTATION_RESULTS (
        annotation_results, // tuple val(meta), path("emma/*")
        sql_config // params.sql_config
    )

    //
    // Updating the SQL database with the LCA and filtered BLAST results
    //

    PUSH_LCA_BLAST_RESULTS (
        SPECIES_VALIDATION.out.full, // tuple path ("lca_combined.${meta.id}.tsv"), path ("blast_combined.${meta.id}.txt"),
        sql_config // params.sql_config
    )

    
    //
    // Subworkflow finishing steps.
    //

    // Collect MultiQC files
    // Need to update this section to include everything
    // ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    ch_versions = ch_versions.mix(PUSH_MTDNA_ASSM_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(SPECIES_VALIDATION.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_MTDNA_ANNOTATION_RESULTS.out.versions.first())
    ch_versions = ch_versions.mix(PUSH_LCA_BLAST_RESULTS.out.versions.first())



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
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
}
