/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Helper functions
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

//mitogenome
include { GETORGANELLE_CONFIG } from '../../modules/nf-core/getorganelle/config/main'
include { GETORGANELLE_FROMREADS } from '../../modules/nf-core/getorganelle/fromreads/main'
include { PUSH_MTDNA_ASSM_RESULTS } from '../../modules/local/upload_results/mtdna/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOMES {

    take:
    fastp_reads // tuple val(meta), path(fastq)
    organelle_type
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

     
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ASSEMBLY USING GET ORGANELLE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    GETORGANELLE_CONFIG (
        organelle_type // val(organelle_type)
    )

    version_ch = GETORGANELLE_CONFIG.out.versions
    .map { versions_file ->
        def content = versions_file.text
        def pattern = /getorganelle:\s*([^\s\n\r]+)/
        def matcher = content =~ pattern
        return matcher ? matcher[0][1].replaceAll(/["']/, '') : "unknown"
    }

    // Could maybe add in a function to make the assembly name to pass into getorganelle so i dont have to make it in the module

    
    combined_input = fastp_reads.combine(GETORGANELLE_CONFIG.out.db).combine(version_ch)

    GETORGANELLE_FROMREADS (
        combined_input // tuple val(meta), path(fastq), val(organelle_type), path(db), val(version)  // getOrganelle has a database and config file
    )

    log_ch = GETORGANELLE_FROMREADS.out.fasta.join(GETORGANELLE_FROMREADS.out.log)
        .map { meta, fasta, log -> 
            def assembly_prefix = fasta.baseName
            [meta, fasta, log, assembly_prefix] }

    PUSH_MTDNA_ASSM_RESULTS (
        log_ch,
        params.sql_config
    )

    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(GETORGANELLE_FROMREADS.out.etc.collect{it[1]})
    ch_versions = ch_versions.mix(GETORGANELLE_CONFIG.out.versions.first())
    ch_versions = ch_versions.mix(GETORGANELLE_FROMREADS.out.versions.first())


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
    mito_assembly = GETORGANELLE_FROMREADS.out.fasta
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]


}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
