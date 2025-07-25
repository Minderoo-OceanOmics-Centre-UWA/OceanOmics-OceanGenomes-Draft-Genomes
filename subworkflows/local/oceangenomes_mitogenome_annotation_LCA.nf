/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//mitogenome
include { DOWNLOAD_BLAST_DB } from '../../modules/local/download_blast_db/main'
include { DOWNLOAD_TAXONKIT_DB } from '../../modules/local/download_taxonkit_db/main'
include { EMMA } from '../../modules/local/EMMA/main'
// include { MITOZ } from '../../modules/local/mitoz/main'
include { BLAST_BLASTN } from '../../modules/nf-core/blast/blastn/main'
// include { BLAST_BLASTP } from '../../modules/nf-core/blast/blastp/main'
include { LCA } from '../../modules/local/LCA/main'
include { SPECIES_VALIDATION } from '../../modules/local/species_validation/main'
include { PUSH_MTDNA_ASSM_RESULTS } from '../../modules/local/upload_results/emma/main'
include { PUSH_LCA_BLAST_RESULTS } from '../../modules/local/upload_results/lca/main'
// Helper functions
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MITOGENOME WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MITOGENOME_ANNOTATION {

    take:
    mito_assembly //  tuple val(meta), path(fasta)
    curated_blast_db // params.curated_blast_db
    sql_config // params.sql_config
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Download taxonomy database
    DOWNLOAD_BLAST_DB(Channel.value("taxdb"))
    // Download taxonkit database
    DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))
    
    // Extracts the assembly name from the fasta file and embeds it into the meta map
    fasta_with_assembly_prefix = mito_assembly
        .map { meta, fasta ->
            def assembly_prefix = fasta.baseName
            def meta_ext = meta + [ assembly_prefix: assembly_prefix ]
            [meta_ext, fasta]
        }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ANNOTATION, EMMA FOR MAIN ANNOTATION AND MITOZ FOR COMPARISON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
               
    EMMA (
        fasta_with_assembly_prefix // tuple val(meta), path(fasta), val(assembly_prefix)
    )

    // // Function to extract assembly name
    def getAnnotationName = { filename ->
        def name = filename.toString().replaceAll(/\.fa$/, '')
        def parts = name.split('\\.', 2)
        return parts.size() > 1 ? parts[1] : name
    }

    // Use mix() to process Co1, 12s and 16s sequences through blast
    combined_sequences = EMMA.out.co1_sequences
        .map { meta, file -> 
            def annotation_name = getAnnotationName(file.name)
            [meta, file, 'CO1', annotation_name] 
        }
        .mix(
            EMMA.out.s12_sequences.map { meta, file -> 
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '12s', annotation_name] 
            },
            EMMA.out.s16_sequences.map { meta, file -> 
                def annotation_name = getAnnotationName(file.name)
                [meta, file, '16s', annotation_name] 
            }
        )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    USING CO1,12s and 16s RUN BLAST TO DETERMINE LCA FOR SPECIES VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    BLAST_BLASTN (
        combined_sequences, // tuple val(meta), path(fasta), val(gene_type), val(annotation_name)
        curated_blast_db,
        DOWNLOAD_BLAST_DB.out.db_files // path(db)
    )

    LCA (
        BLAST_BLASTN.out.filtered, // tuple val(meta), path(blast_filtered), val(gene_type), val(annotation_name)
        DOWNLOAD_TAXONKIT_DB.out.db_files // path(db)
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check the LCA results against nominal species ID and push results to SQL db
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    grouped_lca = LCA.out.lca
        .groupTuple(by: 0, size: 3)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }

    grouped_blast = BLAST_BLASTN.out.validation
        .groupTuple(by: 0, size: 3)
        .map { tuple ->
            def meta = tuple[0]
            def files = tuple[1..-1].flatten()
            [meta, files]
        }
    
    grouped_blast_lca = grouped_blast.join(grouped_lca)
    
    SPECIES_VALIDATION (
        grouped_blast_lca, // tuple val(meta), path(blast_filtered), path(lca_filtered)
        sql_config
    )

    PUSH_MTDNA_ASSM_RESULTS (
        EMMA.out.results, // tuple val(meta), path("emma/*")
        sql_config
    )

    PUSH_LCA_BLAST_RESULTS (
        SPECIES_VALIDATION.out.full, // tuple path ("lca_combined.${meta.id}.tsv"), path ("blast_combined.${meta.id}.txt"),
        sql_config
    )
    // Collect MultiQC files
    ch_multiqc_files = ch_multiqc_files.mix(BLAST_BLASTN.out.summary.collect{it[1]})
    ch_versions = ch_versions.mix(EMMA.out.versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions.first())
    ch_versions = ch_versions.mix(LCA.out.versions.first())



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
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
  


}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
