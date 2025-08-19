//
// Subworkflow with to produce samplesheet for fastq/fastp processes specific to the nf-core/oceangenomes_draftgenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Taxon module, specific to OceanOmics DQL database
include { TAXON                     } from '../../../modules/local/taxon_from_db'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow samplesheetHybrid {

    take:
    repaired_ch
    samplesheet

    main:
    
    // Count items in samplesheet to determine if it's empty
    samplesheet
        .count()
        .set { samplesheet_count }

    // Branch based on whether samplesheet has content
    samplesheet_count
        .branch { count ->
            has_samplesheet: count > 0
            no_samplesheet: count == 0
        }
        .set { decision }

    // Process CSV samplesheet when available
    ch_samplesheet_from_csv = decision.has_samplesheet
        .combine(samplesheet)
        .map { count, samplesheet_content -> samplesheet_content }
        .splitCsv(header: true)
        .map { row ->
            log.info "📋 Using samplesheet CSV for sample: ${row.sample}"
            
            def meta = [
                id   : row.sample,
                run  : row.run,
                date : row.date ?: row.run.tokenize('_')[1],
                prefix : row.prefix,
                taxon : row.taxon,
                class : row.class
            ]

            def reads = [
                file(row.R1),
                file(row.R2)
            ]

            tuple(meta, reads)
        }

    // Process repaired_ch when no samplesheet
    ch_samplesheet_from_repaired = decision.no_samplesheet
        .combine(repaired_ch.ifEmpty { error "❗ No samplesheet provided and no repaired_ch input given. Cannot continue." })
        .map { count, sample_id, repaired_files ->
            log.info "🔍 No samplesheet provided, using repaired output to build samplesheet..."
            log.info "DEBUG Processing sample_id: ${sample_id}"
            log.info "DEBUG repaired_files[0]: ${repaired_files[0]}"
            log.info "DEBUG type: ${repaired_files[0].class}"
            
            def meta_id = sample_id.split('\\.')[0] // Extract meta_id which is the first part of sample name seperated by .
            
            def meta = [
                id: meta_id,
                run: params.run,
                date: params.date,
                prefix: "${meta_id}.ilmn.${params.date}"
            ]

            tuple(meta, repaired_files)
        }

    //
    // MODULE: Retrieve NCBI taxon ID and taxonomic class from OceanOmics SQL database.
    //      Taxon ID and class can be provided in sample sheet if running outside of OceanOmics
    //

    TAXON (
        ch_samplesheet_from_repaired, 
        params.sql_config
    )
    .map { meta, repaired_files, taxon_csv_file ->
        def taxon_row = taxon_csv_file
            .splitCsv(header: true)
            .first()

        def updated_meta = meta + [
            nom_species_id: taxon_row.nominal_species_id,
            taxon_id: taxon_row.taxon_id,
            class   : taxon_row.class
        ]

        // 🐛 DEBUG print to log
        log.info "🔍 Updated meta for ${updated_meta.id}: ${updated_meta}"

        tuple(updated_meta, repaired_files)
    }
    .set { ch_samplesheet_from_repaired_with_taxon }

    // Error code for if no taxon_id found
    ch_samplesheet_from_repaired_with_taxon.map { meta, repaired_files ->
        if (!meta.taxon_id) error "❗ taxon_id not found for sample ${meta.id}"
        tuple(meta, repaired_files)
    }


    // Combine both channels (only one will have data)
    ch_samplesheet = ch_samplesheet_from_csv.mix(ch_samplesheet_from_repaired_with_taxon)

    emit:
    samplesheet = ch_samplesheet
}