import java.nio.file.Path
import java.nio.file.Files
import java.nio.file.StandardOpenOption

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
                nom_species_id: row.nom_species_id,
                taxon_id : row.taxon_id,
                class : row.class
            ]

            def reads = [
                file(row.fastq_1),
                file(row.fastq_2)
            ]

            tuple(meta, reads)
        }

    def csvColumns = [
        'sample',
        'run',
        'date',
        'prefix',
        'nom_species_id',
        'taxon_id',
        'class',
        'fastq_1',
        'fastq_2'
    ]
    def csvHeader = csvColumns.join(',')
    def repairedSamplesheetFile
    def repairedSamplesheetLock = new Object()

    def escapeCsvValue = { value ->
        def str = value ? value.toString() : ''
        def escaped = str.replace('"', '""')
        def needsQuotes = escaped.contains(',') || escaped.contains('"') || escaped.contains('\n') || escaped.contains('\r')
        needsQuotes ? "\"${escaped}\"" : escaped
    }

    def ensureRepairedSamplesheetFile = {
        if (repairedSamplesheetFile) {
            return repairedSamplesheetFile
        }

        def samplesheetFileName = params.samplesheet_prefix ?: 'samplesheet'
        if (!samplesheetFileName.endsWith('.csv')) {
            samplesheetFileName = samplesheetFileName + '.csv'
        }

        Path outDirPath = file(params.outdir ?: '.')
        outDirPath?.toFile()?.mkdirs()
        Path filePath = outDirPath.resolve(samplesheetFileName)

        Files.write(
            filePath,
            (csvHeader + System.lineSeparator()).getBytes('UTF-8'),
            StandardOpenOption.CREATE,
            StandardOpenOption.TRUNCATE_EXISTING
        )
        log.info "📝 Writing repaired samplesheet entries to: ${filePath}"
        repairedSamplesheetFile = filePath
        return repairedSamplesheetFile
    }

    def appendRepairedSamplesheetRow = { meta, repaired_files ->
        def fastq1 = repaired_files && repaired_files.size() > 0 ? repaired_files[0].toString() : ''
        def fastq2 = repaired_files && repaired_files.size() > 1 ? repaired_files[1].toString() : ''

        def row = [
            escapeCsvValue(meta.id),
            escapeCsvValue(meta.run),
            escapeCsvValue(meta.date),
            escapeCsvValue(meta.prefix),
            escapeCsvValue(meta.nom_species_id ?: ''),
            escapeCsvValue(meta.taxon_id ?: ''),
            escapeCsvValue(meta.class ?: ''),
            escapeCsvValue(fastq1),
            escapeCsvValue(fastq2)
        ].join(',') + System.lineSeparator()

        synchronized(repairedSamplesheetLock) {
            Files.write(
                ensureRepairedSamplesheetFile(),
                row.getBytes('UTF-8'),
                StandardOpenOption.CREATE,
                StandardOpenOption.APPEND
            )
        }
    }

    // Process repaired_ch when no samplesheet
    ch_meta_from_repaired = decision.no_samplesheet
        .combine(repaired_ch.ifEmpty { error "❗ No samplesheet provided and no repaired_ch input given. Cannot continue." })
        .map { count, sample_id, repaired_files ->
            log.info "🔍 No samplesheet provided, using repaired output to build samplesheet..."
            log.info "DEBUG Processing sample_id: ${sample_id}"
            log.info "DEBUG repaired_files[0]: ${repaired_files[0]}"
            log.info "DEBUG type: ${repaired_files[0].class}"
            
            def meta_id = sample_id.split('\\.')[0] // Extract meta_id which is the first part of sample name seperated by .
            def date = params.run.tokenize('_')[1] // Extract date from run parameter

            def meta = [
                id: meta_id,
                run: params.run,
                date: date,
                prefix: "${meta_id}.ilmn.${date}"
            ]

            tuple(meta, repaired_files)
        }
    
    TAXON (
        ch_meta_from_repaired,
        params.sql_config
    )
    .map { meta, repaired_files, taxon_csv_file ->
        def taxon_rows = []
        if (Files.exists(taxon_csv_file) && Files.size(taxon_csv_file) > 0) {
            taxon_rows = taxon_csv_file.splitCsv(header: true)
        }

        if (!taxon_rows || taxon_rows.isEmpty()) {
            log.warn "⚠️ No taxon data returned for ${meta.id}; defaulting metadata fields to 'unknown'"
            taxon_rows = [[nominal_species_id: 'unknown', taxon_id: 'unknown', class: 'unknown']]
        }

        def taxon_row = taxon_rows.first()

        def updated_meta = meta + [
            nom_species_id: taxon_row.nominal_species_id ?: 'unknown',
            taxon_id: taxon_row.taxon_id ?: 'unknown',
            class   : taxon_row.class ?: 'unknown'
        ]
        
        appendRepairedSamplesheetRow(updated_meta, repaired_files)
        tuple(updated_meta, repaired_files)
    }
    .set { ch_samplesheet_from_repaired }


    // Combine both channels (only one will have data)
    ch_samplesheet = ch_samplesheet_from_csv.mix(ch_samplesheet_from_repaired)

    emit:
    samplesheet = ch_samplesheet
}
