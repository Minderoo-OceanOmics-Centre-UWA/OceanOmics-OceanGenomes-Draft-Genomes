//
// Subworkflow with to produce samplesheet for fastq/fastp processes specific to the nf-core/oceangenomes_draftgenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow samplesheetHybrid {

    take:
    repaired_ch { optional true }

    main:
    ch_samplesheet = Channel.empty()  // ✅ Safe placeholder declaration
    
    // If providing an sample sheet as csv format it requires the following format: sample,run,R1,R2
    // alternatively you can also provide the date of the sequencing run if your format isnt XXX_$DATE_XXX: sample,run,date,R1,R2
    // R1 and R2 are the file paths to R1 and R2
    if (params.input) {
        log.info "📋 Using samplesheet CSV: ${params.input}"

        ch_samplesheet = Channel
            .fromPath(params.input)
            .splitCsv(header: true, schema: file(params.schema_input))
            .map { row ->
                def meta = [
                    id   : row.sample,
                    run  : row.run,
                    date : row.date ?: row.run.tokenize('_')[1]
                ]

                def reads = [
                    file(row.R1),
                    file(row.R2)
                ]

                tuple(meta, reads)
            }

    } else if (repaired_ch) {
        log.info "🔍 No samplesheet provided, using repaired output to build samplesheet..."
        repaired_ch.view()
        

        ch_samplesheet = repaired_ch
            .map { sample_id, repaired_files ->
                log.info "DEBUG repaired_files[0]: ${repaired_files[0]}"
                log.info "DEBUG type: ${repaired_files[0].class}"
                def run_match = repaired_files[0].getName() =~ /ilmn\.(.+?)\.R1/
                def run_id = run_match ? run_match[0][1] : 'UNKNOWN_RUN'
                def date = run_id.tokenize('_').size() > 1 ? run_id.tokenize('_')[1] : '000000'

                def meta = [
                    id   : sample_id,
                    run  : run_id,
                    date : date
                ]

                tuple(meta, repaired_files)
            }

    } else {
        throw new RuntimeException("❗ No samplesheet provided and no repaired_ch input given. Cannot continue.")
    }

    emit:
    samplesheet = ch_samplesheet
}
