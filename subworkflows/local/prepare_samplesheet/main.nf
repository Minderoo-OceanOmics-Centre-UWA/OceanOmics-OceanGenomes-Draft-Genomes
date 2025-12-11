//
// Subworkflow with functionality specific to the nf-core/oceangenomesmitogenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_SAMPLESHEET      } from '../../../modules/local/create_samplesheet'
include { samplesheetToList         } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_SAMPLESHEET {
    take:
        samplesheet_input
        input_reads
        samplesheet_prefix

    main:
        // println "DEBUG: PREPARE_SAMPLESHEET received input type: ${samplesheet_input?.getClass()} value=${samplesheet_input}"

        // declare channels
        def parsed_from_input_ch = Channel.empty()
        def samplesheet_src_ch   = Channel.empty()
        def create_samplesheet_input_ch = input_reads.collect()
        def default_samplesheet  = file("${params.outdir}/samplesheet/${samplesheet_prefix}_samplesheet.csv")

        if (params.input) {
            log.info "Using provided samplesheet: ${params.input}"
            ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)
        } else if (default_samplesheet.exists()) {
            log.info "Found existing samplesheet, using: ${default_samplesheet}"
            ch_samplesheet = Channel.fromPath(default_samplesheet)
        } else {
            log.info "Creating new samplesheet..."
            CREATE_SAMPLESHEET(
                create_samplesheet_input_ch,
                params.sql_config
            )
            ch_samplesheet = CREATE_SAMPLESHEET.out.samplesheet
        }
        
        ch_samples = ch_samplesheet
            .map { sheet_path ->
                samplesheetToList(sheet_path as String, "${projectDir}/assets/schema_input.json")
            }
        
        
        // ch_samples.view()
        // Normalize and validate
        def samplesheet_ch = ch_samples
            .flatMap { sample_list ->
                sample_list.collect { sample_record ->
                    // sample_record[0] is a meta map like [id:OG1341] so we want to destructure the map.
                    def raw_meta  = sample_record[0]
                    def sample_id = (raw_meta instanceof Map) ? raw_meta.id : raw_meta

                    def meta = [
                        id            : sample_id,
                        run           : sample_record[1],
                        date          : sample_record[2],
                        prefix        : sample_record[3],
                        nom_species_id: sample_record[4],
                        taxon_id      : sample_record[5],
                        class         : sample_record[6],
                    ]
                    def fastq_1 = sample_record[7]
                    def fastq_2 = sample_record[8]
                    [ sample_id, meta, [ fastq_1, fastq_2 ] ]
                }
            }
            .groupTuple()
            .map { sample_id, metas, fastqs -> validateInputSamplesheet(sample_id, metas, fastqs) }
            // .view()

    emit:
        samplesheet = samplesheet_ch
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//

def validateInputSamplesheet(sample_id, metas, fastqs) {
    // metas is a List of meta maps; usually one per sample ID
    def meta = metas[0]

    // OPTIONAL sanity check: all metas for this id should be identical
    // if (metas.unique(false) != [meta]) {
    //     throw new IllegalArgumentException("Inconsistent metadata for sample ${id}: ${metas}")
    // }

    // fastqs is a List of lists, e.g. [ [fq1, fq2], [fq1_lane2, fq2_lane2], ... ]
    def flat_fastqs = fastqs.flatten().findAll { it != null }

    return [ meta, flat_fastqs ]
}