//
// Subworkflow with functionality specific to the nf-core/oceangenomesmitogenomes pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_SAMPLESHEET      } from '../../../modules/local/create_samplesheet'
include { DOWNLOAD_TAXONKIT_DB    } from '../../../modules/local/download_taxonkit_db'
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

            // The NCBI taxdump backs up the curated species table when it has no
            // row for a sample's nominal name. storeDir makes this a cache hit
            // whenever a mitogenome run has already pulled the same dump.
            def taxdump_ch = Channel.value([])
            if (params.taxonkit_db_dir) {
                DOWNLOAD_TAXONKIT_DB(Channel.value("taxdump"))
                taxdump_ch = DOWNLOAD_TAXONKIT_DB.out.db_files
            }
            else {
                log.warn "No --taxonkit_db_dir set: taxon_id/class come from the species table only, so any sample it does not carry will abort the run."
            }

            CREATE_SAMPLESHEET(
                create_samplesheet_input_ch,
                params.sql_config,
                taxdump_ch
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
                // Stop the run before any assembly work when a sample's taxonomy did
                // not resolve. Checked on the whole list so one message names every
                // offending sample, and placed here because the samplesheet has
                // already been written and published by this point: the operator can
                // correct the rows in place and re-run.
                validateTaxonomy(sample_list)

                sample_list.collect { sample_record ->
                    // Every column is annotated with `meta` in schema_input.json, so
                    // nf-schema hands back [meta, fastq_1, fastq_2] with the taxonomy
                    // already keyed by name. Read it by name: the previous positional
                    // form silently reassigned every field if a column moved.
                    def meta      = sample_record[0]
                    def sample_id = meta.id
                    def fastq_1   = sample_record[1]
                    def fastq_2   = sample_record[2]
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

// A samplesheet value that carries no taxonomy: null, blank, or one of the
// placeholders create_samplesheet.py writes when neither the curated `species`
// table nor the NCBI taxdump could resolve the name.
def isUnresolvedTaxon(value) {
    def s = (value == null) ? '' : value.toString().trim()
    return !s || s.equalsIgnoreCase('unknown') || s.equalsIgnoreCase('None')
}

// Abort the run when any sample's taxon_id or class is unresolved.
//
// create_samplesheet.py deliberately writes these rows as `unknown` and exits 0,
// so that the samplesheet and its taxonomy_resolution.tsv are published and can be
// corrected by hand. That makes this the guard that has to stop the run, and it
// has to stop it here: an `unknown` taxon_id reaches FCS-GX as `--tax-id unknown`
// and fails only after MEGAHIT has already burned the SUs, and an `unknown` class
// has no BUSCO lineage. Nothing downstream may start until these are fixed.
def validateTaxonomy(sample_list) {
    def unresolved = []
    sample_list.each { sample_record ->
        def meta = sample_record[0]
        if (!(meta instanceof Map)) return
        def missing = []
        if (isUnresolvedTaxon(meta.taxon_id)) missing << 'taxon_id'
        if (isUnresolvedTaxon(meta.class))    missing << 'class'
        if (!missing) return
        def label = "${meta.id}\tnom_species_id='${meta.nom_species_id ?: ''}'\tmissing ${missing.join(', ')}"
        if (!unresolved.contains(label)) unresolved << label
    }

    if (!unresolved) return

    def sheet = params.outdir ? "${params.outdir}/samplesheet" : "the samplesheet/ directory under --outdir"
    def message = """Taxonomy is unresolved for ${unresolved.size()} sample(s):
  ${unresolved.join('\n  ')}
The samplesheet has been written and published to ${sheet} (see taxonomy_resolution.tsv alongside it), with these rows marked 'unknown'. Nothing has been assembled.
Fix them either way round: correct the nominal_species_id in the sample table and delete the published samplesheet to regenerate it, or edit the taxon_id/class columns in the published samplesheet and re-run."""

    // error() from inside a channel operator surfaces on the console only as
    // "Unexpected error [InvocationTargetException]", with the message buried in
    // .nextflow.log under a stack trace. Print it first so the sample list is
    // actually readable, then abort.
    log.error message
    error "Taxonomy is unresolved for ${unresolved.size()} sample(s); see the message above."
}

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