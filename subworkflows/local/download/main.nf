/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Download file modules
include { BASESPACE              } from '../../../modules/local/basespace'
include { BBMAP_REPAIR           } from '../../../modules/nf-core/bbmap/repair'
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN BASESPACE DATA DOWNLOAD WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNLOAD_READS {

    take:
    run_id // params.run in nextflow.config
    bs_config
    
    main:

    println "Run ID is: ${run_id}"  // ✅ Works fine here
    println "bs_config is: ${bs_config}"  // ✅ Works fine here
    //samplesheet.view { sheet -> "Samplesheet contents: $sheet" }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Using the run_id, downloads the full run from basespace
    //

    BASESPACE (
        run_id, // val run
        bs_config //params.bs_config// path config ??? i dont think this one is needed maybe a carry over from adams script but not used in this code
    )

    //
    // Group fastq reads by meta_id
    //

    reads_by_meta_id = BASESPACE.out.fastqs
        .flatten()
        .map { file -> 
            def matcher = file.name =~ /^([A-Z]+\d+)/
            if (matcher) {
                def meta_id = matcher[0][1]
                return tuple(meta_id, file)
            } else {
                log.info "Skipping non-matching file: ${file.name}"
                return null  // Ignore files not matching OG123 etc.
            }
        }
        .filter { it != null }  // Remove nulls from the channel
        .groupTuple()  // <-- Group by ogid key

    //
    // MODULE: Repair the sequencing reads
    //

    BBMAP_REPAIR (
        reads_by_meta_id, // tuple val(meta_id), path(reads)
        params.interleave// val(interleave)
    )

    //
    // Collect files
    //

    // ch_multiqc_files = ch_multiqc_files.mix(BASESPACE.out.json.collect{it[1]})
    // ch_multiqc_files = ch_multiqc_files.mix(BBMAP_REPAIR.out.log.collect{it})
    ch_versions = ch_versions.mix(BBMAP_REPAIR.out.versions.first())


    //
    // Emit outputs
    //

    emit:
    repaired_reads = BBMAP_REPAIR.out.repaired // channel for fastp - tuple val(ogid), path("*.{R1,R2}.fq.gz")
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]
}

