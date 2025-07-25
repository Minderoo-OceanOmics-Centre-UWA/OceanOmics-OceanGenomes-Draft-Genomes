/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//download files
include { BASESPACE              } from '../modules/local/basespace/main'
include { BBMAP_REPAIR           } from '../modules/nf-core/bbmap/repair/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOWNLOAD {

    take:
    run_id // params.run in nextflow.config
    bs_config
    
    main:

    println "Run ID is: ${run_id}"  // ✅ Works fine here
    println "bs_config is: ${bs_config}"  // ✅ Works fine here
    //samplesheet.view { sheet -> "Samplesheet contents: $sheet" }

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNLOAD THE DATA FROM BASESPACE AND RUN THROUGH FASTP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    BASESPACE (
        run_id, // val run
        bs_config //params.bs_config// path config ??? i dont think this one is needed maybe a carry over from adams script but not used in this code
    )

    // // Run fastqc on all of the fastq files
    // FASTQC (
    //     BASESPACE.out.fastqs.flatten() // tuple val(meta), path(reads)
    // )

    // Group fastq reads by OGID
    reads_by_ogid = BASESPACE.out.fastqs
        .flatten()
        .map { file -> 
            def matcher = file.name =~ /^([A-Z]+\d+)/
            if (matcher) {
                def ogid = matcher[0][1]
                return tuple(ogid, file)
            } else {
                log.info "Skipping non-matching file: ${file.name}"
                return null  // Ignore files not matching OG123 etc.
            }
        }
        .filter { it != null }  // Remove nulls from the channel
        .groupTuple()  // <-- Group by ogid key


    BBMAP_REPAIR (
        reads_by_ogid, // tuple val(ogid), path(reads)
        params.interleave// val(interleave)
    )



    //BBMAP_REPAIR.out.repaired.view()

    repaired_ch = BBMAP_REPAIR.out.repaired
    // Make a new channel from bbmap_repair output or from samplesheet
    samplesheetHybrid(
        repaired_ch
    )


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
    repaired_reads = BBMAP_REPAIR.out.repaired // channel for fastp
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
