/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'

//assembly
include { MERYL_COUNT               } from '../../../modules/nf-core/meryl/count'
include { MERYL_UNIONSUM            } from '../../../modules/nf-core/meryl/unionsum'
include { MERYL_HISTOGRAM           } from '../../../modules/nf-core/meryl/histogram'
include { GENOMESCOPE2              } from '../../../modules/nf-core/genomescope2'
include { MEGAHIT                   } from '../../../modules/nf-core/megahit'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOME_ASSEMBLY {

    take:
    fastp_reads // tuple val(meta), path('*.fastq.gz')
    //samplesheet // channel: samplesheet read in from --input
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    //
    // MODULE: Run Meryl count
    //
    
    MERYL_COUNT (
        fastp_reads, // tuple val(meta), path(reads)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())
    
    //
    // MODULE: Run Meryl unionsum
    //

    MERYL_UNIONSUM (
        MERYL_COUNT.out.meryl_dbs, // tuple val(meta), path(meryl_dbs)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())

    //
    // MODULE: Run Meryl histogram
    //

    MERYL_HISTOGRAM (
        MERYL_UNIONSUM.out.meryl_db, // tuple val(meta), path(meryl_dbs)
        params.kvalue // val kvalue
    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run Genomescope
    //

    GENOMESCOPE2 (
        MERYL_HISTOGRAM.out.hist // tuple val(meta), path(histogram)
    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())

    //
    // MODULE: Run Megahit
    //

    MEGAHIT (
        fastp_reads // tuple val(meta), path(reads1), path(reads2)
    )
    ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())

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
    megahit_assembled_contigs = MEGAHIT.out.contigs // pass to dweconatmination - fcs-gx
    meryl_db = MERYL_UNIONSUM.out.meryl_db  // need to check COUNT output to make sure im passing the right one
    genomescope_summary = GENOMESCOPE2.out.summary // pass into genome QC for size of genome (unique length)
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
