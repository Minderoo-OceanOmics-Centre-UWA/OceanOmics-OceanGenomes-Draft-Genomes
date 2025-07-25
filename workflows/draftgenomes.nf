/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//download files
include { BASESPACE              } from '../modules/local/basespace/main'
include { samplesheetHybrid      } from '../subworkflows/local/samplesheetHybrid'
include { BBMAP_REPAIR           } from '../modules/nf-core/bbmap/repair/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//assembly
include { MERYL_COUNT } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_HISTOGRAM } from '../modules/nf-core/meryl/histogram/main'
include { GENOMESCOPE2 } from '../modules/nf-core/genomescope2/main'
include { MEGAHIT } from '../modules/nf-core/megahit/main'

//decontamination
include { FCSGX_RUNGX } from '../modules/nf-core/fcsgx/rungx/main'
include { FCSGX_CLEANGENOME } from '../modules/nf-core/fcsgx/cleangenome/main'
include { BBMAP_FILTERBYNAME } from '../modules/nf-core/bbmap/filterbyname/main'
include { FCS_FCSADAPTOR } from '../modules/nf-core/fcs/fcsadaptor/main'
include { FCSGX_CLEANGENOME as FSCSGX_CLEANGENOME_ADAPTOR } from '../modules/nf-core/fcsgx/cleangenome/main'
include { TIARA_TIARA } from '../modules/nf-core/tiara/tiara/main'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_TIARA } from '../modules/nf-core/bbmap/filterbyname/main'

//QC
include { BUSCO_BUSCO as BUSCO_ACTI } from '../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_VERT } from '../modules/nf-core/busco/busco/main'
include { BWAMEM2_INDEX } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_ACTI } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_MEM as BWAMEM2_MEM_VERT } from '../modules/nf-core/bwamem2/mem/main'
include { MERQURY_MERQURY } from '../modules/nf-core/merqury/merqury/main'
include { GFASTATS } from '../modules/nf-core/gfastats/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DRAFTGENOMES {

    take:
    //samplesheet // channel: samplesheet read in from --input
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
    
    // View samplesheet structure
    samplesheetHybrid.out.samplesheet.view()

    samplesheet_ch = samplesheetHybrid.out.samplesheet
    // MODULE: Run FastQC
    
    FASTQC (
        samplesheet_ch
    )

    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
   // sample_ch = params.input ?
    //     Channel.fromPath(params.input).map { parse_csv(it) } :
    //     generate_samplesheet_from_fastqs()//



    FASTP (
        samplesheet_ch, // tuple val(meta), path(reads)
        [], // path  adapter_fasta
        [], // val   discard_trimmed_pass
        [], // val   save_trimmed_fail
        [], // val   save_merged
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(FASTP.out.versions.first())


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME ASSEMBLY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//     //
//     // MODULE: Run Meryl
//     //
    
//     MERYL_COUNT (
//         HIFIADAPTERFILT.out.reads, // tuple val(meta), path(reads)
//         params.kvalue // val kvalue
//     )
//     ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())

//     MERYL_UNIONSUM (
//         // tuple val(meta), path(meryl_dbs)
//         // val kvalue
//     )


//     MERYL_HISTOGRAM (
//         MERYL_COUNT.out.meryl_db, // tuple val(meta), path(meryl_dbs)
//         params.kvalue // val kvalue
//     )
//     ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

//     //
//     // MODULE: Run Genomescope
//     //

//     GENOMESCOPE2 (
//         // tuple val(meta), path(histogram)
//     )


//     //
//     // MODULE: Run Megahit
//     //

//     MEGAHIT (
//         // tuple val(meta), path(reads1), path(reads2)
//     )



// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     GENOME DECONTAMINATION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */


//     //
//     // MODULE: Run fcs-gx find contamination
//     //

//     FCSGX_RUNGX (
//         // tuple val(meta), val(taxid), path(fasta)
//         // path gxdb
//         // val ramdisk_path
//     )


//     //
//     // MODULE: Run fcs-gx clean
//     //

//     FCSGX_CLEANGENOME (
//         // tuple val(meta), path(fasta), path(fcsgx_report)
//     )


//     //
//     // MODULE: Run bbmap filter
//     //

//     BBMAP_FILTERBYNAME (
//         // tuple val(meta), path(reads)
//         // val(names_to_filter)
//         // val(output_format)
//         // val(interleaved_output)
//     )



//     //
//     // MODULE: Run fcs adaptor find
//     //

//     FCS_FCSADAPTOR (
//         // tuple val(meta), path(assembly)
//     )

//     //
//     // MODULE: Run fcs-gx clean adaptor
//     //

//     FSCSGX_CLEANGENOME_ADAPTOR (
//         // tuple val(meta), path(fasta), path(fcsgx_report)
//     )


//     //
//     // MODULE: Run Tiara
//     //

//     TIARA_TIARA (
//         // tuple val(meta), path(fasta)
//     )


//     //
//     // MODULE: Run bbmap to filter tiara contamination
//     //

//     BBMAP_FILTERBYNAME_TIARA (
//         // tuple val(meta), path(reads)
//         // val(names_to_filter)
//         // val(output_format)
//         // val(interleaved_output)        
//     )

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     GENOME QUALITY CONTROL STATISTICS
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

//     //
//     // MODULE: Run BUSCO
//     //

//     BUSCO_BUSCO (
//         // tuple val(meta), path(fasta, stageAs:'tmp_input/*')
//         // val mode                              // Required:    One of genome, proteins, or transcriptome
//         // val lineage                           // Required:    lineage for checking against, or "auto/auto_prok/auto_euk" for enabling auto-lineage
//         // path busco_lineages_path              // Recommended: BUSCO lineages file - downloads if not set
//         // path config_file                      // Optional:    BUSCO configuration file
//         // val clean_intermediates               // Optional:    Remove intermediate files
//     )


//     //
//     // MODULE: Run BWA index
//     //

//     BWAMEM2_INDEX (
//         // tuple val(meta), path(fasta)
//     )


//     //
//     // MODULE: Run BWA align
//     //

//     BWAMEM2_MEM (
//         // tuple val(meta), path(reads)
//         // tuple val(meta2), path(index)
//         // tuple val(meta3), path(fasta)
//         // val   sort_bam
//     )


//     //
//     // MODULE: Run Merqury
//     //

//     MERQURY_MERQURY (
//         // tuple val(meta), path(meryl_db), path(assembly)
//     )


//     //
//     // MODULE: Run gfa stats
//     //

//     GFASTATS (
//         // tuple val(meta), path(assembly)
//         // val out_fmt
//         // val genome_size
//         // val target
//         // tuple val(meta2), path(agpfile)
//         // tuple val(meta3), path(include_bed)
//         // tuple val(meta4), path(exclude_bed)
//         // tuple val(meta5), path(instructions)
//     )













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
    fastp_reads = FASTP.out.reads // channel to the mitogenome pipeline
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
