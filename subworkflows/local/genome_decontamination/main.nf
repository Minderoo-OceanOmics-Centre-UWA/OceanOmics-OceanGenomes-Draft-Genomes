/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'

//decontamination
include { FCSGX_RUNGX               } from '../../../modules/nf-core/fcs/fcsgx/rungx'
include { FCSGX_CLEANGENOME         } from '../../../modules/nf-core/fcs/fcsgx/cleangenome'
include { BBMAP_FILTERBYNAME        } from '../../../modules/local/bbmap/filterbyname'
include { BBMAP_REFORMAT            } from '../../../modules/local/bbmap/reformat'
include { FCS_FCSADAPTOR            } from '../../../modules/nf-core/fcs/fcsadaptor'
include { FCSGX_CLEANGENOME as FSCSGX_CLEANGENOME_ADAPTOR   } from '../../../modules/nf-core/fcs/fcsgx/cleangenome'
include { TIARA_TIARA               } from '../../../modules/nf-core/tiara/tiara'
include { BBMAP_FILTERBYNAME as BBMAP_FILTERBYNAME_TIARA    } from '../../../modules/nf-core/bbmap/filterbyname'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOME_DECONTAMINATION {

    take:
    genome_assembly // Results from megahit - tuple val(meta), path("*.v129mh.fasta")
    // may need to zip the fasta for fcs-gx, so will wait to see if it works 
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    fasta_with_assembly_prefix = genome_assembly
    .map { meta, fasta ->
        def meta_ext = meta.containsKey('assembly_prefix') 
            ? meta 
            : meta + [ assembly_prefix: fasta.baseName ]
        [meta_ext, fasta]
    }

    //
    // MODULE: Run fcs-gx find contamination
    //

    FCSGX_RUNGX (
        fasta_with_assembly_prefix, // tuple val(meta), path(fasta)
        params.gxdb, // path gxdb
        params.ramdisk_path, // val ramdisk_path
    )
    ch_versions = ch_versions.mix(FCSGX_RUNGX.out.versions.first())

    filter_contaminants_combined_ch = fasta_with_assembly_prefix
        .join(FCSGX_RUNGX.out.fcsgx_report, by: 0)

    //
    // MODULE: Run fcs-gx clean
    //

    FCSGX_CLEANGENOME (
        filter_contaminants_combined_ch, // tuple val(meta), path(fasta), path(fcsgx_report)
        Channel.value('rc'),// val(cleaned_suffix) - suffix for the cleaned fast output
        Channel.value('contam'),// val(contam_suffix) - suffix for the contaminants output
    )
    ch_versions = ch_versions.mix(FCSGX_CLEANGENOME.out.versions.first())

    fcsgx_combined_ch = FCSGX_RUNGX.out.fcsgx_report
        .join(FCSGX_CLEANGENOME.out.cleaned, by: 0)

    //
    // MODULE: Run bbmap filter
    //

    BBMAP_FILTERBYNAME (
        fcsgx_combined_ch, // tuple val(meta), path(names_to_filter), path(reads)
        Channel.value('fa'), // val(output_format)
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME.out.versions.first())

    //
    // MODULE: Run bbmap reformat
    //

    BBMAP_REFORMAT (
        BBMAP_FILTERBYNAME.out.first_filtered, // tuple val(meta), path("$first_filtered_reads") 
    )
    ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions.first())


    //
    // MODULE: Run fcs adaptor find
    //

    FCS_FCSADAPTOR (
        BBMAP_REFORMAT.out.reads // tuple val(meta), path(reads)
    )
    ch_versions = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    filter_adaptors_combined_ch = BBMAP_REFORMAT.out.reads
        .join(FCS_FCSADAPTOR.out.adaptor_report, by: 0)
    
    //
    // MODULE: Run fcs-gx clean adaptor
    //

    FSCSGX_CLEANGENOME_ADAPTOR (
        filter_adaptors_combined_ch, // tuple val(meta), path(fasta), path(fcsgx_report)
        Channel.value('rmadapt'),// val(cleaned_suffix) - suffix for the cleaned fast output
        Channel.value('adaptor-contam'),// val(contam_suffix) - suffix for the contaminants output
    )
    ch_versions = ch_versions.mix(FSCSGX_CLEANGENOME_ADAPTOR.out.versions.first())

    //
    // MODULE: Run Tiara
    //

    TIARA_TIARA (
        FSCSGX_CLEANGENOME_ADAPTOR.out.cleaned // tuple val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(TIARA_TIARA.out.versions.first())

    tiara_filter_combined_ch = FSCSGX_CLEANGENOME_ADAPTOR.out.cleaned
        .join(TIARA_TIARA.out.contig_removal, by: 0)

    //
    // MODULE: Run bbmap to filter tiara contamination
    //

    BBMAP_FILTERBYNAME_TIARA (
        tiara_filter_combined_ch, // tuple val(meta), path(reads), path(names_to_filter)
        Channel.value('fna'),// val(output_format)
        [],// val(interleaved_output)        
    )
    ch_versions = ch_versions.mix(BBMAP_FILTERBYNAME_TIARA.out.versions.first())

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
    decontamined_assembled_reads = BBMAP_FILTERBYNAME_TIARA.out.reads // channel to the mitogenome pipeline
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]



}

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     THE END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */
