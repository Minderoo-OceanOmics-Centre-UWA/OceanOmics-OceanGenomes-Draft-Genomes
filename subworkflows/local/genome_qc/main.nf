/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML            } from '../../nf-core/utils_nfcore_pipeline'

//QC
include { BUSCO_BUSCO                       } from '../../../modules/nf-core/busco/busco'
include { BWAMEM2_INDEX                     } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                       } from '../../../modules/nf-core/bwamem2/mem'
include { MERQURY_MERQURY                   } from '../../../modules/nf-core/merqury/merqury'
include { GFASTATS                          } from '../../../modules/nf-core/gfastats'


//FUNCTION: Join two [meta, value] channels on a set of keys, then merge metas.
//          By default joins on id, run, date, prefix.
def join_on_keys_and_merge = { ch1, ch2, List keys = ['id','run','date','prefix'] ->

    def keyer = { Map m -> keys.collect { k ->
        if( !m.containsKey(k) )
            throw new IllegalArgumentException("Missing meta key '${k}' in ${m}")
        m[k]
    }}

    def c1 = ch1.map { meta, val -> [ keyer(meta), [meta, val] ] }
    def c2 = ch2.map { meta, val -> [ keyer(meta), [meta, val] ] }

    c1.join(c2).map { keyvals, left, right ->
        def (m1, v1) = left
        def (m2, v2) = right
        // Merge maps; values from m2 override m1 on duplicate keys
        // Return as tuple: [merged_meta, v1, v2]
        return tuple(m1 + m2, v1, v2)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOME_QC {

    take:
    fastp // tuple val(meta), path('*.fastq.gz'),
    assembly // tuple val(meta), path(assembly)
    meryl_db // tuple val(meta), path(meryl_dir)
    genomescope_summary // tuple val(meta), path("${meta.prefix}_summary.txt") 
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME QUALITY CONTROL STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Determining which BUSCO database to use based on meta.class
    ch_with_busco_db = assembly.map { meta, assembly ->
        def busco_db = meta.class == 'Actinopteri' ? 
            params.busco_acti_db : 
            params.busco_vert_db
        
        return [meta, assembly, busco_db]
    }

    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (
        ch_with_busco_db, // tuple val(meta), path(fasta, stageAs:'tmp_input/*'), val(lineage/db)
        Channel.value('genome') // val mode // Required:    One of genome, proteins, or transcriptome
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    //
    // MODULE: Run BWA index
    //

    BWAMEM2_INDEX (
        assembly // tuple val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

    ch_tmp = join_on_keys_and_merge(fastp, BWAMEM2_INDEX.out.index)
    ch_bwamem2_mem_input = ch_tmp.join(assembly, by:0)

    //
    // MODULE: Run BWA align
    //

    BWAMEM2_MEM (
        ch_bwamem2_mem_input // tuple val(meta), path(reads), path(index), path(fasta)
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    //
    // Channel for merqury
    //

    ch_merqury_input = join_on_keys_and_merge(meryl_db, assembly)

    //
    // MODULE: Run Merqury
    //

    MERQURY_MERQURY (
        ch_merqury_input // tuple val(meta), path(meryl_db), path(assembly)
    )
    ch_versions = ch_versions.mix(MERQURY_MERQURY.out.versions.first())

    //
    // Channel for gfa stats
    //

    ch_gfastats_input = join_on_keys_and_merge(assembly, genomescope_summary)

    //
    // MODULE: Run gfa stats
    //

    GFASTATS (
        ch_gfastats_input, // tuple val(meta), path(assembly), path(genomescope_summary)
        "fa", // val out_fmt
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())

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
    results = GFASTATS.out.assembly_summary
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_collated_versions              // channel: [ path(versions.yml) ]
}

