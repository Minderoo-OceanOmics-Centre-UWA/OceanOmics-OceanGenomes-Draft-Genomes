/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//QC
include { BUSCO_BUSCO                       } from '../../../modules/nf-core/busco/busco'
include { EXTRACT_BUSCO_SEQUENCES           } from '../../../modules/local/extract_busco_sequences'
include { BWAMEM2_INDEX                     } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                       } from '../../../modules/nf-core/bwamem2/mem'
include { MERQURY_MERQURY                   } from '../../../modules/nf-core/merqury/merqury'
include { GFASTATS                          } from '../../../modules/nf-core/gfastats'


//FUNCTION: Join multiple [meta, value] channels on a set of keys, then merge metas.
//          By default joins on id, run, date, prefix.
def join_on_keys_and_merge = { channels, List keys = ['id','run','date','prefix'] ->
    // println "DEBUG: Starting join with ${channels.size()} channels"
    
    def keyer = { Map m -> 
        def keyVals = keys.collect { k ->
            if( !m.containsKey(k) )
                throw new IllegalArgumentException("Missing meta key '${k}' in ${m}")
            m[k]
        }
        // println "DEBUG: Generated key for ${m.id}: ${keyVals}"
        return keyVals
    }
    
    // Convert all channels to [key, [meta, value]] format
    def keyedChannels = channels.collect { ch ->
        ch.map { meta, val -> 
            def keyVals = keyer(meta)
            [ keyVals, [meta, val] ] 
        }
        // .view { "DEBUG: Keyed channel entry: ${it[0]} -> meta.id: ${it[1][0].id}" }
    }
    
    // Join all channels sequentially
    def joined = keyedChannels[0]
        // .view { "DEBUG: First channel entry: ${it}" }
    for (int i = 1; i < keyedChannels.size(); i++) {
        // println "DEBUG: About to join with channel ${i+1}"
        joined = joined.join(keyedChannels[i])
        // joined = joined.view { "DEBUG: After joining with channel ${i+1}: ${it}" }
    }
    
    // Merge metas and extract values
    joined.map { tuple ->
        // println "DEBUG: Final mapping tuple: ${tuple}"
        def keyvals = tuple[0]
        def metaVals = tuple[1..-1]
        
        def merged_meta = [:]
        def values = []
        
        metaVals.each { metaVal ->
            def (meta, val) = metaVal
            merged_meta = merged_meta + meta
            values << val
        }
        
        def result = [merged_meta] + values
        // println "DEBUG: Final result: ${result[0].id}"
        return result
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
    genomescope_summary // tuple val(meta), path("${meta.prefix}_summary.txt") 
    assembly // tuple val(meta), path(assembly)
    meryl_db // tuple val(meta), path(meryl_dir)

    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_merqury_results = Channel.empty()
    ch_gfastats_results = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME QUALITY CONTROL STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    
    // Determining which BUSCO database to use based on meta.class
    ch_with_busco_db = assembly.map { meta, file ->
        def busco_db = meta.class == 'Actinopteri' ? 
            params.busco_acti_db : 
            params.busco_vert_db
        
        return [meta, file, busco_db]
    }

    
    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (
        ch_with_busco_db, // tuple val(meta), path(fasta, stageAs:'tmp_input/*'), val(lineage/db)
        Channel.value('genome') // val mode // Required:    One of genome, proteins, or transcriptome
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

    // Channel for extract busco sequences
    ch_extract_busco_input = join_on_keys_and_merge([BUSCO_BUSCO.out.full_table, assembly])
    // Returns: tuple val(meta), path(busco_table), path(genome_fasta)


    //
    // MODULE: Extract the coding sequences from the genome for all the busco results
    //

    EXTRACT_BUSCO_SEQUENCES(ch_extract_busco_input)


    //
    // MODULE: Run BWA index
    //

    BWAMEM2_INDEX (
        assembly // tuple val(meta), path(fasta)
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())

    
    // Channel for BWA mem
    ch_bwamem2_mem_input = join_on_keys_and_merge([fastp, BWAMEM2_INDEX.out.index, assembly])
    // Returns: [merged_meta, fastp_files, index, assembly]

    //
    // MODULE: Run BWA align
    //

    BWAMEM2_MEM (
        ch_bwamem2_mem_input // tuple val(meta), path(reads), path(index), path(fasta)
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())


    // Channel for merqury
    ch_merqury_input = join_on_keys_and_merge([meryl_db, assembly])
    
    //
    // MODULE: Run Merqury
    //

    MERQURY_MERQURY (
        ch_merqury_input // tuple val(meta), path(meryl_db), path(assembly)
    )
    ch_versions = ch_versions.mix(MERQURY_MERQURY.out.versions.first())
    ch_merqury_results = MERQURY_MERQURY.out.stats.join(MERQURY_MERQURY.out.assembly_qv, by:0) // channel: tuple val(meta), path("*.completeness.stats")


    // Channel for gfa stats
    ch_gfastats_input = join_on_keys_and_merge([assembly, genomescope_summary])


    //
    // MODULE: Run gfa stats
    //

    GFASTATS (
        ch_gfastats_input, // tuple val(meta), path(assembly), path(genomescope_summary)
        "fa", // val out_fmt
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())
    ch_gfastats_results = GFASTATS.out.assembly_summary // channel: tuple val(meta), path("*.assembly_summary")
    
    //
    // Collect files
    //

    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_SEQUENCES.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    ch_versions = ch_versions.mix(MERQURY_MERQURY.out.versions.first())
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())


    emit:
    busco_short_summary = BUSCO_BUSCO.out.short_summaries_json // channel: tuple val(meta), path('*.busco.short_summary.txt')
    merqury_results = ch_merqury_results // channel: tuple val(meta), path("*.completeness.stats"), path("${prefix}.qv")
    gfastats_results = ch_gfastats_results // channel: tuple val(meta), path("*.assembly_summary")
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]
}

