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
include { GFASTATS                          } from '../../../modules/local/gfastats'
include { SEQKIT_STATS                      } from '../../../modules/nf-core/seqkit/stats'
include { RECHECK_GENOME_SIZE               } from '../../../modules/local/coverage/recheck'
include { buscoDbForClass                   } from '../utils_nfcore_oceangenomes_draftgenomes_pipeline'
include { MEASURE_KMER_COVERAGE             } from '../../../modules/local/coverage/kmer_depth'
include { SUMMARISE_KMER_COVERAGE           } from '../../../modules/local/coverage/kmer_depth_summarise'
include { GENOMESCOPE_RESEED                } from '../../../modules/local/genomescope_reseed'


//FUNCTION: Join multiple [meta, v1, v2, ...] channels on a set of keys, then merge metas.
//          By default joins on id, run, date, prefix.
//
//          A channel may emit any number of values after the meta, not just one. Every
//          value is concatenated positionally into the result, in the order the channels
//          are listed, so the downstream process input must be declared in that same
//          order. MEASURE_KMER_COVERAGE.out.bedcov emits [meta, bedcov, status] and an
//          earlier two-argument version of the inner map threw MissingMethodException on
//          it, taking the whole run down.
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
    
    // Convert all channels to [key, [meta, values]] format. Takes the emission as a single
    // argument rather than destructuring it, so a channel emitting more than one value
    // after the meta keys the same way as one emitting a single value.
    def keyedChannels = channels.collect { ch ->
        ch.map { entry ->
            def items = entry instanceof List ? entry : [entry]
            def meta = items[0]
            def vals = items.size() > 1 ? items[1..-1] : []
            def keyVals = keyer(meta)
            [ keyVals, [meta, vals] ]
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
            def (meta, vals) = metaVal
            merged_meta = merged_meta + meta
            values.addAll(vals)
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
    coverage_summary // tuple val(meta), path("*_coverage_summary.json") -- provisional
    fastp_json // tuple val(meta), path("*.fastp.json")
    genomescope_model // tuple val(meta), path("*_model.txt")
    meryl_hist // tuple val(meta), path("*.meryl.hist")

    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_inputs = Channel.empty()
    ch_sample_multiqc_inputs = Channel.empty()
    ch_merqury_results = Channel.empty()
    ch_gfastats_results = Channel.empty()
    ch_seqkit_stats_results = Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME QUALITY CONTROL STATISTICS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    
    // The class-to-lineage ladder lives in the pipeline utils subworkflow so the
    // precomputed-BUSCO branch in draftgenomes.nf reaches the same answer when it has to
    // choose between result sets from two different lineages.
    ch_with_busco_db = assembly.map { meta, file -> [meta, file, buscoDbForClass(meta)] }

    
    //
    // MODULE: Run BUSCO
    //

    BUSCO_BUSCO (
        ch_with_busco_db, // tuple val(meta), path(fasta, stageAs:'tmp_input/*'), val(lineage/db)
        Channel.value('genome') // val mode // Required:    One of genome, proteins, or transcriptome
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    ch_multiqc_inputs = ch_multiqc_inputs.mix(BUSCO_BUSCO.out.short_summaries_txt)

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
    // MODULE: Run Seqkit stats
    //

    SEQKIT_STATS (
        assembly // tuple val(meta), path(assembly)
    )
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())
    ch_seqkit_stats_results = SEQKIT_STATS.out.stats // channel: tuple val(meta), path("*.seqkit_stats.tsv")

    //
    // MODULE: Measure haploid k-mer coverage from read depth
    //
    // An independent estimate of the quantity GenomeScope fits, from the alignments rather
    // than the k-mer histogram. It needs the BAM and the BUSCO table, so it can only run
    // here -- GENOME_ASSEMBLY has neither.
    //
    // Split across two containers: samtools has no python, python has no samtools.
    MEASURE_KMER_COVERAGE (
        join_on_keys_and_merge([BWAMEM2_MEM.out.bam, BUSCO_BUSCO.out.full_table])
    )
    ch_versions = ch_versions.mix(MEASURE_KMER_COVERAGE.out.versions.first())

    ch_kmer_depth_summarise = join_on_keys_and_merge(
        [MEASURE_KMER_COVERAGE.out.bedcov, MEASURE_KMER_COVERAGE.out.contig_depth,
         fastp_json, ch_seqkit_stats_results, genomescope_model])

    SUMMARISE_KMER_COVERAGE (
        ch_kmer_depth_summarise // tuple val(meta), path(bedcov), path(status), path(contig_depth), path(fastp_json), path(seqkit_stats), path(model)
    )
    ch_versions = ch_versions.mix(SUMMARISE_KMER_COVERAGE.out.versions.first())
    ch_kmer_depth = SUMMARISE_KMER_COVERAGE.out.kmer_depth

    //
    // MODULE: Refit GenomeScope seeded with the measured coverage, where the first failed
    //
    // Only runs for samples the gate flagged, and only where the measured lambda is high
    // enough that there is something to fit. It keeps the reseeded fit only if that fit
    // beats the original on its own terms; otherwise the original stands.
    //
    ch_reseed_input = join_on_keys_and_merge(
        [meryl_hist, coverage_summary, ch_kmer_depth, genomescope_summary, genomescope_model])

    GENOMESCOPE_RESEED (
        ch_reseed_input // tuple val(meta), path(hist), path(coverage_json), path(kmer_depth_json), path(summary), path(model)
    )
    ch_versions = ch_versions.mix(GENOMESCOPE_RESEED.out.versions.first())
    ch_reseed_decision = GENOMESCOPE_RESEED.out.decision

    //
    // MODULE: Cross-check the GenomeScope genome size against the assembly
    //
    // The coverage summary written during assembly is provisional: at that point MEGAHIT
    // had not run, so the one reliability check that needs an assembly could not be made.
    // This republishes it as the final verdict.
    //
    ch_recheck_input = join_on_keys_and_merge(
        [coverage_summary, ch_seqkit_stats_results, ch_kmer_depth])

    RECHECK_GENOME_SIZE (
        ch_recheck_input // tuple val(meta), path(coverage_json), path(seqkit_stats), path(kmer_depth_json)
    )
    ch_versions = ch_versions.mix(RECHECK_GENOME_SIZE.out.versions.first())
    ch_coverage_summary_final = RECHECK_GENOME_SIZE.out.coverage_json


    //
    // Collect files
    //

    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(ch_multiqc_inputs)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(BUSCO_BUSCO.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(BWAMEM2_INDEX.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(BWAMEM2_MEM.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MERQURY_MERQURY.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(GFASTATS.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(SEQKIT_STATS.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(MEASURE_KMER_COVERAGE.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(SUMMARISE_KMER_COVERAGE.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(GENOMESCOPE_RESEED.out.tool_params)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(RECHECK_GENOME_SIZE.out.multiqc)
    ch_sample_multiqc_inputs = ch_sample_multiqc_inputs.mix(RECHECK_GENOME_SIZE.out.tool_params)
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_inputs.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BUSCO_BUSCO.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BWAMEM2_INDEX.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BWAMEM2_MEM.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MERQURY_MERQURY.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GFASTATS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MEASURE_KMER_COVERAGE.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(SUMMARISE_KMER_COVERAGE.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(GENOMESCOPE_RESEED.out.tool_params.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(RECHECK_GENOME_SIZE.out.tool_params.collect { it[1] })
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    ch_versions = ch_versions.mix(EXTRACT_BUSCO_SEQUENCES.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    ch_versions = ch_versions.mix(MERQURY_MERQURY.out.versions.first())
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())


    emit:
    busco_short_summary = BUSCO_BUSCO.out.short_summaries_json // channel: tuple val(meta), path('*.busco.short_summary.txt')
    merqury_results = ch_merqury_results // channel: tuple val(meta), path("*.completeness.stats"), path("${prefix}.qv")
    gfastats_results = ch_gfastats_results // channel: tuple val(meta), path("*.assembly_summary")
    seqkit_stats_results = ch_seqkit_stats_results // channel: tuple val(meta), path("*.seqkit_stats.tsv")
    coverage_summary = ch_coverage_summary_final // channel: tuple val(meta), path("*_coverage_summary.json")
    kmer_depth = ch_kmer_depth // channel: tuple val(meta), path("*_kmer_depth.json")
    reseed_decision = ch_reseed_decision // channel: tuple val(meta), path("*_reseed_decision.json")
    multiqc_files = ch_multiqc_files             // channel: [ path(multiqc_files) ]
    multiqc_inputs = ch_sample_multiqc_inputs    // channel: [ tuple(meta), path(multiqc_file) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]
}
