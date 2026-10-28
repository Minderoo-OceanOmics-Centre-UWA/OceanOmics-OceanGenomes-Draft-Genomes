#!/bin/bash
# Check every sample under $rundir for the files the pipeline should have published,
# before the final backup moves them. Writes _FileCheck.tsv: one row per expected file
# (patterns are relative to the sample directory, <OG> standing for the sample ID) and
# one column per sample, each cell "Found <path>" or "Missing".

# Load in the configfile
. ../configfile.txt

## Use this variable to override configfile settings if needed, otherwise comment them out
# rundir=/scratch/pawsey1348/tpeirce/NOVA_251215_AD/draftgenomes
rm -rf FileCheck
rm -f _FileCheck.tsv

mkdir -p FileCheck

# Expected files, relative to the sample directory, with OG standing for the sample ID.
# Each entry is ONE row of the report however many files its glob matches. The previous
# version expanded the globs into the list itself, so a pattern matching two files (as
# coverage/*.json does since the depth summary was added) added a row for that sample
# only and shifted every row below it against the labels.
expected=(
    # Read QC
    'fastp'
    'fastp/OG.ilmn.*.R1.fastq.gz'
    'fastp/OG.ilmn.*.R2.fastq.gz'
    'fastp/*.fastp.json'
    'fastp/OG.ilmn.*.fastp.html'
    'fastp/OG.ilmn.*.fastp.log'
    'fastp/fastqc'
    'fastp/fastqc/OG.ilmn.*.R1_fastqc.html'
    'fastp/fastqc/OG.ilmn.*.R1_fastqc.zip'
    'fastp/fastqc/OG.ilmn.*.R2_fastqc.html'
    'fastp/fastqc/OG.ilmn.*.R2_fastqc.zip'
    # Read-level decontamination (absent when run with --skip_kraken2_decontamination)
    'kraken2'
    'kraken2/OG.ilmn.*.kraken2.report.txt'
    'kraken2/OG.ilmn.*.kraken2_retention.txt'
    # Genome size and coverage verdict
    'coverage'
    'coverage/OG.ilmn.*_coverage_summary.json'
    'coverage/OG.ilmn.*_sequencing_coverage.txt'
    'coverage/OG.ilmn.*_kmer_depth.json'
    # K-mers
    'kmers'
    'kmers/OG.ilmn.*.meryl.hist'
    'kmers/OG.ilmn.*.meryl'
    'kmers/OG.ilmn.*.meryl/merylIndex'
    'kmers/OG.ilmn.*genomescope'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_summary.txt'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_model.txt'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_linear_plot.png'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_log_plot.png'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_transformed_linear_plot.png'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_genomescope2*_transformed_log_plot.png'
    'kmers/OG.ilmn.*genomescope/OG.ilmn.*_reseed_decision.json'
    'kmers/OG.ilmn.*.v129mh.merqury.qv'
    'kmers/*merqury.completeness.stats'
    'kmers/OG.ilmn.*.v129mh.merqury.OG.ilmn.*.v129mh.fna.qv'
    'kmers/OG.ilmn.*.v129mh.merqury.OG.ilmn.*.v129mh.fna.spectra-cn.fl.png'
    'kmers/OG.ilmn.*.v129mh.merqury.OG.ilmn.*.v129mh.fna.spectra-cn.ln.png'
    'kmers/OG.ilmn.*.v129mh.merqury.OG.ilmn.*.v129mh.fna.spectra-cn.st.png'
    'kmers/OG.ilmn.*.v129mh.merqury.spectra-asm.fl.png'
    'kmers/OG.ilmn.*.v129mh.merqury.spectra-asm.ln.png'
    'kmers/OG.ilmn.*.v129mh.merqury.spectra-asm.st.png'
    # Assembly and decontamination
    'assemblies'
    'assemblies/genome'
    'assemblies/genome/OG.ilmn.*.v129mh.fasta'
    'assemblies/genome/OG.ilmn.*.megahit_checkpoint.txt'
    'assemblies/genome/OG.ilmn.*.rmadapt.fasta'
    'assemblies/genome/OG.ilmn.*.adaptor-contam.fasta'
    'assemblies/genome/OG.ilmn.*.v129mh.fna'
    'assemblies/genome/NCBI'
    'assemblies/genome/NCBI/*.fcs_gx_report.txt'
    'assemblies/genome/NCBI/*.taxonomy.rpt'
    'assemblies/genome/NCBI/*.summary.txt'
    'assemblies/genome/NCBI/OG.ilmn.*.v129mh.rc.fasta'
    'assemblies/genome/NCBI/OG*.contam.fasta'
    'assemblies/genome/NCBI/OG.ilmn.*.filter_report.txt'
    'assemblies/genome/NCBI/OG.ilmn.*.review_scaffolds.txt'
    'assemblies/genome/NCBI/OG.ilmn.*.contig_count_500bp.txt'
    'assemblies/genome/NCBI/OG.ilmn.*.v129mh.fa'
    'assemblies/genome/NCBI/adaptor'
    'assemblies/genome/NCBI/adaptor/*cleaned_sequences.fa.gz'
    'assemblies/genome/NCBI/adaptor/*fcs_adaptor.log'
    'assemblies/genome/NCBI/adaptor/*fcs_adaptor_report.txt'
    'assemblies/genome/NCBI/adaptor/*pipeline_args.yaml'
    'assemblies/genome/NCBI/adaptor/*skipped_trims.jsonl'
    'assemblies/genome/tiara'
    'assemblies/genome/tiara/OG.ilmn.*.tiara.txt'
    'assemblies/genome/tiara/log_OG.ilmn.*.tiara.txt'
    'assemblies/genome/tiara/OG.ilmn.*.tiara.contig_removal.txt'
    'assemblies/genome/tiara/OG.ilmn.*.tiara_filter_summary.txt'
    # Genome QC
    'assemblies/genome/busco'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.short_summary.txt'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.short_summary.json'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.full_table.tsv'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.missing_busco_list.tsv'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.busco_sequences.tar.gz'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.busco_sequences.tar.gz.md5'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/busco.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/bbtools_err.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/bbtools_out.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/hmmsearch_err.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/hmmsearch_out.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/miniprot_align*_err.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/miniprot_align*out.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/miniprot_index*err.log'
    'assemblies/genome/busco/OG.ilmn.*.v129mh.busco.*.logs/miniprot_index*out.log'
    'assemblies/genome/busco/busco_sequences/*.fasta'
    'assemblies/genome/busco/busco_sequences/*.bed'
    'assemblies/genome/bwa'
    'assemblies/genome/bwa/OG.sorted.bam'
    'assemblies/genome/bwa/OG-sn_results.tsv'
    'assemblies/genome/gfastats'
    'assemblies/genome/gfastats/*assembly_summary'
    'assemblies/genome/seqkit'
    'assemblies/genome/seqkit/OG.ilmn.*.v129mh.seqkit_stats.tsv'
    # Per-sample MultiQC report (copied to S3 by the final backup)
    'multiqc'
    'multiqc/OG_multiqc_report.html'
)
# The meryl database is 64 shards, 0x000000 to 0x111111, each with Data and Index.
for shard in {0,1}{0,1}{0,1}{0,1}{0,1}{0,1}; do
    expected+=("kmers/OG.ilmn.*.meryl/0x${shard}.merylData" "kmers/OG.ilmn.*.meryl/0x${shard}.merylIndex")
done

# Row labels, written once so every sample column lines up with them.
{
    echo "Files"
    for pattern in "${expected[@]}"; do
        echo "<OG>/${pattern//OG/<OG>}"
    done
} > FileCheck/00_labels.tsv

for OGdir in "$rundir"/*; do
    [ -d "$OGdir" ] || continue
    OG=$(basename "$OGdir")
    echo "Checking $OG"

    {
        echo "$OG"
        for pattern in "${expected[@]}"; do
            matches=()
            # OG in a pattern is the sample ID; the remaining * are globs.
            while IFS= read -r match; do
                matches+=("$match")
            done < <(compgen -G "$OGdir/${pattern//OG/$OG}")

            if [ "${#matches[@]}" -gt 0 ]; then
                # Touch each match so an unreadable or offline file shows up here and
                # not part way through the backup.
                for match in "${matches[@]}"; do
                    if [ -f "$match" ]; then
                        head -c 1 "$match" >/dev/null 2>&1 || true
                    elif [ -d "$match" ]; then
                        ls "$match" >/dev/null 2>&1 || true
                    fi
                done
                echo "Found ${matches[0]}"
            else
                echo "Missing"
            fi
        done
    } > "FileCheck/$OG.files.tsv"
done

paste FileCheck/00_labels.tsv $(ls FileCheck/*.files.tsv) > _FileCheck.tsv
echo "Wrote _FileCheck.tsv ($(( $(wc -l < FileCheck/00_labels.tsv) - 1 )) expected files per sample)"
missing=$(tail -n +2 _FileCheck.tsv | grep -c "Missing" || true)
echo "Rows with at least one missing file: $missing"
