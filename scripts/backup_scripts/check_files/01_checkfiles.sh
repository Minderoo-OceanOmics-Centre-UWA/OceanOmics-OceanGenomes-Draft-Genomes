#!/bin/bash
# Load in the configfile
. ../configfile.txt

# everytime you run this new you need to remove the FileCheck directory if it exists and the _FileCheck_report.tsv file
## Use this variable to override configfile settings if needed, otherwise comment them out
# rundir=/scratch/pawsey0964/tpeirce/NOVA_251215_AD/draftgenomes
rm -r FileCheck
rm _FileCheck.tsv

mkdir -p FileCheck

for OGdir in $rundir/*; do
# Set OG and DATE variables (replace with actual values)
    OG=$(basename $OGdir)
    TAX=vert #$(cat $results/taxon.txt | grep -w $OG | awk -F'\t' '{print substr($2, 1, 4)}')
    echo *
    # List of expected files and directories
    file_list=(
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.busco_sequences.tar.gz
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.busco_sequences.tar.gz.md5
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.full_table.tsv
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/bbtools_err.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/bbtools_out.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/busco.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/hmmsearch_err.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/hmmsearch_out.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/miniprot_align*_err.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/miniprot_align*out.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/miniprot_index*err.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs/miniprot_index*out.log
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.logs
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.missing_busco_list.tsv
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.short_summary.json
        $OGdir/assemblies/genome/busco/$OG.ilmn.*.v129mh.busco.*.short_summary.txt
        $OGdir/assemblies/genome/busco/busco_sequences/*.bed
        $OGdir/assemblies/genome/busco/busco_sequences/*.fasta
        $OGdir/assemblies/genome/busco
        $OGdir/assemblies/genome/bwa/$OG.sorted.bam
        $OGdir/assemblies/genome/bwa/$OG-sn_results.tsv
        $OGdir/assemblies/genome/bwa
        $OGdir/assemblies/genome/NCBI/adaptor/*cleaned_sequences.fa.gz
        $OGdir/assemblies/genome/NCBI/adaptor/*fcs_adaptor.log
        $OGdir/assemblies/genome/NCBI/adaptor/*fcs_adaptor_report.txt
        $OGdir/assemblies/genome/NCBI/adaptor/*pipeline_args.yaml
        $OGdir/assemblies/genome/NCBI/adaptor/*skipped_trims.jsonl
        $OGdir/assemblies/genome/NCBI/adaptor
        $OGdir/assemblies/genome/NCBI/$OG*.contam.fasta
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.*.contig_count_500bp.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.*.filter_report.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.*.review_scaffolds_1kb.txt
        $OGdir/assemblies/genome/NCBI/*.fcs_gx_report.txt
        $OGdir/assemblies/genome/NCBI/*.taxonomy.rpt
        $OGdir/assemblies/genome/NCBI/*.summary.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.*.v129mh.rc.fasta
        $OGdir/assemblies/genome/NCBI
        $OGdir/assemblies/genome/$OG.ilmn.*.adaptor-contam.fasta
        $OGdir/assemblies/genome/$OG.ilmn.*.rmadapt.fasta
        $OGdir/assemblies/genome/$OG.ilmn.*.v129mh.fasta
        $OGdir/assemblies/genome/$OG.ilmn.*.v129mh.fna
        $OGdir/assemblies/genome/tiara/log_$OG.ilmn.*.tiara.txt
        $OGdir/assemblies/genome/tiara/$OG.ilmn.*.tiara.contig_removal.txt
        $OGdir/assemblies/genome/tiara/$OG.ilmn.*.tiara.txt
        $OGdir/assemblies/genome/tiara/$OG.ilmn.*.tiara_filter_summary.txt
        $OGdir/assemblies/genome/tiara
        $OGdir/assemblies/genome/gfastats/*assembly_summary
        $OGdir/assemblies/genome/gfastats/*.fa
        $OGdir/assemblies/genome/gfastats
        $OGdir/assemblies/genome
        $OGdir/assemblies
        $OGdir/fastp/*.fastp.json
        $OGdir/fastp/$OG.ilmn.*.fastp.html
        $OGdir/fastp/$OG.ilmn.*.fastp.log
        $OGdir/fastp/$OG.ilmn.*.R1.fastq.gz
        $OGdir/fastp/$OG.ilmn.*.R2.fastq.gz
        $OGdir/fastp/fastqc/$OG.ilmn.*.R1_fastqc.html
        $OGdir/fastp/fastqc/$OG.ilmn.*.R1_fastqc.zip
        $OGdir/fastp/fastqc/$OG.ilmn.*.R2_fastqc.html
        $OGdir/fastp/fastqc/$OG.ilmn.*.R2_fastqc.zip
        $OGdir/fastp/fastqc
        $OGdir/fastp
        $OGdir/coverage/*.json
        $OGdir/coverage/*.txt
        $OGdir/coverage
        $OGdir/kmers/$OG.ilmn.*genomescope/$OG.ilmn.*linear_plot.png
        $OGdir/kmers/$OG.ilmn.*genomescope/$OG.ilmn.*log_plot.png
        $OGdir/kmers/$OG.ilmn.*genomescope/$OG.ilmn.*model.txt
        $OGdir/kmers/$OG.ilmn.*genomescope/$OG.ilmn.*summary.txt
        $OGdir/kmers/$OG.ilmn.*genomescope/*transformed_linear_plot.png
        $OGdir/kmers/$OG.ilmn.*genomescope/$OG.ilmn.*transformed_log_plot.png
        $OGdir/kmers/$OG.ilmn.*genomescope
        $OGdir/kmers/$OG.ilmn.*.meryl/merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl
        $OGdir/kmers/$OG.ilmn.*.meryl.hist
        $OGdir/kmers/*merqury.completeness.stats
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.$OG.ilmn.*.v129mh.fna.qv
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.$OG.ilmn.*.v129mh.fna.spectra-cn.fl.png
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.$OG.ilmn.*.v129mh.fna.spectra-cn.ln.png
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.$OG.ilmn.*.v129mh.fna.spectra-cn.st.png
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.qv
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.spectra-asm.fl.png
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.spectra-asm.ln.png
        $OGdir/kmers/$OG.ilmn.*.v129mh.merqury.spectra-asm.st.png
        $OGdir/kmers
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x000111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x001111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x010111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x011111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x100111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x101111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x110111.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111000.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111000.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111001.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111001.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111010.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111010.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111011.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111011.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111100.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111100.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111101.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111101.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111110.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111110.merylIndex
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111111.merylData
        $OGdir/kmers/$OG.ilmn.*.meryl/0x111111.merylIndex
    )

    # Iterate over the list and check if each file or directory exists
    echo -e $OG > FileCheck/$OG.files.tsv
    for file in "${file_list[@]}"; do
        
        if [ -e "$file" ]; then
            echo -e "$file Found"
            if [ -f "$file" ]; then
                head -c 1 "$file" >/dev/null 2>&1 || true
            elif [ -d "$file" ]; then
                ls "$file" >/dev/null 2>&1 || true
            fi
        else
            echo -e "$file Missing"
        fi
    done >> FileCheck/$OG.files.tsv
done



paste file_list.tsv FileCheck/* > _FileCheck.tsv