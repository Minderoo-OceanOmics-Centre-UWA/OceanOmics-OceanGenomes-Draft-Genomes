process EXTRACT_BUSCO_SEQUENCES {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
      'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    tuple val(meta), path(busco_table), path(genome_fasta)

    output:
    tuple val(meta), path("*.busco_sequences.fasta"), emit: fasta
    tuple val(meta), path("*.busco_coordinates.bed"),  emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Build BED (0-based) and put the full desired FASTA header in column 4.
    # BUSCO table columns (tab): 1=busco_id 2=status 3=sequence 4=start 5=end 6=strand 7=score 8=length 9=url 10=description
    awk -v MID='${meta.id}' 'BEGIN{FS=OFS="\\t"}
        /^#/ || NF==0 || \$1==\"Busco id\" { next }     # skip comments/header/blank
        \$2==\"Missing\" { next }                       # skip Missing
        NF<7 { next }                                # need at least first 7 fields

        {
            seq = \$3
            s1  = \$4 + 0
            e1  = \$5 + 0
            st  = \$6
            sc  = \$7

            # Clamp any 0/negative coordinates to 1 (BUSCO sometimes emits 0)
            if (s1 < 1) s1 = 1
            if (e1 < 1) e1 = 1

            # Put in genomic order (BED ignores strand for coords)
            start1 = (s1 <= e1) ? s1 : e1
            end1   = (s1 <= e1) ? e1 : s1

            # Convert to BED: 0-based start, 1-based end
            s0 = start1 - 1
            e0 = end1

            # Ensure end > start (at least 1 bp)
            if (e0 <= s0) e0 = s0 + 1

            len  = (NF>=8 ? \$8 : "NA")
            desc = (NF>=10 ? \$10 : (NF>=9 ? \$9 : ""))
            gsub(/\\t/, " ", desc)

            header = MID "|" \$1 "|" \$2 "|" seq "|" \$4 "-" \$5 "|" st "|score=" sc "|length=" len "|" desc

            # BED6: chrom start end name score strand
            print seq, s0, e0, header, sc, st
        }' ${busco_table} > ${prefix}.busco_coordinates.bed

    if [ ! -s ${prefix}.busco_coordinates.bed ]; then
        echo "No complete or fragmented BUSCO genes found to extract"
        : > ${prefix}.busco_sequences.fasta
    else
        # Use nameOnly→ header = BED col 4; -s respects strand
        bedtools getfasta \\
        -fi ${genome_fasta} \\
        -bed ${prefix}.busco_coordinates.bed \\
        -s -nameOnly \\
        -fo ${prefix}.busco_sequences.fasta
    fi

    cat > versions.yml <<-YML
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    YML
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.assembly_prefix}"

    """
    #!/usr/bin/env bash
    set -euo pipefail

    # Create empty output files
    : > ${prefix}.busco_coordinates.bed
    : > ${prefix}.busco_sequences.fasta

    cat > versions.yml <<-YML
    "${task.process}":
        bedtools: "2.31.1"
    YML
    """
}