# BUSCO Sequence Extraction Module

A Nextflow DSL2 module to extract BUSCO gene sequences from genome assemblies and create concatenated FASTA files with comprehensive metadata headers.

## Description

This module processes BUSCO `full_table.tsv` output files to:

1. Parse BUSCO gene coordinates and metadata
2. Extract corresponding sequences from genome FASTA files
3. Create a concatenated FASTA file with detailed headers containing all BUSCO information
4. Generate a BED file with BUSCO gene coordinates

## Usage

### Basic Usage

```nextflow
include { EXTRACT_BUSCO_SEQUENCES } from './extract_busco_sequences.nf'

workflow {
    input_ch = Channel.of([
        [ id: 'sample1' ],
        file('path/to/full_table.tsv'),
        file('path/to/genome.fasta')
    ])
    
    EXTRACT_BUSCO_SEQUENCES(input_ch)
}
```

### With Samplesheet

Create a CSV file (`samples.csv`) with the following structure:
```csv
sample_id,busco_table,genome_fasta
sample1,/path/to/sample1_full_table.tsv,/path/to/sample1_genome.fasta
sample2,/path/to/sample2_full_table.tsv,/path/to/sample2_genome.fasta
```

Then use it in your workflow:
```nextflow
params.input = 'samples.csv'

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header:true, sep:',')
        .map { row ->
            def meta = [id: row.sample_id]
            [meta, file(row.busco_table), file(row.genome_fasta)]
        }
        .set { samples_ch }
    
    EXTRACT_BUSCO_SEQUENCES(samples_ch)
}
```

## Input

- **meta** (map): Groovy Map containing sample information (e.g., `[id: 'sample1']`)
- **busco_table** (file): BUSCO `full_table.tsv` output file
- **genome_fasta** (file): Genome assembly FASTA file

## Output

- **fasta** (file): Concatenated FASTA file with extracted BUSCO sequences
- **bed** (file): BED file with BUSCO gene coordinates  
- **versions** (file): Software versions used

## Output Format

### FASTA Headers
Each sequence in the output FASTA file has a comprehensive header with the following format:
```
>sample_id|busco_id|status|sequence|start-end|strand|score=X.X|length=XXX|description
```

Example:
```
>sample1|903at7742|Complete|k141_1235172|10178-48080|-|score=7477.2|length=3867|HEPN domain
ATGCGATCGATCG...
>sample1|2178at7742|Complete|k141_966436|12918-63580|-|score=4418.5|length=2543|protein SZT2
CGATCGATCGATC...
```

### BED File
The BED file contains coordinates in standard format (0-based) with the following columns:
1. Chromosome/Contig name
2. Start position (0-based)
3. End position
4. BUSCO ID
5. Score
6. Strand

## Features

- **Handles all BUSCO statuses**: Processes Complete and Fragmented BUSCOs (skips Missing ones)
- **Strand-aware extraction**: Correctly extracts sequences on both forward and reverse strands
- **Comprehensive metadata**: Includes all BUSCO table information in sequence headers
- **Error handling**: Gracefully handles empty tables or missing sequences
- **Coordinate correction**: Automatically handles reverse strand coordinate issues

## Dependencies

- **bedtools** (≥2.31.1): For sequence extraction
- **seqtk** (≥1.4): Alternative sequence tools
- **python** (≥3.9): For header processing

## Configuration

You can customize the process behavior using Nextflow configuration:

```nextflow
process {
    withName: 'EXTRACT_BUSCO_SEQUENCES' {
        cpus = 4
        memory = '8.GB'
        time = '2.h'
        
        publishDir = [
            path: "${params.outdir}/busco_sequences",
            mode: 'copy'
        ]
    }
}
```

## Running the Module

```bash
# With conda
nextflow run example_workflow.nf -profile conda --input samples.csv

# With docker
nextflow run example_workflow.nf -profile docker --input samples.csv

# With singularity
nextflow run example_workflow.nf -profile singularity --input samples.csv
```

## Example Output Structure

```
results/
├── busco_sequences/
│   ├── sample1.busco_sequences.fasta
│   ├── sample1.busco_coordinates.bed
│   ├── sample2.busco_sequences.fasta
│   └── sample2.busco_coordinates.bed
└── concatenated/
    └── all_samples_busco_sequences.fasta  # If using collectFile
```

## Notes

- Only Complete and Fragmented BUSCO genes are extracted (Missing ones are skipped)
- The module automatically handles coordinate system differences (1-based to 0-based conversion)
- Sequences are extracted with proper strand orientation
- Headers contain all original BUSCO metadata for downstream analysis