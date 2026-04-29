// Calculate theoretical sequencing coverage from reads and estimated genome size
process CALCULATE_SEQUENCING_COVERAGE {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fastp_json), path(genomescope_summary)

    output:
    tuple val(meta), path("*_sequencing_coverage.txt"), emit: coverage_report
    path("*_coverage_summary.json"), emit: coverage_json
    tuple val(meta), path("*_coverage_summary_mqc.yaml"), emit: multiqc
    tuple val(meta), path("20_calculate_sequencing_coverage.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
def prefix = task.ext.prefix ?: "${meta.prefix}"
"""
#!/usr/bin/env python3

import json
import re
import platform
import html

# Read FastP JSON report to get sequencing statistics
with open("${fastp_json}", 'r') as f:
    fastp_data = json.load(f)

# Extract read statistics from FastP
total_reads_before = fastp_data['summary']['before_filtering']['total_reads']
total_bases_before = fastp_data['summary']['before_filtering']['total_bases']
total_reads_after = fastp_data['summary']['after_filtering']['total_reads']
total_bases_after = fastp_data['summary']['after_filtering']['total_bases']

# Read GenomeScope summary to extract estimated genome size
estimated_genome_size = 0
genome_size_min = 0
genome_size_max = 0
heterozygosity_min = 0
heterozygosity_max = 0
model_fit_min = 0
model_fit_max = 0

with open("${genomescope_summary}", 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('Genome Haploid Length'):
            # Extract number from format like "Genome Haploid Length    4,123,456 bp"
            size_matches = re.findall(r'([0-9,]+) bp', line)
            if size_matches:
                genome_size_min = int(size_matches[0].replace(',', ''))
                genome_size_max = int(size_matches[-1].replace(',', ''))
                estimated_genome_size = genome_size_max
        elif re.match(r"^Heterozyg(?:osity|ous)", line, flags=re.I):
            # Handles "Heterozygosity" or "Heterozygous (ab)"
            het_matches = re.findall(r'([0-9.]+)%', line)
            if het_matches:
                heterozygosity_min = float(het_matches[0])
                heterozygosity_max = float(het_matches[-1])
        elif line.startswith('Model Fit'):
            fit_matches = re.findall(r'([0-9.]+)%', line)
            if fit_matches:
                model_fit_min = float(fit_matches[0])
                model_fit_max = float(fit_matches[-1])

# Calculate coverage metrics
coverage_before = total_bases_before / estimated_genome_size if estimated_genome_size > 0 else 0
coverage_after = total_bases_after / estimated_genome_size if estimated_genome_size > 0 else 0

# Calculate other useful metrics
bases_removed = total_bases_before - total_bases_after
reads_removed = total_reads_before - total_reads_after
filtering_efficiency = (bases_removed / total_bases_before * 100) if total_bases_before > 0 else 0
data_retention = total_bases_after / total_bases_before * 100 if total_bases_before > 0 else 0

if coverage_after >= 50:
    coverage_status = "EXCELLENT"
    coverage_recommendation = "Suitable for high-quality genome assembly and variant calling."
elif coverage_after >= 30:
    coverage_status = "GOOD"
    coverage_recommendation = "Suitable for genome assembly and most downstream analyses."
elif coverage_after >= 20:
    coverage_status = "ADEQUATE"
    coverage_recommendation = "Suitable for basic genome assembly."
elif coverage_after >= 10:
    coverage_status = "LOW"
    coverage_recommendation = "Draft assembly may be possible, but additional sequencing may help."
else:
    coverage_status = "INSUFFICIENT"
    coverage_recommendation = "Additional sequencing is strongly recommended."

# Write detailed report
with open("${prefix}_sequencing_coverage.txt", 'w') as f:
    f.write("SEQUENCING COVERAGE ANALYSIS\\n")
    f.write("=" * 50 + "\\n\\n")
    f.write(f"Sample ID: ${meta.id}\\n\\n")
    
    f.write("GENOME SIZE ESTIMATION (GenomeScope):\\n")
    f.write(f"  Estimated genome size: {estimated_genome_size:,} bp\\n")
    if genome_size_min and genome_size_max:
        f.write(f"  Genome size range: {genome_size_min:,} - {genome_size_max:,} bp\\n")
    f.write(f"  Estimated heterozygosity: {heterozygosity_min:.4f}% - {heterozygosity_max:.4f}%\\n")
    if model_fit_min or model_fit_max:
        mf_min = f"{model_fit_min:.4f}%" if model_fit_min else "NA"
        mf_max = f"{model_fit_max:.4f}%" if model_fit_max else "NA"
        f.write(f"  Model fit: {mf_min} - {mf_max}\\n")
    f.write("\\n")
    
    f.write("SEQUENCING STATISTICS (FastP):\\n")
    f.write("  Before filtering:\\n")
    f.write(f"    Total reads: {total_reads_before:,}\\n")
    f.write(f"    Total bases: {total_bases_before:,} bp\\n")
    f.write(f"    Theoretical coverage: {coverage_before:.1f}x\\n\\n")
    
    f.write("  After filtering:\\n")
    f.write(f"    Total reads: {total_reads_after:,}\\n")
    f.write(f"    Total bases: {total_bases_after:,} bp\\n")
    f.write(f"    Theoretical coverage: {coverage_after:.1f}x\\n\\n")
    
    f.write("FILTERING SUMMARY:\\n")
    f.write(f"  Reads removed: {reads_removed:,} ({reads_removed/total_reads_before*100:.1f}%)\\n")
    f.write(f"  Bases removed: {bases_removed:,} ({filtering_efficiency:.1f}%)\\n")
    f.write(f"  Data retention: {data_retention:.1f}%\\n\\n")
    
    f.write("COVERAGE ASSESSMENT:\\n")
    f.write(f"  Status: {coverage_status} coverage ({coverage_after:.1f}x)\\n")
    f.write(f"  Recommendation: {coverage_recommendation}\\n")

# Create JSON summary for downstream processes
summary_data = {
    "sample_id": "${meta.id}",
    "estimated_genome_size": estimated_genome_size,
    "genome_size_min": genome_size_min,
    "genome_size_max": genome_size_max,
    "model_fit_min": model_fit_min,
    "model_fit_max": model_fit_max,
    "heterozygosity_min": heterozygosity_min,
    "heterozygosity_max": heterozygosity_max,
    "total_bases_after_filtering": total_bases_after,
    "total_reads_after_filtering": total_reads_after,
    "theoretical_coverage": coverage_after,
    "coverage_before_filtering": coverage_before,
    "filtering_efficiency": filtering_efficiency,
    "data_retention_percent": data_retention,
    "coverage_status": coverage_status,
    "coverage_recommendation": coverage_recommendation
}

with open("${prefix}_coverage_summary.json", 'w') as f:
    json.dump(summary_data, f, indent=2)

coverage_rows = [
    ("Sample ID", "${meta.id}"),
    ("Estimated genome size", f"{estimated_genome_size:,} bp" if estimated_genome_size else "NA"),
    ("Genome size range", f"{genome_size_min:,} - {genome_size_max:,} bp" if genome_size_min and genome_size_max else "NA"),
    ("Heterozygosity", f"{heterozygosity_min:.4f}% - {heterozygosity_max:.4f}%"),
    ("Model fit", f"{model_fit_min:.4f}% - {model_fit_max:.4f}%" if model_fit_min or model_fit_max else "NA"),
    ("Coverage before filtering", f"{coverage_before:.1f}x"),
    ("Coverage after filtering", f"{coverage_after:.1f}x"),
    ("Reads after filtering", f"{total_reads_after:,}"),
    ("Bases after filtering", f"{total_bases_after:,} bp"),
    ("Filtering efficiency", f"{filtering_efficiency:.1f}%"),
    ("Data retention", f"{data_retention:.1f}%"),
    ("Coverage assessment", f"{coverage_status} ({coverage_after:.1f}x)"),
    ("Recommendation", coverage_recommendation),
]

coverage_html = ["<table class=\\"table table-condensed\\">", "<tbody>"]
for label, value in coverage_rows:
    coverage_html.append(
        f"<tr><th>{html.escape(str(label))}</th><td>{html.escape(str(value))}</td></tr>"
    )
coverage_html.extend(["</tbody>", "</table>"])

mqc_lines = [
    "id: 'nf-core-oceangenomesdraftgenomes-coverage-summary'",
    "description: 'Theoretical sequencing coverage derived from FastP and GenomeScope2.'",
    "section_name: 'Coverage Summary'",
    "plot_type: 'html'",
    "data: |",
]
mqc_lines.extend([f"    {line}" for line in coverage_html])

with open("${prefix}_coverage_summary_mqc.yaml", 'w') as f:
    f.write("\\n".join(mqc_lines) + "\\n")

tool_params_html = (
    "<tr><td>Calculate Sequencing Coverage</td>"
    "<td><samp>inline Python coverage summary using FastP JSON and GenomeScope2 summary inputs</samp></td>"
    "<td>Computes theoretical pre/post-filtering coverage and writes ${prefix}_coverage_summary.json.</td></tr>"
)

with open("20_calculate_sequencing_coverage.tool_params_mqcrow.html", "w") as f:
    f.write(tool_params_html + "\\n")

python_ver = platform.python_version()

with open("versions.yml", "w") as vf:
    # Indentation is YAML-significant; keep the two spaces before "python:"
    vf.write(f'"${task.process}":\\n')
    vf.write(f"  python: {python_ver}\\n")
"""
}
