// Calculate theoretical sequencing coverage from reads and estimated genome size
process CALCULATE_SEQUENCING_COVERAGE {
    tag "$meta.id"
    label 'process_low'
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(fastp_json), path(genomescope_summary), path(genomescope_model), path(meryl_hist)

    output:
    tuple val(meta), path("*_sequencing_coverage.txt"), emit: coverage_report
    path("*_coverage_summary.json"), emit: coverage_json
    tuple val(meta), path("*_coverage_summary.json"), emit: coverage_json_meta
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
import platform
import html
import subprocess

# Read FastP JSON report to get sequencing statistics
with open("${fastp_json}", 'r') as f:
    fastp_data = json.load(f)

# Extract read statistics from FastP
total_reads_before = fastp_data['summary']['before_filtering']['total_reads']
total_bases_before = fastp_data['summary']['before_filtering']['total_bases']
total_reads_after = fastp_data['summary']['after_filtering']['total_reads']
total_bases_after = fastp_data['summary']['after_filtering']['total_bases']

# --- GenomeScope reliability gate -------------------------------------------
#
# Every metric below this point divides by the estimated genome size, so a bad
# GenomeScope fit silently turns into a confident-looking coverage verdict. The gate
# itself lives in bin/genomescope_reliability.py so it can be unit tested and re-run
# offline over published results; see that file for why each check exists.
#
# No assembly cross-check here: this process runs before MEGAHIT, so there is no
# assembly to compare against yet. RECHECK_GENOME_SIZE adds that flag during QC and
# publishes the authoritative summary. What this writes is provisional.
gate = json.loads(subprocess.check_output([
    "genomescope_reliability.py",
    "--summary", "${genomescope_summary}",
    "--model", "${genomescope_model}",
    "--histogram", "${meryl_hist}",
    "--min-model-fit", "${params.genomescope_min_model_fit}",
    "--max-bias", "${params.genomescope_max_model_bias}",
    "--min-peak-kmercov-ratio", "${params.genomescope_min_peak_kmercov_ratio}",
    "--max-peak-kmercov-ratio", "${params.genomescope_max_peak_kmercov_ratio}",
    "--peak-search-max", "${params.genomescope_peak_search_max}",
], text=True))

estimated_genome_size = gate["estimated_genome_size"] or 0
genome_size_flags = gate["genome_size_flags"]
genome_size_reliable = gate["genome_size_reliable"]
# GenomeScope's two Model Fit values are not bounds: the first is the share of ALL
# retained k-mers the model accounts for, the second the fit within the model's own
# region. See bin/genomescope_reliability.py.
model_fit_allkmers, model_fit_full = gate["model_fit_allkmers"], gate["model_fit_full"]
heterozygosity_min = gate["heterozygosity_min"]
heterozygosity_max = gate["heterozygosity_max"]
peak_coverage, kmercov = gate["kmer_peak_coverage"], gate["kmercov"]

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

if not genome_size_reliable:
    # Keep the computed numbers -- they are still the best available -- but do not
    # present them as a coverage grade.
    coverage_status = "UNRELIABLE_GENOME_SIZE_ESTIMATE"
    coverage_recommendation = (
        "GenomeScope could not fit this sample (" + ", ".join(genome_size_flags) + "), so the "
        "genome size estimate and every coverage figure derived from it are unreliable. "
        "Assess coverage against the decontaminated assembly size instead."
    )


def pct(value):
    return "NA" if value is None else f"{value:.4f}%"


# Write detailed report
with open("${prefix}_sequencing_coverage.txt", 'w') as f:
    f.write("SEQUENCING COVERAGE ANALYSIS\\n")
    f.write("=" * 50 + "\\n\\n")
    f.write(f"Sample ID: ${meta.id}\\n\\n")

    f.write("GENOME SIZE ESTIMATION (GenomeScope):\\n")
    f.write(f"  Estimated genome size: {estimated_genome_size:,} bp\\n")
    if gate["genome_size_min"] and gate["genome_size_max"]:
        f.write(f"  Genome size range: {gate['genome_size_min']:,} - {gate['genome_size_max']:,} bp\\n")
    f.write(f"  Estimated heterozygosity: {pct(heterozygosity_min)} - {pct(heterozygosity_max)}\\n")
    f.write(f"  Model fit (full model): {pct(model_fit_full)}\\n")
    f.write(f"  K-mers modelled (all k-mers): {pct(model_fit_allkmers)}\\n")
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

    f.write(f"  K-mer peak coverage: {peak_coverage}x\\n" if peak_coverage else "  K-mer peak coverage: none detected\\n")
    f.write(f"  Fitted k-mer coverage (GenomeScope kmercov): {kmercov:.1f}x\\n" if kmercov else "  Fitted k-mer coverage: NA\\n")
    f.write(f"  Genome size estimate reliable: {'yes' if genome_size_reliable else 'no'}\\n")
    if genome_size_flags:
        f.write(f"  Genome size warnings: {', '.join(genome_size_flags)}\\n")
    f.write("\\n")

    f.write("COVERAGE ASSESSMENT:\\n")
    f.write(f"  Status: {coverage_status} coverage ({coverage_after:.1f}x)\\n")
    f.write(f"  Recommendation: {coverage_recommendation}\\n")

# Create JSON summary for downstream processes
summary_data = {
    "sample_id": "${meta.id}",
    "estimated_genome_size": estimated_genome_size,
    "genome_size_min": gate["genome_size_min"],
    "genome_size_max": gate["genome_size_max"],
    "model_fit_allkmers": model_fit_allkmers,
    "model_fit_full": model_fit_full,
    # Deprecated aliases, kept for one release. model_fit_min was never a minimum.
    "model_fit_min": model_fit_allkmers,
    "model_fit_max": model_fit_full,
    "heterozygosity_min": heterozygosity_min,
    "heterozygosity_max": heterozygosity_max,
    "total_bases_after_filtering": total_bases_after,
    "total_reads_after_filtering": total_reads_after,
    "theoretical_coverage": coverage_after,
    "coverage_before_filtering": coverage_before,
    "filtering_efficiency": filtering_efficiency,
    "data_retention_percent": data_retention,
    "coverage_status": coverage_status,
    "coverage_recommendation": coverage_recommendation,
    "genome_size_reliable": genome_size_reliable,
    "genome_size_flags": "; ".join(genome_size_flags),
    "kmer_peak_coverage": peak_coverage,
    # One name for this number. RECHECK_GENOME_SIZE, the MultiQC row and
    # draft_genome_stats.py all read "kmercov"; fitted_kmer_coverage is a deprecated alias
    # kept for one release so existing readers of published summaries keep working.
    "kmercov": kmercov,
    "fitted_kmer_coverage": kmercov,
    "assembly_cross_check": "pending"
}

with open("${prefix}_coverage_summary.json", 'w') as f:
    json.dump(summary_data, f, indent=2)

coverage_rows = [
    ("Sample ID", "${meta.id}"),
    ("Estimated genome size", f"{estimated_genome_size:,} bp" if estimated_genome_size else "NA"),
    ("Genome size range",
     f"{gate['genome_size_min']:,} - {gate['genome_size_max']:,} bp"
     if gate["genome_size_min"] and gate["genome_size_max"] else "NA"),
    ("Heterozygosity", f"{pct(heterozygosity_min)} - {pct(heterozygosity_max)}"),
    ("Model fit (full model)", pct(model_fit_full)),
    ("K-mers modelled (all k-mers)", pct(model_fit_allkmers)),
    ("Coverage before filtering", f"{coverage_before:.1f}x"),
    ("Coverage after filtering", f"{coverage_after:.1f}x"),
    ("Reads after filtering", f"{total_reads_after:,}"),
    ("Bases after filtering", f"{total_bases_after:,} bp"),
    ("Filtering efficiency", f"{filtering_efficiency:.1f}%"),
    ("Data retention", f"{data_retention:.1f}%"),
    ("K-mer peak coverage", f"{peak_coverage}x" if peak_coverage else "none detected"),
    ("Fitted k-mer coverage", f"{kmercov:.1f}x" if kmercov else "NA"),
    ("Genome size estimate reliable", "yes" if genome_size_reliable else "no"),
    ("Genome size warnings", ", ".join(genome_size_flags) if genome_size_flags else "none"),
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
    "<td><samp>genomescope_reliability.py --min-model-fit ${params.genomescope_min_model_fit} "
    "--max-bias ${params.genomescope_max_model_bias} "
    "--peak-search-max ${params.genomescope_peak_search_max}</samp></td>"
    "<td>Computes theoretical pre/post-filtering coverage and gates it on whether the "
    "GenomeScope fit can be trusted. Writes ${prefix}_coverage_summary.json.</td></tr>"
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
