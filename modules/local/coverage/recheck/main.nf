// Re-check the GenomeScope genome size against the decontaminated assembly.
//
// This is the one reliability check that cannot be made from the k-mer histogram alone.
// A GenomeScope model can converge cleanly onto something that is not the genome: OG3037
// fitted at 86% and reported 209 Mb against a 1.26 Gb assembly, and no fit statistic
// objected. Comparing the estimate to an assembly built from the same reads does object.
//
// It lives here, and not in CALCULATE_SEQUENCING_COVERAGE, because that process runs
// before MEGAHIT -- at that point there is no assembly to compare against. The summary
// this writes supersedes the provisional one written there.
process RECHECK_GENOME_SIZE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    // The provisional summary is staged into provisional/ rather than the task root, and
    // the name it is staged under is the name this process writes its own output as.
    // Staged in the root, that collision cost the whole NOVA_260909_LA re-run: Nextflow
    // excludes input files from output matching, so `*_coverage_summary.json` matched
    // nothing and every task failed with MissingFileException under errorStrategy ignore
    // (silently -- exit 0, empty stderr), while the write itself resolved through the
    // input symlink and overwrote CALCULATE_SEQUENCING_COVERAGE's own work-dir copy.
    tuple val(meta), path(coverage_json, stageAs: 'provisional/*'), path(seqkit_stats), path(kmer_depth_json)

    output:
    tuple val(meta), path("*_coverage_summary.json"), emit: coverage_json
    tuple val(meta), path("*_genome_size_recheck_mqc.yaml"), emit: multiqc
    tuple val(meta), path("28_recheck_genome_size.tool_params_mqcrow.html"), emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

script:
def prefix = task.ext.prefix ?: "${meta.prefix}"
"""
#!/usr/bin/env python3

import html
import json
import platform
import subprocess

merged = json.loads(subprocess.check_output([
    "genomescope_reliability.py",
    "--recheck-assembly",
    "--coverage-json", "${coverage_json}",
    "--seqkit-stats", "${seqkit_stats}",
    "--kmer-depth-json", "${kmer_depth_json}",
    "--min-assembly-ratio", "${params.genomescope_min_assembly_ratio}",
    "--max-assembly-ratio", "${params.genomescope_max_assembly_ratio}",
    "--max-lambda-ratio", "${params.kmer_depth_max_ratio}",
    "--min-seedable-lambda", "${params.kmer_depth_min_seedable}",
], text=True))

with open("${prefix}_coverage_summary.json", "w") as f:
    json.dump(merged, f, indent=2)

size = merged.get("estimated_genome_size") or 0
assembly_size = merged.get("assembly_size") or 0
lam = merged.get("lambda_depth")
kmercov = merged.get("kmercov")
host_ratio = merged.get("busco_depth_over_assembly_depth")
host_size = merged.get("host_assembly_size") or 0
host_fraction = merged.get("host_assembly_fraction")
basis = merged.get("assembly_ratio_basis") or "total"
denominator = merged.get("assembly_ratio_denominator") or 0
rows = [
    ("Estimated genome size", f"{size:,} bp" if size else "NA"),
    ("Decontaminated assembly size", f"{assembly_size:,} bp" if assembly_size else "NA"),
    # How much of this assembly is the animal. A lower bound: a symbiont sitting at the
    # host's own depth counts as host, and a high-copy host repeat counts as symbiont.
    ("Host assembly size (lower bound)", f"{host_size:,} bp" if host_size else "NA"),
    ("Host fraction of assembly (lower bound)",
     f"{host_fraction * 100:.1f}%" if host_fraction is not None else "NA"),
    ("Depth partition", merged.get("partition_status") or "NA"),
    # Named, because GenomeScope estimates the host genome and the denominator is the
    # host mass wherever the partition succeeded, not the whole assembly.
    (f"Estimate / {basis} assembly",
     f"{size / denominator:.2f}" if size and denominator else "NA"),
    ("Fitted k-mer coverage (GenomeScope)", f"{kmercov:.1f}x" if kmercov else "NA"),
    ("Measured k-mer coverage (BUSCO read depth)", f"{lam:.1f}x" if lam else "NA"),
    ("Fitted / measured", f"{kmercov / lam:.2f}" if kmercov and lam else "NA"),
    ("Host single-copy depth / assembly-wide depth",
     f"{host_ratio:.2f}" if host_ratio is not None else "NA"),
    ("Genome size estimate reliable", "yes" if merged.get("genome_size_reliable") else "no"),
    ("Genome size warnings", merged.get("genome_size_flags") or "none"),
    ("Coverage assessment", merged.get("coverage_status", "NA")),
]

table = ["<table class=\\"table table-condensed\\">", "<tbody>"]
for label, value in rows:
    table.append(f"<tr><th>{html.escape(str(label))}</th><td>{html.escape(str(value))}</td></tr>")
table.extend(["</tbody>", "</table>"])

mqc_lines = [
    "id: 'nf-core-oceangenomesdraftgenomes-genome-size-recheck'",
    "description: 'GenomeScope genome size estimate checked against the decontaminated assembly size.'",
    "section_name: 'Genome Size Cross-check'",
    "plot_type: 'html'",
    "data: |",
]
mqc_lines.extend([f"    {line}" for line in table])

with open("${prefix}_genome_size_recheck_mqc.yaml", "w") as f:
    f.write("\\n".join(mqc_lines) + "\\n")

tool_params_html = (
    "<tr><td>Recheck Genome Size</td>"
    "<td><samp>genomescope_reliability.py --recheck-assembly "
    "--min-assembly-ratio ${params.genomescope_min_assembly_ratio} "
    "--max-assembly-ratio ${params.genomescope_max_assembly_ratio}</samp></td>"
    "<td>Flags a GenomeScope size estimate that disagrees with the assembly built from "
    "the same reads, or a fitted k-mer coverage that disagrees with read depth over "
    "single-copy BUSCO genes, and republishes ${prefix}_coverage_summary.json as the "
    "final verdict.</td></tr>"
)

with open("28_recheck_genome_size.tool_params_mqcrow.html", "w") as f:
    f.write(tool_params_html + "\\n")

with open("versions.yml", "w") as vf:
    # Indentation is YAML-significant; keep the two spaces before "python:"
    vf.write(f'"${task.process}":\\n')
    vf.write(f"  python: {platform.python_version()}\\n")
"""
}
