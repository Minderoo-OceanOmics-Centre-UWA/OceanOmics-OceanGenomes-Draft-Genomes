process COMPILE_JSON_TO_CSV {
    tag "Compiling JSON files to CSV"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    path json_files, stageAs: "input_jsons/*"

    output:
    path "*.csv"       , emit: csv
    path "*_mqc.yaml"  , emit: multiqc
    path "21_compile_json_to_csv.tool_params_mqcrow.html", emit: tool_params
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: 'compiled_results'
    """
    compile_json_to_csv.py \\
        --input_dir input_jsons \\
        --output ${prefix}.csv \\
        $args

    python - <<'PY'
import csv
import html
from pathlib import Path

csv_path = Path("${prefix}.csv")
rows = list(csv.DictReader(csv_path.open()))

columns = [
    ("sample_id", "Sample"),
    ("estimated_genome_size", "Genome Size (bp)"),
    ("coverage_before_filtering", "Coverage Before (x)"),
    ("theoretical_coverage", "Coverage After (x)"),
    ("data_retention_percent", "Data Retention (%)"),
    ("filtering_efficiency", "Bases Removed (%)"),
    ("heterozygosity_max", "Het. Max (%)"),
    ("model_fit_max", "Model Fit Max (%)"),
    ("coverage_status", "Status"),
]

def fmt_value(key, value):
    if value in (None, ""):
        return "NA"
    if key == "estimated_genome_size":
        try:
            return f"{int(float(value)):,}"
        except ValueError:
            return str(value)
    if key == "sample_id":
        return str(value)
    if key == "coverage_status":
        return str(value)
    try:
        return f"{float(value):.2f}"
    except ValueError:
        return str(value)

table_lines = [
    "<p>Compiled coverage metrics derived from FastP and GenomeScope2 for all samples.</p>",
    "<table class=\\"table table-condensed\\">",
    "<thead><tr>" + "".join(f"<th>{html.escape(label)}</th>" for _, label in columns) + "</tr></thead>",
    "<tbody>",
]

for row in rows:
    cells = []
    for key, _label in columns:
        cells.append(f"<td>{html.escape(fmt_value(key, row.get(key)))}</td>")
    table_lines.append("<tr>" + "".join(cells) + "</tr>")

table_lines.extend(["</tbody>", "</table>"])

mqc_lines = [
    "id: 'nf-core-oceangenomesdraftgenomes-coverage-summary'",
    "description: 'Compiled theoretical sequencing coverage derived from FastP and GenomeScope2.'",
    "section_name: 'Coverage Summary'",
    "plot_type: 'html'",
    "data: |",
]
mqc_lines.extend([f"    {line}" for line in table_lines])

Path("${prefix}_mqc.yaml").write_text("\\n".join(mqc_lines) + "\\n")
Path("21_compile_json_to_csv.tool_params_mqcrow.html").write_text(
    "<tr><td>Compile JSON to CSV</td><td><samp>compile_json_to_csv.py --input_dir input_jsons --output ${prefix}.csv</samp></td><td>Aggregates per-sample coverage JSON files into the Coverage Summary table.</td></tr>\\n"
)
PY

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version | sed 's/Python //g')" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: 'compiled_results'
    """
    touch ${prefix}.csv
    touch ${prefix}_mqc.yaml
    cat <<-END_TOOL_PARAMS > 21_compile_json_to_csv.tool_params_mqcrow.html
    <tr><td>Compile JSON to CSV</td><td><samp>compile_json_to_csv.py --input_dir input_jsons --output ${prefix}.csv</samp></td><td>Aggregates per-sample coverage JSON files into the Coverage Summary table.</td></tr>
    END_TOOL_PARAMS

    printf '"%s":\\n    python: %s\\n' \\
        "${task.process}" \\
        "\$(python --version | sed 's/Python //g')" > versions.yml
    """
}
