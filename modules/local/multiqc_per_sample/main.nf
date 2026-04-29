process MULTIQC_PER_SAMPLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.27--pyhdfd78af_0' :
        'biocontainers/multiqc:1.27--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*")
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(replace_names)
    path(sample_names)

    output:
    tuple val(meta), path("*multiqc_report.html"), emit: report
    tuple val(meta), path("*_data")              , emit: data
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    """
    python <<'PY'
from collections import OrderedDict
from pathlib import Path
import html
import re

import yaml


def to_builtin(obj):
    if isinstance(obj, dict):
        return {key: to_builtin(value) for key, value in obj.items()}
    return obj


def simplify_group(group_name):
    group_name = str(group_name).strip().strip("'").strip('"')
    return group_name.rsplit(':', 1)[-1]


def load_version_file(path):
    try:
        data = yaml.safe_load(path.read_text()) or {}
    except Exception:
        return OrderedDict()

    parsed = OrderedDict()
    if not isinstance(data, dict):
        return parsed

    for group_name, tools in data.items():
        simple_group = simplify_group(group_name)
        if not isinstance(tools, dict):
            continue

        parsed[simple_group] = OrderedDict(
            (str(tool_name), str(version).strip())
            for tool_name, version in tools.items()
        )

    return parsed


for pattern in ('*_mqc_versions.yml', '*_mqc_versions.yaml'):
    for path in Path('.').rglob(pattern):
        path.unlink()

version_files = []
for pattern in ('software_versions.yml', 'software_versions.yaml', 'versions.yml', 'versions.yaml'):
    version_files.extend(
        sorted(
            path for path in Path('.').rglob(pattern)
            if '_mqc_versions.' not in path.name
        )
    )

version_files = sorted(
    version_files,
    key=lambda path: (0 if path.name.startswith('software_versions') else 1, str(path)),
)

versions_by_group = OrderedDict()
for version_file in version_files:
    for group_name, tools in load_version_file(version_file).items():
        merged_tools = versions_by_group.setdefault(group_name, OrderedDict())
        for tool_name, version in tools.items():
            merged_tools.setdefault(tool_name, version)

if versions_by_group:
    Path('software_versions_mqc_versions.yml').write_text(
        yaml.safe_dump(to_builtin(versions_by_group), sort_keys=False),
        encoding='utf-8',
    )

tool_to_group = {
    'bbmaprepair': 'BBMAP_REPAIR',
    'fastp': 'FASTP',
    'fastqc': 'FASTQC',
    'merylcount': 'MERYL_COUNT',
    'merylunionsum': 'MERYL_UNIONSUM',
    'merylhistogram': 'MERYL_HISTOGRAM',
    'genomescope2': 'GENOMESCOPE2',
    'calculatesequencingcoverage': 'CALCULATE_SEQUENCING_COVERAGE',
    'compilejsontocsv': 'COMPILE_JSON_TO_CSV',
    'megahit': 'MEGAHIT',
    'fcsgxrungx': 'FCSGX_RUNGX',
    'fcsgxcleangenome': 'FCSGX_CLEANGENOME',
    'bbmapfilterbyname': 'BBMAP_FILTERBYNAME',
    'fcsadaptor': 'FCS_FCSADAPTOR',
    'fcsgxcleangenomeadaptor': 'FSCSGX_CLEANGENOME_ADAPTOR',
    'tiara': 'TIARA_TIARA',
    'bbmapfilterbynametiara': 'BBMAP_FILTERBYNAME_TIARA',
    'busco': 'BUSCO_BUSCO',
    'bwamem2index': 'BWAMEM2_INDEX',
    'bwamem2mem': 'BWAMEM2_MEM',
    'merqury': 'MERQURY_MERQURY',
    'gfastats': 'GFASTATS',
}

tool_order = {
    'BBMAP_REPAIR': 5,
    'FASTP': 10,
    'FASTQC': 20,
    'MERYL_COUNT': 30,
    'MERYL_UNIONSUM': 40,
    'MERYL_HISTOGRAM': 50,
    'GENOMESCOPE2': 60,
    'CALCULATE_SEQUENCING_COVERAGE': 70,
    'COMPILE_JSON_TO_CSV': 75,
    'MEGAHIT': 80,
    'FCSGX_RUNGX': 90,
    'FCSGX_CLEANGENOME': 100,
    'BBMAP_FILTERBYNAME': 110,
    'FCS_FCSADAPTOR': 120,
    'FSCSGX_CLEANGENOME_ADAPTOR': 130,
    'TIARA_TIARA': 140,
    'BBMAP_FILTERBYNAME_TIARA': 150,
    'BUSCO_BUSCO': 160,
    'BWAMEM2_INDEX': 170,
    'BWAMEM2_MEM': 180,
    'MERQURY_MERQURY': 190,
    'GFASTATS': 200,
}

row_pattern = re.compile(
    r'<tr><td>(?P<tool>.*?)</td><td>(?P<params>.*?)</td><td>(?P<notes>.*?)</td></tr>',
    re.IGNORECASE | re.DOTALL,
)

tool_rows = []
for row_file in sorted(Path('.').rglob('*.tool_params_mqcrow.html')):
    match = row_pattern.search(row_file.read_text().strip())
    if not match:
        continue

    tool_html = match.group('tool').strip()
    params_html = match.group('params').strip()
    notes_html = match.group('notes').strip()
    tool_key = re.sub(r'[^a-z0-9]+', '', html.unescape(tool_html).lower())

    version_html = 'NA'
    group_name = tool_to_group.get(tool_key)
    if group_name and group_name in versions_by_group:
        group_versions = versions_by_group[group_name]
        if len(group_versions) == 1:
            version_html = f"<samp>{html.escape(next(iter(group_versions.values())))}</samp>"
        else:
            formatted_versions = [
                f"{html.escape(tool_name)} {html.escape(version)}"
                for tool_name, version in group_versions.items()
            ]
            version_html = f"<samp>{'; '.join(formatted_versions)}</samp>"

    tool_rows.append(
        (
            tool_order.get(group_name, 999),
            html.unescape(tool_html).lower(),
            f"<tr><td>{tool_html}</td><td>{version_html}</td><td>{params_html}</td><td>{notes_html}</td></tr>",
        )
    )

if tool_rows:
    tool_params_lines = [
        "id: 'nf-core-oceangenomesdraftgenomes-sample-tool-parameters'",
        "description: 'Exact tool parameters captured from the module command scripts for this sample.'",
        "section_name: 'Tool Parameters Used'",
        "plot_type: 'html'",
        "data: |",
        '    <table class="table table-condensed">',
        '    <thead><tr><th>Tool</th><th>Version</th><th>Effective Parameters</th><th>Notes</th></tr></thead>',
        "    <tbody>",
    ]
    tool_params_lines.extend(
        f"    {row}"
        for _order, _tool_name, row in sorted(tool_rows, key=lambda item: (item[0], item[1]))
    )
    tool_params_lines.extend(["    </tbody>", "    </table>"])
    newline = chr(10)
    Path('sample_tool_parameters_mqc.yaml').write_text(newline.join(tool_params_lines) + newline, encoding='utf-8')
PY

    multiqc \\
        --force \\
        $args \\
        $config \\
        $prefix \\
        $extra_config \\
        $logo \\
        $replace \\
        $samples \\
        .

    printf '"%s":\\n    multiqc: %s\\n' \\
        "${task.process}" \\
        "\$(multiqc --version | sed -e 's/multiqc, version //g')" > versions.yml
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_report.html

    printf '"%s":\\n    multiqc: %s\\n' \\
        "${task.process}" \\
        "\$(multiqc --version | sed -e 's/multiqc, version //g')" > versions.yml
    """
}
