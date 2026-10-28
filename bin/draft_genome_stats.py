#!/usr/bin/env python3
"""Shared parsers and PostgreSQL helpers for draft-genome statistics."""

from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Iterable, Optional, Union

PathLike = Union[str, Path]


FAMILY_COLUMNS = {
    "fastp": (
        "passed_filter_reads", "low_quality_reads", "too_many_n_reads",
        "too_short_reads", "too_long_reads", "raw_total_reads",
        "raw_total_bases", "raw_q20_bases", "raw_q30_bases",
        "raw_q20_rate", "raw_q30_rate", "raw_read1_mean_length",
        "raw_read2_mean_length", "raw_gc_content", "total_reads",
        "total_bases", "q20_bases", "q30_bases", "q20_rate", "q30_rate",
        "read1_mean_length", "read2_mean_length", "gc_content", "mach",
        "initial",
    ),
    "assembly": (
        "homozygosity", "heterozygosity", "genomesize", "repeatsize",
        "uniquesize", "modelfit", "errorrate",
        # Both values GenomeScope reports, not just the right-hand column. See
        # parse_genomescope for why modelfit_allkmers is not called modelfit_min.
        "homozygosity_min", "heterozygosity_min", "genomesize_min",
        "repeatsize_min", "uniquesize_min", "modelfit_allkmers",
        # The verdict, not just the numbers. A stored genome size nobody can audit is how
        # a 48% fit got published as "EXCELLENT" in the first place.
        "genome_size_reliable", "genome_size_flags", "kmercov", "lambda_depth",
        # How big the host genome is, and how much of the assembly is the animal, from
        # per-contig read depth rather than a k-mer peak. Both are lower bounds: a
        # symbiont at the host's own depth counts as host, a high-copy host repeat does
        # not. The only answer available below about 10x host coverage.
        "host_assembly_size", "host_assembly_fraction",
        # The same host mass with haplotype-split pairs counted once, and whether the
        # partition's BUSCO anchor could be trusted at all. A low-confidence partition is
        # one whose host size should not be quoted.
        "host_haploid_size", "partition_confidence",
        # Why a partition is low confidence, and how many BUSCO genes it and the depth rest
        # on. busco_read_depth is the per-base read depth over those genes: the sequencing
        # depth ENA's COVERAGE wants. lambda_depth is the same depth halved and moved into
        # k-mer space, for checking GenomeScope, and must not be sent as a coverage.
        "partition_note", "busco_regions", "busco_read_depth",
        # Which GenomeScope fit the numbers above are: the original, or the reseeded refit
        # that GENOMESCOPE_RESEED keeps only when it wins.
        "genomescope_fit",
    ),
    "decontamination": (
        "num_contigs_mitochondrion", "num_contigs_plastid",
        "num_contigs_prokarya", "bp_mitochondrion", "bp_plastid",
        "bp_prokarya", "num_contigs_exclude", "bp_exclude",
        "num_contigs_trim", "bp_trim", "num_contigs_review", "bp_review",
        "num_contigs",
    ),
    "busco": (
        "complete", "single_copy", "multi_copy", "fragmented", "missing",
        "n_markers", "domain", "number_of_scaffolds", "number_of_contigs",
        "total_length", "percent_gaps", "scaffold_n50", "contigs_n50",
        "internal_stop_codon_count", "internal_stop_codon_percent",
    ),
    "merqury": (
        "unique_k_mers_assembly", "k_mers_total", "qv", "error",
        "k_mer_set", "solid_k_mers", "total_k_mers", "completeness",
    ),
    "gfastats": (
        "gfa_num_contigs", "gfa_contig_n50", "gfa_num_scaffolds",
        "gfa_scaffold_n50", "gfa_largest_scaffold",
        "gfa_total_scaffold_length", "gfa_gc_content_percent",
    ),
}

FLOAT_COLUMNS = {
    "raw_q20_rate", "raw_q30_rate", "raw_gc_content", "q20_rate",
    "q30_rate", "gc_content", "homozygosity", "heterozygosity",
    "modelfit", "errorrate", "modelfit_allkmers", "homozygosity_min",
    "heterozygosity_min", "lambda_depth", "kmercov", "host_assembly_fraction",
    "busco_read_depth",
    "complete", "single_copy", "multi_copy",
    "fragmented", "missing", "percent_gaps", "scaffold_n50",
    "internal_stop_codon_percent", "qv", "error", "completeness",
    "gfa_gc_content_percent",
}

# These GenomeScope percentage columns were NUMERIC(5,2) in draft_genomes, which rounded
# on insert: 0.381807% and 0.421044% both became ~0.4, so the low-heterozygosity fishes
# were indistinguishable. They are NUMERIC(8,4) now. The set is kept because rows written
# before that migration are still rounded to two places, and the verification tolerance
# has to allow for it until they are backfilled again.
TWO_DECIMAL_COLUMNS = {
    "homozygosity", "heterozygosity", "modelfit", "errorrate",
    "homozygosity_min", "heterozygosity_min", "modelfit_allkmers",
}

BUSCO_KEYMAP = {
    "complete": ("complete_percentage", "complete", "c"),
    "single_copy": ("single_copy_percentage", "single_percentage", "single_copy", "single"),
    "multi_copy": ("multi_copy_percentage", "duplicated_percentage", "multi_copy", "duplicated"),
    "fragmented": ("fragmented_percentage", "fragmented", "f"),
    "missing": ("missing_percentage", "missing", "m"),
    "n_markers": ("n_markers", "number_of_buscos", "total_buscos", "lineage_dataset_size"),
    "domain": ("domain", "lineage", "lineage_dataset", "lineage_name"),
    "number_of_scaffolds": ("number_of_scaffolds", "scaffolds"),
    "number_of_contigs": ("number_of_contigs", "contigs"),
    "total_length": ("total_length", "assembly_size", "genome_size_bp"),
    "percent_gaps": ("percent_gaps", "gap_percent", "n_percent"),
    "scaffold_n50": ("scaffold_n50", "n50_scaffold", "scaffold_n50_bp"),
    "contigs_n50": ("contigs_n50", "n50_contig", "contig_n50", "contig_n50_bp"),
    "internal_stop_codon_count": ("internal_stop_codon_count",),
    "internal_stop_codon_percent": ("internal_stop_codon_percent",),
}


def load_kv_config(path: PathLike) -> dict[str, str]:
    config: dict[str, str] = {}
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            config[key.strip().lower()] = value.strip()
    return config


def db_params(path: PathLike) -> dict[str, Any]:
    config = load_kv_config(path)
    aliases = {
        "dbname": ("dbname", "db_name"), "user": ("user", "db_user"),
        "password": ("password", "db_password"), "host": ("host", "db_host"),
        "port": ("port", "db_port"), "sslmode": ("sslmode", "db_sslmode"),
    }
    values = {name: next((config[k] for k in keys if config.get(k)), None)
              for name, keys in aliases.items()}
    missing = [key for key in ("dbname", "user", "password", "host") if not values[key]]
    if missing:
        raise ValueError(f"missing database config keys: {', '.join(missing)}")
    params: dict[str, Any] = {
        "dbname": values["dbname"], "user": values["user"],
        "password": values["password"], "host": values["host"],
        "port": int(values["port"] or "5432"), "connect_timeout": 10,
    }
    if values["sslmode"]:
        params["sslmode"] = values["sslmode"]
    return params


def connect(path: PathLike):
    try:
        import psycopg2
    except ImportError as exc:
        raise RuntimeError("psycopg2 is required for database access") from exc
    return psycopg2.connect(**db_params(path))


def _int(value: Any, default: Optional[int] = None) -> Optional[int]:
    if value is None or str(value).strip().lower() in {"", "na", "n/a", "none", "null", "-", "n/e", "ne"}:
        return default
    try:
        return int(float(str(value).replace(",", "").strip()))
    except (TypeError, ValueError):
        return default


def _float(value: Any, default: Optional[float] = None) -> Optional[float]:
    if value is None or str(value).strip().lower() in {"", "na", "n/a", "none", "null", "-", "n/e", "ne"}:
        return default
    try:
        parsed = float(str(value).replace(",", "").replace("%", "").strip())
        return parsed if math.isfinite(parsed) else default
    except (TypeError, ValueError):
        return default


def ids_from_dot_name(value: PathLike) -> tuple[Optional[str], Optional[str]]:
    parts = os.path.basename(str(value)).split(".")
    return (parts[0], parts[2].split("_", 1)[0]) if len(parts) >= 3 else (None, None)


def _record(og_id: Any, seq_date: Any, values: dict[str, Any]) -> dict[str, Any]:
    if not og_id or not seq_date:
        raise ValueError("could not infer og_id and seq_date")
    return {"og_id": str(og_id), "seq_date": str(seq_date), **values}


def parse_fastp(path: PathLike) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    parts = Path(path).name.split(".")
    if len(parts) < 3:
        raise ValueError(f"{path}: expected OGID.<tech>.MACHINE_DATE_INITIAL.*.json")
    run_parts = parts[2].split("_", 2)
    if len(run_parts) != 3 or not all(run_parts):
        raise ValueError(f"{path}: third dot token must be MACHINE_DATE_INITIAL")
    mach, seq_date, initial = run_parts
    filtering = data.get("filtering_result") or {}
    before = (data.get("summary") or {}).get("before_filtering") or {}
    after = (data.get("summary") or {}).get("after_filtering") or {}
    required_summary = {"total_reads", "total_bases", "q20_bases", "q30_bases",
                        "q20_rate", "q30_rate", "read1_mean_length",
                        "read2_mean_length", "gc_content"}
    for label, section in (("before_filtering", before), ("after_filtering", after)):
        missing = sorted(required_summary - set(section))
        if missing:
            raise ValueError(f"{path}: fastp {label} missing fields: {', '.join(missing)}")
    if "passed_filter_reads" not in filtering:
        raise ValueError(f"{path}: fastp filtering_result missing passed_filter_reads")
    values = {
        "passed_filter_reads": _int(filtering.get("passed_filter_reads"), 0),
        "low_quality_reads": _int(filtering.get("low_quality_reads"), 0),
        "too_many_n_reads": _int(filtering.get("too_many_n_reads", filtering.get("too_many_N_reads")), 0),
        "too_short_reads": _int(filtering.get("too_short_reads"), 0),
        "too_long_reads": _int(filtering.get("too_long_reads"), 0),
        "mach": mach, "initial": initial,
    }
    for prefix, section in (("raw_", before), ("", after)):
        for key in ("total_reads", "total_bases", "q20_bases", "q30_bases", "read1_mean_length", "read2_mean_length"):
            values[prefix + key] = _int(section.get(key), 0)
        for key in ("q20_rate", "q30_rate", "gc_content"):
            values[prefix + key] = _float(section.get(key), 0.0)
    return _record(parts[0], seq_date, values)


def parse_genomescope(path: PathLike) -> dict[str, Any]:
    """Both values GenomeScope reports for each property, not just the right-hand one.

    The summary prints two columns under a "min   max" header, and until now only the
    second reached the database. For the length and heterozygosity rows the pair really is
    a lower and upper bound, so half of an interval was being dropped -- OG2949 spans
    331-1113 Mb and only 1113 was stored.

    Model Fit is different and worse. GenomeScope's own R source emits two unrelated
    statistics there:

        sprintf(format_column_2, percentage_format(model_fit_allscore[1])),
        sprintf(format_column_3, percentage_format(model_fit_fullscore[1]))

    allscore is the share of every retained k-mer the model accounts for; fullscore is the
    fit within the model's own region. They are not bounds, which is why they can invert
    (OG2647 82/25). Storing only fullscore meant a sample could be recorded as an 84% fit
    while the model explained 17% of its k-mers -- OG2644 exactly. The columns are named
    `modelfit` (fullscore, unchanged) and `modelfit_allkmers` rather than min/max, so the
    schema does not inherit GenomeScope's mislabel.
    """
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    # The name prefix first, then the input file, taking the first that names a sample.
    # GENOMESCOPE_RESEED used to write "name prefix = reseed", which names nothing; its
    # input file line still does (OG2949.ilmn.260909.meryl.hist).
    og_id = seq_date = None
    for wanted in ("name prefix", "input file"):
        for line in lines:
            if line.lower().startswith(wanted):
                match = re.search(r"=\s*(.+)$", line)
                if match:
                    og_id, seq_date = ids_from_dot_name(match.group(1).strip())
                    break
        if og_id and seq_date:
            break

    # label prefix -> (column for the first value, column for the second)
    targets = {
        "homozygous": ("homozygosity_min", "homozygosity"),
        "heterozygous": ("heterozygosity_min", "heterozygosity"),
        "genome haploid length": ("genomesize_min", "genomesize"),
        "genome repeat length": ("repeatsize_min", "repeatsize"),
        "genome unique length": ("uniquesize_min", "uniquesize"),
        "model fit": ("modelfit_allkmers", "modelfit"),
        # GenomeScope prints the same error rate twice; one column is enough.
        "read error rate": (None, "errorrate"),
    }
    values: dict[str, Any] = {
        column: None for pair in targets.values() for column in pair if column
    }
    for line in lines:
        columns = re.split(r"\s{2,}", line.strip())
        if len(columns) < 3:
            continue
        label = columns[0].lower()
        for prefix, (first, second) in targets.items():
            if label.startswith(prefix):
                if first:
                    values[first] = _float(columns[-2].replace("bp", ""))
                values[second] = _float(columns[-1].replace("bp", ""))
                break

    # A diverged model writes "Inf bp" for a bound, which _float turns into None. That
    # used to fail the whole record and the sample got no GenomeScope row at all (OG3043).
    # Store the finite side and leave the other NULL -- the reliability flags carry the
    # reason, and a partial row beats no row.
    required = ("genomesize", "genomesize_min")
    if all(values.get(key) is None for key in required):
        raise ValueError(f"{path}: no finite genome size in either column")
    for key in ("genomesize", "genomesize_min", "repeatsize", "repeatsize_min",
                "uniquesize", "uniquesize_min"):
        if values.get(key) is not None:
            values[key] = int(round(values[key]))
    return _record(og_id, seq_date, values)


def parse_coverage_summary(path: PathLike) -> dict[str, Any]:
    """The reliability verdict from a published *_coverage_summary.json.

    Stored alongside the GenomeScope numbers because the numbers on their own cannot be
    audited later: `genomesize` reads the same whether the fit was trustworthy or not.
    `lambda_depth` is the coverage measured from read depth over single-copy BUSCO genes,
    so a row carries both what GenomeScope fitted and what the alignments say.
    """
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    og_id, seq_date = ids_from_dot_name(Path(path).name)
    values = {
        "genome_size_reliable": data.get("genome_size_reliable"),
        "genome_size_flags": data.get("genome_size_flags") or None,
        # One canonical name, read in one place. `fitted_kmer_coverage` is the legacy
        # spelling: summaries published before the rename carry only that, and the
        # --skip_genome_assembly tiers read exactly such files.
        "kmercov": _float(data.get("kmercov") or data.get("fitted_kmer_coverage")),
        "lambda_depth": _float(data.get("lambda_depth")),
        "host_assembly_size": _int(data.get("host_assembly_size")),
        "host_assembly_fraction": _float(data.get("host_assembly_fraction")),
        "host_haploid_size": _int(data.get("host_haploid_size")),
        "partition_confidence": data.get("partition_confidence") or None,
        "partition_note": data.get("partition_note") or None,
        "busco_read_depth": _float(data.get("busco_mean_depth")),
        "genomescope_fit": data.get("genomescope_fit") or None,
    }
    # Summaries written before RECHECK_GENOME_SIZE carried it have no busco_regions; leave
    # the column to *_kmer_depth.json then rather than writing a NULL over it.
    if "busco_regions" in data:
        values["busco_regions"] = _int(data.get("busco_regions"))
    return _record(og_id, seq_date, values)


# The columns that come from read depth over the assembly, and nothing else. A depth-only
# upload writes these and leaves the GenomeScope numbers and the verdict alone.
DEPTH_COLUMNS = (
    "lambda_depth", "busco_read_depth", "busco_regions",
    "host_assembly_size", "host_assembly_fraction", "host_haploid_size",
    "partition_confidence", "partition_note",
)


def parse_kmer_depth(path: PathLike) -> dict[str, Any]:
    """Read depth and the host partition from a published *_kmer_depth.json.

    This is where SUMMARISE_KMER_COVERAGE writes them first; the coverage summary copies
    them. It is the only source for a depth-only re-run, which has no new coverage
    summary, and the only source of busco_regions for summaries that predate it. A sample
    whose assembly gave too few BUSCO regions has a row here with every depth null and a
    status saying so, and that is uploaded as nulls, not skipped.
    """
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    og_id, seq_date = ids_from_dot_name(Path(path).name)
    values = {
        "lambda_depth": _float(data.get("lambda_depth")),
        "busco_read_depth": _float(data.get("busco_mean_depth")),
        "busco_regions": _int(data.get("busco_regions")),
        "host_assembly_size": _int(data.get("host_assembly_size")),
        "host_assembly_fraction": _float(data.get("host_assembly_fraction")),
        "host_haploid_size": _int(data.get("host_haploid_size")),
        "partition_confidence": data.get("partition_confidence") or None,
        "partition_note": data.get("partition_note") or None,
    }
    return _record(og_id, seq_date, values)


def parse_tiara(path: PathLike) -> dict[str, Any]:
    aliases = {"mitochondria": "mitochondrion", "mito": "mitochondrion",
               "chloroplast": "plastid", "prokaryote": "prokarya",
               "bacteria": "prokarya", "bacterial": "prokarya"}
    totals = {kind: [0, 0] for kind in ("mitochondrion", "plastid", "prokarya")}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        normalized = {re.sub(r"[^a-z0-9_]", "_", (name or "").lower().strip()): name
                      for name in (reader.fieldnames or [])}
        if not {"category", "num_contigs", "bp"}.issubset(normalized):
            raise ValueError(f"{path}: Tiara TSV requires Category, num_contigs and bp")
        for row in reader:
            category = str(row[normalized["category"]]).strip().lower()
            category = aliases.get(category, category)
            if category in totals:
                count = _int(row[normalized["num_contigs"]])
                bp = _int(row[normalized["bp"]])
                if count is None or bp is None:
                    raise ValueError(f"{path}: invalid Tiara numeric value for {category}")
                totals[category][0] += count
                totals[category][1] += bp
    og_id, seq_date = ids_from_dot_name(path)
    values = {}
    for category, (count, bp) in totals.items():
        values[f"num_contigs_{category}"] = count
        values[f"bp_{category}"] = bp
    return _record(og_id, seq_date, values)


def parse_filter_report(path: PathLike) -> dict[str, Any]:
    values = {"num_contigs_exclude": 0, "bp_exclude": 0,
              "num_contigs_trim": 0, "bp_trim": 0,
              "num_contigs_review": 0, "bp_review": 0}
    matched = 0
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        match = re.match(r"^(EXCLUDE|TRIM|REVIEW)\s+(\d+)\s+(\d+)$", line.strip(), re.I)
        if match:
            matched += 1
            tag = match.group(1).lower()
            values[f"num_contigs_{tag}"] = int(match.group(2))
            values[f"bp_{tag}"] = int(match.group(3))
    if not matched:
        raise ValueError(f"{path}: no EXCLUDE/TRIM/REVIEW rows found")
    og_id, seq_date = ids_from_dot_name(path)
    return _record(og_id, seq_date, values)


def parse_lt500(path: PathLike) -> dict[str, Any]:
    count = None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        match = re.search(r"(\d+)\s*$", line.strip())
        if match:
            count = int(match.group(1))
            break
    if count is None:
        raise ValueError(f"{path}: could not parse contigs <500bp count")
    og_id, seq_date = ids_from_dot_name(path)
    return _record(og_id, seq_date, {"num_contigs": count})


def parse_decontamination(tiara: PathLike, filter_report: PathLike,
                          contigs: PathLike) -> dict[str, Any]:
    records = [parse_tiara(tiara), parse_filter_report(filter_report), parse_lt500(contigs)]
    keys = {(record["og_id"], record["seq_date"]) for record in records}
    if len(keys) != 1:
        raise ValueError(f"decontamination inputs have mixed identifiers: {sorted(keys)}")
    combined = dict(records[0])
    for record in records[1:]:
        combined.update({key: value for key, value in record.items() if key not in ("og_id", "seq_date")})
    return combined


def _norm_key(value: str) -> str:
    return re.sub(r"[^a-z0-9_]", "_", value.lower().strip())


def parse_busco(path: PathLike) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    merged = {}
    for section in ("metrics", "results"):
        if isinstance(data.get(section), dict):
            merged.update(data[section])
    normalized = {_norm_key(key): value for key, value in merged.items()}
    values = {}
    for column, candidates in BUSCO_KEYMAP.items():
        value = next((normalized[key] for key in candidates if normalized.get(key) is not None), None)
        if column == "domain":
            values[column] = str(value).strip() if value is not None and str(value).strip() else None
        elif column in FLOAT_COLUMNS:
            values[column] = _float(value)
        else:
            values[column] = _int(value)
    if values["complete"] is None or values["n_markers"] is None:
        raise ValueError(f"{path}: BUSCO JSON is missing core complete/n_markers metrics")
    source = ((data.get("parameters") or {}).get("in") or
              (data.get("parameters") or {}).get("out"))
    og_id, seq_date = ids_from_dot_name(source or path)
    return _record(og_id, seq_date, values)


def _read_merqury(path: PathLike, kind: str) -> dict[str, Any]:
    rows = list(csv.reader(Path(path).read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    rows = [row for row in rows if row]
    if not rows:
        raise ValueError(f"{path}: empty Merqury file")
    metric_columns = (["unique_k_mers_assembly", "k_mers_total", "qv", "error"]
                      if kind == "qv" else ["k_mer_set", "solid_k_mers", "total_k_mers", "completeness"])
    normalized_header = [_norm_key(value) for value in rows[0]]
    if set(metric_columns).issubset(normalized_header) and ({"sample"}.issubset(normalized_header) or
                                                             {"og_id", "seq_date"}.issubset(normalized_header)):
        names, data_rows = normalized_header, rows[1:]
    else:
        if len(rows[0]) == 5:
            names = ["sample", *metric_columns]
        elif len(rows[0]) == 6:
            names = ["og_id", "seq_date", *metric_columns]
        else:
            raise ValueError(f"{path}: expected 5 or 6 Merqury columns")
        data_rows = rows
    if len(data_rows) != 1:
        raise ValueError(f"{path}: expected exactly one Merqury data row, found {len(data_rows)}")
    row = dict(zip(names, data_rows[0]))
    if row.get("sample"):
        og_id, seq_date = ids_from_dot_name(row["sample"])
    else:
        og_id, seq_date = row.get("og_id"), row.get("seq_date")
    values = {}
    for column in metric_columns:
        if column == "k_mer_set":
            values[column] = row.get(column) or None
        elif column in FLOAT_COLUMNS:
            values[column] = _float(row.get(column))
        else:
            values[column] = _int(row.get(column))
    missing = [column for column, value in values.items() if value is None]
    if missing:
        raise ValueError(f"{path}: invalid or missing Merqury fields: {', '.join(missing)}")
    return _record(og_id, seq_date, values)


def parse_merqury(completeness: PathLike, qv: PathLike) -> dict[str, Any]:
    comp_record, qv_record = _read_merqury(completeness, "completeness"), _read_merqury(qv, "qv")
    if (comp_record["og_id"], comp_record["seq_date"]) != (qv_record["og_id"], qv_record["seq_date"]):
        raise ValueError("Merqury completeness and QV inputs have mixed identifiers")
    return {**comp_record, **{key: value for key, value in qv_record.items()
                              if key not in ("og_id", "seq_date")}}


def parse_gfastats(path: PathLike) -> dict[str, Any]:
    values = {column: None for column in FAMILY_COLUMNS["gfastats"]}
    patterns = {
        "gfa_num_scaffolds": r"^\s*#\s*scaffolds\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_total_scaffold_length": r"^\s*total\s+scaffold\s+length\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_scaffold_n50": r"^\s*scaffold\s+N50\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_largest_scaffold": r"^\s*largest\s+scaffold\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_num_contigs": r"^\s*#\s*contigs\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_contig_n50": r"^\s*contig\s+N50\s*:\s*(-?\d+(?:\.\d+)?)\s*$",
        "gfa_gc_content_percent": r"^\s*GC\s*content.*?:\s*([0-9]+(?:\.[0-9]+)?)",
    }
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        for column, pattern in patterns.items():
            match = re.match(pattern, line, re.I)
            if match:
                values[column] = (_float(match.group(1)) if column in FLOAT_COLUMNS else _int(match.group(1)))
                break
    missing = [key for key, value in values.items() if value is None]
    if missing:
        raise ValueError(f"{path}: missing gfastats fields: {', '.join(missing)}")
    og_id, seq_date = ids_from_dot_name(path)
    return _record(og_id, seq_date, values)


def assert_identity(record: dict[str, Any], og_id: str, seq_date: str, family: str) -> None:
    actual = (record.get("og_id"), record.get("seq_date"))
    if actual != (og_id, seq_date):
        raise ValueError(f"{family}: parsed {actual[0]}/{actual[1]}, expected {og_id}/{seq_date}")


# Phase R (docs/assembly_relations_plan.md §4): the draft is one row in draft_genome and each
# process writes only its own table, keyed on draft_id. FAMILY_COLUMNS above keeps the legacy
# column names, which is what parsers, verification and the draft_genomes compatibility view
# use; this map says where each family's columns are stored. Table names are unqualified on
# purpose: before the R.6 merge a test session reaches the v3 tables with
# search_path = v3, public, and after it they are in public.
#
# Columns not listed here live on draft_genome itself (mach, initial). A legacy name stored
# under a different name maps to (table, stored_name).
FAMILY_TABLES: dict[str, dict[str, Union[str, tuple[str, str]]]] = {
    "fastp": {column: "draft_fastp" for column in FAMILY_COLUMNS["fastp"]
              if column not in ("mach", "initial")},
    "assembly": {column: "draft_genomescope" for column in FAMILY_COLUMNS["assembly"]},
    "decontamination": {
        **{column: "draft_tiara" for column in (
            "num_contigs_mitochondrion", "num_contigs_plastid", "num_contigs_prokarya",
            "bp_mitochondrion", "bp_plastid", "bp_prokarya")},
        **{column: "draft_contamination" for column in (
            "num_contigs_exclude", "bp_exclude", "num_contigs_trim", "bp_trim",
            "num_contigs_review", "bp_review")},
        # The count of contigs under 500 bp (*contig_count_500bp.txt), not a total.
        "num_contigs": ("draft_contamination", "num_contigs_lt_500bp"),
    },
    "busco": {column: "draft_busco" for column in FAMILY_COLUMNS["busco"]},
    "merqury": {column: "draft_merqury" for column in FAMILY_COLUMNS["merqury"]},
    "gfastats": {column: "draft_gfastats" for column in FAMILY_COLUMNS["gfastats"]},
}

# mach -> sequencing.instrument, for the run link. iSEQ QC runs are never a draft's source.
MACH_INSTRUMENT = {"NOVA": "NovaSeq", "NEXT": "NextSeq"}


def _table_columns(family: str, record: dict[str, Any]) -> dict[str, list[tuple[str, str]]]:
    """{table: [(legacy_name, stored_name), ...]} for the family columns this record carries."""
    tables: dict[str, list[tuple[str, str]]] = {}
    for column, target in FAMILY_TABLES[family].items():
        if column not in record:
            continue
        table, stored = target if isinstance(target, tuple) else (target, column)
        tables.setdefault(table, []).append((column, stored))
    return tables


def link_illumina_runs(cursor, draft_id: int, record: dict[str, Any]) -> int:
    """Link a draft to the Illumina sequencing rows it was built from; return how many.

    By (og_id, seq_date), Illumina DNA libraries only (RNA libraries share the technology), with
    the instrument matching mach, as the Phase R backfill does (migration 055). No
    match while the workbook lags is a warning. A match that exists only on non-Illumina
    sequencing (the 15 Hi-C "drafts" migration 049 removed) is an error: a draft genome is an
    Illumina assembly.
    """
    instrument = MACH_INSTRUMENT.get(str(record.get("mach") or ""))
    cursor.execute(
        "INSERT INTO draft_genome_run (draft_id, sequencing_id) "
        "SELECT %(draft_id)s, s.sequencing_id FROM v_sequencing_browse s "
        "WHERE s.og_id = %(og_id)s AND s.seq_date = %(seq_date)s "
        "AND s.technology = 'Illumina' AND s.library_type = 'illumina' "
        "AND s.instrument = %(instrument)s "
        "ON CONFLICT DO NOTHING",
        {"draft_id": draft_id, "og_id": record["og_id"], "seq_date": record["seq_date"],
         "instrument": instrument})
    cursor.execute("SELECT count(*) FROM draft_genome_run WHERE draft_id = %s", (draft_id,))
    linked = cursor.fetchone()[0]
    if not linked:
        cursor.execute(
            "SELECT string_agg(DISTINCT coalesce(technology, '?') || ' ' || coalesce(library_type, '?'), "
            "', ') FROM v_sequencing_browse "
            "WHERE og_id = %s AND seq_date = %s", (record["og_id"], record["seq_date"]))
        others = cursor.fetchone()[0]
        if others:
            raise ValueError(f"{record['og_id']}/{record['seq_date']}: only {others} sequencing "
                             "matches this draft; a draft genome must come from an Illumina DNA "
                             "library on the run its mach names")
        print(f"WARNING: {record['og_id']}/{record['seq_date']}: no sequencing row yet; the run "
              "link will be added by a later upload once the workbook has it")
    return linked


def upsert_record(cursor, family: str, record: dict[str, Any]) -> None:
    columns = tuple(column for column in FAMILY_COLUMNS[family] if column in record)
    if not columns:
        raise ValueError(f"{family}: record contains no recognized metric columns")
    identity = {"og_id": record.get("og_id"), "seq_date": record.get("seq_date"),
                "mach": record.get("mach"), "initial": record.get("initial")}
    if family == "fastp":
        cursor.execute(
            "INSERT INTO draft_genome (og_id, seq_date, mach, initial) "
            "VALUES (%(og_id)s, %(seq_date)s, %(mach)s, %(initial)s) "
            "ON CONFLICT (og_id, seq_date) DO UPDATE SET mach = EXCLUDED.mach, "
            "initial = EXCLUDED.initial RETURNING draft_id", identity)
    else:
        cursor.execute(
            "INSERT INTO draft_genome (og_id, seq_date) VALUES (%(og_id)s, %(seq_date)s) "
            "ON CONFLICT (og_id, seq_date) DO UPDATE SET og_id = EXCLUDED.og_id "
            "RETURNING draft_id", identity)
    draft_id = cursor.fetchone()[0]

    for table, pairs in _table_columns(family, record).items():
        values = {stored: record.get(column) for column, stored in pairs}
        # A process row exists only where the process produced a value (each table has a
        # CHECK num_nonnulls(...) > 0), so an all-empty record writes nothing.
        if all(value is None for value in values.values()):
            continue
        names = ("draft_id", *values)
        assignments = ", ".join(f"{name} = EXCLUDED.{name}" for name in values)
        sql = (f"INSERT INTO {table} ({', '.join(names)}) VALUES "
               f"({', '.join('%(' + name + ')s' for name in names)}) "
               f"ON CONFLICT (draft_id) DO UPDATE SET {assignments}, loaded_at = now()")
        cursor.execute(sql, {"draft_id": draft_id, **values})

    if family == "fastp":
        link_illumina_runs(cursor, draft_id, record)


def upload_records(config_path: PathLike,
                   records: Iterable[tuple[str, dict[str, Any]]]) -> None:
    """Upload one or more records in one transaction for compatibility CLIs."""
    connection = connect(config_path)
    try:
        with connection.cursor() as cursor:
            for family, record in records:
                upsert_record(cursor, family, record)
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def fetch_database_record(cursor, og_id: str, seq_date: str) -> Optional[dict[str, Any]]:
    columns = tuple(dict.fromkeys(column for family in FAMILY_COLUMNS.values() for column in family))
    cursor.execute(f"SELECT {', '.join(columns)} FROM draft_genomes WHERE og_id = %s AND seq_date = %s",
                   (og_id, seq_date))
    row = cursor.fetchone()
    return None if row is None else dict(zip(columns, row))


def values_equal(expected: Any, actual: Any, column: str,
                 rel_tol: float = 1e-9, abs_tol: float = 1e-9) -> bool:
    if expected is None or actual is None:
        return expected is None and actual is None
    if column in FLOAT_COLUMNS:
        column_abs_tol = 0.005000001 if column in TWO_DECIMAL_COLUMNS else abs_tol
        return math.isclose(float(expected), float(actual), rel_tol=rel_tol, abs_tol=column_abs_tol)
    return expected == actual


def verify_families(expected: dict[str, dict[str, Any]], actual: Optional[dict[str, Any]]) -> list[dict[str, str]]:
    results = []
    for family, record in expected.items():
        mismatches = []
        for column in FAMILY_COLUMNS[family]:
            got = None if actual is None else actual.get(column)
            if not values_equal(record.get(column), got, column):
                mismatches.append(f"{column}: expected={record.get(column)!r}, actual={got!r}")
        results.append({"family": family, "status": "PASS" if not mismatches else "FAIL",
                        "details": "; ".join(mismatches)})
    return results


def merge_identity(records: Iterable[dict[str, Any]]) -> tuple[str, str]:
    keys = {(str(record["og_id"]), str(record["seq_date"])) for record in records}
    if len(keys) != 1:
        raise ValueError(f"records have mixed identifiers: {sorted(keys)}")
    return next(iter(keys))
