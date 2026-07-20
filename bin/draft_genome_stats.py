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
    "modelfit", "errorrate", "complete", "single_copy", "multi_copy",
    "fragmented", "missing", "percent_gaps", "scaffold_n50",
    "internal_stop_codon_percent", "qv", "error", "completeness",
    "gfa_gc_content_percent",
}

# These GenomeScope percentage columns are NUMERIC(..., 2) in draft_genomes.
# PostgreSQL therefore rounds source values to two decimal places on insert.
TWO_DECIMAL_COLUMNS = {
    "homozygosity", "heterozygosity", "modelfit", "errorrate",
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
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    source = None
    for wanted in ("name prefix", "input file"):
        for line in lines:
            if line.lower().startswith(wanted):
                match = re.search(r"=\s*(.+)$", line)
                if match:
                    source = match.group(1).strip()
                    break
        if source:
            break
    og_id, seq_date = ids_from_dot_name(source or "")
    targets = {
        "homozygous": ("homozygosity", _float),
        "heterozygous": ("heterozygosity", _float),
        "genome haploid length": ("genomesize", _float),
        "genome repeat length": ("repeatsize", _float),
        "genome unique length": ("uniquesize", _float),
        "model fit": ("modelfit", _float), "read error rate": ("errorrate", _float),
    }
    values = {column: None for column, _ in targets.values()}
    for line in lines:
        columns = re.split(r"\s{2,}", line.strip())
        if len(columns) < 3:
            continue
        label = columns[0].lower()
        for prefix, (column, parser) in targets.items():
            if label.startswith(prefix):
                values[column] = parser(columns[-1].replace("bp", ""))
                break
    missing = [key for key, value in values.items() if value is None]
    if missing:
        raise ValueError(f"{path}: missing GenomeScope fields: {', '.join(missing)}")
    for key in ("genomesize", "repeatsize", "uniquesize"):
        values[key] = int(round(values[key]))
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


def upsert_record(cursor, family: str, record: dict[str, Any]) -> None:
    columns = tuple(column for column in FAMILY_COLUMNS[family] if column in record)
    if not columns:
        raise ValueError(f"{family}: record contains no recognized metric columns")
    names = ("og_id", "seq_date", *columns)
    assignments = ", ".join(f"{name} = EXCLUDED.{name}" for name in columns)
    sql = (f"INSERT INTO draft_genomes ({', '.join(names)}) VALUES "
           f"({', '.join('%(' + name + ')s' for name in names)}) "
           f"ON CONFLICT (og_id, seq_date) DO UPDATE SET {assignments}")
    cursor.execute(sql, {name: record.get(name) for name in names})


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
