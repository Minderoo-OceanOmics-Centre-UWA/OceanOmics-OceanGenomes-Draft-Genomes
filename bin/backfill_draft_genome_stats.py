#!/usr/bin/env python3
"""Validate and backfill draft-genome statistics from a mounted archive."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional

from draft_genome_stats import (
    FAMILY_COLUMNS,
    _read_merqury,
    assert_identity,
    connect,
    db_params,
    fetch_database_record,
    parse_busco,
    parse_decontamination,
    parse_fastp,
    parse_genomescope,
    parse_gfastats,
    upsert_record,
    verify_families,
)


COMPONENTS: dict[str, tuple[str, Callable[[Path], dict[str, Any]]]] = {
    "fastp": ("**/*.fastp.json", parse_fastp),
    "genomescope": ("kmers/*genomescope/*summary.txt", parse_genomescope),
    "tiara": ("assemblies/genome/tiara/*tiara_filter_summary.txt", lambda path: _parse_decontam_part(path, "tiara")),
    "filter_report": ("assemblies/genome/NCBI/*filter_report.txt", lambda path: _parse_decontam_part(path, "filter")),
    "contigs_lt500": ("assemblies/genome/NCBI/*contig_count_500bp.txt", lambda path: _parse_decontam_part(path, "contigs")),
    "busco": ("assemblies/genome/busco/*short_summary.json", parse_busco),
    "merqury_completeness": ("kmers/*completeness.stats", lambda path: _read_merqury(path, "completeness")),
    "merqury_qv": ("kmers/*.merqury.qv", lambda path: _read_merqury(path, "qv")),
    "gfastats": ("assemblies/genome/gfastats/*assembly_summary*", parse_gfastats),
}

INVENTORY_FIELDS = ("og_id", "seq_date", *COMPONENTS.keys(), "status", "details")
UPLOAD_FIELDS = ("og_id", "seq_date", "status", "details")
VERIFY_FIELDS = ("og_id", "seq_date", "family", "status", "details")


def _parse_decontam_part(path: Path, kind: str) -> dict[str, Any]:
    from draft_genome_stats import parse_filter_report, parse_lt500, parse_tiara
    return {"tiara": parse_tiara, "filter": parse_filter_report, "contigs": parse_lt500}[kind](path)


@dataclass
class Sample:
    og_id: str
    seq_date: str
    paths: dict[str, Path] = field(default_factory=dict)
    records: dict[str, dict[str, Any]] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)

    @property
    def valid(self) -> bool:
        return not self.errors


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate or atomically backfill all draft-genome statistics from mounted Acacia."
    )
    parser.add_argument("--archive-root", required=True, type=Path,
                        help="Mounted genomes.v2 directory containing one directory per OG ID")
    parser.add_argument("--fastp-root", type=Path,
                        help="Optional mounted S3 draft-genomes root containing per-OG fastp directories; defaults to archive root")
    parser.add_argument("--manifest", required=True, type=Path,
                        help="Tab-separated file with og_id and seq_date columns")
    parser.add_argument("--db-config", required=True, type=Path,
                        help="PostgreSQL key=value connection config")
    parser.add_argument("--report-dir", required=True, type=Path,
                        help="Directory for inventory, expected, upload and verification reports")
    parser.add_argument("--apply", action="store_true",
                        help="Write validated records to PostgreSQL; omitted means validation only")
    return parser.parse_args(argv)


def read_manifest(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        raise ValueError(f"manifest is not a readable file: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"og_id", "seq_date"}.issubset(reader.fieldnames):
            raise ValueError("manifest must be TSV with og_id and seq_date headers")
        targets = []
        seen = set()
        for line_number, row in enumerate(reader, start=2):
            og_id = (row.get("og_id") or "").strip()
            seq_date = (row.get("seq_date") or "").strip()
            if not og_id or not seq_date:
                raise ValueError(f"manifest line {line_number}: og_id and seq_date are required")
            if not __import__("re").fullmatch(r"OG\d+[A-Za-z]?(?:-\d+)?(?:_HIC[A-Za-z]*\d*)?", og_id):
                raise ValueError(f"manifest line {line_number}: invalid OG ID {og_id!r}")
            if not __import__("re").fullmatch(r"\d{6}", seq_date):
                raise ValueError(f"manifest line {line_number}: seq_date must be six digits")
            key = (og_id, seq_date)
            if key in seen:
                raise ValueError(f"manifest line {line_number}: duplicate target {og_id}/{seq_date}")
            seen.add(key)
            targets.append(key)
    if not targets:
        raise ValueError("manifest contains no samples")
    return targets


def _select_component(sample: Sample, root: Path, component: str) -> None:
    pattern, parser = COMPONENTS[component]
    candidates = sorted(root.glob(pattern))
    matches: list[tuple[Path, dict[str, Any]]] = []
    parse_errors: list[str] = []
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            record = parser(candidate)
            if (record.get("og_id"), record.get("seq_date")) == (sample.og_id, sample.seq_date):
                matches.append((candidate, record))
        except Exception as exc:
            if sample.og_id in candidate.name and sample.seq_date in candidate.name:
                parse_errors.append(f"{candidate}: {exc}")
    if len(matches) == 1 and not parse_errors:
        sample.paths[component], sample.records[component] = matches[0]
    elif len(matches) > 1:
        sample.errors.append(f"{component}: ambiguous matches: " + ", ".join(str(item[0]) for item in matches))
    elif parse_errors:
        sample.errors.append(f"{component}: malformed candidate(s): " + " | ".join(parse_errors))
    else:
        sample.errors.append(f"{component}: no file matching {pattern!r} for {sample.og_id}/{sample.seq_date}")


def discover_sample(archive_root: Path, og_id: str, seq_date: str,
                    fastp_root: Optional[Path] = None) -> Sample:
    sample = Sample(og_id, seq_date)
    sample_root = archive_root / og_id
    if not sample_root.is_dir():
        sample.errors.append(f"sample directory is missing or unreadable: {sample_root}")
        return sample
    for component in COMPONENTS:
        component_root = (fastp_root / og_id) if component == "fastp" and fastp_root else sample_root
        _select_component(sample, component_root, component)
    if not sample.valid:
        return sample
    try:
        decontamination = parse_decontamination(
            sample.paths["tiara"], sample.paths["filter_report"], sample.paths["contigs_lt500"]
        )
        merqury = {**sample.records["merqury_completeness"],
                   **{key: value for key, value in sample.records["merqury_qv"].items()
                      if key not in ("og_id", "seq_date")}}
        families = {
            "fastp": sample.records["fastp"],
            "assembly": sample.records["genomescope"],
            "decontamination": decontamination,
            "busco": sample.records["busco"],
            "merqury": merqury,
            "gfastats": sample.records["gfastats"],
        }
        for family, record in families.items():
            assert_identity(record, og_id, seq_date, family)
            missing = [column for column in FAMILY_COLUMNS[family] if column not in record]
            if missing:
                raise ValueError(f"{family}: parser omitted columns: {', '.join(missing)}")
        sample.records = families
    except Exception as exc:
        sample.errors.append(str(exc))
    return sample


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_validation_reports(report_dir: Path, samples: list[Sample]) -> None:
    inventory = []
    for sample in samples:
        row: dict[str, Any] = {"og_id": sample.og_id, "seq_date": sample.seq_date}
        row.update({name: str(sample.paths.get(name, "")) for name in COMPONENTS})
        row.update({"status": "PASS" if sample.valid else "FAIL", "details": " | ".join(sample.errors)})
        inventory.append(row)
    _write_tsv(report_dir / "inventory.tsv", INVENTORY_FIELDS, inventory)
    with (report_dir / "expected.jsonl").open("w", encoding="utf-8") as handle:
        for sample in samples:
            if sample.valid:
                handle.write(json.dumps({"og_id": sample.og_id, "seq_date": sample.seq_date,
                                         "families": sample.records}, sort_keys=True) + "\n")


def runtime_smoke_check() -> None:
    missing = []
    for module in ("psycopg2", "pandas", "numpy"):
        try:
            importlib.import_module(module)
        except ImportError:
            missing.append(module)
    if missing:
        raise RuntimeError("uploader environment is missing required modules: " + ", ".join(missing))


def apply_samples(samples: list[Sample], config: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    runtime_smoke_check()
    connection = connect(config)
    uploads: list[dict[str, str]] = []
    verifications: list[dict[str, str]] = []
    try:
        with connection.cursor() as cursor:
            cursor.execute("SELECT 1")
            if cursor.fetchone() != (1,):
                raise RuntimeError("database connectivity smoke check returned an unexpected result")
        connection.rollback()
        for sample in samples:
            sample_verification: list[dict[str, str]] = []
            try:
                with connection.cursor() as cursor:
                    for family in FAMILY_COLUMNS:
                        upsert_record(cursor, family, sample.records[family])
                    actual = fetch_database_record(cursor, sample.og_id, sample.seq_date)
                    sample_verification = verify_families(sample.records, actual)
                    failures = [row for row in sample_verification if row["status"] != "PASS"]
                    if failures:
                        raise RuntimeError("source-to-database verification failed before commit")
                connection.commit()
                uploads.append({"og_id": sample.og_id, "seq_date": sample.seq_date,
                                "status": "COMMITTED", "details": "all six families verified"})
            except Exception as exc:
                connection.rollback()
                uploads.append({"og_id": sample.og_id, "seq_date": sample.seq_date,
                                "status": "ROLLED_BACK", "details": str(exc)})
                if not sample_verification:
                    sample_verification = [{"family": family, "status": "NOT_VERIFIED",
                                            "details": "transaction failed before verification"}
                                           for family in FAMILY_COLUMNS]
            for row in sample_verification:
                verifications.append({"og_id": sample.og_id, "seq_date": sample.seq_date, **row})
    finally:
        connection.close()
    return uploads, verifications


def write_summary(path: Path, samples: list[Sample], apply: bool,
                  uploads: list[dict[str, str]], verifications: list[dict[str, str]]) -> None:
    invalid = [f"{sample.og_id}/{sample.seq_date}" for sample in samples if not sample.valid]
    failed_uploads = [f"{row['og_id']}/{row['seq_date']}" for row in uploads if row["status"] != "COMMITTED"]
    failed_verify = [f"{row['og_id']}/{row['seq_date']}:{row['family']}"
                     for row in verifications if row["status"] != "PASS"]
    lines = [
        f"mode: {'apply' if apply else 'validation-only'}",
        f"manifest_samples: {len(samples)}",
        f"validation_passed: {sum(sample.valid for sample in samples)}",
        f"validation_failed: {len(invalid)}",
        f"transactions_committed: {sum(row['status'] == 'COMMITTED' for row in uploads)}",
        f"transactions_failed: {len(failed_uploads)}",
        f"verification_failures: {len(failed_verify)}",
        "safe_to_rerun: " + (", ".join(sorted(set(invalid + failed_uploads))) or "none"),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    try:
        archive_root = args.archive_root.resolve()
        if not archive_root.is_dir() or not os.access(archive_root, os.R_OK | os.X_OK):
            raise ValueError(f"archive root is missing or unreadable on this node: {archive_root}")
        fastp_root = args.fastp_root.resolve() if args.fastp_root else archive_root
        if not fastp_root.is_dir() or not os.access(fastp_root, os.R_OK | os.X_OK):
            raise ValueError(f"fastp root is missing or unreadable on this node: {fastp_root}")
        if not args.db_config.is_file():
            raise ValueError(f"database config is missing: {args.db_config}")
        db_params(args.db_config)  # Validate syntax without connecting in validation-only mode.
        targets = read_manifest(args.manifest)
        args.report_dir.mkdir(parents=True, exist_ok=True)
        samples = [discover_sample(archive_root, og_id, seq_date, fastp_root)
                   for og_id, seq_date in targets]
        write_validation_reports(args.report_dir, samples)
        uploads: list[dict[str, str]] = []
        verifications: list[dict[str, str]] = []
        if all(sample.valid for sample in samples) and args.apply:
            try:
                uploads, verifications = apply_samples(samples, args.db_config)
            except Exception as exc:
                details = f"database/environment preflight failed: {exc}"
                uploads = [{"og_id": sample.og_id, "seq_date": sample.seq_date,
                            "status": "ROLLED_BACK", "details": details}
                           for sample in samples]
                verifications = [
                    {"og_id": sample.og_id, "seq_date": sample.seq_date,
                     "family": family, "status": "NOT_VERIFIED", "details": details}
                    for sample in samples for family in FAMILY_COLUMNS
                ]
        _write_tsv(args.report_dir / "upload.tsv", UPLOAD_FIELDS, uploads)
        _write_tsv(args.report_dir / "verification.tsv", VERIFY_FIELDS, verifications)
        write_summary(args.report_dir / "summary.txt", samples, args.apply, uploads, verifications)
        if not all(sample.valid for sample in samples):
            print(f"Validation failed; see {args.report_dir / 'inventory.tsv'}", file=sys.stderr)
            return 1
        if args.apply and (any(row["status"] != "COMMITTED" for row in uploads) or
                           any(row["status"] != "PASS" for row in verifications)):
            print(f"Backfill failed; see {args.report_dir}", file=sys.stderr)
            return 1
        print(f"{'Backfill complete' if args.apply else 'Validation complete'}: {args.report_dir}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
