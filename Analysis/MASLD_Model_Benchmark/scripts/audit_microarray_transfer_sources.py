#!/usr/bin/env python3
"""Acquire official GEO metadata for independent MASLD microarray transfer cohorts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import urllib.request


SOURCES = {
    "GSE49541": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49541/soft/GSE49541_family.soft.gz",
    "GSE83452": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE83nnn/GSE83452/soft/GSE83452_family.soft.gz",
}
EXPECTED_SAMPLE_COUNTS = {"GSE49541": 72, "GSE83452": 231}
MAX_DOWNLOAD_BYTES = 256 * 1024 * 1024


class MicroarraySourceAuditError(RuntimeError):
    """Raised when an official GEO source violates the frozen census."""


def fetch(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "masld-bench-microarray-source-audit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        declared = response.headers.get("Content-Length")
        payload = response.read(MAX_DOWNLOAD_BYTES + 1)
        if len(payload) > MAX_DOWNLOAD_BYTES:
            raise MicroarraySourceAuditError(f"official SOFT exceeds size cap: {url}")
        if declared is not None and int(declared) != len(payload):
            raise MicroarraySourceAuditError(f"Content-Length differs: {url}")
        headers = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    try:
        gzip.decompress(payload)
    except (OSError, EOFError) as error:
        raise MicroarraySourceAuditError(f"official SOFT is not valid gzip: {url}") from error
    return payload, headers


def parse_soft(payload: bytes) -> tuple[dict[str, list[str]], list[dict[str, object]]]:
    try:
        text = gzip.decompress(payload).decode("utf-8", errors="strict")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise MicroarraySourceAuditError("SOFT is not strict UTF-8 gzip text") from error
    series: dict[str, list[str]] = {}
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    list_fields = {
        "!Sample_description": "description",
        "!Sample_characteristics_ch1": "characteristics",
        "!Sample_data_processing": "data_processing",
        "!Sample_relation": "relations",
    }
    scalar_fields = {
        "!Sample_title": "title",
        "!Sample_source_name_ch1": "source_name",
        "!Sample_organism_ch1": "organism",
        "!Sample_platform_id": "platform_id",
        "!Sample_library_strategy": "library_strategy",
    }
    for line in text.splitlines():
        if line.startswith("!Series_") and " = " in line:
            key, value = line.split(" = ", 1)
            series.setdefault(key, []).append(value.strip())
            continue
        if line.startswith("^SAMPLE = "):
            current = {
                "accession": line.split(" = ", 1)[1].strip(),
                "title": "",
                "source_name": "",
                "organism": "",
                "platform_id": "",
                "library_strategy": "",
                "description": [],
                "characteristics": [],
                "data_processing": [],
                "relations": [],
                "supplementary_files": [],
            }
            samples.append(current)
            continue
        if current is None or not line.startswith("!Sample_") or " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        value = value.strip()
        if key in scalar_fields:
            current[scalar_fields[key]] = value
        elif key in list_fields:
            current[list_fields[key]].append(value)  # type: ignore[union-attr]
        elif key.startswith("!Sample_supplementary_file"):
            current["supplementary_files"].append(value)  # type: ignore[union-attr]
    if not samples:
        raise MicroarraySourceAuditError("SOFT has no sample entities")
    if any(
        not row["accession"]
        or not row["title"]
        or not row["organism"]
        or not row["platform_id"]
        for row in samples
    ):
        raise MicroarraySourceAuditError("sample identity, organism, or platform is absent")
    if len({str(row["accession"]) for row in samples}) != len(samples):
        raise MicroarraySourceAuditError("sample accessions are not unique")
    return series, samples


def characteristic_key_counts(samples: list[dict[str, object]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in samples:
        for value in row["characteristics"]:  # type: ignore[union-attr]
            key = str(value).split(":", 1)[0].strip() if ":" in str(value) else "__unkeyed__"
            counts[key] += 1
    return dict(sorted(counts.items()))


def write_samples(path: Path, series_id: str, rows: list[dict[str, object]]) -> None:
    fields = (
        "series",
        "accession",
        "title",
        "source_name",
        "organism",
        "platform_id",
        "library_strategy",
        "description_json",
        "characteristics_json",
        "data_processing_json",
        "relations_json",
        "supplementary_files_json",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "series": series_id,
                    "accession": row["accession"],
                    "title": row["title"],
                    "source_name": row["source_name"],
                    "organism": row["organism"],
                    "platform_id": row["platform_id"],
                    "library_strategy": row["library_strategy"],
                    "description_json": json.dumps(row["description"], sort_keys=True),
                    "characteristics_json": json.dumps(row["characteristics"], sort_keys=True),
                    "data_processing_json": json.dumps(row["data_processing"], sort_keys=True),
                    "relations_json": json.dumps(row["relations"], sort_keys=True),
                    "supplementary_files_json": json.dumps(
                        row["supplementary_files"], sort_keys=True
                    ),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=False)
    raw = arguments.output / "raw"
    sample_root = arguments.output / "samples"
    raw.mkdir()
    sample_root.mkdir()
    audits: dict[str, object] = {}
    receipts: dict[str, object] = {}
    for series_id, url in SOURCES.items():
        payload, headers = fetch(url)
        path = raw / f"{series_id}_family.soft.gz"
        path.write_bytes(payload)
        series, samples = parse_soft(payload)
        expected = EXPECTED_SAMPLE_COUNTS[series_id]
        if len(samples) != expected:
            raise MicroarraySourceAuditError(
                f"{series_id} sample count differs: {len(samples)} != {expected}"
            )
        if {str(row["organism"]) for row in samples} != {"Homo sapiens"}:
            raise MicroarraySourceAuditError(f"{series_id} includes non-human samples")
        write_samples(sample_root / f"{series_id}.tsv", series_id, samples)
        audits[series_id] = {
            "sample_records": len(samples),
            "unique_titles": len({str(row["title"]) for row in samples}),
            "platform_ids": dict(
                sorted(Counter(str(row["platform_id"]) for row in samples).items())
            ),
            "source_names": dict(
                sorted(Counter(str(row["source_name"]) for row in samples).items())
            ),
            "characteristic_key_record_counts": characteristic_key_counts(samples),
            "series_title": series.get("!Series_title", []),
            "series_summary": series.get("!Series_summary", []),
            "overall_design": series.get("!Series_overall_design", []),
            "pubmed_ids": series.get("!Series_pubmed_id", []),
            "relations": series.get("!Series_relation", []),
            "supplementary_files": series.get("!Series_supplementary_file", []),
            "participant_join_authoritative": False,
            "phenotype_mapping_authoritative": False,
            "model_training_activated": False,
        }
        receipts[series_id] = {
            "url": url,
            "bytes": len(payload),
            "sha256": sha256(payload).hexdigest(),
            "http_headers": headers,
        }
    summary = {
        "schema_version": "masld-bench-microarray-transfer-source-audit-v1",
        "status": "pass_source_metadata_only",
        "cohorts": audits,
        "total_geo_sample_records": sum(EXPECTED_SAMPLE_COUNTS.values()),
        "biological_unit_not_yet_frozen": True,
        "participant_join_authoritative": False,
        "phenotype_mapping_authoritative": False,
        "rights_audit_complete": False,
        "raw_array_acquisition_complete": False,
        "model_training_activated": False,
        "next_step": "cohort_specific_participant_repeat_histology_rights_raw_CEL_and_platform_probe_audit",
    }
    (arguments.output / "source_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (arguments.output / "source_receipt.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
