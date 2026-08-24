#!/usr/bin/env python3
"""Acquire official GEO SOFT metadata for independent bulk expansion cohorts."""

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
    "GSE268273": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268273/soft/GSE268273_family.soft.gz",
    "GSE260666": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE260nnn/GSE260666/soft/GSE260666_family.soft.gz",
    "GSE274114": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274114/soft/GSE274114_family.soft.gz",
}
EXPECTED_SAMPLE_COUNTS = {"GSE268273": 109, "GSE260666": 16, "GSE274114": 39}
MAX_DOWNLOAD_BYTES = 64 * 1024 * 1024


class BulkSourceAuditError(RuntimeError):
    """Raised when official GEO metadata violates the source census."""


def fetch(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "masld-bench-bulk-source-audit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        declared = response.headers.get("Content-Length")
        payload = response.read(MAX_DOWNLOAD_BYTES + 1)
        if len(payload) > MAX_DOWNLOAD_BYTES:
            raise BulkSourceAuditError(f"source exceeds size cap: {url}")
        if declared is not None and int(declared) != len(payload):
            raise BulkSourceAuditError(f"Content-Length differs: {url}")
        receipt = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    try:
        gzip.decompress(payload)
    except (OSError, EOFError) as error:
        raise BulkSourceAuditError(f"source is not valid gzip: {url}") from error
    return payload, receipt


def parse_soft(payload: bytes) -> tuple[dict[str, list[str]], list[dict[str, object]]]:
    try:
        text = gzip.decompress(payload).decode("utf-8", errors="strict")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise BulkSourceAuditError("SOFT is not strict UTF-8 gzip text") from error
    series: dict[str, list[str]] = {}
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
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
        if key == "!Sample_title":
            current["title"] = value
        elif key == "!Sample_source_name_ch1":
            current["source_name"] = value
        elif key == "!Sample_organism_ch1":
            current["organism"] = value
        elif key == "!Sample_description":
            current["description"].append(value)  # type: ignore[union-attr]
        elif key == "!Sample_characteristics_ch1":
            current["characteristics"].append(value)  # type: ignore[union-attr]
        elif key == "!Sample_data_processing":
            current["data_processing"].append(value)  # type: ignore[union-attr]
        elif key == "!Sample_relation":
            current["relations"].append(value)  # type: ignore[union-attr]
        elif key.startswith("!Sample_supplementary_file"):
            current["supplementary_files"].append(value)  # type: ignore[union-attr]
    if not samples:
        raise BulkSourceAuditError("SOFT has no sample entities")
    if any(not row["accession"] or not row["title"] for row in samples):
        raise BulkSourceAuditError("sample accession or title is absent")
    if len({str(row["accession"]) for row in samples}) != len(samples):
        raise BulkSourceAuditError("sample accessions are not unique")
    return series, samples


def write_samples(path: Path, series_id: str, rows: list[dict[str, object]]) -> None:
    fields = (
        "series",
        "accession",
        "title",
        "source_name",
        "organism",
        "description_json",
        "characteristics_json",
        "data_processing_json",
        "relations_json",
        "supplementary_files_json",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "series": series_id,
                    "accession": row["accession"],
                    "title": row["title"],
                    "source_name": row["source_name"],
                    "organism": row["organism"],
                    "description_json": json.dumps(row["description"], sort_keys=True),
                    "characteristics_json": json.dumps(
                        row["characteristics"], sort_keys=True
                    ),
                    "data_processing_json": json.dumps(
                        row["data_processing"], sort_keys=True
                    ),
                    "relations_json": json.dumps(row["relations"], sort_keys=True),
                    "supplementary_files_json": json.dumps(
                        row["supplementary_files"], sort_keys=True
                    ),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    raw = args.output / "raw"
    sample_dir = args.output / "samples"
    raw.mkdir()
    sample_dir.mkdir()
    audits: dict[str, object] = {}
    receipts: dict[str, object] = {}
    for series_id, url in SOURCES.items():
        payload, headers = fetch(url)
        path = raw / f"{series_id}_family.soft.gz"
        path.write_bytes(payload)
        series, samples = parse_soft(payload)
        expected = EXPECTED_SAMPLE_COUNTS[series_id]
        if len(samples) != expected:
            raise BulkSourceAuditError(
                f"{series_id} sample count differs: {len(samples)} != {expected}"
            )
        if {str(row["organism"]) for row in samples} != {"Homo sapiens"}:
            raise BulkSourceAuditError(f"{series_id} includes non-human samples")
        write_samples(sample_dir / f"{series_id}.tsv", series_id, samples)
        audits[series_id] = {
            "sample_records": len(samples),
            "unique_titles": len({str(row["title"]) for row in samples}),
            "source_names": dict(
                sorted(Counter(str(row["source_name"]) for row in samples).items())
            ),
            "series_title": series.get("!Series_title", []),
            "series_summary": series.get("!Series_summary", []),
            "overall_design": series.get("!Series_overall_design", []),
            "pubmed_ids": series.get("!Series_pubmed_id", []),
            "relations": series.get("!Series_relation", []),
            "supplementary_files": series.get("!Series_supplementary_file", []),
            "candidate_record_is_active_dataset": False,
            "participant_join_authoritative": False,
            "model_training_activated": False,
        }
        receipts[series_id] = {
            "url": url,
            "bytes": len(payload),
            "sha256": sha256(payload).hexdigest(),
            "http_headers": headers,
        }
    summary = {
        "schema_version": "masld-bench-bulk-expansion-source-audit-v1",
        "status": "pass",
        "cohorts": audits,
        "total_sample_records": sum(EXPECTED_SAMPLE_COUNTS.values()),
        "candidate_record_is_active_dataset": False,
        "model_training_activated": False,
        "next_step": "authoritative_participant_group_matrix_rights_and_overlap_audit_per_cohort",
    }
    (args.output / "source_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.output / "source_receipt.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
