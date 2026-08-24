#!/usr/bin/env python3
"""Acquire and audit public GEO SOFT metadata for new multimodal cohorts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request


SOURCES = {
    "GSE269412": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269412/soft/GSE269412_family.soft.gz",
    "GSE267119": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267119/soft/GSE267119_family.soft.gz",
    "GSE105127": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105127/soft/GSE105127_family.soft.gz",
}
MAX_DOWNLOAD_BYTES = 50 * 1024 * 1024
ZONE_TITLE = re.compile(r"^(?P<donor>[0-9]+)_(?P<zone>CV|IZ|PP)_(?P<assay>RNA|RRBS)$")


class SourceAuditError(RuntimeError):
    """Raised when source metadata violates its frozen structural contract."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def fetch(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "masld-bench-metadata-audit/1.0"},
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        header_length = response.headers.get("Content-Length")
        payload = response.read(MAX_DOWNLOAD_BYTES + 1)
        if len(payload) > MAX_DOWNLOAD_BYTES:
            raise SourceAuditError(f"source exceeds {MAX_DOWNLOAD_BYTES} bytes: {url}")
        if header_length is not None and int(header_length) != len(payload):
            raise SourceAuditError(
                f"Content-Length mismatch for {url}: {header_length} != {len(payload)}"
            )
        headers = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    try:
        gzip.decompress(payload)
    except (OSError, EOFError) as error:
        raise SourceAuditError(f"invalid gzip payload from {url}") from error
    return payload, headers


def parse_soft(payload: bytes) -> list[dict[str, object]]:
    try:
        text = gzip.decompress(payload).decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise SourceAuditError("SOFT payload is not strict UTF-8 gzip text") from error
    records: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE = "):
            current = {
                "accession": line.split(" = ", 1)[1].strip(),
                "title": "",
                "source_name": "",
                "characteristics": [],
                "data_processing": [],
                "supplementary_files": [],
            }
            records.append(current)
            continue
        if current is None or not line.startswith("!Sample_") or " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        value = value.strip()
        if key == "!Sample_title":
            current["title"] = value
        elif key == "!Sample_source_name_ch1":
            current["source_name"] = value
        elif key == "!Sample_characteristics_ch1":
            current["characteristics"].append(value)  # type: ignore[union-attr]
        elif key == "!Sample_data_processing":
            current["data_processing"].append(value)  # type: ignore[union-attr]
        elif key.startswith("!Sample_supplementary_file"):
            current["supplementary_files"].append(value)  # type: ignore[union-attr]
    if not records:
        raise SourceAuditError("SOFT payload contains no SAMPLE entities")
    if any(not record["accession"] or not record["title"] for record in records):
        raise SourceAuditError("SOFT sample is missing accession or title")
    accessions = [str(record["accession"]) for record in records]
    if len(accessions) != len(set(accessions)):
        raise SourceAuditError("SOFT payload contains duplicate sample accessions")
    return records


def candidate_participant_id(title: str) -> str:
    return title.split("_", 1)[0].strip().upper()


def audit_znf469(
    rna: list[dict[str, object]], h3k27ac: list[dict[str, object]]
) -> dict[str, object]:
    rna_ids = [candidate_participant_id(str(record["title"])) for record in rna]
    h3_ids = [candidate_participant_id(str(record["title"])) for record in h3k27ac]
    rna_counts = Counter(rna_ids)
    h3_counts = Counter(h3_ids)
    shared = sorted(set(rna_counts) & set(h3_counts))
    exact = [key for key in shared if rna_counts[key] == h3_counts[key] == 1]
    return {
        "schema_version": "masld-bench-gse267145-join-audit-v1",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "source_paper_participants": 108,
        "rna_sample_records": len(rna),
        "h3k27ac_sample_records": len(h3k27ac),
        "rna_candidate_ids": len(rna_counts),
        "h3k27ac_candidate_ids": len(h3_counts),
        "shared_candidate_ids": len(shared),
        "exact_one_to_one_candidate_ids": len(exact),
        "rna_duplicate_candidate_ids": {
            key: value for key, value in sorted(rna_counts.items()) if value > 1
        },
        "h3k27ac_duplicate_candidate_ids": {
            key: value for key, value in sorted(h3_counts.items()) if value > 1
        },
        "rna_only_candidate_ids": sorted(set(rna_counts) - set(h3_counts)),
        "h3k27ac_only_candidate_ids": sorted(set(h3_counts) - set(rna_counts)),
        "candidate_join_method": "uppercase_title_prefix_before_first_underscore",
        "candidate_join_is_authoritative": False,
        "admission_ready": False,
        "remaining_blockers": [
            "authoritative participant and repeated-library join",
            "rights and redistribution audit",
            "source-native label and missingness audit",
            "RNA and H3K27ac reference/QC crosswalks",
        ],
    }


def audit_zonated(records: list[dict[str, object]]) -> dict[str, object]:
    parsed: list[tuple[str, str, str]] = []
    unmatched: list[str] = []
    for record in records:
        title = str(record["title"])
        match = ZONE_TITLE.fullmatch(title)
        if match is None:
            unmatched.append(title)
        else:
            parsed.append((match["donor"], match["zone"], match["assay"]))
    counts = Counter(donor for donor, _, _ in parsed)
    combinations = Counter(parsed)
    expected_pairs = {
        (zone, assay) for zone in ("CV", "IZ", "PP") for assay in ("RNA", "RRBS")
    }
    complete_donors = []
    for donor in sorted(counts):
        observed = {(zone, assay) for d, zone, assay in parsed if d == donor}
        if observed == expected_pairs and all(
            combinations[(donor, zone, assay)] == 1 for zone, assay in expected_pairs
        ):
            complete_donors.append(donor)
    if len(records) != 114 or unmatched or len(counts) != 19 or len(complete_donors) != 19:
        raise SourceAuditError(
            "GSE105127 no longer satisfies the frozen 19 donor x 3 zone x 2 assay topology"
        )
    return {
        "schema_version": "masld-bench-gse105127-topology-audit-v1",
        "cohort_family_id": "gse105127_zonated_rna_rrbs",
        "sample_records": len(records),
        "candidate_donors": len(counts),
        "complete_candidate_donors": len(complete_donors),
        "zones": ["CV", "IZ", "PP"],
        "assays": ["RNA", "RRBS"],
        "records_per_donor": sorted(set(counts.values())),
        "title_topology_pass": True,
        "title_join_is_authoritative": False,
        "admission_ready": False,
        "remaining_blockers": [
            "authoritative same-zone aliquot topology",
            "rights and redistribution audit",
            "coverage-aware RRBS QC and likelihood",
            "RNA annotation and hg19-to-GRCh38p14 crosswalk",
        ],
    }


def write_samples(path: Path, series: str, records: list[dict[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "series",
                "accession",
                "title",
                "source_name",
                "candidate_participant_id",
                "characteristics_json",
                "data_processing_json",
                "supplementary_files_json",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(
                {
                    "series": series,
                    "accession": record["accession"],
                    "title": record["title"],
                    "source_name": record["source_name"],
                    "candidate_participant_id": candidate_participant_id(
                        str(record["title"])
                    ),
                    "characteristics_json": json.dumps(
                        record["characteristics"], sort_keys=True, separators=(",", ":")
                    ),
                    "data_processing_json": json.dumps(
                        record["data_processing"], sort_keys=True, separators=(",", ":")
                    ),
                    "supplementary_files_json": json.dumps(
                        record["supplementary_files"],
                        sort_keys=True,
                        separators=(",", ":"),
                    ),
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SourceAuditError(f"refusing to overwrite output: {output}")
    (output / "raw").mkdir(parents=True)
    (output / "samples").mkdir()
    receipts: dict[str, object] = {}
    samples: dict[str, list[dict[str, object]]] = {}
    for series, url in SOURCES.items():
        payload, headers = fetch(url)
        target = output / "raw" / f"{series}_family.soft.gz"
        with target.open("xb") as handle:
            handle.write(payload)
        records = parse_soft(payload)
        samples[series] = records
        write_samples(output / "samples" / f"{series}.tsv", series, records)
        receipts[series] = {
            "url": url,
            "bytes": len(payload),
            "sha256": sha256_bytes(payload),
            "headers": headers,
            "sample_records": len(records),
        }
    znf = audit_znf469(samples["GSE269412"], samples["GSE267119"])
    zonated = audit_zonated(samples["GSE105127"])
    (output / "gse267145_join_audit.json").write_text(
        json.dumps(znf, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    (output / "gse105127_topology_audit.json").write_text(
        json.dumps(zonated, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    (output / "source_receipt.json").write_text(
        json.dumps(
            {
                "schema_version": "masld-bench-geo-source-receipt-v1",
                "sources": receipts,
                "all_sources_admission_ready": False,
            },
            sort_keys=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "sample_records": {
                    key: len(value) for key, value in sorted(samples.items())
                },
                "gse267145_shared_candidate_ids": znf["shared_candidate_ids"],
                "gse105127_complete_candidate_donors": zonated[
                    "complete_candidate_donors"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
