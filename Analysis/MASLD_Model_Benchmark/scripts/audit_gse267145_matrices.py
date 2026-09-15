#!/usr/bin/env python3
"""Freeze and structurally audit the public GSE267145 cohort count matrices."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shlex
import urllib.request

from masld_bench.artifacts import verify_frozen_tree
from scripts.audit_geo_multimodal_sources import candidate_participant_id


MATRICES = {
    "GSE269412_RNA": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269412/suppl/GSE269412_DRX0FLO_rnaseq_count_mat.txt.gz",
    "GSE267119_H3K27ac": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267119/suppl/GSE267119_liver_biopsies_h3k27ac_cutrun_cnts.txt.gz",
}
MAX_DOWNLOAD_BYTES = 100 * 1024 * 1024


class MatrixAuditError(RuntimeError):
    """Raised when a public matrix does not meet its structural safety requirements."""


def fetch(url: str) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(
        url, headers={"User-Agent": "masld-bench-matrix-audit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=180) as response:
        declared = response.headers.get("Content-Length")
        payload = response.read(MAX_DOWNLOAD_BYTES + 1)
        if len(payload) > MAX_DOWNLOAD_BYTES:
            raise MatrixAuditError(f"matrix exceeds {MAX_DOWNLOAD_BYTES} bytes")
        if declared is not None and int(declared) != len(payload):
            raise MatrixAuditError(
                f"Content-Length mismatch: {declared} != {len(payload)}"
            )
        headers = {
            key.lower(): value
            for key in ("Content-Length", "Last-Modified", "ETag")
            if (value := response.headers.get(key)) is not None
        }
    return payload, headers


def _clean_header(value: str) -> str:
    return value.strip().strip('"').strip("'")


def _split_fields(line: str, delimiter: str) -> list[str]:
    if delimiter == "tab":
        return line.rstrip("\r\n").split("\t")
    if delimiter == "shell_whitespace":
        try:
            return shlex.split(line, comments=False, posix=True)
        except ValueError as error:
            raise MatrixAuditError("malformed quoted whitespace-delimited row") from error
    raise MatrixAuditError(f"unsupported delimiter: {delimiter}")


def audit_matrix(
    payload: bytes,
    *,
    expected_candidate_ids: set[str],
) -> dict[str, object]:
    try:
        reader = io.TextIOWrapper(
            gzip.GzipFile(fileobj=io.BytesIO(payload), mode="rb"),
            encoding="utf-8",
            errors="strict",
            newline="",
        )
        header_line = reader.readline()
    except (OSError, EOFError, UnicodeDecodeError) as error:
        raise MatrixAuditError("matrix is not strict UTF-8 gzip text") from error
    if not header_line:
        raise MatrixAuditError("matrix is empty")
    delimiter = "tab" if "\t" in header_line else "shell_whitespace"
    header = [_clean_header(value) for value in _split_fields(header_line, delimiter)]
    if len(header) < 2:
        raise MatrixAuditError("matrix header is empty or has no detectable delimiter")
    candidate_by_column = {
        index: candidate_participant_id(value)
        for index, value in enumerate(header)
        if candidate_participant_id(value) in expected_candidate_ids
    }
    matched = list(candidate_by_column.values())
    matched_counts = Counter(matched)
    duplicate_sample_ids = {
        key: value for key, value in sorted(matched_counts.items()) if value > 1
    }
    row_count = 0
    width_counts: Counter[int] = Counter()
    feature_ids: set[str] = set()
    duplicate_feature_rows = 0
    numeric_values = 0
    nonzero_values = 0
    integer_values = 0
    minimum = math.inf
    maximum = -math.inf
    row_offset: int | None = None
    for line_number, line in enumerate(reader, start=2):
        values = _split_fields(line, delimiter)
        if len(values) == 1 and not values[0]:
            continue
        row_count += 1
        width_counts[len(values)] += 1
        if row_offset is None:
            if len(values) == len(header):
                row_offset = 0
            elif (
                len(values) == len(header) + 1
                and len(candidate_by_column) == len(header)
            ):
                row_offset = 1
            else:
                row_offset = 0
        expected_width = len(header) + row_offset
        if len(values) != expected_width:
            continue
        sample_value_columns = {
            header_index + row_offset for header_index in candidate_by_column
        }
        non_sample = [
            values[index]
            for index in range(len(values))
            if index not in sample_value_columns
        ]
        feature_id = "\t".join(non_sample)
        if feature_id in feature_ids:
            duplicate_feature_rows += 1
        else:
            feature_ids.add(feature_id)
        for index in sample_value_columns:
            try:
                value = float(values[index])
            except ValueError as error:
                raise MatrixAuditError(
                    f"non-numeric sample value at line {line_number}, column {index + 1}"
                ) from error
            if not math.isfinite(value) or value < 0:
                raise MatrixAuditError(
                    f"invalid nonnegative count at line {line_number}, column {index + 1}"
                )
            numeric_values += 1
            nonzero_values += int(value != 0)
            integer_values += int(value.is_integer())
            minimum = min(minimum, value)
            maximum = max(maximum, value)
    if row_count == 0:
        raise MatrixAuditError("matrix contains no feature rows")
    if row_offset is None:
        raise MatrixAuditError("matrix row offset could not be determined")
    consistent_width = width_counts == Counter({len(header) + row_offset: row_count})
    return {
        "delimiter": delimiter,
        "header_columns": len(header),
        "implicit_row_name_column": row_offset == 1,
        "feature_rows": row_count,
        "consistent_row_width": consistent_width,
        "row_width_counts": {str(key): value for key, value in sorted(width_counts.items())},
        "matched_candidate_sample_columns": len(matched),
        "matched_candidate_ids": sorted(set(matched)),
        "missing_expected_candidate_ids": sorted(expected_candidate_ids - set(matched)),
        "unexpected_candidate_ids": sorted(set(matched) - expected_candidate_ids),
        "duplicate_sample_candidate_ids": duplicate_sample_ids,
        "duplicate_feature_rows": duplicate_feature_rows,
        "numeric_sample_values": numeric_values,
        "nonzero_sample_values": nonzero_values,
        "all_sample_values_nonnegative_integers": numeric_values == integer_values,
        "minimum_sample_value": None if numeric_values == 0 else minimum,
        "maximum_sample_value": None if numeric_values == 0 else maximum,
        "matrix_sample_join_is_authoritative": False,
    }


def read_metadata_ids(metadata_audit: Path, series: str) -> set[str]:
    path = metadata_audit / "samples" / f"{series}.tsv"
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            index = header.index("candidate_participant_id")
        except ValueError as error:
            raise MatrixAuditError(f"missing candidate participant column: {path}") from error
        values = {
            line.rstrip("\n").split("\t")[index]
            for line in handle
            if line.rstrip("\n")
        }
    if not values:
        raise MatrixAuditError(f"no candidate participant IDs: {path}")
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-audit", type=Path, required=True)
    parser.add_argument("--metadata-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = args.metadata_audit.resolve(strict=True)
    if len(args.metadata_artifacts_sha256) != 64:
        raise MatrixAuditError("metadata ARTIFACTS SHA-256 is malformed")
    observed_metadata = hashlib.sha256((metadata / "ARTIFACTS.json").read_bytes()).hexdigest()
    if observed_metadata != args.metadata_artifacts_sha256:
        raise MatrixAuditError("metadata ARTIFACTS SHA-256 mismatch")
    verify_frozen_tree(metadata)
    expected = {
        "GSE269412_RNA": read_metadata_ids(metadata, "GSE269412"),
        "GSE267119_H3K27ac": read_metadata_ids(metadata, "GSE267119"),
    }
    output = args.output.resolve()
    if output.exists():
        raise MatrixAuditError(f"refusing to overwrite output: {output}")
    (output / "raw").mkdir(parents=True)
    reports: dict[str, object] = {}
    receipts: dict[str, object] = {}
    for matrix_id, url in MATRICES.items():
        payload, headers = fetch(url)
        path = output / "raw" / f"{matrix_id}.txt.gz"
        with path.open("xb") as handle:
            handle.write(payload)
        report = audit_matrix(payload, expected_candidate_ids=expected[matrix_id])
        reports[matrix_id] = report
        receipts[matrix_id] = {
            "url": url,
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "headers": headers,
        }
    rna_ids = set(reports["GSE269412_RNA"]["matched_candidate_ids"])  # type: ignore[index]
    h3_ids = set(reports["GSE267119_H3K27ac"]["matched_candidate_ids"])  # type: ignore[index]
    audit = {
        "schema_version": "masld-bench-gse267145-matrix-audit-v1",
        "cohort_family_id": "gse267145_znf469_human_liver",
        "metadata_audit_artifacts_sha256": args.metadata_artifacts_sha256,
        "matrices": reports,
        "paired_matrix_candidate_ids": len(rna_ids & h3_ids),
        "h3k27ac_matrix_ids_missing_from_rna": sorted(h3_ids - rna_ids),
        "matrix_join_is_authoritative": False,
        "admission_ready": False,
        "remaining_blockers": [
            "authoritative participant and source-paper cohort join",
            "rights and redistribution audit",
            "source-native labels, QC, and reference crosswalks",
            "train-only split and normalization contract",
        ],
    }
    (output / "matrix_audit.json").write_text(
        json.dumps(audit, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    (output / "source_receipt.json").write_text(
        json.dumps(
            {
                "schema_version": "masld-bench-gse267145-matrix-receipt-v1",
                "sources": receipts,
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
                "paired_matrix_candidate_ids": audit["paired_matrix_candidate_ids"],
                "rna_features": reports["GSE269412_RNA"]["feature_rows"],  # type: ignore[index]
                "h3k27ac_features": reports["GSE267119_H3K27ac"]["feature_rows"],  # type: ignore[index]
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
