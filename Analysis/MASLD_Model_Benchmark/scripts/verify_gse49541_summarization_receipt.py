#!/usr/bin/env python3
"""Adjudicate the GSE49541 single-array summarization, fail-closed.

The R stage emits matrices and a receipt.  This module decides whether that
output actually satisfies the frozen requirements before anything is published:
every array summarized on its own, both axes intact, no hard QC failure, no
sample removed by a descriptive flag, and no label, series matrix, all-sample
RMA, or across-array normalization anywhere in the path.

The matrices themselves are re-read, so a receipt cannot claim a shape the files
do not have.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping

EXPECTED_ARRAYS = 72
EXPECTED_PLATFORM_FEATURES = 54_675
EXPECTED_CDF = "HG-U133_Plus_2"
REQUIRED_FALSE_FLAGS = (
    "all_sample_RMA_run",
    "across_array_quantile_normalization_run",
    "GEO_series_matrix_read",
    "labels_read",
    "model_training_activated",
    "sealed_outcomes_read",
    "automatic_exclusion_from_descriptive_flag",
    "cross_array_QC_visible_to_training",
)
EXPECTED_PACKAGE_VERSIONS = {
    "affy": "1.88.0",
    "affyio": "1.80.0",
    "frma": "1.62.0",
    "hgu133plus2cdf": "2.18.0",
    "hgu133plus2frmavecs": "1.5.0",
}


class SummarizationError(RuntimeError):
    """Raised when the summarization does not satisfy the frozen requirements."""


def _scalar(value: Any) -> Any:
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def read_matrix_axes(path: Path, id_column: str) -> tuple[list[str], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[0] != id_column:
            raise SummarizationError(f"{path.name} first column is not {id_column}")
        columns = header[1:]
        rows = []
        for record in reader:
            if len(record) != len(header):
                raise SummarizationError(f"{path.name} is ragged at row {len(rows) + 1}")
            rows.append(record[0])
    if len(set(rows)) != len(rows):
        raise SummarizationError(f"{path.name} repeats a feature identifier")
    if len(set(columns)) != len(columns):
        raise SummarizationError(f"{path.name} repeats a sample identifier")
    return rows, columns


def verify_summarization(
    receipt: Mapping[str, Any], *, platform_matrix: Path, gene_matrix: Path
) -> dict[str, Any]:
    if receipt.get("schema_version") != "masld-bench-gse49541-gpl570-single-array-summary-v1":
        raise SummarizationError("summarization receipt schema differs")
    if receipt.get("series") != "GSE49541" or receipt.get("platform_id") != "GPL570":
        raise SummarizationError("summarization receipt names a different source")
    if receipt.get("method") != "frma_with_exact_hgu133plus2frmavecs":
        raise SummarizationError("summarization used a method the contract does not admit")
    if _scalar(receipt.get("arrays_read_per_frma_call")) != 1:
        raise SummarizationError("a summarization call read more than one array")
    if _scalar(receipt.get("label_blind")) is not True:
        raise SummarizationError("summarization is not marked label blind")
    if _scalar(receipt.get("samples_excluded")) != 0:
        raise SummarizationError("a sample was excluded from a label-blind stage")

    for flag in REQUIRED_FALSE_FLAGS:
        if _scalar(receipt.get(flag)) is not False:
            raise SummarizationError(f"summarization receipt does not clear {flag}")

    versions = {
        key: _scalar(value) for key, value in (receipt.get("package_versions") or {}).items()
    }
    if versions != EXPECTED_PACKAGE_VERSIONS:
        raise SummarizationError(f"summarization package versions differ: {versions!r}")

    failures = {
        key: int(_scalar(value)) for key, value in (receipt.get("hard_failures") or {}).items()
    }
    if not failures:
        raise SummarizationError("summarization receipt carries no hard-failure tally")
    offending = sorted(key for key, value in failures.items() if value != 0)
    if offending:
        raise SummarizationError(f"summarization hit a hard QC failure: {offending}")

    platform_rows, platform_columns = read_matrix_axes(platform_matrix, "platform_feature_id")
    gene_rows, gene_columns = read_matrix_axes(gene_matrix, "ensembl_gene_id")
    if len(platform_rows) != EXPECTED_PLATFORM_FEATURES:
        raise SummarizationError("platform feature matrix axis differs from GPL570")
    if len(platform_columns) != EXPECTED_ARRAYS or len(gene_columns) != EXPECTED_ARRAYS:
        raise SummarizationError(f"a matrix does not carry {EXPECTED_ARRAYS} arrays")
    if platform_columns != gene_columns:
        raise SummarizationError("platform and gene matrices disagree on the sample axis")
    if len(gene_rows) < 1:
        raise SummarizationError("gene matrix retained no gene")
    if _scalar(receipt.get("platform_features")) != len(platform_rows):
        raise SummarizationError("receipt platform axis differs from the written matrix")
    if _scalar(receipt.get("gene_features")) != len(gene_rows):
        raise SummarizationError("receipt gene axis differs from the written matrix")
    if _scalar(receipt.get("arrays_summarized")) != len(platform_columns):
        raise SummarizationError("receipt array count differs from the written matrix")

    per_array = receipt.get("per_array") or []
    if len(per_array) != EXPECTED_ARRAYS:
        raise SummarizationError(f"per-array QC does not cover {EXPECTED_ARRAYS} arrays")
    accessions = [_scalar(record.get("sample_accession")) for record in per_array]
    if sorted(accessions) != sorted(platform_columns):
        raise SummarizationError("per-array QC does not match the matrix sample axis")
    for record in per_array:
        accession = _scalar(record.get("sample_accession"))
        if _scalar(record.get("arrays_read_in_this_call")) != 1:
            raise SummarizationError(f"{accession} was summarized alongside another array")
        if _scalar(record.get("cdf_name")) != EXPECTED_CDF:
            raise SummarizationError(f"{accession} is not a {EXPECTED_CDF} array")
        if _scalar(record.get("array_type_matches_GPL570")) is not True:
            raise SummarizationError(f"{accession} array type does not match GPL570")
        if _scalar(record.get("summary_finite")) != EXPECTED_PLATFORM_FEATURES:
            raise SummarizationError(f"{accession} summary carries a nonfinite value")
        if _scalar(record.get("excluded")) is not False:
            raise SummarizationError(f"{accession} was excluded")
        if _scalar(record.get("descriptive_flag_only")) is not True:
            raise SummarizationError(f"{accession} QC is not marked descriptive")

    return {
        "schema_version": "masld-bench-gse49541-summarization-verdict-v1",
        "status": "pass_label_blind_single_array_summarization",
        "series": "GSE49541",
        "platform_id": "GPL570",
        "cohort_family_id": "gse31803_gse49541_fibrosis_array",
        "method": "frma_with_exact_hgu133plus2frmavecs",
        "arrays_summarized": len(platform_columns),
        "arrays_read_per_frma_call": 1,
        "platform_features": len(platform_rows),
        "gene_features": len(gene_rows),
        "eligible_platform_features": _scalar(receipt.get("eligible_platform_features")),
        "preprocessing_object_sha256": _scalar(receipt.get("preprocessing_object_sha256")),
        "platform_feature_matrix_sha256": _scalar(receipt.get("platform_feature_matrix_sha256")),
        "gene_matrix_sha256": _scalar(receipt.get("gene_matrix_sha256")),
        "package_versions": versions,
        "hard_failures": failures,
        "samples_excluded": 0,
        "automatic_exclusion_from_descriptive_flag": False,
        "cross_array_QC_visible_to_training": False,
        "label_blind": True,
        "all_sample_RMA_run": False,
        "across_array_quantile_normalization_run": False,
        "GEO_series_matrix_read": False,
        "labels_read": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--platform-matrix", type=Path, required=True)
    parser.add_argument("--gene-matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    verdict = verify_summarization(
        json.loads(arguments.receipt.read_text(encoding="utf-8")),
        platform_matrix=arguments.platform_matrix,
        gene_matrix=arguments.gene_matrix,
    )
    arguments.output.write_text(
        json.dumps(verdict, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
