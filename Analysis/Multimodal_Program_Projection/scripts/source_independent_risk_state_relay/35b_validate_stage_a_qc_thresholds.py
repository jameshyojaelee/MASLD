#!/usr/bin/env python3
"""Independently validate the pre-QC Stage-A threshold freeze."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def main() -> None:
    seal = json.loads(
        (CANDIDATE_ROOT / "STAGE_A_QC_THRESHOLDS_FROZEN.json").read_text(encoding="utf-8")
    )
    if seal.get("status") != "stage_a_qc_thresholds_frozen_before_raw_qc":
        raise RuntimeError("Invalid Stage-A threshold-freeze status")
    for field in [
        "raw_qc_accessed_before_freeze", "scientific_condition_labels_opened",
        "scientific_outcomes_inspected",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Threshold freeze improperly sets {field}")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Threshold-freeze output hash mismatch: {name}")
    manifest = read_tsv(CANDIDATE_ROOT / "qc_threshold_input_manifest.tsv")
    if {row["role"] for row in manifest} != {
        "threshold_input", "threshold_freeze_signoff", "qc_unblinding_contract_seal"
    }:
        raise RuntimeError("Threshold-freeze input universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Threshold-freeze source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Threshold-freeze source hash drift: {path}")
    thresholds = read_tsv(CANDIDATE_ROOT / "frozen_blinded_assay_qc_thresholds.tsv")
    if len(thresholds) != int(seal["n_thresholds"]):
        raise RuntimeError("Threshold count drift")
    if len({(row["assay_id"], row["metric_name"]) for row in thresholds}) != len(thresholds):
        raise RuntimeError("Frozen thresholds are duplicated")
    if any(
        yes(row["raw_qc_accessed_before_freeze"])
        or not yes(row["review_concordant"])
        or row["reviewer_1"] == row["reviewer_2"]
        for row in thresholds
    ):
        raise RuntimeError("Frozen threshold review state drift")
    print(
        "STAGE_A_QC_THRESHOLD_VALIDATION_PASS "
        f"assays={seal['n_assays']} thresholds={len(thresholds)} outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
