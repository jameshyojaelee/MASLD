#!/usr/bin/env python3
"""Freeze assay-specific Stage-A QC thresholds before raw QC is inspected."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, atomic_write_json, read_tsv, sha256_file,
    write_tsv,
)


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def exact_rows(input_root: Path, contract_root: Path, basename: str) -> list[dict[str, str]]:
    path = input_root / basename
    template = contract_root / basename.replace(".tsv", "_template.tsv")
    if not path.is_file() or not template.is_file() or header(path) != header(template):
        raise RuntimeError(f"QC-threshold schema/source failure: {basename}")
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"QC-threshold input is empty: {basename}")
    return rows


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage-A QC thresholds: {CANDIDATE_ROOT}")
    contract_root = candidate_source("PLAN45_STAGE_A_QC_CONTRACT_ROOT")
    input_root = candidate_source("PLAN45_STAGE_A_QC_THRESHOLD_INPUT_ROOT")
    contract_seal_path = contract_root / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json"
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    if contract_seal.get("status") != "sealed_outcome_blind_stage_a_qc_unblinding_contract":
        raise RuntimeError("Invalid Stage-A QC contract dependency")
    for name, expected in contract_seal["output_sha256"].items():
        path = contract_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"QC-contract source drift: {path}")

    thresholds_path = input_root / "blinded_assay_qc_thresholds.tsv"
    thresholds = exact_rows(input_root, contract_root, thresholds_path.name)
    signoffs = exact_rows(input_root, contract_root, "threshold_freeze_signoff.tsv")
    if len(signoffs) != 1:
        raise RuntimeError("Threshold freeze requires exactly one signoff row")
    signoff = signoffs[0]
    if (
        Path(signoff["thresholds_path"]).as_posix()
        != thresholds_path.relative_to(PROJECT_ROOT).as_posix()
        or signoff["thresholds_sha256"] != sha256_file(thresholds_path)
    ):
        raise RuntimeError("Threshold signoff is not bound to the exact threshold table")
    if (
        not signoff["reviewer_1"] or not signoff["reviewer_2"]
        or signoff["reviewer_1"] == signoff["reviewer_2"]
        or not yes(signoff["review_concordant"])
        or yes(signoff["raw_qc_accessed_before_freeze"])
        or yes(signoff["scientific_condition_inspected"])
        or not signoff["signed_utc"]
    ):
        raise RuntimeError("Threshold freeze signoff violates the blinded pre-QC firewall")

    seen: set[tuple[str, str]] = set()
    for row in thresholds:
        key = (row["assay_id"], row["metric_name"])
        if key in seen or not all(key):
            raise RuntimeError("QC threshold identifiers are empty or duplicated")
        seen.add(key)
        if row["operator"] not in {"ge", "le", "eq"}:
            raise RuntimeError("QC threshold operator is invalid")
        value = float(row["threshold"])
        if not math.isfinite(value) or not row["metric_unit"] or not row["threshold_source"]:
            raise RuntimeError("QC threshold value/unit/source is invalid")
        if (
            yes(row["raw_qc_accessed_before_freeze"])
            or not row["frozen_utc"]
            or not yes(row["review_concordant"])
            or row["reviewer_1"] == row["reviewer_2"]
        ):
            raise RuntimeError("A QC threshold lacks valid independent pre-QC review")

    CANDIDATE_ROOT.mkdir(parents=True)
    frozen_path = CANDIDATE_ROOT / "frozen_blinded_assay_qc_thresholds.tsv"
    write_tsv(frozen_path, thresholds, header(thresholds_path))
    signoff_path = CANDIDATE_ROOT / "frozen_threshold_freeze_signoff.tsv"
    write_tsv(signoff_path, signoffs, header(input_root / "threshold_freeze_signoff.tsv"))
    manifest_rows = []
    for role, path in [
        ("threshold_input", thresholds_path),
        ("threshold_freeze_signoff", input_root / "threshold_freeze_signoff.tsv"),
        ("qc_unblinding_contract_seal", contract_seal_path),
    ]:
        manifest_rows.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest_path = CANDIDATE_ROOT / "qc_threshold_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, ["role", "source_path", "size_bytes", "sha256"])
    outputs = [frozen_path, signoff_path, manifest_path]
    payload = {
        "status": "stage_a_qc_thresholds_frozen_before_raw_qc",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_assays": len({row["assay_id"] for row in thresholds}),
        "n_thresholds": len(thresholds),
        "raw_qc_accessed_before_freeze": False,
        "scientific_condition_labels_opened": False,
        "scientific_outcomes_inspected": False,
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_QC_THRESHOLDS_FROZEN.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
