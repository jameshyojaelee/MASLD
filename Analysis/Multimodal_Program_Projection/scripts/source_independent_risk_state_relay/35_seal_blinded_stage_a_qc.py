#!/usr/bin/env python3
"""Seal all Stage-A technical QC decisions before condition unblinding."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, atomic_write_json, read_tsv, sha256_file,
    write_tsv,
)


INPUT_FILES = [
    "blinded_sample_file_manifest.tsv",
    "blinded_sample_qc_metrics.tsv",
    "target_independent_guide_response.tsv",
    "blinded_qc_signoff.tsv",
]


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def number(value: object, field: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Non-finite {field}: {value!r}")
    return result


def bh_adjust(p_values: list[float]) -> list[float]:
    if any(not 0 <= value <= 1 for value in p_values):
        raise RuntimeError("BH input p-value outside [0,1]")
    n = len(p_values)
    order = sorted(range(n), key=lambda index: (p_values[index], index))
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[index] * n / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def metric_pass(value: float, operator: str, threshold: float) -> bool:
    if operator == "ge":
        return value >= threshold
    if operator == "le":
        return value <= threshold
    if operator == "eq":
        return value == threshold
    raise RuntimeError(f"Unsupported QC operator: {operator}")


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
        raise RuntimeError(f"Blinded-QC schema/source failure: {basename}")
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"Blinded-QC input is empty: {basename}")
    return rows


def validate_output_hashes(root: Path, seal: dict[str, object]) -> None:
    for name, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / str(name)
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Sealed blinded-QC dependency drift: {path}")


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite blinded Stage-A QC: {CANDIDATE_ROOT}")
    execution_root = candidate_source("PLAN45_STAGE_A_EXECUTION_ROOT")
    contract_root = candidate_source("PLAN45_STAGE_A_QC_CONTRACT_ROOT")
    threshold_root = candidate_source("PLAN45_STAGE_A_QC_THRESHOLD_ROOT")
    input_root = candidate_source("PLAN45_STAGE_A_QC_INPUT_ROOT")
    execution_seal_path = execution_root / "FROZEN_STAGE_A_EXECUTION.json"
    contract_seal_path = contract_root / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json"
    threshold_seal_path = threshold_root / "STAGE_A_QC_THRESHOLDS_FROZEN.json"
    execution_seal = json.loads(execution_seal_path.read_text(encoding="utf-8"))
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    threshold_seal = json.loads(threshold_seal_path.read_text(encoding="utf-8"))
    if execution_seal.get("status") != "stage_a_execution_frozen_ready_for_blinded_outcome_generation":
        raise RuntimeError("Blinded QC requires a validated actual Stage-A execution freeze")
    if execution_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Execution freeze reports scientific outcome access")
    if contract_seal.get("status") != "sealed_outcome_blind_stage_a_qc_unblinding_contract":
        raise RuntimeError("Invalid blinded-QC contract dependency")
    if contract_seal.get("scientific_condition_labels_opened") is not False:
        raise RuntimeError("QC contract reports condition-label access")
    if threshold_seal.get("status") != "stage_a_qc_thresholds_frozen_before_raw_qc":
        raise RuntimeError("Blinded QC requires a separate pre-QC threshold freeze")
    if threshold_seal.get("raw_qc_accessed_before_freeze") is not False:
        raise RuntimeError("Threshold release reports raw-QC access before freezing")
    validate_output_hashes(execution_root, execution_seal)
    validate_output_hashes(contract_root, contract_seal)
    validate_output_hashes(threshold_root, threshold_seal)

    contract_manifest = read_tsv(
        contract_root / "stage_a_qc_unblinding_contract_input_manifest.tsv"
    )
    contract_execution = next(
        row for row in contract_manifest if row["role"] == "execution_contract_seal"
    )
    execution_manifest = read_tsv(execution_root / "stage_a_execution_input_manifest.tsv")
    actual_execution_contract = next(
        row for row in execution_manifest
        if row["role"] == "stage_a_execution_contract_seal"
    )
    if (
        actual_execution_contract["source_path"] != contract_execution["source_path"]
        or actual_execution_contract["sha256"] != contract_execution["sha256"]
    ):
        raise RuntimeError("Actual execution freeze was built under a different execution contract")
    threshold_manifest = read_tsv(threshold_root / "qc_threshold_input_manifest.tsv")
    threshold_contract = next(
        row for row in threshold_manifest if row["role"] == "qc_unblinding_contract_seal"
    )
    if (
        threshold_contract["source_path"] != str(contract_seal_path.relative_to(PROJECT_ROOT))
        or threshold_contract["sha256"] != sha256_file(contract_seal_path)
    ):
        raise RuntimeError("Threshold release was built under a different QC contract")

    inputs = {name: exact_rows(input_root, contract_root, name) for name in INPUT_FILES}
    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    expected_sample = {
        row["sample_id"]: (row["blinded_label"], row["assay_id"])
        for row in randomization
    }
    if len(expected_sample) != len(randomization):
        raise RuntimeError("Execution freeze has duplicate sample IDs")

    files = inputs["blinded_sample_file_manifest.tsv"]
    if {row["sample_id"] for row in files} != set(expected_sample) or len(files) != len(expected_sample):
        raise RuntimeError("Blinded file manifest does not cover exactly the randomized samples")
    for row in files:
        if (row["blinded_label"], row["assay_id"]) != expected_sample[row["sample_id"]]:
            raise RuntimeError("Blinded file identity differs from the execution freeze")
        if not yes(row["file_readable"]) or yes(row["scientific_condition_inspected"]):
            raise RuntimeError("Blinded file audit opened a condition or retained an unreadable file")
        if not yes(row["review_concordant"]) or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Blinded file audit lacks independent review")
        path = (PROJECT_ROOT / row["raw_source_path"]).resolve()
        if PROJECT_ROOT.resolve() not in path.parents or not path.is_file():
            raise RuntimeError(f"Blinded raw source absent/outside project: {path}")
        if path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Blinded raw source drift: {path}")

    thresholds = read_tsv(threshold_root / "frozen_blinded_assay_qc_thresholds.tsv")
    if not thresholds:
        raise RuntimeError("Frozen QC threshold release is empty")
    threshold_by_key: dict[tuple[str, str], dict[str, str]] = {}
    assay_ids = set(expected[1] for expected in expected_sample.values())
    for row in thresholds:
        key = (row["assay_id"], row["metric_name"])
        if key in threshold_by_key or row["assay_id"] not in assay_ids:
            raise RuntimeError("QC threshold is duplicated or names an unknown assay")
        if row["operator"] not in {"ge", "le", "eq"}:
            raise RuntimeError("QC threshold operator is invalid")
        number(row["threshold"], "QC threshold")
        if yes(row["raw_qc_accessed_before_freeze"]) or not row["frozen_utc"]:
            raise RuntimeError("QC threshold was not frozen before raw-QC access")
        if not yes(row["review_concordant"]) or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("QC threshold lacks independent review")
        threshold_by_key[key] = row
    if {row["assay_id"] for row in thresholds} != assay_ids:
        raise RuntimeError("Every assay must have at least one frozen QC threshold")

    metrics = inputs["blinded_sample_qc_metrics.tsv"]
    metric_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in metrics:
        key = (row["sample_id"], row["metric_name"])
        if key in metric_by_key or row["sample_id"] not in expected_sample:
            raise RuntimeError("Blinded metric is duplicated or names an unknown sample")
        if (row["blinded_label"], row["assay_id"]) != expected_sample[row["sample_id"]]:
            raise RuntimeError("Blinded metric identity differs from the execution freeze")
        threshold = threshold_by_key.get((row["assay_id"], row["metric_name"]))
        if threshold is None or row["metric_unit"] != threshold["metric_unit"]:
            raise RuntimeError("Blinded metric lacks a matching frozen threshold/unit")
        if yes(row["scientific_condition_inspected"]):
            raise RuntimeError("Blinded metric audit opened a condition label")
        if not yes(row["review_concordant"]) or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Blinded metric lacks independent review")
        metric_path = (PROJECT_ROOT / row["metric_source_path"]).resolve()
        if PROJECT_ROOT.resolve() not in metric_path.parents or not metric_path.is_file() or sha256_file(metric_path) != row["metric_source_sha256"]:
            raise RuntimeError("Blinded metric source provenance failed")
        number(row["metric_value"], "QC metric")
        metric_by_key[key] = row
    expected_metrics = {
        (sample_id, metric_name)
        for sample_id, (_, assay_id) in expected_sample.items()
        for threshold_assay, metric_name in threshold_by_key
        if threshold_assay == assay_id
    }
    if set(metric_by_key) != expected_metrics:
        raise RuntimeError(
            f"Blinded metric universe drift: missing={len(expected_metrics-set(metric_by_key))} extra={len(set(metric_by_key)-expected_metrics)}"
        )

    status_rows = []
    for sample_id in sorted(expected_sample):
        label, assay_id = expected_sample[sample_id]
        failures = []
        sample_metrics = [
            row for (observed_sample, _), row in metric_by_key.items()
            if observed_sample == sample_id
        ]
        for row in sample_metrics:
            threshold = threshold_by_key[(assay_id, row["metric_name"])]
            if not metric_pass(
                number(row["metric_value"], "QC metric"), threshold["operator"],
                number(threshold["threshold"], "QC threshold"),
            ):
                failures.append(row["metric_name"])
        status_rows.append({
            "sample_id": sample_id,
            "blinded_label": label,
            "assay_id": assay_id,
            "n_qc_metrics": len(sample_metrics),
            "technical_qc_pass": str(not failures).lower(),
            "preunblinding_exclusion_status": (
                "included_after_blinded_qc" if not failures else "excluded_before_unblinding"
            ),
            "preunblinding_exclusion_reason": (
                "" if not failures else "failed_frozen_metrics:" + ";".join(sorted(failures))
            ),
            "scientific_condition_inspected": "false",
        })

    guide_rows = inputs["target_independent_guide_response.tsv"]
    if len(guide_rows) < 500 or len({row["gene_id"] for row in guide_rows}) != len(guide_rows):
        raise RuntimeError("Guide-response QC lacks a full, unique observable gene family")
    p_values = [number(row["p_value"], "guide-response p") for row in guide_rows]
    adjusted = bh_adjust(p_values)
    guide_outputs = []
    q_threshold = float(contract_seal["guide_response_q_threshold"])
    effect_threshold = float(contract_seal["guide_response_abs_effect_threshold"])
    for row, q_value in zip(guide_rows, adjusted):
        if row["contrast"] != "NTC_vs_untransduced" or yes(row["condition_labels_opened"]):
            raise RuntimeError("Guide-response QC is not target-independent/blinded")
        if int(row["n_backgrounds"]) < 3:
            raise RuntimeError("Guide-response QC has fewer than three backgrounds")
        if not yes(row["review_concordant"]) or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Guide-response QC lacks independent review")
        if abs(number(row["q_value"], "guide-response q") - q_value) > 1e-8:
            raise RuntimeError("Guide-response BH q-value does not rederive")
        effect = number(row["effect"], "guide-response effect")
        excluded = q_value < q_threshold and abs(effect) >= effect_threshold
        if yes(row["exclusion_gene"]) != excluded:
            raise RuntimeError("Guide-response exclusion flag does not rederive")
        guide_outputs.append({
            "gene_id": row["gene_id"],
            "contrast": row["contrast"],
            "effect": row["effect"],
            "p_value": row["p_value"],
            "q_value": f"{q_value:.17g}",
            "exclude_from_all_stage_a_axes": str(excluded).lower(),
        })

    signoffs = inputs["blinded_qc_signoff.tsv"]
    required_scopes = {"files", "thresholds", "metrics", "guide_response"}
    if {row["review_scope"] for row in signoffs} != required_scopes or len(signoffs) != len(required_scopes):
        raise RuntimeError("Blinded-QC signoff universe drift")
    if any(
        row["reviewer_1"] == row["reviewer_2"] or not yes(row["review_concordant"])
        or yes(row["scientific_condition_inspected"]) or not row["signed_utc"]
        for row in signoffs
    ):
        raise RuntimeError("Blinded-QC signoff is invalid")

    retention = []
    status_by_assay: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in status_rows:
        status_by_assay[row["assay_id"]].append(row)
    for assay_id in sorted(status_by_assay):
        rows = status_by_assay[assay_id]
        retention.append({
            "assay_id": assay_id,
            "n_randomized_samples": len(rows),
            "n_included_after_blinded_qc": sum(yes(row["technical_qc_pass"]) for row in rows),
            "n_excluded_before_unblinding": sum(not yes(row["technical_qc_pass"]) for row in rows),
        })

    input_manifest = []
    for name in INPUT_FILES:
        path = input_root / name
        input_manifest.append({
            "role": name.removesuffix(".tsv"),
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    for role, path in [
        ("stage_a_execution_seal", execution_seal_path),
        ("qc_unblinding_contract_seal", contract_seal_path),
        ("qc_threshold_seal", threshold_seal_path),
    ]:
        input_manifest.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })

    outputs: list[Path] = []
    status_path = CANDIDATE_ROOT / "blinded_sample_qc_status.tsv"
    write_tsv(status_path, status_rows, list(status_rows[0]))
    outputs.append(status_path)
    guide_path = CANDIDATE_ROOT / "guide_response_exclusion_genes.tsv"
    write_tsv(guide_path, guide_outputs, list(guide_outputs[0]))
    outputs.append(guide_path)
    retention_path = CANDIDATE_ROOT / "blinded_qc_design_retention.tsv"
    write_tsv(retention_path, retention, list(retention[0]))
    outputs.append(retention_path)
    manifest_path = CANDIDATE_ROOT / "blinded_qc_input_manifest.tsv"
    write_tsv(manifest_path, input_manifest, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)
    signoff_copy = CANDIDATE_ROOT / "frozen_blinded_qc_signoff.tsv"
    write_tsv(signoff_copy, signoffs, header(input_root / "blinded_qc_signoff.tsv"))
    outputs.append(signoff_copy)

    n_pass = sum(yes(row["technical_qc_pass"]) for row in status_rows)
    payload = {
        "status": (
            "blinded_stage_a_qc_sealed"
            if n_pass else "blinded_stage_a_qc_sealed_no_samples_pass"
        ),
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_randomized_samples": len(status_rows),
        "n_samples_included": n_pass,
        "n_samples_excluded": len(status_rows) - n_pass,
        "n_guide_response_genes_tested": len(guide_outputs),
        "n_guide_response_genes_excluded": sum(yes(row["exclude_from_all_stage_a_axes"]) for row in guide_outputs),
        "sample_exclusions_frozen": True,
        "guide_response_exclusions_frozen": True,
        "scientific_condition_labels_opened": False,
        "scientific_outcomes_inspected": False,
        "unblinding_authorized": False,
        "stage_b_design_frozen": False,
        "next_gate": "data-manager unblinding key validation after this seal is final",
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "BLINDED_STAGE_A_QC_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
