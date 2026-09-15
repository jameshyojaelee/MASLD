#!/usr/bin/env python3
"""Freeze LS-GKM exposed-development evaluation choices before evaluator access."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


class EvaluationLockError(RuntimeError):
    """Raised when an evaluation authority or separation differs."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EvaluationLockError(f"JSON object differs: {path}")
    return value


def project_file(root: Path, relative_text: str, expected: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise EvaluationLockError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if path.is_symlink() or not path.is_file() or sha256(path) != expected:
        raise EvaluationLockError(f"file authority differs: {relative_text}")
    return path


def bind_manifest(root: Path, authority: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    path = project_file(root, authority["artifacts_path"], authority["artifacts_sha256"])
    manifest = load_json(path)
    roster = {
        row["path"]: row
        for row in manifest.get("artifacts", [])
        if isinstance(row, dict) and isinstance(row.get("path"), str)
    }
    if manifest.get("schema_version") != "masld-bench-artifacts-v1" or len(roster) != len(manifest.get("artifacts", [])):
        raise EvaluationLockError("artifact manifest differs")
    return manifest, roster


def require_record(roster: Mapping[str, Mapping[str, Any]], member: str, expected: str) -> None:
    record = roster.get(member)
    if record is None or record.get("sha256") != expected:
        raise EvaluationLockError(f"artifact member record differs: {member}")


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version")
        != "masld-bench-lsgkm-gse281364-dinucleotide-development-evaluation-input-v1"
        or config.get("status") != "freeze_development_evaluator_after_raw_prediction_lock"
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
    ):
        raise EvaluationLockError("evaluation identity differs")
    candidates = [
        (row.get("model_id"), row.get("head_id"))
        for row in config.get("candidates", [])
        if isinstance(row, dict)
    ]
    if candidates != [
        ("available_simple_controls", "allele_identity_ridge"),
        ("lsgkm_exact_dinucleotide_null", "direct_gkmsvm"),
        ("lsgkm_exact_dinucleotide_null", "deltasvm"),
    ]:
        raise EvaluationLockError("candidate roster differs")
    row = config["row_contract"]
    if (
        row.get("elements") != 1033
        or row.get("long_range_blocks") != 239
        or row.get("contexts") != ["HepG2_control", "HepG2_PAOA"]
        or row.get("fixed_seeds") != [1103, 2909, 4721, 6673, 8111]
        or row.get("rows_per_candidate") != 10330
        or row.get("allele_effect_sign") != "ALT_minus_REF"
    ):
        raise EvaluationLockError("row contract differs")
    metric = config["metric_contract"]
    uncertainty = config["uncertainty"]
    if (
        metric.get("primary")
        != "tanh_mean_Fisher_z_Spearman_across_two_assay_contexts"
        or metric.get("reference_candidate")
        != "available_simple_controls/allele_identity_ridge"
        or metric.get("rmse_or_calibration_for_uncalibrated_lsgkm") is not False
        or metric.get("readouts_count_as_independent_model_families") is not False
        or uncertainty.get("method") != "paired_long_range_block_bootstrap"
        or uncertainty.get("resamples") != 10000
        or uncertainty.get("seed") != 20260825
        or uncertainty.get("minimum_valid_resamples") != 9500
        or uncertainty.get("minimum_positive_gain_seeds") != 4
    ):
        raise EvaluationLockError("metric or uncertainty contract differs")
    firewall = config["stage_firewall"]
    if any(
        firewall.get(field) is not False
        for field in (
            "lock_builder_prediction_values_open_authorized",
            "lock_builder_control_prediction_values_open_authorized",
            "lock_builder_outcome_access_authorized",
            "model_fit_or_calibration_authorized",
            "candidate_selection_or_shortlist_authorized",
            "threshold_change_authorized",
            "sealed_asset_access_authorized",
            "external_claim_authorized",
            "champion_claim_authorized",
        )
    ) or any(
        firewall.get(field) is not True
        for field in (
            "evaluator_prediction_values_open_authorized",
            "evaluator_control_prediction_values_open_authorized",
            "evaluator_outcome_access_authorized",
            "evaluator_metric_calculation_authorized",
        )
    ):
        raise EvaluationLockError("stage firewall differs")


def freeze_lock(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    config = load_json(config_path)
    validate_config(config)
    if output.exists():
        raise EvaluationLockError("lock output exists")
    output.mkdir(parents=True, mode=0o750)
    for name in (
        "raw_prediction_authority",
        "independent_aggregation_audit",
        "row_authority",
    ):
        _manifest, roster = bind_manifest(root, config[name])
        require_record(roster, config[name]["member"], config[name]["member_sha256"])
    for name in ("outcome_authority", "shared_control_authority"):
        _manifest, roster = bind_manifest(root, config[name])
        require_record(roster, config[name]["member"], config[name]["member_sha256"])
    incident = config["protocol_incident"]
    incident_record = load_json(project_file(root, incident["path"], incident["sha256"]))
    if (
        incident_record.get("status")
        != "contained_narrow_control_value_exposure_no_label_or_metric_leakage"
        or incident_record.get("required_downstream_handling", {}).get(
            "freeze_lsgkm_candidate_prediction_hashes_before_control_or_outcome_evaluation"
        )
        is not True
    ):
        raise EvaluationLockError("protocol incident disposition differs")
    result = {
        "schema_version": "masld-bench-lsgkm-gse281364-development-evaluation-lock-v1",
        "status": "locked_before_evaluator_control_and_outcome_access",
        "config_sha256": sha256(config_path),
        "raw_prediction_artifacts_sha256": config["raw_prediction_authority"]["artifacts_sha256"],
        "raw_prediction_member_sha256": config["raw_prediction_authority"]["member_sha256"],
        "aggregation_audit_artifacts_sha256": config["independent_aggregation_audit"]["artifacts_sha256"],
        "candidates": config["candidates"],
        "row_contract": config["row_contract"],
        "metric_contract": config["metric_contract"],
        "uncertainty": config["uncertainty"],
        "reporting_limits": config["reporting_limits"],
        "protocol_incident_id": incident_record["incident_id"],
        "prediction_values_opened": False,
        "control_prediction_values_opened": False,
        "outcomes_opened": False,
        "metrics_calculated": False,
        "sealed_assets_opened": False,
    }
    (output / "selection_lock.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze_lock(
        args.project_root.resolve(strict=True),
        args.config.resolve(strict=True),
        args.output,
    )


if __name__ == "__main__":
    main()
