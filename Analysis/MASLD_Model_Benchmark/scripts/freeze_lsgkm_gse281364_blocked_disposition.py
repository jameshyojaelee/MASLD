#!/usr/bin/env python3
"""Audit and freeze the terminal LS-GKM matched-input activation failure."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


class BlockedDispositionError(RuntimeError):
    """Raised when the failed activation evidence differs from its selection record."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise BlockedDispositionError(f"JSON object differs: {path}")
    return value


def resolve_project_path(root: Path, relative_text: str) -> Path:
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts:
        raise BlockedDispositionError("unsafe project-relative path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    return path


def build_disposition(root: Path, config_path: Path) -> dict[str, Any]:
    config = load_json(config_path)
    if (
        config.get("schema_version")
        != "masld-bench-lsgkm-gse281364-blocked-disposition-input-v1"
        or config.get("status") != "freeze_terminal_blocked_disposition"
    ):
        raise BlockedDispositionError("blocked-disposition config differs")

    activation_binding = config["activation_config"]
    activation_config = resolve_project_path(root, activation_binding["path"])
    if file_sha256(activation_config) != activation_binding["sha256"]:
        raise BlockedDispositionError("activation config differs")
    activation = load_json(activation_config)
    if any(
        activation.get(field) is not False
        for field in (
            "outcome_access_authorized",
            "sealed_asset_access_authorized",
            "production_training_authorized",
            "production_prediction_authorized",
        )
    ):
        raise BlockedDispositionError("activation action firewall differs")

    evidence_binding = config["failed_evidence"]
    evidence_root = resolve_project_path(root, evidence_binding["path"])
    expected_files = dict(evidence_binding["files"])
    observed_files = {
        path.relative_to(evidence_root).as_posix()
        for path in evidence_root.rglob("*")
        if path.is_file()
    }
    if observed_files != set(expected_files):
        raise BlockedDispositionError("failed evidence file roster differs")
    for relative, expected_hash in expected_files.items():
        if file_sha256(evidence_root / relative) != expected_hash:
            raise BlockedDispositionError(f"failed evidence hash differs: {relative}")
    forbidden_suffixes = (".fa", ".fasta", ".fa.gz", ".fasta.gz", ".model", ".weights")
    if any(relative.endswith(forbidden_suffixes) or "prediction" in relative.lower() for relative in observed_files):
        raise BlockedDispositionError("forbidden production artifact exists")

    failure = load_json(evidence_root / "activation/contract/matching_gate_failure.json")
    expected = config["expected_failed_fit"]
    gate = failure["gate"]
    exact_values = {
        "split_id": failure["split_id"],
        "seed": failure["seed"],
        "genomic_test_fold": failure["genomic_test_fold"],
        "raw_positive_windows": failure["raw_positive_windows"],
        "eligible_unique_positive_windows": failure["eligible_unique_positive_windows"],
        "raw_training_candidate_windows": failure["raw_training_candidate_windows"],
        "eligible_unique_training_candidate_windows": failure["eligible_unique_training_candidate_windows"],
        "matched_pairs": gate["matched_pairs"],
        "minimum_matched_pairs": gate["minimum_matched_pairs"],
        "minimum_positive_coverage": gate["minimum_positive_coverage"],
    }
    if exact_values != {key: expected[key] for key in exact_values} or not math.isclose(
        float(gate["positive_coverage"]),
        float(expected["positive_coverage"]),
        rel_tol=0.0,
        abs_tol=1e-15,
    ):
        raise BlockedDispositionError("failed fit census differs")
    if (
        failure.get("status") != "terminally_blocked_prespecified_matching_gate"
        or gate.get("pair_count_gate_passed") is not False
        or gate.get("positive_coverage_gate_passed") is not False
        or gate.get("passed") is not False
        or failure.get("input_fastas_written") != 0
        or failure.get("production_fits_executed") != 0
        or failure.get("production_predictions_generated") != 0
        or failure.get("outcomes_read") is not False
        or failure.get("reporter_counts_read") is not False
        or failure.get("sealed_assets_read") is not False
    ):
        raise BlockedDispositionError("terminal failure firewall differs")

    design = config["production_design"]
    if (
        design.get("planned_shared_fits") != 25
        or design.get("evaluated_fit_units_before_terminal_stop") != 1
        or design.get("remaining_fit_units_not_materialized") != 24
        or design.get("production_submission_authorized") is not False
    ):
        raise BlockedDispositionError("blocked production design differs")
    terminal_rule = config["terminal_rule"]
    if (
        terminal_rule.get("thresholds_may_be_relaxed_in_this_campaign") is not False
        or terminal_rule.get("matching_tolerances_may_be_relaxed_in_this_campaign") is not False
        or terminal_rule.get("future_independent_redesign_requires_new_prespecification_and_campaign_id") is not True
    ):
        raise BlockedDispositionError("terminal redesign rule differs")

    return {
        "schema_version": "masld-bench-lsgkm-gse281364-terminal-disposition-v1",
        "status": "terminally_blocked_300bp_matched_negative_design",
        "family_id": "shared_lsgkm_gkmsvm_deltasvm",
        "activation_config_sha256": file_sha256(activation_config),
        "failed_evidence_file_sha256": expected_files,
        "first_failed_fit": failure,
        "production_design": design,
        "admitted_authorities": config["admitted_authorities"],
        "terminal_rule": terminal_rule,
        "claim_limits": config["claim_limits"],
        "input_fastas_written": 0,
        "production_fits_executed": 0,
        "production_predictions_generated": 0,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_assets_read": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise BlockedDispositionError("refusing to overwrite disposition")
    disposition = build_disposition(
        arguments.project_root.resolve(strict=True),
        arguments.config.resolve(strict=True),
    )
    arguments.output.write_text(
        json.dumps(disposition, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(disposition, sort_keys=True))


if __name__ == "__main__":
    main()
