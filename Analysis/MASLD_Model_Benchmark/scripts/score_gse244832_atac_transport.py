#!/usr/bin/env python3
"""Score committed GSE244832 observed-ATAC predictions without condition labels."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from typing import Any, Callable, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree
from scripts.commit_gse244832_atac_transport_predictions import load_registration
from scripts.score_gse281367_atac_transport import (
    SCHEMA,
    ATACTransportScoringError,
    _digest,
    _load_json,
    deviance_per_insertion,
    spearman_or_zero,
)


DATASET_ID = "gse244832"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
DONOR_COUNT = 18
AXIS_CLASS = "gse244832_reference_guarded_atac_exchange_axis"
COMMIT_CLASS = "gse244832_atac_transport_prediction_commit"


def score_arrays(
    observed: Any,
    predicted: Any,
    missing: Any,
    *,
    block_ids: list[str],
    output_family: str,
) -> list[dict[str, Any]]:
    import numpy as np

    truth = np.asarray(observed)
    score = np.asarray(predicted)
    states = np.asarray(missing)
    profile = output_family.endswith("profile")
    if truth.shape != score.shape or states.shape != truth.shape[:3]:
        raise ATACTransportScoringError("scoring axes differ")
    if truth.shape[:3] != (DONOR_COUNT, 4, 16000) or len(block_ids) != 16000:
        raise ATACTransportScoringError("target transport geometry differs")
    rows: list[dict[str, Any]] = []
    positions_by_block = {
        block: np.asarray([index for index, value in enumerate(block_ids) if value == block])
        for block in sorted(set(block_ids))
    }
    for donor in range(DONOR_COUNT):
        for lineage in range(4):
            for block, positions in positions_by_block.items():
                block_states = states[donor, lineage, positions]
                keep = block_states == 0
                if not keep.any():
                    continue
                block_truth = truth[donor, lineage, positions][keep]
                block_score = score[donor, lineage, positions][keep]
                if profile:
                    count_truth = block_truth.sum(axis=1)
                    count_score = block_score.sum(axis=1)
                else:
                    count_truth, count_score = block_truth, block_score
                if block_truth.sum() <= 0:
                    continue
                rows.append(
                    {
                        "donor_index": donor,
                        "lineage_id": LINEAGES[lineage],
                        "block_id": block,
                        "windows_scored": int(keep.sum()),
                        "coverage_fraction": float(keep.mean()),
                        "deviance_per_insertion": deviance_per_insertion(block_truth, block_score),
                        "count_spearman": spearman_or_zero(count_truth, count_score),
                    }
                )
    if not rows:
        raise ATACTransportScoringError("no eligible donor/block score was produced")
    return rows


def validate_commit(
    commit_root: Path,
    expected_sha256: str,
    registration_artifacts_sha256: str,
) -> dict[str, Any]:
    commit_root = commit_root.resolve(strict=True)
    manifest = verify_frozen_tree(commit_root)
    if (
        _digest(commit_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != COMMIT_CLASS
        or manifest.get("metadata", {}).get("sequence_execution_authorized") is not False
        or manifest.get("metadata", {}).get("execution_registration_artifacts_sha256")
        != registration_artifacts_sha256
        or manifest.get("metadata", {}).get("status") != "passed"
    ):
        raise ATACTransportScoringError("prediction commit authority differs")
    receipt = _load_json(commit_root / "commit.json")
    if (
        receipt.get("schema_version") != "masld-bench-atac-transport-prediction-commit-v1"
        or receipt.get("status") != "committed_before_outcome_open"
        or receipt.get("input_regime") != "observed_atac"
        or receipt.get("prediction_layout") != "donor_lineage"
        or receipt.get("rna_assay_read") is not False
        or receipt.get("development_outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
        or receipt.get("sequence_execution_authorized") is not False
        or receipt.get("execution_registration_artifacts_sha256")
        != registration_artifacts_sha256
    ):
        raise ATACTransportScoringError("prediction was not blindly committed")
    return receipt


def _read_windows(axis_root: Path, target_role: str) -> tuple[list[int], list[str]]:
    with (axis_root / "windows.tsv").open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "window_index", "window_id", "role", "contig", "start", "end"
        ):
            raise ATACTransportScoringError("label-free window axis differs")
        rows = [row for row in reader if row["role"] == target_role]
    if len(rows) != 16000:
        raise ATACTransportScoringError("target-role window count differs")
    return [int(row["window_index"]) for row in rows], [row["contig"] for row in rows]


def load_outcomes(
    outcome_root: Path,
    outcome_artifacts_sha256: str,
    axis_root: Path,
    output_family: str,
    target_role: str,
) -> tuple[Any, Any, list[str]]:
    """Open biological arrays only after validate_commit has returned."""

    import numpy as np

    outcome_root = outcome_root.resolve(strict=True)
    manifest = verify_frozen_tree(outcome_root)
    if (
        _digest(outcome_root / "ARTIFACTS.json") != outcome_artifacts_sha256
        or manifest.get("metadata", {}).get("artifact_class") != "atac_transport_evaluator_outcomes"
    ):
        raise ATACTransportScoringError("outcome authority differs")
    positions, blocks = _read_windows(axis_root, target_role)
    source = outcome_root / "evaluator_outcomes/gse244832"
    if output_family.endswith("profile"):
        observed = np.load(
            source / "unit_fragment_profile.uint32.npy", mmap_mode="r", allow_pickle=False
        )
    else:
        observed = np.load(
            source / "unit_fragment_total.uint32.npy", mmap_mode="r", allow_pickle=False
        )
    eligible = np.load(
        source / "eligible_min_50_cells.bool.npy", mmap_mode="r", allow_pickle=False
    )
    if observed.shape[:3] != (DONOR_COUNT, 4, 32000) or eligible.shape != (DONOR_COUNT, 4):
        raise ATACTransportScoringError("outcome array geometry differs")
    return observed[:, :, positions], eligible, blocks


def _write_tsv(path: Path, rows: list[Mapping[str, Any]]) -> None:
    fields = (
        "donor_index",
        "lineage_id",
        "block_id",
        "windows_scored",
        "coverage_fraction",
        "deviance_per_insertion",
        "count_spearman",
    )
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def score_committed(
    *,
    commit_root: Path,
    commit_artifacts_sha256: str,
    prediction_root: Path,
    outcome_root: Path,
    outcome_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
    registration_root: Path,
    registration_artifacts_sha256: str,
    output: Path,
    outcome_loader: Callable[..., tuple[Any, Any, list[str]]] = load_outcomes,
) -> dict[str, Any]:
    import numpy as np

    if output.exists():
        raise ATACTransportScoringError("refusing to overwrite scoring output")
    commit = validate_commit(
        commit_root,
        commit_artifacts_sha256,
        registration_artifacts_sha256,
    )
    registration, unit_states = load_registration(
        registration_root,
        registration_artifacts_sha256,
        axis_artifacts_sha256,
    )
    axis_manifest = verify_frozen_tree(axis_root)
    axis_summary = _load_json(axis_root / "exchange_summary.json")
    if (
        _digest(axis_root / "ARTIFACTS.json") != axis_artifacts_sha256
        or commit.get("axis_artifacts_sha256") != axis_artifacts_sha256
        or axis_manifest.get("metadata", {}).get("artifact_class") != AXIS_CLASS
        or axis_manifest.get("metadata", {}).get("sequence_execution_authorized") is not False
        or axis_summary.get("non_observed_atac_executable_rows") != 0
        or commit.get("execution_registration_artifacts_sha256")
        != registration_artifacts_sha256
        or commit.get("model_id") not in registration.get("allowed_model_ids", ())
    ):
        raise ATACTransportScoringError("axis/commit binding differs")
    verify_frozen_tree(prediction_root)
    prediction_manifest = _load_json(prediction_root / "prediction_bundle.json")
    if _digest(prediction_root / "ARTIFACTS.json") != commit.get("prediction_artifacts_sha256"):
        raise ATACTransportScoringError("prediction tree differs after commit")
    arrays = prediction_manifest["arrays"]
    predicted = np.load(
        prediction_root / arrays["predictions"]["path"], mmap_mode="r", allow_pickle=False
    )
    missing = np.load(
        prediction_root / arrays["missing_state"]["path"], mmap_mode="r", allow_pickle=False
    )
    target_role = "test" if commit["rotation_id"] == "valid_context_test_target" else "valid"
    observed, eligible, blocks = outcome_loader(
        outcome_root,
        outcome_artifacts_sha256,
        axis_root,
        commit["output_family"],
        target_role,
    )
    if predicted.shape != observed.shape or missing.shape != observed.shape[:3]:
        raise ATACTransportScoringError("prediction/outcome geometry differs")
    if not np.array_equal(np.asarray(eligible), np.asarray(unit_states) == 0):
        raise ATACTransportScoringError("registered and evaluator eligibility differ")
    rows = score_arrays(
        observed,
        predicted,
        missing,
        block_ids=blocks,
        output_family=commit["output_family"],
    )
    output.mkdir(parents=True, mode=0o750)
    _write_tsv(output / "donor_lineage_block_metrics.tsv", rows)
    by_lineage: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        by_lineage[str(row["lineage_id"])].append(float(row["deviance_per_insertion"]))
    summary = {
        "schema_version": SCHEMA,
        "status": "condition_blind_development_transport_scored",
        "model_id": commit["model_id"],
        "input_regime": commit["input_regime"],
        "output_family": commit["output_family"],
        "rotation_id": commit["rotation_id"],
        "commit_artifacts_sha256": commit_artifacts_sha256,
        "execution_registration_artifacts_sha256": registration_artifacts_sha256,
        "metric_rows": len(rows),
        "mean_deviance_by_lineage": {
            key: sum(values) / len(values) for key, values in sorted(by_lineage.items())
        },
        "condition_labels_read": False,
        "condition_labels_used": False,
        "phenotype_values_read": False,
        "rna_assay_read": False,
        "sequence_execution_authorized": False,
        "gse244832_champion_claim_allowed": False,
        "family_native_outputs_ranked_separately": True,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    freeze_tree(
        output,
        {
            "artifact_class": "gse244832_atac_transport_condition_blind_score",
            "model_id": commit["model_id"],
            "commit_artifacts_sha256": commit_artifacts_sha256,
            "execution_registration_artifacts_sha256": registration_artifacts_sha256,
            "condition_labels_used": False,
            "rna_assay_read": False,
            "sequence_execution_authorized": False,
            "champion_claim_allowed": False,
            "status": "passed",
        },
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--commit-root", required=True, type=Path)
    parser.add_argument("--commit-artifacts-sha256", required=True)
    parser.add_argument("--prediction-root", required=True, type=Path)
    parser.add_argument("--outcome-root", required=True, type=Path)
    parser.add_argument("--outcome-artifacts-sha256", required=True)
    parser.add_argument("--axis-root", required=True, type=Path)
    parser.add_argument("--axis-artifacts-sha256", required=True)
    parser.add_argument("--registration-root", required=True, type=Path)
    parser.add_argument("--registration-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(json.dumps(score_committed(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
