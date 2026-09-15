#!/usr/bin/env python3
"""Evaluate a frozen scBasset validation-fold roster after prediction record."""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy import sparse

from scripts import evaluate_scbasset_valid as base


EXPECTED_FOLDS = (1, 2, 3, 4)
EXPECTED_SEED = 20260824
LOCK_SCHEMA = "masld-bench-scbasset-crossfold-prediction-lock-v1"


class ScBassetCrossfoldEvaluationError(ValueError):
    """Raised when the prediction record or crossed-fold requirement differs."""


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ScBassetCrossfoldEvaluationError(f"missing or linked {label}: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScBassetCrossfoldEvaluationError(f"invalid {label}: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetCrossfoldEvaluationError(f"{label} must be an object")
    return value


def _verify_artifact_root(path: str, expected_sha256: str, label: str) -> Path:
    root = Path(path)
    if not root.is_absolute():
        raise ScBassetCrossfoldEvaluationError(f"{label} path must be absolute")
    root = root.resolve(strict=True)
    if base.sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise ScBassetCrossfoldEvaluationError(f"{label} ARTIFACTS hash differs")
    return root


def load_prediction_lock(
    lock_root: Path, lock_artifacts_sha256: str
) -> list[dict[str, Any]]:
    """Verify all predictions and inputs before any outcome authority is opened."""

    # The Python 3.11 driver fully verifies every frozen tree before entering
    # the Python 3.10 scientific runtime.  The evaluator independently checks
    # the reviewed ARTIFACTS bindings and never imports the control plane.
    root = lock_root.resolve(strict=True)
    if base.sha256_file(root / "ARTIFACTS.json") != lock_artifacts_sha256:
        raise ScBassetCrossfoldEvaluationError("prediction-lock ARTIFACTS hash differs")
    lock = _read_json(root / "prediction_lock.json", "prediction lock")
    if (
        lock.get("schema_version") != LOCK_SCHEMA
        or lock.get("status") != "pass"
        or lock.get("seed") != EXPECTED_SEED
        or lock.get("folds") != list(EXPECTED_FOLDS)
        or lock.get("outcomes_read") is not False
    ):
        raise ScBassetCrossfoldEvaluationError("prediction-lock header differs")
    records = lock.get("records")
    if not isinstance(records, list) or len(records) != len(EXPECTED_FOLDS):
        raise ScBassetCrossfoldEvaluationError("prediction-lock roster differs")

    validated: list[dict[str, Any]] = []
    raw_donors: set[str] = set()
    for expected_fold, record in zip(EXPECTED_FOLDS, records):
        if not isinstance(record, Mapping) or int(record.get("outer_fold", -1)) != expected_fold:
            raise ScBassetCrossfoldEvaluationError("prediction-lock fold order differs")
        input_root = _verify_artifact_root(
            str(record.get("input_root", "")),
            str(record.get("input_artifacts_sha256", "")),
            f"fold {expected_fold} input",
        )
        prediction_root = _verify_artifact_root(
            str(record.get("prediction_root", "")),
            str(record.get("prediction_artifacts_sha256", "")),
            f"fold {expected_fold} prediction",
        )
        summary = _read_json(input_root / "inputs/summary.json", "input summary")
        receipt = _read_json(
            prediction_root / "predictions/prediction_receipt.json",
            "prediction receipt",
        )
        artifact_manifest = _read_json(
            prediction_root / "ARTIFACTS.json", "prediction ARTIFACTS"
        )
        split_id = f"donor{expected_fold}_genomic{expected_fold}"
        valid_fold = (expected_fold + 1) % 5
        if (
            summary.get("status") != "pass"
            or summary.get("split_id") != split_id
            or summary.get("donor_test_fold") != expected_fold
            or summary.get("donor_valid_fold") != valid_fold
            or summary.get("held_donor_atac_exported") is not False
            or summary.get("genomic_test_atac_exported") is not False
            or receipt.get("status") != "pass"
            or receipt.get("role") != "valid"
            or receipt.get("split_id") != split_id
            or receipt.get("held_donor_atac_used") is not False
            or receipt.get("held_cell_embedding_available") is not False
            or artifact_manifest.get("metadata", {}).get("seed") != EXPECTED_SEED
        ):
            raise ScBassetCrossfoldEvaluationError(
                f"fold {expected_fold} input/prediction contract differs"
            )
        donors = base.read_tsv(
            input_root / "inputs/held_valid_donors.tsv", base.DONOR_FIELDS
        )
        donor_ids = [row["donor_id"] for row in donors]
        if (
            len(donor_ids) != int(receipt.get("held_donors", -1))
            or len(set(donor_ids)) != len(donor_ids)
            or any(int(row["outer_fold"]) != valid_fold for row in donors)
            or any(row["evaluation_role"] != "valid" for row in donors)
            or raw_donors.intersection(donor_ids)
        ):
            raise ScBassetCrossfoldEvaluationError(
                f"fold {expected_fold} held-donor roster differs"
            )
        raw_donors.update(donor_ids)
        validated.append(
            {
                **dict(record),
                "input_root": input_root,
                "prediction_root": prediction_root,
                "split_id": split_id,
                "valid_fold": valid_fold,
                "summary": summary,
                "receipt": receipt,
                "donors": donors,
            }
        )
    return validated


def _load_dynamic_bundle(
    prediction_root: Path, expected_rows: int, expected_donors: int
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    bundle_path = prediction_root / "predictions/bundle/prediction_bundle.json"
    bundle = _read_json(bundle_path, "prediction bundle")
    metadata = bundle.get("metadata", {})
    if (
        bundle.get("task_id") != "rna_conditioned_atac"
        or bundle.get("model_id") != "scbasset"
        or bundle.get("n_predictions") != expected_rows
        or metadata.get("evaluation_role") != "valid"
        or metadata.get("held_donor_count") != expected_donors
        or metadata.get("region_count") != 16_000
        or metadata.get("held_atac_input_exposed") is not False
        or metadata.get("observed_atac_exported") is not False
        or metadata.get("rna_input_exposed") is not False
    ):
        raise ScBassetCrossfoldEvaluationError("prediction bundle contract differs")
    record = bundle.get("standardized_table", {})
    table = bundle_path.parent / str(record.get("path", ""))
    if (
        not table.is_file()
        or table.is_symlink()
        or table.stat().st_size != record.get("size_bytes")
        or base.sha256_file(table) != record.get("sha256")
    ):
        raise ScBassetCrossfoldEvaluationError("prediction table differs")
    rows = base.read_tsv(table, base.PREDICTION_FIELDS)
    if len(rows) != expected_rows:
        raise ScBassetCrossfoldEvaluationError("prediction row count differs")
    return bundle, rows


def evaluate_fold(
    record: Mapping[str, Any],
    *,
    donor_bigwigs: Path,
    split_root: Path,
    workers: int,
    output: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    """Evaluate one already fixed validation fold against evaluator-only ATAC."""

    fold = int(record["outer_fold"])
    valid_fold = int(record["valid_fold"])
    input_root = Path(record["input_root"])
    prediction_root = Path(record["prediction_root"])
    summary = record["summary"]
    donors = record["donors"]
    donor_ids = [row["donor_id"] for row in donors]
    n_donors = len(donor_ids)
    expected_rows = n_donors * len(base.LINEAGES) * 16_000
    bundle, prediction_rows = _load_dynamic_bundle(
        prediction_root, expected_rows, n_donors
    )
    namespace = str(bundle["unit_id_namespace"])
    inputs = input_root / "inputs"
    regions = base.read_tsv(inputs / "valid_regions.tsv", base.REGION_FIELDS)
    cells = base.read_tsv(inputs / "training_cells.tsv", base.CELL_FIELDS)
    training_cells = int(summary["training_cells"])
    if (
        len(regions) != 16_000
        or len(cells) != training_cells
        or len(donors) != int(record["receipt"]["held_donors"])
        or any(row["role"] != "valid" for row in regions)
    ):
        raise ScBassetCrossfoldEvaluationError(f"fold {fold} validation axes differ")

    donor_index = {value: index for index, value in enumerate(donor_ids)}
    lineage_index = {value: index for index, value in enumerate(base.LINEAGES)}
    region_index = {row["region_id"]: index for index, row in enumerate(regions)}
    split_windows = {
        row["window_id"]: row
        for row in base.read_tsv(
            split_root / "ccre_evaluation_windows.tsv", base.CCRE_FIELDS
        )
        if int(row["genomic_fold"]) == valid_fold
    }
    if set(split_windows) != set(region_index):
        raise ScBassetCrossfoldEvaluationError(
            f"fold {fold} validation cCREs differ from split"
        )
    for row in regions:
        window = split_windows[row["region_id"]]
        midpoint = (int(window["output_start"]) + int(window["output_end"])) // 2
        if (
            row["block_id"] != window["contig"]
            or int(row["sequence_start"]) != midpoint - 672
            or int(row["sequence_end"]) != midpoint + 672
        ):
            raise ScBassetCrossfoldEvaluationError(
                f"fold {fold} sequence/output window join differs"
            )

    expected: dict[str, tuple[int, int, int]] = {}
    for donor in donor_ids:
        for lineage in base.LINEAGES:
            for region in regions:
                row_hash = base.join_hash(
                    namespace,
                    "row",
                    f"{donor}\0{lineage}\0{region['region_id']}",
                )
                expected[row_hash] = (
                    donor_index[donor],
                    lineage_index[lineage],
                    region_index[region["region_id"]],
                )
    candidate = np.full((n_donors, 5, 16_000), np.nan, dtype=np.float64)
    for row in prediction_rows:
        position = expected.pop(row["row_hash"], None)
        if position is None:
            raise ScBassetCrossfoldEvaluationError("prediction row identity differs")
        donor = donor_ids[position[0]]
        lineage = base.LINEAGES[position[1]]
        region = regions[position[2]]
        if (
            row["donor_hash"] != base.join_hash(namespace, "unit", donor)
            or row["block_hash"]
            != base.join_hash(namespace, "block", region["block_id"])
            or row["stratum"] != lineage
        ):
            raise ScBassetCrossfoldEvaluationError("prediction join fields differ")
        value = float(row["predicted"])
        if not math.isfinite(value) or value <= 0:
            raise ScBassetCrossfoldEvaluationError("prediction value differs")
        candidate[position] = value
    if expected or np.any(~np.isfinite(candidate)):
        raise ScBassetCrossfoldEvaluationError("prediction universe is incomplete")

    training = sparse.load_npz(inputs / "m_valid.npz").tocsr()
    if training.shape != (16_000, training_cells):
        raise ScBassetCrossfoldEvaluationError("training-only validation axes differ")
    training_labels = np.asarray([row["lineage"] for row in cells])
    lineage_mean = np.empty((5, 16_000), dtype=np.float64)
    for lineage, index in lineage_index.items():
        positions = np.flatnonzero(training_labels == lineage)
        if positions.size == 0:
            raise ScBassetCrossfoldEvaluationError("training lineage is empty")
        lineage_mean[index] = np.asarray(training[:, positions].mean(axis=1)).ravel()
    global_mean = np.asarray(training.mean(axis=1)).ravel()
    lineage_mean += 1.0e-8
    global_mean += 1.0e-8
    lineage_mean /= lineage_mean.sum(axis=1, keepdims=True)
    global_mean /= global_mean.sum()

    manifest = base.read_tsv(
        donor_bigwigs / "bigwig_manifest.tsv", base.BIGWIG_FIELDS
    )
    manifest_index = {
        (row["donor_id"], row["lineage_id"]): row
        for row in manifest
        if row["donor_id"] in donor_index and row["lineage_id"] in lineage_index
    }
    if len(manifest_index) != n_donors * len(base.LINEAGES):
        raise ScBassetCrossfoldEvaluationError("held-valid bigWig census differs")
    windows = [
        (
            split_windows[row["region_id"]]["contig"],
            int(split_windows[row["region_id"]]["output_start"]),
            int(split_windows[row["region_id"]]["output_end"]),
        )
        for row in regions
    ]
    extraction = [
        {
            **row,
            "lineage_id": lineage,
            "path": str(donor_bigwigs / row["path"]),
            "windows": windows,
        }
        for (_, lineage), row in sorted(manifest_index.items())
    ]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        extracted = list(executor.map(base._extract_bigwig, extraction))
    observed = np.zeros((n_donors, 5, 16_000), dtype=np.uint32)
    for item in extracted:
        observed[
            donor_index[item["donor_id"]], lineage_index[item["lineage"]]
        ] = item["counts"]
    if np.any(observed.sum(axis=2) <= 0):
        raise ScBassetCrossfoldEvaluationError("held-valid outcome is empty")

    block_positions = {
        block: np.asarray(
            [index for index, region in enumerate(regions) if region["block_id"] == block],
            dtype=np.int64,
        )
        for block in sorted({row["block_id"] for row in regions})
    }
    unit_rows: list[dict[str, Any]] = []
    full_rows: list[dict[str, Any]] = []
    for donor_offset, donor in enumerate(donor_ids):
        for lineage_offset, lineage in enumerate(base.LINEAGES):
            truth = observed[donor_offset, lineage_offset].astype(np.float64)
            vectors = {
                "scbasset": candidate[donor_offset, lineage_offset],
                "training_lineage_mean": lineage_mean[lineage_offset],
                "training_global_mean": global_mean,
            }
            for model_id, vector in vectors.items():
                full_rows.append(
                    {
                        "outer_fold": fold,
                        "model_id": model_id,
                        "donor_hash": base.join_hash(namespace, "unit", donor),
                        "stratum": lineage,
                        "peak_auprc": format(base.average_precision(truth, vector), ".17g"),
                        "profile_spearman": format(
                            base.spearman_or_zero(truth, vector), ".17g"
                        ),
                    }
                )
                for block, positions in block_positions.items():
                    block_truth = truth[positions]
                    if block_truth.sum() <= 0:
                        continue
                    unit_rows.append(
                        {
                            "outer_fold": fold,
                            "model_id": model_id,
                            "donor_hash": base.join_hash(namespace, "unit", donor),
                            "block_hash": base.join_hash(namespace, "block", block),
                            "stratum": lineage,
                            "observed_insertions": int(block_truth.sum()),
                            "regions": len(positions),
                            "deviance_per_insertion": format(
                                base.deviance_per_insertion(
                                    block_truth, vector[positions]
                                ),
                                ".17g",
                            ),
                        }
                    )
    output.mkdir(parents=True, mode=0o750)
    base.write_tsv(
        output / "donor_lineage_block_metrics.tsv",
        (
            "outer_fold",
            "model_id",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed_insertions",
            "regions",
            "deviance_per_insertion",
        ),
        unit_rows,
    )
    base.write_tsv(
        output / "donor_lineage_secondary_metrics.tsv",
        (
            "outer_fold",
            "model_id",
            "donor_hash",
            "stratum",
            "peak_auprc",
            "profile_spearman",
        ),
        full_rows,
    )
    fold_summary = summarize_unit_rows(unit_rows)
    fold_summary.update(
        {
            "outer_fold": fold,
            "valid_fold": valid_fold,
            "split_id": record["split_id"],
            "n_donors": n_donors,
            "prediction_artifacts_sha256": record[
                "prediction_artifacts_sha256"
            ],
            "input_artifacts_sha256": record["input_artifacts_sha256"],
        }
    )
    (output / "evaluation.json").write_text(
        json.dumps(fold_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return fold_summary, unit_rows, full_rows


def summarize_unit_rows(rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    models = ("scbasset", "training_lineage_mean", "training_global_mean")
    if not rows or {str(row["model_id"]) for row in rows} != set(models):
        raise ScBassetCrossfoldEvaluationError("unit-metric model roster differs")
    model_mean = {
        model: sum(
            float(row["deviance_per_insertion"])
            for row in rows
            if row["model_id"] == model
        )
        / sum(row["model_id"] == model for row in rows)
        for model in models
    }
    strongest = min(
        ("training_lineage_mean", "training_global_mean"),
        key=lambda model: (model_mean[model], model),
    )
    relative = (
        model_mean[strongest] - model_mean["scbasset"]
    ) / model_mean[strongest]
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["model_id"]), str(row["stratum"]))].append(
            float(row["deviance_per_insertion"])
        )
    lineage_gain = {}
    for lineage in base.LINEAGES:
        baseline = grouped[(strongest, lineage)]
        candidate = grouped[("scbasset", lineage)]
        if not baseline or not candidate:
            raise ScBassetCrossfoldEvaluationError("lineage metric roster differs")
        baseline_mean = sum(baseline) / len(baseline)
        candidate_mean = sum(candidate) / len(candidate)
        lineage_gain[lineage] = (baseline_mean - candidate_mean) / baseline_mean
    return {
        "mean_deviance_per_insertion": model_mean,
        "strongest_training_only_baseline": strongest,
        "scbasset_relative_deviance_reduction": relative,
        "lineage_relative_deviance_reduction": lineage_gain,
        "development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_lineages": 4,
            "maximum_allowed_lineage_worsening": -0.02,
            "overall_threshold_passed": relative >= 0.05,
            "improved_lineages": sum(value > 0 for value in lineage_gain.values()),
            "no_lineage_worse_than_threshold": min(lineage_gain.values()) >= -0.02,
        },
    }


def evaluate(
    *,
    prediction_lock: Path,
    prediction_lock_artifacts_sha256: str,
    donor_bigwigs: Path,
    donor_bigwigs_artifacts_sha256: str,
    split_root: Path,
    split_artifacts_sha256: str,
    output: Path,
    workers: int,
) -> dict[str, Any]:
    if output.exists() or not 1 <= workers <= 16:
        raise ScBassetCrossfoldEvaluationError("output or worker contract differs")
    records = load_prediction_lock(
        prediction_lock, prediction_lock_artifacts_sha256
    )
    # Outcome authorities are resolved only after every prediction/input verifies.
    donor_bigwigs = _verify_artifact_root(
        donor_bigwigs.as_posix(),
        donor_bigwigs_artifacts_sha256,
        "donor bigWigs",
    )
    split_root = _verify_artifact_root(
        split_root.as_posix(), split_artifacts_sha256, "sequence split"
    )
    output.mkdir(parents=True, mode=0o750)
    all_unit_rows: list[dict[str, Any]] = []
    all_full_rows: list[dict[str, Any]] = []
    fold_summaries = []
    for record in records:
        summary, unit_rows, full_rows = evaluate_fold(
            record,
            donor_bigwigs=donor_bigwigs,
            split_root=split_root,
            workers=workers,
            output=output / "folds" / f"fold{record['outer_fold']}",
        )
        fold_summaries.append(summary)
        all_unit_rows.extend(unit_rows)
        all_full_rows.extend(full_rows)
    base.write_tsv(
        output / "donor_lineage_block_metrics.tsv",
        (
            "outer_fold",
            "model_id",
            "donor_hash",
            "block_hash",
            "stratum",
            "observed_insertions",
            "regions",
            "deviance_per_insertion",
        ),
        all_unit_rows,
    )
    base.write_tsv(
        output / "donor_lineage_secondary_metrics.tsv",
        (
            "outer_fold",
            "model_id",
            "donor_hash",
            "stratum",
            "peak_auprc",
            "profile_spearman",
        ),
        all_full_rows,
    )
    result = {
        "schema_version": "masld-bench-scbasset-crossfold-valid-evaluation-v1",
        "status": "pass",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "evaluation_role": "valid",
        "seed": EXPECTED_SEED,
        "outer_folds": list(EXPECTED_FOLDS),
        "missing_outer_fold": 0,
        "n_donors": sum(int(item["n_donors"]) for item in fold_summaries),
        "n_lineages": len(base.LINEAGES),
        "n_regions_per_fold": 16_000,
        "one_seed_four_crossed_fold_screen": True,
        "folds_are_not_independent_biological_replicates": True,
        "promotion_gate_finalized": False,
        "champion_claim_allowed": False,
        "biological_unit": "donor",
        "cells_used_as_independent_replicates": False,
        "prediction_lock_verified_before_outcomes": True,
        "held_valid_atac_read_by_evaluator_only": True,
        "held_valid_atac_exposed_to_model": False,
        "test_regions_read": False,
        "test_atac_read": False,
        "test_predictions_read": False,
        "external_evaluation": False,
        "prediction_lock_artifacts_sha256": prediction_lock_artifacts_sha256,
        "donor_bigwigs_artifacts_sha256": donor_bigwigs_artifacts_sha256,
        "split_artifacts_sha256": split_artifacts_sha256,
        "fold_summaries": fold_summaries,
        **summarize_unit_rows(all_unit_rows),
    }
    (output / "evaluation.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-lock", type=Path, required=True)
    parser.add_argument("--prediction-lock-artifacts-sha256", required=True)
    parser.add_argument("--donor-bigwigs", type=Path, required=True)
    parser.add_argument("--donor-bigwigs-artifacts-sha256", required=True)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--split-artifacts-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    arguments = parser.parse_args()
    print(json.dumps(evaluate(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
