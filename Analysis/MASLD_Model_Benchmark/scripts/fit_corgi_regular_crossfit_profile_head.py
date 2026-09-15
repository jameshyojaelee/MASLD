#!/usr/bin/env python3
"""Fit training-only Corgi profile heads and export held-fold predictions."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from scripts.evaluate_corgi_regional_smoke import (
    CONTEXT_ARMS,
    PREDICTION_FIELDS,
    SCORE_NAME,
    _artifact_member,
    _artifact_metadata,
    _bigwig_sums,
    _load_outcome_paths,
    multinomial_deviance_per_insertion,
    read_tsv,
    validate_prediction_tables,
    verify_complete_marker,
    verify_selected_member,
)


SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-fit-v1"
CONTRACT_SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-preflight-v1"
ALPHA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
PRIMARY_CONTEXT_ARM = "actual_released_rank_masked"
EPSILON_TOTAL_MASS = 1.0e-6
FOLDS = tuple(range(5))


class CorgiHeadFitError(RuntimeError):
    """Raised when the cross-fitted head or outcome separation differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def smooth_profile(values: Any, epsilon_total_mass: float = EPSILON_TOTAL_MASS) -> np.ndarray:
    vector = np.asarray(values, dtype=np.float64)
    if (
        vector.ndim != 1
        or vector.size < 3
        or np.any(~np.isfinite(vector))
        or np.any(vector < 0)
        or not 0.0 < epsilon_total_mass < 0.01
    ):
        raise CorgiHeadFitError("profile vector differs")
    total = float(vector.sum())
    probability = (
        np.full(vector.shape, 1.0 / vector.size, dtype=np.float64)
        if total <= 0
        else vector / total
    )
    return (1.0 - epsilon_total_mass) * probability + epsilon_total_mass / vector.size


def geometric_profile(corgi: Any, baseline: Any, alpha: float) -> np.ndarray:
    if alpha not in ALPHA_GRID:
        raise CorgiHeadFitError("profile-head alpha differs")
    left = smooth_profile(corgi)
    right = smooth_profile(baseline)
    if left.shape != right.shape:
        raise CorgiHeadFitError("Corgi and baseline profile shapes differ")
    log_score = alpha * np.log(left) + (1.0 - alpha) * np.log(right)
    score = np.exp(log_score - float(np.max(log_score)))
    probability = score / float(score.sum())
    if np.any(~np.isfinite(probability)) or np.any(probability <= 0):
        raise CorgiHeadFitError("profile-head probability differs")
    return probability


def select_alpha(losses: Mapping[float, Sequence[float]]) -> tuple[float, list[dict[str, float | int]]]:
    if set(losses) != set(ALPHA_GRID):
        raise CorgiHeadFitError("alpha loss grid differs")
    rows: list[dict[str, float | int]] = []
    for alpha in ALPHA_GRID:
        values = np.asarray(losses[alpha], dtype=np.float64)
        if values.ndim != 1 or values.size < 2 or np.any(~np.isfinite(values)):
            raise CorgiHeadFitError("alpha training losses differ")
        rows.append(
            {
                "alpha": alpha,
                "units": int(values.size),
                "mean_fit_loss": float(np.mean(values)),
                "standard_error": float(np.std(values, ddof=1) / math.sqrt(values.size)),
            }
        )
    minimum = min(rows, key=lambda row: (float(row["mean_fit_loss"]), float(row["alpha"])))
    threshold = float(minimum["mean_fit_loss"]) + float(minimum["standard_error"])
    eligible = [row for row in rows if float(row["mean_fit_loss"]) <= threshold + 1.0e-15]
    selected = min(float(row["alpha"]) for row in eligible)
    for row in rows:
        row["one_se_threshold"] = threshold
        row["selected"] = int(float(row["alpha"]) == selected)
    return selected, rows


def validate_contract_receipt(value: Mapping[str, Any]) -> list[dict[str, Any]]:
    if (
        value.get("schema_version") != CONTRACT_SCHEMA_VERSION
        or value.get("status") != "pass_outcome_free_crossfit_head_contract"
        or value.get("model_id") != "corgi_regular"
        or value.get("task_id") != "rna_conditioned_atac"
        or value.get("dataset_id") != "gse296875"
        or value.get("lineage_id") != "hepatocyte"
        or value.get("adaptation_rung") != "new_assay_head"
        or value.get("donor_and_genomic_fold_crossfit") is not True
        or value.get("source_signal_arrays_read") is not False
        or value.get("development_outcomes_read") is not False
        or value.get("test_or_sealed_features_or_outcomes_read") is not False
        or value.get("model_fit_performed") is not False
        or value.get("evaluation_performed") is not False
        or value.get("promotion_gate_evaluated") is not False
        or value.get("champion_or_external_claim_allowed") is not False
        or value.get("global_census_modified") is not False
    ):
        raise CorgiHeadFitError("head-contract receipt differs")
    rows = value.get("folds")
    if not isinstance(rows, list) or len(rows) != len(FOLDS):
        raise CorgiHeadFitError("head-contract fold count differs")
    normalized: list[dict[str, Any]] = []
    for row in rows:
        evaluation = int(row["evaluation_fold"])
        fit_valid = [int(item) for item in row["fit_valid_folds"]]
        expected_valid = [fold for fold in FOLDS if fold != evaluation]
        expected_base = [(fold - 1) % len(FOLDS) for fold in expected_valid]
        if (
            evaluation not in FOLDS
            or int(row["evaluation_base_outer_fold"]) != (evaluation - 1) % len(FOLDS)
            or fit_valid != expected_valid
            or [int(item) for item in row["fit_base_outer_folds"]] != expected_base
        ):
            raise CorgiHeadFitError("head-contract fold firewall differs")
        normalized.append(dict(row))
    if {int(row["evaluation_fold"]) for row in normalized} != set(FOLDS):
        raise CorgiHeadFitError("head-contract folds are incomplete")
    return sorted(normalized, key=lambda row: int(row["evaluation_fold"]))


def _load_base_fold(predictions: Path, outer_fold: int) -> dict[str, Any]:
    fold_root = predictions / f"fold{outer_fold}"
    for filename in (
        "receipt.json",
        "prediction_records.tsv",
        "window_records.tsv",
        "regional_predictions.npz",
    ):
        relative = f"fold{outer_fold}/{filename}"
        member = _artifact_member(predictions, relative)
        verify_selected_member(predictions, relative, str(member["sha256"]))
    receipt = json.loads((fold_root / "receipt.json").read_text(encoding="utf-8"))
    record_fields, records = read_tsv(fold_root / "prediction_records.tsv")
    window_fields, windows = read_tsv(fold_root / "window_records.tsv")
    if record_fields != PREDICTION_FIELDS:
        raise CorgiHeadFitError("prediction record fields differ")
    validate_prediction_tables(records, windows, receipt)
    with np.load(fold_root / "regional_predictions.npz", allow_pickle=False) as archive:
        if SCORE_NAME not in archive.files:
            raise CorgiHeadFitError("Corgi score is absent")
        scores = np.asarray(archive[SCORE_NAME], dtype=np.float64)
    if scores.shape != (len(records), len(windows)) or np.any(~np.isfinite(scores)) or np.any(scores < 0):
        raise CorgiHeadFitError("Corgi score matrix differs")
    by_contig: dict[str, list[int]] = {}
    for index, window in enumerate(windows):
        by_contig.setdefault(window["contig"], []).append(index)
    return {
        "receipt": receipt,
        "record_fields": record_fields,
        "records": records,
        "window_fields": window_fields,
        "windows": windows,
        "scores": scores,
        "by_contig": by_contig,
    }


def _fit_one_head(
    *,
    evaluation_fold: int,
    fit_base_outer_folds: Sequence[int],
    base_folds: Mapping[int, Mapping[str, Any]],
    donor_paths: Mapping[tuple[str, int], Path],
    training_paths: Mapping[int, Path],
) -> tuple[float, list[dict[str, float | int]], int, set[str]]:
    losses = {alpha: [] for alpha in ALPHA_GRID}
    fit_donor_ids: set[str] = set()
    for outer_fold in fit_base_outer_folds:
        fold = base_folds[outer_fold]
        valid_fold = int(fold["receipt"]["valid_donor_fold"])
        if valid_fold == evaluation_fold or int(fold["receipt"]["valid_genomic_fold"]) == evaluation_fold:
            raise CorgiHeadFitError("held donor or genomic fold entered head fit")
        baseline = np.asarray(_bigwig_sums(training_paths[outer_fold], fold["windows"]), dtype=np.float64)
        for row_index, record in enumerate(fold["records"]):
            if record["context_arm"] != PRIMARY_CONTEXT_ARM:
                continue
            donor = record["donor_id"]
            path = donor_paths.get((donor, valid_fold))
            if path is None:
                raise CorgiHeadFitError("fit donor ATAC is unavailable")
            observed = np.asarray(_bigwig_sums(path, fold["windows"]), dtype=np.float64)
            fit_donor_ids.add(donor)
            for indices in fold["by_contig"].values():
                if len(indices) < 3 or float(observed[indices].sum()) <= 0:
                    continue
                for alpha in ALPHA_GRID:
                    prediction = geometric_profile(
                        fold["scores"][row_index, indices], baseline[indices], alpha
                    )
                    losses[alpha].append(
                        multinomial_deviance_per_insertion(observed[indices], prediction)
                    )
    selected, rows = select_alpha(losses)
    return selected, rows, len(losses[ALPHA_GRID[0]]), fit_donor_ids


def _apply_one_head(
    *, fold: Mapping[str, Any], baseline_path: Path, selected_alpha: float
) -> np.ndarray:
    baseline = np.asarray(_bigwig_sums(baseline_path, fold["windows"]), dtype=np.float64)
    output = np.full(fold["scores"].shape, np.nan, dtype=np.float32)
    for row_index in range(len(fold["records"])):
        for indices in fold["by_contig"].values():
            output[row_index, indices] = geometric_profile(
                fold["scores"][row_index, indices], baseline[indices], selected_alpha
            ).astype(np.float32)
    if np.any(~np.isfinite(output)) or np.any(output <= 0):
        raise CorgiHeadFitError("held-fold profile prediction differs")
    return output


def fit(
    *,
    contract_preflight: Path,
    predictions: Path,
    donor_bigwigs: Path,
    training_bigwigs: Path,
    output: Path,
) -> Mapping[str, Any]:
    if output.exists():
        raise CorgiHeadFitError("refusing to overwrite profile-head fit")
    for root in (contract_preflight, predictions, donor_bigwigs, training_bigwigs):
        verify_complete_marker(root)
    contract_metadata = _artifact_metadata(contract_preflight)
    if (
        contract_metadata.get("artifact_class")
        != "corgi_regular_crossfit_profile_head_contract_preflight"
        or contract_metadata.get("development_outcomes_read") is not False
        or contract_metadata.get("evaluation_performed") is not False
        or contract_metadata.get("global_census_modified") is not False
    ):
        raise CorgiHeadFitError("contract preflight metadata differs")
    contract_member = _artifact_member(contract_preflight, "contract/contract_receipt.json")
    contract_path = verify_selected_member(
        contract_preflight, "contract/contract_receipt.json", str(contract_member["sha256"])
    )
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    folds = validate_contract_receipt(contract)
    prediction_metadata = _artifact_metadata(predictions)
    donor_metadata = _artifact_metadata(donor_bigwigs)
    training_metadata = _artifact_metadata(training_bigwigs)
    if (
        prediction_metadata.get("artifact_class")
        != "Corgi_outcome_aligned_tile_smoke_bundle_v3"
        or prediction_metadata.get("outcomes_read") is not False
        or prediction_metadata.get("outcome_role") != "valid"
        or donor_metadata.get("artifact_class") != "gse296875_deduplicated_tn5_bigwigs"
        or donor_metadata.get("biological_unit") != "donor"
        or training_metadata.get("artifact_class") != "gse296875_training_fold_pseudobulk"
        or training_metadata.get("biological_outer_unit") != "donor"
    ):
        raise CorgiHeadFitError("fit source metadata differs")
    donor_paths, training_paths = _load_outcome_paths(donor_bigwigs, training_bigwigs)
    base_folds = {outer: _load_base_fold(predictions, outer) for outer in FOLDS}
    output.mkdir(parents=True, exist_ok=False)
    fit_rows: list[dict[str, object]] = []
    fold_receipts: list[dict[str, object]] = []
    all_fit_donor_hashes: set[str] = set()
    for contract_fold in folds:
        evaluation_fold = int(contract_fold["evaluation_fold"])
        evaluation_base = int(contract_fold["evaluation_base_outer_fold"])
        fit_base = [int(value) for value in contract_fold["fit_base_outer_folds"]]
        selected, alpha_rows, training_units, fit_donors = _fit_one_head(
            evaluation_fold=evaluation_fold,
            fit_base_outer_folds=fit_base,
            base_folds=base_folds,
            donor_paths=donor_paths,
            training_paths=training_paths,
        )
        fold = base_folds[evaluation_base]
        if (
            int(fold["receipt"]["valid_donor_fold"]) != evaluation_fold
            or int(fold["receipt"]["valid_genomic_fold"]) != evaluation_fold
        ):
            raise CorgiHeadFitError("held prediction fold differs")
        held_donors = {row["donor_id"] for row in fold["records"]}
        if held_donors & fit_donors:
            raise CorgiHeadFitError("held donor entered profile-head fit")
        profile = _apply_one_head(
            fold=fold,
            baseline_path=training_paths[evaluation_base],
            selected_alpha=selected,
        )
        fold_output = output / f"fold{evaluation_fold}"
        fold_output.mkdir()
        np.savez_compressed(fold_output / "profile_predictions.npz", profile_probability=profile)
        _write_tsv(fold_output / "prediction_records.tsv", fold["record_fields"], fold["records"])
        _write_tsv(fold_output / "window_records.tsv", fold["window_fields"], fold["windows"])
        for row in alpha_rows:
            fit_rows.append({"evaluation_fold": evaluation_fold, **row})
        all_fit_donor_hashes.update(
            sha256(f"gse296875-corgi-head-fit-v1\0{donor}".encode()).hexdigest()
            for donor in fit_donors
        )
        fold_receipt = {
            "schema_version": SCHEMA_VERSION,
            "status": "pass_training_only_head_fit_and_held_prediction_export",
            "evaluation_fold": evaluation_fold,
            "evaluation_base_outer_fold": evaluation_base,
            "fit_valid_folds": contract_fold["fit_valid_folds"],
            "fit_base_outer_folds": fit_base,
            "selected_alpha": selected,
            "training_units": training_units,
            "training_donors": len(fit_donors),
            "held_donors": len(held_donors),
            "held_donor_overlap_with_fit": 0,
            "held_genomic_fold_overlap_with_fit": 0,
            "context_arms": list(CONTEXT_ARMS),
            "same_alpha_applied_to_all_context_arms": True,
            "held_fold_atac_signal_values_used_during_fit": False,
            "held_fold_atac_signal_values_used_during_export": False,
            "histology_or_disease_labels_read": False,
            "test_or_sealed_features_or_outcomes_read": False,
            "evaluation_metrics_computed": False,
            "profile_predictions_sha256": digest(fold_output / "profile_predictions.npz"),
        }
        (fold_output / "receipt.json").write_text(
            json.dumps(fold_receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        fold_receipts.append(fold_receipt)
    _write_tsv(
        output / "fit_loss_grid.tsv",
        (
            "evaluation_fold",
            "alpha",
            "units",
            "mean_fit_loss",
            "standard_error",
            "one_se_threshold",
            "selected",
        ),
        fit_rows,
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "status": "pass_training_only_crossfit_head_fit_and_prediction_export",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "contract_preflight_artifacts_sha256": digest(contract_preflight / "ARTIFACTS.json"),
        "base_predictions_artifacts_sha256": digest(predictions / "ARTIFACTS.json"),
        "donor_bigwigs_artifacts_sha256": digest(donor_bigwigs / "ARTIFACTS.json"),
        "training_bigwigs_artifacts_sha256": digest(training_bigwigs / "ARTIFACTS.json"),
        "folds": fold_receipts,
        "fit_donor_hashes": sorted(all_fit_donor_hashes),
        "development_atac_outcomes_read_for_training_only_fit": True,
        "held_fold_atac_signal_values_used_during_each_head_fit": False,
        "artifact_hash_verification_may_read_bound_file_bytes": True,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "benchmark_metrics_computed": False,
        "evaluation_performed": False,
        "promotion_gate_evaluated": False,
        "champion_or_external_claim_allowed": False,
        "global_census_modified": False,
    }
    (output / "fit_receipt.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract-preflight", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--donor-bigwigs", type=Path, required=True)
    parser.add_argument("--training-bigwigs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = fit(
        contract_preflight=args.contract_preflight.resolve(strict=True),
        predictions=args.predictions.resolve(strict=True),
        donor_bigwigs=args.donor_bigwigs.resolve(strict=True),
        training_bigwigs=args.training_bigwigs.resolve(strict=True),
        output=args.output.resolve(strict=False),
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
