#!/usr/bin/env python3
"""Freeze the primary-only GSE267145 RNA finalist development audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Mapping, Sequence

import numpy as np


PRIMARY_METRIC = "stage3_macro_f1"
CALIBRATION_METRIC = "stage3_multiclass_brier"
EXPECTED_SEEDS = (1701, 1709, 1721, 1723, 1733)
EXPECTED_OUTER_FOLDS = tuple(range(5))


class RNAFinalistSelectionError(RuntimeError):
    """Raised when the frozen development selection cannot be reproduced."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RNAFinalistSelectionError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _finite(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise RNAFinalistSelectionError("metric value is not numeric") from error
    if not math.isfinite(result):
        raise RNAFinalistSelectionError("metric value is not finite")
    return result


def select_rna_finalist(
    candidate_rows: Sequence[Mapping[str, Any]],
    candidate_model_ids: Sequence[str],
) -> dict[str, Any]:
    expected = tuple(candidate_model_ids)
    if len(expected) != 4 or len(set(expected)) != 4:
        raise RNAFinalistSelectionError("RNA candidate universe differs")
    by_model: dict[str, Mapping[str, Any]] = {}
    for row in candidate_rows:
        model_id = str(row.get("model_id", ""))
        if model_id in by_model:
            raise RNAFinalistSelectionError("RNA candidate metric is duplicated")
        by_model[model_id] = row
    if set(by_model) != set(expected):
        raise RNAFinalistSelectionError("RNA candidate metric universe differs")
    primary = {model: _finite(by_model[model][PRIMARY_METRIC]) for model in expected}
    best_value = max(primary.values())
    primary_winners = [model for model in expected if primary[model] == best_value]
    if len(primary_winners) != 1:
        raise RNAFinalistSelectionError("primary stage3 macro-F1 has no unique RNA winner")
    selected = primary_winners[0]
    calibration = {
        model: _finite(by_model[model][CALIBRATION_METRIC]) for model in expected
    }
    best_calibration = min(calibration.values())
    calibration_winners = [
        model for model in expected if calibration[model] == best_calibration
    ]
    if calibration_winners != [selected]:
        raise RNAFinalistSelectionError(
            "primary RNA winner is not also the unique best calibrated candidate"
        )
    return {
        "selected_model_id": selected,
        "primary_metric_id": PRIMARY_METRIC,
        "primary_estimate": primary[selected],
        "calibration_metric_id": CALIBRATION_METRIC,
        "calibration_estimate": calibration[selected],
        "unique_primary_winner": True,
        "unique_calibration_winner": True,
    }


def rederive_one_standard_error(
    audit: Mapping[str, Any],
) -> dict[str, Any]:
    candidates = audit.get("candidates")
    selected = audit.get("selected")
    if not isinstance(candidates, list) or not isinstance(selected, Mapping):
        raise RNAFinalistSelectionError("stage3 candidate audit schema differs")
    checked: list[dict[str, Any]] = []
    for candidate in candidates:
        if not isinstance(candidate, Mapping):
            raise RNAFinalistSelectionError("stage3 candidate record differs")
        scores = candidate.get("fold_scores")
        if (
            not isinstance(scores, list)
            or len(scores) != 4
            or candidate.get("failure_reason") is not None
            or candidate.get("outer_training_refit_valid") is not True
        ):
            continue
        values = [_finite(value) for value in scores]
        record = dict(candidate)
        record["mean_score"] = mean(values)
        record["standard_error"] = stdev(values) / math.sqrt(4)
        checked.append(record)
    if not checked:
        raise RNAFinalistSelectionError("stage3 candidate audit has no valid candidate")
    best = min(checked, key=lambda row: (-row["mean_score"], row["candidate_id"]))
    threshold = best["mean_score"] - best["standard_error"]
    eligible = [row for row in checked if row["mean_score"] >= threshold]
    rederived = min(
        eligible,
        key=lambda row: (tuple(row["complexity"]), row["candidate_id"]),
    )
    for key, value in (
        ("candidate_id", rederived["candidate_id"]),
        ("best_mean_score", best["mean_score"]),
        ("best_standard_error", best["standard_error"]),
        ("one_standard_error_threshold", threshold),
        ("one_standard_error_candidate_count", len(eligible)),
    ):
        observed = selected.get(key)
        if isinstance(value, float):
            if not math.isclose(_finite(observed), value, rel_tol=0.0, abs_tol=1e-12):
                raise RNAFinalistSelectionError(f"one-standard-error field differs: {key}")
        elif observed != value:
            raise RNAFinalistSelectionError(f"one-standard-error field differs: {key}")
    return {
        "selected_candidate_id": rederived["candidate_id"],
        "candidate_count": len(candidates),
        "valid_candidate_count": len(checked),
        "one_standard_error_candidate_count": len(eligible),
        "one_standard_error_threshold": threshold,
        "best_mean_score": best["mean_score"],
        "best_standard_error": best["standard_error"],
        "selected_mean_score": rederived["mean_score"],
    }


def _check_hash(path: Path, expected: str, label: str) -> None:
    if sha256_file(path) != expected:
        raise RNAFinalistSelectionError(f"{label} SHA-256 differs")


def audit_selection(
    *,
    benchmark_root: Path,
    contract_path: Path,
    contract_sha256: str,
    output: Path,
) -> dict[str, Any]:
    from masld_bench.artifacts import verify_frozen_tree

    if output.exists():
        raise RNAFinalistSelectionError(f"refusing to overwrite selection audit: {output}")
    _check_hash(contract_path, contract_sha256, "selection contract")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if (
        contract.get("status") != "registered_development_selection_audit_pending"
        or contract.get("selection_rule", {}).get("secondary_endpoint_use") != "forbidden"
        or contract.get("claim_boundary", {}).get("champion_claim_allowed") is not False
    ):
        raise RNAFinalistSelectionError("selection contract differs")
    sources = contract["source_artifacts"]
    production = benchmark_root / sources["production_root"]
    scoring = benchmark_root / sources["scoring_root"]
    for root, expected, label in (
        (production, sources["production_artifacts_sha256"], "production"),
        (scoring, sources["scoring_artifacts_sha256"], "scoring"),
    ):
        _check_hash(root / "ARTIFACTS.json", expected, label)
        verify_frozen_tree(root)
    _check_hash(
        production / "aggregate/ARTIFACTS.json",
        sources["aggregate_artifacts_sha256"],
        "production aggregate",
    )
    _check_hash(
        scoring / "scores/ARTIFACTS.json",
        sources["scores_artifacts_sha256"],
        "scores",
    )
    metrics_path = scoring / "scores/standardized_metrics.tsv"
    _check_hash(metrics_path, sources["standardized_metrics_sha256"], "metrics table")
    _check_hash(
        benchmark_root / sources["task_contract"],
        sources["task_contract_sha256"],
        "task contract",
    )
    surface_path = benchmark_root / sources["benchmark_surface"]
    _check_hash(surface_path, sources["benchmark_surface_sha256"], "benchmark surface")
    surface = json.loads(surface_path.read_text(encoding="utf-8"))
    candidate_models = tuple(contract["candidate_model_ids"])
    rna_baselines = surface["lanes"]["rna_only"]["baselines"]
    if set(candidate_models) != set(rna_baselines):
        raise RNAFinalistSelectionError("surface RNA candidate universe differs")
    if (
        surface["independent_evaluator"]["primary"] != [PRIMARY_METRIC]
        or surface["independent_evaluator"]["calibration"] != [CALIBRATION_METRIC]
        or surface["resampling"]["inner_hyperparameter_selection"]["rule"]
        != "one_standard_error"
        or rna_baselines["rna_hvg_pca_linear_svm"]["probability_calibration"]
        != "inner_oof_multinomial_logistic_calibrator"
    ):
        raise RNAFinalistSelectionError("registered primary, calibration, or one-SE rule differs")

    _, metric_rows = read_tsv(metrics_path)
    if len(metric_rows) != 495:
        raise RNAFinalistSelectionError("standardized metric census differs")
    decision_rows = [
        row
        for row in metric_rows
        if row["model_id"] in candidate_models
        and row["endpoint_id"] == "stage3"
        and row["metric_id"] in {PRIMARY_METRIC, CALIBRATION_METRIC}
        and row["applicability_state"] == "observed"
    ]
    if len(decision_rows) != 2 * len(candidate_models):
        raise RNAFinalistSelectionError("RNA decision metric census differs")
    by_model: dict[str, dict[str, Any]] = {
        model: {"model_id": model} for model in candidate_models
    }
    candidate_metric_rows: list[dict[str, Any]] = []
    metric_lookup: dict[tuple[str, str], Mapping[str, str]] = {}
    for row in decision_rows:
        key = (row["model_id"], row["metric_id"])
        if key in metric_lookup:
            raise RNAFinalistSelectionError("RNA decision metric is duplicated")
        metric_lookup[key] = row
        by_model[row["model_id"]][row["metric_id"]] = row["estimate"]
    selected = select_rna_finalist(list(by_model.values()), candidate_models)
    if selected["selected_model_id"] != contract["expected_finalist_if_audit_passes"]:
        raise RNAFinalistSelectionError("RNA finalist differs from registered expectation")
    for model in candidate_models:
        primary = metric_lookup[(model, PRIMARY_METRIC)]
        calibration = metric_lookup[(model, CALIBRATION_METRIC)]
        if (
            primary["endpoint_role"] != "primary"
            or primary["direction"] != "maximize"
            or calibration["endpoint_role"] != "calibration"
            or calibration["direction"] != "minimize"
            or primary["confirmatory_inference_allowed"] != "false"
            or calibration["confirmatory_inference_allowed"] != "false"
        ):
            raise RNAFinalistSelectionError("RNA decision metric contract differs")
        candidate_metric_rows.append(
            {
                "model_id": model,
                "stage3_macro_f1": primary["estimate"],
                "stage3_macro_f1_ci95_low": primary["ci95_low"],
                "stage3_macro_f1_ci95_high": primary["ci95_high"],
                "stage3_multiclass_brier": calibration["estimate"],
                "stage3_multiclass_brier_ci95_low": calibration["ci95_low"],
                "stage3_multiclass_brier_ci95_high": calibration["ci95_high"],
                "selected_rna_finalist": str(model == selected["selected_model_id"]).lower(),
            }
        )

    one_se_rows: list[dict[str, Any]] = []
    for outer_fold in EXPECTED_OUTER_FOLDS:
        for seed in EXPECTED_SEEDS:
            unit = production / f"units/outer_{outer_fold}/seed_{seed}"
            verify_frozen_tree(unit)
            audit_root = unit / f"fit/fit_receipts/outer_{outer_fold}/seed_{seed}"
            for model in candidate_models:
                audit_path = audit_root / f"candidate_audit--{model}--stage3.json"
                audit = json.loads(audit_path.read_text(encoding="utf-8"))
                if audit.get("model_id") != model or audit.get("endpoint_id") != "stage3":
                    raise RNAFinalistSelectionError("stage3 candidate audit identity differs")
                rederived = rederive_one_standard_error(audit)
                one_se_rows.append(
                    {
                        "outer_fold": outer_fold,
                        "model_seed": seed,
                        "model_id": model,
                        **rederived,
                        "rederived_match": "true",
                    }
                )
    if len(one_se_rows) != 100:
        raise RNAFinalistSelectionError("stage3 one-SE audit census differs")

    output.mkdir(parents=True)
    write_tsv(
        output / "rna_candidate_primary_calibration.tsv",
        tuple(candidate_metric_rows[0]),
        candidate_metric_rows,
    )
    write_tsv(
        output / "stage3_hyperparameter_one_se_audit.tsv",
        tuple(one_se_rows[0]),
        one_se_rows,
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-rna-finalist-selection-audit-v1",
        "status": "passed_frozen_development_rna_finalist_selection",
        "selection_id": contract["selection_id"],
        "selected_model_id": selected["selected_model_id"],
        "candidate_model_count": len(candidate_models),
        "candidate_model_ids": list(candidate_models),
        "primary_metric_id": PRIMARY_METRIC,
        "primary_estimate": selected["primary_estimate"],
        "primary_unique_best": True,
        "calibration_metric_id": CALIBRATION_METRIC,
        "calibration_estimate": selected["calibration_estimate"],
        "selected_model_unique_best_calibrated": True,
        "stage3_hyperparameter_audits_rederived": len(one_se_rows),
        "stage3_hyperparameter_one_standard_error_match": True,
        "stage3_calibration": "inner_oof_multinomial_logistic_calibrator",
        "decision_metric_rows": len(decision_rows),
        "secondary_endpoint_values_used_for_selection": False,
        "source_stage5_used_for_selection": False,
        "recorded_sex_used_for_selection": False,
        "external_features_or_labels_accessed": False,
        "unit_of_inference": "participant",
        "model_seeds_are_biological_replicates": False,
        "full_gse267145_fit_gate_open": True,
        "rna_finalist_role": "internal_development_selection_only",
        "champion_claim_allowed": False,
        "external_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "clinical_claim_allowed": False,
        "confirmatory_inference_allowed": False,
        "source_artifacts": sources,
    }
    with (output / "selection_receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-root", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = audit_selection(
        benchmark_root=arguments.benchmark_root,
        contract_path=arguments.contract,
        contract_sha256=arguments.contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
