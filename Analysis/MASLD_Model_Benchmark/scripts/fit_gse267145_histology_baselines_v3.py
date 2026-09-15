#!/usr/bin/env python3
"""Run the v2 fitter with the prespecified v3 fibrosis-group fallback."""

from __future__ import annotations

from collections import Counter
import json
import math
from pathlib import Path
from statistics import mean, stdev
import sys
from typing import Any, Callable, Mapping, Sequence

import numpy as np

import scripts.fit_gse267145_histology_baselines as v2


REVISION_ID = "model-training-069-production-v3"
FALLBACK_ID = "training_only_fibrosis_group3_class_prior_v1"
CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "config/evaluation/gse267145_histology_production_v3.json"
)
CONTRACT_SHA256 = "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f"
FROZEN_C_VALUES = (0.01, 0.1, 1.0, 10.0)
FROZEN_L1_RATIOS = (0.0, 0.5, 1.0)
CONVERGENCE_FAILURE = (
    "HistologyBaselineError:elastic-net logistic candidate did not converge"
)


def _valid_candidates(
    candidates: Sequence[Mapping[str, Any]],
    outputs: Mapping[str, tuple[np.ndarray, np.ndarray]],
) -> list[Mapping[str, Any]]:
    valid: list[Mapping[str, Any]] = []
    for candidate in candidates:
        scores = np.asarray(candidate.get("fold_scores", []), dtype=float)
        if (
            scores.shape == (4,)
            and np.all(np.isfinite(scores))
            and candidate.get("failure_reason") is None
            and candidate.get("outer_training_refit_valid") is True
            and candidate.get("candidate_id") in outputs
        ):
            valid.append(candidate)
    return valid


def _training_only_prior_fallback(
    *,
    candidates: Sequence[Mapping[str, Any]],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    groups: np.ndarray,
) -> dict[str, Any]:
    if len(candidates) != len(FROZEN_C_VALUES) * len(FROZEN_L1_RATIOS):
        raise v2.HistologyBaselineError(
            "fibrosis-group fallback requires the exact 12-candidate grid"
        )
    observed_grid = [
        (
            float(candidate.get("parameters", {}).get("c")),
            float(candidate.get("parameters", {}).get("l1_ratio")),
        )
        for candidate in candidates
    ]
    expected_grid = [
        (c_value, l1_ratio)
        for c_value in FROZEN_C_VALUES
        for l1_ratio in FROZEN_L1_RATIOS
    ]
    if observed_grid != expected_grid:
        raise v2.HistologyBaselineError(
            "fibrosis-group fallback candidate grid differs"
        )
    failures = Counter(str(candidate.get("failure_reason")) for candidate in candidates)
    if failures != Counter({CONVERGENCE_FAILURE: len(candidates)}):
        raise v2.HistologyBaselineError(
            "fibrosis-group fallback requires convergence-only failures"
        )

    oof = np.full((len(outer_training_indices), 3), np.nan, dtype=np.float64)
    inner_counts: dict[str, list[int]] = {}
    fold_scores: list[float] = []
    for inner_fold in range(4):
        fit_local = np.flatnonzero(inner_assignment != inner_fold)
        val_local = np.flatnonzero(inner_assignment == inner_fold)
        counts = np.bincount(groups[fit_local], minlength=3).astype(np.float64)
        if counts.shape != (3,) or np.any(counts <= 0) or counts.sum() <= 0:
            raise v2.HistologyBaselineError(
                "fibrosis-group fallback inner training lacks a class"
            )
        probability = counts / counts.sum()
        oof[val_local] = probability
        inner_counts[str(inner_fold)] = [int(value) for value in counts]
        fold_scores.append(
            v2.macro_f1(
                groups[val_local],
                np.argmax(oof[val_local], axis=1),
                range(3),
            )
        )
    if not np.all(np.isfinite(oof)):
        raise v2.HistologyBaselineError(
            "fibrosis-group fallback OOF probabilities are incomplete"
        )
    outer_counts = np.bincount(groups, minlength=3).astype(np.float64)
    if outer_counts.shape != (3,) or np.any(outer_counts <= 0):
        raise v2.HistologyBaselineError(
            "fibrosis-group fallback outer training lacks a class"
        )
    outer_probability = outer_counts / outer_counts.sum()
    test = np.repeat(outer_probability[None, :], len(outer_test_indices), axis=0)
    standard_error = stdev(fold_scores) / math.sqrt(4)
    selected = {
        "candidate_id": v2._candidate_id({"fallback_id": FALLBACK_ID}),
        "parameters": {"fallback_id": FALLBACK_ID},
        "fold_scores": fold_scores,
        "mean_score": mean(fold_scores),
        "standard_error": standard_error,
        "complexity": ["deterministic_training_only_class_prior"],
        "failure_reason": None,
        "outer_training_refit_valid": True,
        "candidate_count": len(candidates),
        "valid_candidate_count": 0,
        "candidate_failure_count": len(candidates),
        "candidate_failure_reasons": {CONVERGENCE_FAILURE: len(candidates)},
        "fallback_triggered": True,
        "fallback_id": FALLBACK_ID,
        "fallback_scope": "secondary_fibrosis_group3_only",
        "fallback_activation_condition": (
            "zero valid candidates after the exact 12-candidate frozen SAGA grid; "
            "all failures are convergence failures"
        ),
        "inner_training_class_counts": inner_counts,
        "outer_training_class_counts": [int(value) for value in outer_counts],
        "outer_training_class_probabilities": [
            float(value) for value in outer_probability
        ],
        "one_standard_error_applied": False,
        "one_standard_error_not_applicable_no_fallback_hyperparameters": True,
        "outer_test_features_used": False,
        "outer_test_outcomes_used": False,
        "recorded_sex_used": False,
        "source_stage5_used": False,
        "randomness_used": False,
        "seed_lineage": {
            "fallback": "none_deterministic_partition_class_counts",
            "selected_predictions_reuse_exact_training_only_fallback": True,
        },
    }
    return {"oof": oof, "test": test, "selected": selected}


def tune_group_classifier_v3(
    *,
    cache: v2.RepresentationCache,
    representation: Mapping[str, Any],
    outer_training_indices: np.ndarray,
    outer_test_indices: np.ndarray,
    inner_assignment: np.ndarray,
    groups: np.ndarray,
    c_values: Sequence[float],
    l1_ratios: Sequence[float],
    seed: int,
    audit_callback: Callable[[str, Mapping[str, Any]], None] | None = None,
) -> dict[str, Any]:
    """Keep the v2 grid and use priors only after exact convergence-only exhaustion."""

    candidates: list[dict[str, Any]] = []
    candidate_outputs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for c_value in c_values:
        for l1_ratio in l1_ratios:
            scores: list[float] = []
            nonzero: list[int] = []
            failure_reason: str | None = None
            oof = np.full((len(outer_training_indices), 3), np.nan, dtype=np.float64)
            for inner_fold in range(4):
                fit_local = np.flatnonzero(inner_assignment != inner_fold)
                val_local = np.flatnonzero(inner_assignment == inner_fold)
                if set(int(value) for value in groups[fit_local]) != {0, 1, 2}:
                    raise v2.HistologyBaselineError(
                        "inner fibrosis-group training lacks a class"
                    )
                try:
                    rep = cache.get(
                        f"inner{inner_fold}",
                        outer_training_indices[fit_local],
                        outer_training_indices[val_local],
                        representation["feature_request"],
                    )
                    components = min(
                        representation["pca_components"], rep["max_components"]
                    )
                    model = v2._logistic(
                        rep["fit"][:, :components],
                        groups[fit_local],
                        c_value=float(c_value),
                        l1_ratio=float(l1_ratio),
                        seed=seed + inner_fold,
                    )
                    predicted = model.predict(rep["evaluation"][:, :components])
                    oof[val_local] = model.predict_proba(
                        rep["evaluation"][:, :components]
                    )
                    scores.append(
                        v2.macro_f1(groups[val_local], predicted, range(3))
                    )
                    nonzero.append(int(np.count_nonzero(model.coef_)))
                except (v2.HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
                    break
            parameters = {"c": float(c_value), "l1_ratio": float(l1_ratio)}
            candidate_id = v2._candidate_id(parameters)
            test: np.ndarray | None = None
            if failure_reason is None:
                try:
                    final_rep = cache.get(
                        "outer_final",
                        outer_training_indices,
                        outer_test_indices,
                        representation["feature_request"],
                    )
                    components = min(
                        representation["pca_components"],
                        final_rep["max_components"],
                    )
                    final = v2._logistic(
                        final_rep["fit"][:, :components],
                        groups,
                        c_value=float(c_value),
                        l1_ratio=float(l1_ratio),
                        seed=seed + 1000,
                    )
                    test = final.predict_proba(
                        final_rep["evaluation"][:, :components]
                    )
                except (v2.HistologyBaselineError, ValueError) as error:
                    failure_reason = f"{type(error).__name__}:{error}"
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "parameters": parameters,
                    "fold_scores": scores if failure_reason is None else [],
                    "failure_reason": failure_reason,
                    "outer_training_refit_valid": failure_reason is None,
                    "complexity": [
                        representation["pca_components"],
                        mean(nonzero) if nonzero else 10**9,
                        float(c_value),
                        float(l1_ratio),
                    ],
                }
            )
            if failure_reason is None and test is not None:
                candidate_outputs[candidate_id] = (oof, np.asarray(test))

    valid = _valid_candidates(candidates, candidate_outputs)
    if valid:
        selected = v2.select_with_failure_audit(candidates, maximize=True)
        selected["seed_lineage"] = {
            "candidate_inner": "base_seed + inner_fold",
            "candidate_outer_training_refit": "base_seed + 1000",
            "selected_predictions_reuse_exact_candidate_fits": True,
        }
        selected["fallback_triggered"] = False
        selected["fallback_id"] = FALLBACK_ID
        result = {
            "oof": candidate_outputs[selected["candidate_id"]][0],
            "test": candidate_outputs[selected["candidate_id"]][1],
            "selected": selected,
        }
    else:
        result = _training_only_prior_fallback(
            candidates=candidates,
            outer_training_indices=outer_training_indices,
            outer_test_indices=outer_test_indices,
            inner_assignment=inner_assignment,
            groups=groups,
        )
    if audit_callback is not None:
        audit_callback(
            "fibrosis_group3",
            {
                "selected": result["selected"],
                "candidates": v2._candidate_audit(candidates),
                "fallback_contract": {
                    "revision_id": REVISION_ID,
                    "fallback_id": FALLBACK_ID,
                    "eligible_endpoint": "fibrosis_group3",
                    "candidate_grid_count": 12,
                    "activation_requires_convergence_only_failures": True,
                    "outer_test_outcomes_used": False,
                },
            },
        )
    return result


_selection_receipt_v2 = v2._selection_receipt
_fit_preflight_v2 = v2.fit_preflight


def _selection_receipt_v3(
    stage: Mapping[str, Any], secondary: Mapping[str, Any]
) -> dict[str, Any]:
    receipt = _selection_receipt_v2(stage, secondary)
    group = receipt["fibrosis_group3"]
    triggered = group.get("fallback_triggered") is True
    receipt["one_standard_error_applied_to_every_selection"] = not triggered
    receipt["one_standard_error_applied_to_every_valid_frozen_grid_selection"] = True
    receipt["fibrosis_group3_fallback_policy"] = FALLBACK_ID
    receipt["fibrosis_group3_fallback_triggered"] = triggered
    receipt["primary_endpoint_selection_changed_by_v3"] = False
    receipt["other_secondary_endpoint_selection_changed_by_v3"] = False
    return receipt


def _fallback_summary(output: Path, outer_fold: int, seed: int) -> dict[str, Any]:
    selection_path = (
        output
        / "fit_receipts"
        / f"outer_{outer_fold}"
        / f"seed_{seed}"
        / "selection_receipts.json"
    )
    selections = json.loads(selection_path.read_text(encoding="utf-8"))
    triggered = sorted(
        model_id
        for model_id, receipt in selections.items()
        if isinstance(receipt, dict)
        and receipt.get("fibrosis_group3_fallback_triggered") is True
    )
    return {
        "fallback_id": FALLBACK_ID,
        "activation_count": len(triggered),
        "activated_model_ids": triggered,
        "eligible_endpoint": "fibrosis_group3",
        "eligible_endpoint_role": "secondary_only",
        "outer_test_outcomes_used": False,
        "metrics_calculated": False,
    }


def fit_preflight_v3(
    *,
    molecular: Path,
    folds: Path,
    fit_views: Path,
    surface_path: Path,
    output: Path,
    outer_fold: int,
    seed: int,
) -> dict[str, Any]:
    receipt = _fit_preflight_v2(
        molecular=molecular,
        folds=folds,
        fit_views=fit_views,
        surface_path=surface_path,
        output=output,
        outer_fold=outer_fold,
        seed=seed,
    )
    receipt["schema_version"] = (
        "masld-bench-gse267145-histology-baseline-fit-preflight-v3"
    )
    receipt["revision_id"] = REVISION_ID
    receipt["fibrosis_group3_fallback"] = _fallback_summary(
        output, outer_fold, seed
    )
    receipt["primary_endpoint_code_path_changed"] = False
    receipt["other_secondary_endpoint_code_paths_changed"] = False
    receipt["outer_test_outcomes_read"] = False
    receipt["outer_test_metrics_calculated"] = False
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def _verify_contract() -> None:
    if (
        not CONTRACT_PATH.is_file()
        or CONTRACT_PATH.is_symlink()
        or v2.sha256_file(CONTRACT_PATH) != CONTRACT_SHA256
    ):
        raise v2.HistologyBaselineError("v3 fallback contract SHA-256 differs")
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    fallback = contract.get("fallback_contract", {})
    if (
        contract.get("revision_id") != REVISION_ID
        or fallback.get("fallback_id") != FALLBACK_ID
        or fallback.get("eligible_endpoint") != "fibrosis_group3"
        or contract.get("scoring_authorized") is not False
    ):
        raise v2.HistologyBaselineError("v3 fallback contract differs")


def main() -> int:
    _verify_contract()
    if "--preflight-outer-fold" not in sys.argv or "--preflight-seed" not in sys.argv:
        raise v2.HistologyBaselineError(
            "v3 fitter must run as one frozen outer-fold/seed unit"
        )
    v2.tune_group_classifier = tune_group_classifier_v3
    v2._selection_receipt = _selection_receipt_v3
    v2.fit_preflight = fit_preflight_v3
    return v2.main()


if __name__ == "__main__":
    raise SystemExit(main())
