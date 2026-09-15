#!/usr/bin/env python3
"""Independently score prespecified GSE281364 OOF sequence-head predictions."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import pearsonr, rankdata, spearmanr

from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.freeze_gse281364_open_sequence_taskspec import (
    BASELINES,
    CONTEXTS,
    HEADS,
    MISSING_BASELINES,
    MODELS,
    PREDICTION_FIELDS,
    SEEDS,
    file_sha256,
    load_config,
    validate_config,
)


SCHEMA = "masld-bench-gse281364-open-sequence-development-evaluation-v1"
ROW_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
)
CANDIDATES = tuple(
    [("available_simple_controls", head) for head in BASELINES]
    + [(model, head) for model in MODELS for head in HEADS]
)


class HeadEvaluationError(RuntimeError):
    """Raised when prediction identities or independent metrics differ."""


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HeadEvaluationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise HeadEvaluationError(f"{label} must be a JSON object")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    if not rows:
        raise HeadEvaluationError(f"empty TSV: {path}")
    return rows


def read_gzip_predictions(path: Path) -> list[dict[str, str]]:
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != PREDICTION_FIELDS:
                raise HeadEvaluationError("prediction header differs")
            return [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise HeadEvaluationError(f"cannot read predictions: {error}") from error


def finite_number(value: str, *, label: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise HeadEvaluationError(f"{label} is not numeric") from error
    if not math.isfinite(number):
        raise HeadEvaluationError(f"{label} is not finite")
    return number


def fisher_macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = [float(np.clip(value, -0.999999, 0.999999)) for value in values]
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def metric_values(observed: np.ndarray, prediction: np.ndarray) -> dict[str, float | None]:
    if observed.shape != prediction.shape or observed.ndim != 1 or observed.size < 3:
        raise HeadEvaluationError("metric vectors differ")
    residual = observed - prediction
    result: dict[str, float | None] = {
        "spearman": None,
        "pearson": None,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "r2": None,
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    if denominator > 0:
        result["r2"] = float(1.0 - np.sum(residual**2) / denominator)
    if not np.all(prediction == prediction[0]) and not np.all(observed == observed[0]):
        result["spearman"] = float(spearmanr(observed, prediction).statistic)
        result["pearson"] = float(pearsonr(observed, prediction).statistic)
        slope, intercept = np.polyfit(prediction, observed, 1)
        result["calibration_intercept"] = float(intercept)
        result["calibration_slope"] = float(slope)
    return result


def fast_spearman(observed: np.ndarray, prediction: np.ndarray) -> float | None:
    if observed.size < 3 or np.all(observed == observed[0]) or np.all(prediction == prediction[0]):
        return None
    left = rankdata(observed)
    right = rankdata(prediction)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def bootstrap_primary(
    *,
    observed: Mapping[str, np.ndarray],
    predictions: Mapping[tuple[str, str], Mapping[str, np.ndarray]],
    block_ids: np.ndarray,
    resamples: int,
    seed: int,
) -> dict[tuple[str, str], np.ndarray]:
    unique_blocks = sorted(set(block_ids.tolist()))
    by_block = {
        block: np.flatnonzero(block_ids == block) for block in unique_blocks
    }
    rng = np.random.default_rng(seed)
    output = {
        candidate: np.full(resamples, np.nan, dtype=np.float64)
        for candidate in predictions
    }
    for bootstrap_index in range(resamples):
        sampled = rng.integers(0, len(unique_blocks), size=len(unique_blocks))
        indices = np.concatenate([by_block[unique_blocks[index]] for index in sampled])
        for candidate, by_context in predictions.items():
            values = [
                fast_spearman(observed[context][indices], by_context[context][indices])
                for context in CONTEXTS
            ]
            macro = fisher_macro(values)
            if macro is not None:
                output[candidate][bootstrap_index] = macro
    return output


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise HeadEvaluationError(f"cannot write empty table: {path}")
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def display(value: float | None) -> str:
    return "not_applicable_constant_prediction" if value is None else format(value, ".17g")


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    config = load_config(arguments.config.resolve(strict=True))
    try:
        validate_config(config)
    except Exception as error:
        raise HeadEvaluationError(str(error)) from error
    if arguments.output.exists():
        raise HeadEvaluationError("evaluation output exists")
    fit_root = arguments.fit_root.resolve(strict=True)
    fit_receipt = load_json(fit_root / "receipt.json", label="fit receipt")
    if (
        fit_receipt.get("status") != "pass_prespecified_exposed_development_oof_fit"
        or fit_receipt.get("candidate_count") != 12
        or fit_receipt.get("prediction_rows") != 123960
        or fit_receipt.get("task_spec_artifacts_sha256") != arguments.task_spec_sha256
        or fit_receipt.get("metrics_calculated")
        or fit_receipt.get("experimental_replicates_used_as_independent_donors")
        or fit_receipt.get("sealed_assets_read")
        or fit_receipt.get("stack_fit")
        or fit_receipt.get("residual_correlation_calculated")
        or fit_receipt.get("conditional_model_built_or_fit")
    ):
        raise HeadEvaluationError("fit receipt or firewall differs")
    root = arguments.root.resolve(strict=True)
    row_root = (root / config["row_authority"]["tree_path"]).resolve(strict=True)
    outcome_root = (root / config["outcome_authority"]["tree_path"]).resolve(strict=True)
    row_path = row_root / config["row_authority"]["row_member"]
    outcome_path = outcome_root / config["outcome_authority"]["outcome_member"]
    if file_sha256(row_path) != config["row_authority"]["row_sha256"]:
        raise HeadEvaluationError("row-universe member differs")
    if file_sha256(outcome_path) != config["outcome_authority"]["outcome_sha256"]:
        raise HeadEvaluationError("outcome member differs")
    row_rows = read_tsv(row_path)
    if tuple(row_rows[0]) != ROW_FIELDS or len(row_rows) != 10330:
        raise HeadEvaluationError("row-universe schema differs")
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    if len(row_lookup) != 10330:
        raise HeadEvaluationError("row-universe identity differs")
    elements = sorted({row["element_id"] for row in row_rows})
    element_index = {element: index for index, element in enumerate(elements)}
    element_metadata = {}
    for row in row_rows:
        metadata = (
            row["source_locus_group_id"],
            row["long_range_block_id"],
            row["outer_fold"],
            row["unit_hash"],
            row["block_hash"],
        )
        if element_metadata.setdefault(row["element_id"], metadata) != metadata:
            raise HeadEvaluationError("element metadata differs")
    outcomes = load_outcomes(outcome_path, set(elements))
    observed = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements])
        for context in CONTEXTS
    }
    block_ids = np.asarray([element_metadata[element][1] for element in elements])
    if len(set(block_ids.tolist())) != 239:
        raise HeadEvaluationError("bootstrap block census differs")

    rows = read_gzip_predictions(fit_root / "oof_predictions.tsv.gz")
    if len(rows) != 123960:
        raise HeadEvaluationError("prediction row count differs")
    values: dict[tuple[str, str, int, str], np.ndarray] = {
        (model, head, seed, context): np.full(1033, np.nan, dtype=np.float64)
        for model, head in CANDIDATES
        for seed in SEEDS
        for context in CONTEXTS
    }
    identities = set()
    for index, row in enumerate(rows):
        seed = int(row["seed"])
        candidate = (row["model_id"], row["head_id"])
        identity = (*candidate, seed, row["row_hash"])
        if candidate not in CANDIDATES or seed not in SEEDS or identity in identities:
            raise HeadEvaluationError("prediction candidate, seed, or identity differs")
        identities.add(identity)
        authority = row_lookup.get((seed, row["row_hash"]))
        if authority is None or any(
            row[field] != authority[field]
            for field in (
                "unit_hash",
                "block_hash",
                "stratum",
                "outer_fold",
                "study_id",
                "assay_context_id",
                "element_id",
                "source_locus_group_id",
                "long_range_block_id",
            )
        ):
            raise HeadEvaluationError("prediction row metadata differs")
        if (
            row["experimental_replicates"] != "4"
            or row["biological_donors"] != "0"
            or row["outcome_role"] != "exposed_development_MPRA_only"
        ):
            raise HeadEvaluationError("prediction replication or role differs")
        context = row["assay_context_id"]
        element = element_index[row["element_id"]]
        expected_observed = float(observed[context][element])
        if abs(finite_number(row["observed"], label="observed") - expected_observed) > 1.0e-12:
            raise HeadEvaluationError("fit output observed value differs from independent aggregation")
        array = values[(candidate[0], candidate[1], seed, context)]
        if math.isfinite(array[element]):
            raise HeadEvaluationError("duplicate prediction element")
        array[element] = finite_number(row["prediction"], label="prediction")
    if len(identities) != 123960 or any(not np.isfinite(array).all() for array in values.values()):
        raise HeadEvaluationError("prediction coverage is incomplete")

    per_seed_rows = []
    per_seed_primary: dict[tuple[str, str, int], float | None] = {}
    ensemble_predictions: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    ensemble_rows = []
    point_primary: dict[tuple[str, str], float | None] = {}
    for candidate in CANDIDATES:
        model, head = candidate
        for seed in SEEDS:
            context_metrics = []
            for context in CONTEXTS:
                metric = metric_values(observed[context], values[(model, head, seed, context)])
                context_metrics.append(metric)
                per_seed_rows.append(
                    {
                        "model_id": model,
                        "head_id": head,
                        "seed": seed,
                        "assay_context_id": context,
                        "elements": 1033,
                        **{name: display(value) for name, value in metric.items()},
                    }
                )
            primary = fisher_macro([metric["spearman"] for metric in context_metrics])
            per_seed_primary[(model, head, seed)] = primary
            per_seed_rows.append(
                {
                    "model_id": model,
                    "head_id": head,
                    "seed": seed,
                    "assay_context_id": "macro",
                    "elements": 2066,
                    "spearman": display(primary),
                    "pearson": "not_applicable_macro",
                    "rmse": "not_applicable_macro",
                    "mae": "not_applicable_macro",
                    "r2": "not_applicable_macro",
                    "calibration_intercept": "not_applicable_macro",
                    "calibration_slope": "not_applicable_macro",
                }
            )
        by_context = {
            context: np.mean(
                np.vstack([values[(model, head, seed, context)] for seed in SEEDS]),
                axis=0,
            )
            for context in CONTEXTS
        }
        ensemble_predictions[candidate] = by_context
        context_metrics = []
        for context in CONTEXTS:
            metric = metric_values(observed[context], by_context[context])
            context_metrics.append(metric)
            ensemble_rows.append(
                {
                    "model_id": model,
                    "head_id": head,
                    "assay_context_id": context,
                    "elements": 1033,
                    **{name: display(value) for name, value in metric.items()},
                }
            )
        primary = fisher_macro([metric["spearman"] for metric in context_metrics])
        point_primary[candidate] = primary
        ensemble_rows.append(
            {
                "model_id": model,
                "head_id": head,
                "assay_context_id": "macro",
                "elements": 2066,
                "spearman": display(primary),
                "pearson": "not_applicable_macro",
                "rmse": "not_applicable_macro",
                "mae": "not_applicable_macro",
                "r2": "not_applicable_macro",
                "calibration_intercept": "not_applicable_macro",
                "calibration_slope": "not_applicable_macro",
            }
        )
    baseline_candidates = [
        candidate
        for candidate in CANDIDATES
        if candidate[0] == "available_simple_controls" and point_primary[candidate] is not None
    ]
    if not baseline_candidates:
        raise HeadEvaluationError("no nonconstant available simple control")
    strongest_baseline = max(
        baseline_candidates,
        key=lambda candidate: (float(point_primary[candidate]), candidate[1]),
    )
    bootstrap = bootstrap_primary(
        observed=observed,
        predictions=ensemble_predictions,
        block_ids=block_ids,
        resamples=int(config["uncertainty"]["resamples"]),
        seed=int(config["uncertainty"]["seed"]),
    )
    baseline_bootstrap = bootstrap[strongest_baseline]
    bootstrap_rows = []
    summaries = []
    intervals: dict[tuple[str, str], tuple[float, float, float, float, float]] = {}
    for candidate in CANDIDATES:
        samples = bootstrap[candidate]
        valid = np.isfinite(samples)
        point = point_primary[candidate]
        if point is None or int(valid.sum()) < 9500:
            low = high = standard_error = math.nan
        else:
            low, high = (float(value) for value in np.quantile(samples[valid], (0.025, 0.975)))
            standard_error = float(np.std(samples[valid], ddof=1))
        gain = None if point is None else float(point - float(point_primary[strongest_baseline]))
        gain_samples = samples - baseline_bootstrap
        gain_valid = np.isfinite(gain_samples)
        gain_low = gain_high = math.nan
        if gain is not None and int(gain_valid.sum()) >= 9500:
            gain_low, gain_high = (
                float(value) for value in np.quantile(gain_samples[gain_valid], (0.025, 0.975))
            )
        intervals[candidate] = (low, high, standard_error, gain_low, gain_high)
        bootstrap_rows.append(
            {
                "model_id": candidate[0],
                "head_id": candidate[1],
                "bootstrap_unit": "long_range_block_id",
                "bootstrap_resamples": 10000,
                "valid_primary_resamples": int(valid.sum()),
                "primary_ci_low": display(None if math.isnan(low) else low),
                "primary_ci_high": display(None if math.isnan(high) else high),
                "primary_bootstrap_se": display(None if math.isnan(standard_error) else standard_error),
                "gain_vs_strongest_available_control": display(gain),
                "gain_ci_low": display(None if math.isnan(gain_low) else gain_low),
                "gain_ci_high": display(None if math.isnan(gain_high) else gain_high),
            }
        )
    finite_candidates = [candidate for candidate in CANDIDATES if point_primary[candidate] is not None]
    leader = max(finite_candidates, key=lambda candidate: (float(point_primary[candidate]), candidate))
    leader_se = intervals[leader][2]
    if not math.isfinite(leader_se):
        raise HeadEvaluationError("leader bootstrap standard error is unavailable")
    one_se_threshold = float(point_primary[leader]) - leader_se
    one_se_rows = []
    for candidate in CANDIDATES:
        point = point_primary[candidate]
        positive_seeds = 0
        if point is not None:
            for seed in SEEDS:
                candidate_seed = per_seed_primary[(*candidate, seed)]
                baseline_seed = per_seed_primary[(*strongest_baseline, seed)]
                if (
                    candidate_seed is not None
                    and baseline_seed is not None
                    and candidate_seed > baseline_seed
                ):
                    positive_seeds += 1
        within_one_se = point is not None and point >= one_se_threshold
        control = candidate[0] == "available_simple_controls"
        partial_one_se = within_one_se and (control or positive_seeds >= 4)
        low, high, standard_error, gain_low, gain_high = intervals[candidate]
        summary = {
            "model_id": candidate[0],
            "head_id": candidate[1],
            "primary_metric": "tanh_mean_Fisher_z_Spearman_across_two_assay_contexts",
            "primary_score": display(point),
            "primary_ci_low": display(None if math.isnan(low) else low),
            "primary_ci_high": display(None if math.isnan(high) else high),
            "primary_bootstrap_se": display(None if math.isnan(standard_error) else standard_error),
            "strongest_available_control_model_id": strongest_baseline[0],
            "strongest_available_control_head_id": strongest_baseline[1],
            "gain_vs_strongest_available_control": display(
                None if point is None else float(point - float(point_primary[strongest_baseline]))
            ),
            "gain_ci_low": display(None if math.isnan(gain_low) else gain_low),
            "gain_ci_high": display(None if math.isnan(gain_high) else gain_high),
            "positive_gain_seeds": positive_seeds,
            "within_one_standard_error_of_partial_leader": str(within_one_se).lower(),
            "partial_campaign_one_se_screen": str(partial_one_se).lower(),
            "development_shortlist": "false",
            "shortlist_blocked": "true",
            "finalist_claim": "false",
            "complementarity_eligible": "false",
            "champion_claim": "false",
            "confirmatory_p_value": "not_reported_development_screen",
            "multiplicity": "none_no_confirmatory_family",
        }
        summaries.append(summary)
        if partial_one_se:
            one_se_rows.append(summary)
    summaries.sort(
        key=lambda row: (
            -float(row["primary_score"])
            if row["primary_score"] != "not_applicable_constant_prediction"
            else math.inf,
            row["model_id"],
            row["head_id"],
        )
    )
    if not one_se_rows:
        raise HeadEvaluationError("partial-campaign one-standard-error screen is empty")
    arguments.output.mkdir(parents=True, mode=0o750)
    write_tsv(arguments.output / "per_seed_metrics.tsv", per_seed_rows)
    write_tsv(arguments.output / "ensemble_metrics.tsv", ensemble_rows)
    write_tsv(arguments.output / "bootstrap_intervals.tsv", bootstrap_rows)
    write_tsv(arguments.output / "candidate_summary.tsv", summaries)
    write_tsv(arguments.output / "partial_campaign_one_se_set.tsv", one_se_rows)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_independent_exposed_development_evaluation",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "candidate_count": 12,
        "prediction_rows": 123960,
        "selected_elements": 1033,
        "contexts": list(CONTEXTS),
        "fixed_seeds": list(SEEDS),
        "outer_folds": 5,
        "long_range_blocks": 239,
        "bootstrap_resamples": 10000,
        "bootstrap_unit": "long_range_block_id",
        "strongest_available_simple_control": {
            "model_id": strongest_baseline[0],
            "head_id": strongest_baseline[1],
            "primary_score": point_primary[strongest_baseline],
        },
        "partial_campaign_leader": {
            "model_id": leader[0],
            "head_id": leader[1],
            "primary_score": point_primary[leader],
        },
        "one_standard_error_threshold": one_se_threshold,
        "partial_campaign_one_se_count": len(one_se_rows),
        "development_shortlist_count": 0,
        "campaign_scope": "partial_common_head_campaign",
        "mandatory_baselines_complete": False,
        "mandatory_baselines_missing": list(MISSING_BASELINES),
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
        "conditional_trigger_blocked": True,
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "outcome_role": "exposed_development_MPRA_only",
        "confirmatory_family": "none_exposed_development_screen",
        "confirmatory_p_values_reported": False,
        "champion_claim": False,
        "external_evaluation": False,
        "sealed_assets_read": False,
        "stack_fit": False,
        "residual_correlation_calculated": False,
        "conditional_model_built_or_fit": False,
        "fit_implementation_imported": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--fit-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    evaluate(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
