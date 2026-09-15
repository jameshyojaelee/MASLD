#!/usr/bin/env python3
"""Independently evaluate frozen DNA-language OOF heads by model and head."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.stats import pearsonr, rankdata, spearmanr

from masld_bench.artifacts import verify_frozen_tree
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    CANDIDATES,
    CONTEXTS,
    MISSING_BASELINES,
    MODELS,
    PREDICTION_FIELDS,
    SEEDS,
    file_sha256,
    load_config,
)


SCHEMA = "masld-bench-gse281364-dna-language-seeded-head-evaluation-v1"
ROW_FIELDS = PREDICTION_FIELDS[:11]


class DNASeededEvaluationError(RuntimeError):
    """Raised when frozen identities or independent metrics differ."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DNASeededEvaluationError(f"invalid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise DNASeededEvaluationError(f"JSON object required: {path}")
    return value


def _read_tsv(path: Path, fields: Sequence[str], *, compressed: bool = False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else Path.open
    if compressed:
        handle_context = opener(path, "rt", encoding="utf-8", newline="")
    else:
        handle_context = opener(path, "r", encoding="utf-8", newline="")
    with handle_context as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise DNASeededEvaluationError(f"TSV schema differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise DNASeededEvaluationError(f"empty TSV: {path}")
    return rows


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _number(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise DNASeededEvaluationError("prediction is not numeric") from error
    if not math.isfinite(result):
        raise DNASeededEvaluationError("prediction is not finite")
    return result


def _display(value: float | None) -> str:
    return "not_applicable_constant_prediction" if value is None else format(value, ".17g")


def _macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = np.clip(np.asarray(values), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> dict[str, float | None]:
    residual = observed - predicted
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    values: dict[str, float | None] = {
        "spearman": None,
        "pearson": None,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "r2": None if denominator == 0 else float(1.0 - np.sum(residual**2) / denominator),
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    if not np.all(predicted == predicted[0]) and not np.all(observed == observed[0]):
        values["spearman"] = float(spearmanr(observed, predicted).statistic)
        values["pearson"] = float(pearsonr(observed, predicted).statistic)
        slope, intercept = np.polyfit(predicted, observed, 1)
        values["calibration_intercept"] = float(intercept)
        values["calibration_slope"] = float(slope)
    return values


def _fast_spearman(observed: np.ndarray, predicted: np.ndarray) -> float | None:
    if np.all(observed == observed[0]) or np.all(predicted == predicted[0]):
        return None
    left = rankdata(observed)
    right = rankdata(predicted)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def _bootstrap(
    observed: Mapping[str, np.ndarray],
    predictions: Mapping[tuple[str, str], Mapping[str, np.ndarray]],
    blocks: np.ndarray,
    *,
    resamples: int,
    seed: int,
) -> dict[tuple[str, str], np.ndarray]:
    unique = sorted(set(blocks.tolist()))
    by_block = {block: np.flatnonzero(blocks == block) for block in unique}
    output = {candidate: np.full(resamples, np.nan) for candidate in predictions}
    rng = np.random.default_rng(seed)
    for iteration in range(resamples):
        sampled = rng.integers(0, len(unique), size=len(unique))
        indices = np.concatenate([by_block[unique[index]] for index in sampled])
        for candidate, contexts in predictions.items():
            value = _macro([_fast_spearman(observed[context][indices], contexts[context][indices]) for context in CONTEXTS])
            if value is not None:
                output[candidate][iteration] = value
    return output


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config.resolve(strict=True))
    if arguments.output.exists():
        raise DNASeededEvaluationError("evaluation output exists")
    task_tree = arguments.task_spec_tree.resolve(strict=True)
    fit_tree = arguments.fit_tree.resolve(strict=True)
    if file_sha256(task_tree / "ARTIFACTS.json") != arguments.task_spec_sha256:
        raise DNASeededEvaluationError("TaskSpec identity differs")
    if file_sha256(fit_tree / "ARTIFACTS.json") != arguments.fit_sha256:
        raise DNASeededEvaluationError("fit identity differs")
    verify_frozen_tree(task_tree)
    verify_frozen_tree(fit_tree)
    fit = _load_json(fit_tree / "receipt.json")
    if (
        fit.get("status") != "pass_complete_dna_language_seeded_oof_fit"
        or fit.get("candidate_count") != len(CANDIDATES)
        or fit.get("prediction_rows") != len(CANDIDATES) * 10330
        or fit.get("task_spec_artifacts_sha256") != arguments.task_spec_sha256
        or fit.get("metrics_calculated")
        or fit.get("features_imputed")
        or fit.get("obsolete_6kb_folds_reused")
        or fit.get("models_ranked")
        or fit.get("shortlist_created")
        or fit.get("champion_claim")
        or fit.get("sealed_assets_read")
    ):
        raise DNASeededEvaluationError("fit firewall differs")
    row_binding = config["row_authority"]
    outcome_binding = config["outcome_authority"]
    row_tree = (root / row_binding["tree_path"]).resolve(strict=True)
    outcome_tree = (root / outcome_binding["tree_path"]).resolve(strict=True)
    if (
        file_sha256(row_tree / "ARTIFACTS.json") != row_binding["artifacts_sha256"]
        or file_sha256(outcome_tree / "ARTIFACTS.json") != outcome_binding["artifacts_sha256"]
        or file_sha256(row_tree / row_binding["member"]) != row_binding["member_sha256"]
        or file_sha256(outcome_tree / outcome_binding["member"]) != outcome_binding["member_sha256"]
    ):
        raise DNASeededEvaluationError("evaluation authority differs")
    verify_frozen_tree(row_tree)
    verify_frozen_tree(outcome_tree)
    row_rows = _read_tsv(row_tree / row_binding["member"], ROW_FIELDS)
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    if len(row_rows) != 10330 or len(row_lookup) != 10330:
        raise DNASeededEvaluationError("row identities differ")
    elements = sorted({row["element_id"] for row in row_rows})
    index = {element: position for position, element in enumerate(elements)}
    metadata: dict[str, str] = {}
    for row in row_rows:
        if metadata.setdefault(row["element_id"], row["long_range_block_id"]) != row["long_range_block_id"]:
            raise DNASeededEvaluationError("block metadata differs")
    blocks = np.asarray([metadata[element] for element in elements])
    if len(elements) != 1033 or len(set(blocks.tolist())) != 239:
        raise DNASeededEvaluationError("element/block denominator differs")
    outcome_data = load_outcomes(outcome_tree / outcome_binding["member"], set(elements))
    observed = {context: np.asarray([outcome_data[(element, context)]["mean"] for element in elements]) for context in CONTEXTS}

    prediction_rows = _read_tsv(fit_tree / "oof_predictions.tsv.gz", PREDICTION_FIELDS, compressed=True)
    values = {
        (model, head, seed, context): np.full(1033, np.nan)
        for model, head in CANDIDATES for seed in SEEDS for context in CONTEXTS
    }
    identities = set()
    for row in prediction_rows:
        candidate = (row["model_id"], row["head_id"])
        seed = int(row["seed"])
        identity = (*candidate, seed, row["row_hash"])
        authority = row_lookup.get((seed, row["row_hash"]))
        if candidate not in CANDIDATES or seed not in SEEDS or identity in identities or authority is None:
            raise DNASeededEvaluationError("prediction identity differs")
        identities.add(identity)
        if any(row[field] != authority[field] for field in ROW_FIELDS):
            raise DNASeededEvaluationError("prediction metadata differs")
        if row["experimental_replicates"] != "4" or row["biological_donors"] != "0" or row["outcome_role"] != "exposed_development_MPRA_only":
            raise DNASeededEvaluationError("prediction topology differs")
        slot = values[(*candidate, seed, row["assay_context_id"])]
        position = index[row["element_id"]]
        if math.isfinite(slot[position]):
            raise DNASeededEvaluationError("duplicate prediction")
        slot[position] = _number(row["prediction"])
    expected = len(CANDIDATES) * 10330
    if len(prediction_rows) != expected or len(identities) != expected or any(not np.isfinite(array).all() for array in values.values()):
        raise DNASeededEvaluationError("prediction coverage differs")

    model_info = {row["model_id"]: row for row in config["models"]}
    per_seed_rows = []
    per_seed_primary: dict[tuple[str, str, int], float | None] = {}
    ensemble: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    ensemble_rows = []
    points: dict[tuple[str, str], float | None] = {}
    for candidate in CANDIDATES:
        model, head = candidate
        family = "task_native_control" if model == "available_simple_controls" else "dna_language"
        feature_space = head if model == "available_simple_controls" else model_info[model]["family_native_input"] + "/" + head
        for seed in SEEDS:
            context_metrics = []
            for context in CONTEXTS:
                metrics = _metrics(observed[context], values[(model, head, seed, context)])
                context_metrics.append(metrics)
                per_seed_rows.append({"family": family, "feature_space": feature_space, "model_id": model, "head_id": head, "seed": seed, "assay_context_id": context, "elements": 1033, **{name: _display(value) for name, value in metrics.items()}})
            primary = _macro([item["spearman"] for item in context_metrics])
            per_seed_primary[(*candidate, seed)] = primary
            per_seed_rows.append({"family": family, "feature_space": feature_space, "model_id": model, "head_id": head, "seed": seed, "assay_context_id": "macro", "elements": 2066, "spearman": _display(primary)})
        by_context = {context: np.mean(np.vstack([values[(model, head, seed, context)] for seed in SEEDS]), axis=0) for context in CONTEXTS}
        ensemble[candidate] = by_context
        context_metrics = []
        for context in CONTEXTS:
            metrics = _metrics(observed[context], by_context[context])
            context_metrics.append(metrics)
            ensemble_rows.append({"family": family, "feature_space": feature_space, "model_id": model, "head_id": head, "assay_context_id": context, "elements": 1033, **{name: _display(value) for name, value in metrics.items()}})
        points[candidate] = _macro([item["spearman"] for item in context_metrics])
        ensemble_rows.append({"family": family, "feature_space": feature_space, "model_id": model, "head_id": head, "assay_context_id": "macro", "elements": 2066, "spearman": _display(points[candidate])})

    controls = [candidate for candidate in CANDIDATES if candidate[0] == "available_simple_controls" and points[candidate] is not None]
    reference = max(controls, key=lambda candidate: (float(points[candidate]), candidate[1]))
    samples = _bootstrap(observed, ensemble, blocks, resamples=int(config["uncertainty"]["resamples"]), seed=int(config["uncertainty"]["seed"]))
    reference_samples = samples[reference]
    summary_rows = []
    interval_rows = []
    for candidate in CANDIDATES:
        current = samples[candidate]
        valid = np.isfinite(current)
        if points[candidate] is None or int(valid.sum()) < 9500:
            low = high = standard_error = None
        else:
            low, high = (float(value) for value in np.quantile(current[valid], (0.025, 0.975)))
            standard_error = float(np.std(current[valid], ddof=1))
        gain = None if points[candidate] is None else float(points[candidate] - points[reference])
        gain_samples = current - reference_samples
        gain_valid = np.isfinite(gain_samples)
        if gain is None or int(gain_valid.sum()) < 9500:
            gain_low = gain_high = None
        else:
            gain_low, gain_high = (float(value) for value in np.quantile(gain_samples[gain_valid], (0.025, 0.975)))
        positive = sum(
            per_seed_primary[(*candidate, seed)] is not None
            and per_seed_primary[(*reference, seed)] is not None
            and per_seed_primary[(*candidate, seed)] > per_seed_primary[(*reference, seed)]
            for seed in SEEDS
        )
        model, head = candidate
        common = {
            "family": "task_native_control" if model == "available_simple_controls" else "dna_language",
            "model_id": model,
            "head_id": head,
            "family_native_input": "not_applicable_control" if model == "available_simple_controls" else model_info[model]["family_native_input"],
            "restricted": "false" if model == "available_simple_controls" else str(bool(model_info[model]["restricted"])).lower(),
            "primary_score": _display(points[candidate]),
            "primary_ci_low": _display(low),
            "primary_ci_high": _display(high),
            "primary_bootstrap_se": _display(standard_error),
            "available_control_reference": f"{reference[0]}/{reference[1]}",
            "gain_vs_available_control": _display(gain),
            "gain_ci_low": _display(gain_low),
            "gain_ci_high": _display(gain_high),
            "positive_gain_seeds": positive,
            "bootstrap_resamples": int(config["uncertainty"]["resamples"]),
            "bootstrap_unit": "long_range_block_id",
        }
        interval_rows.append(common)
        summary_rows.append({**common, "cross_family_rank": "not_reported", "mandatory_baselines_complete": "false", "shortlist": "false", "champion_claim": "false", "sealed_claim": "false", "confirmatory_p_value": "not_reported_development_screen"})
    arguments.output.mkdir(parents=True, mode=0o750)
    _write_tsv(arguments.output / "per_seed_metrics.tsv", per_seed_rows)
    _write_tsv(arguments.output / "ensemble_metrics.tsv", ensemble_rows)
    _write_tsv(arguments.output / "paired_block_bootstrap.tsv", interval_rows)
    _write_tsv(arguments.output / "model_head_summary.tsv", summary_rows)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_independent_dna_language_seeded_development_evaluation",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "models": list(MODELS),
        "candidate_count": len(CANDIDATES),
        "prediction_rows": expected,
        "elements": 1033,
        "long_range_blocks": 239,
        "fixed_seeds": list(SEEDS),
        "available_control_reference": {"model_id": reference[0], "head_id": reference[1], "primary_score": points[reference]},
        "bootstrap_resamples": int(config["uncertainty"]["resamples"]),
        "bootstrap_unit": "long_range_block_id",
        "models_reported_separately": True,
        "cross_family_ranking": False,
        "overall_leader_reported": False,
        "features_imputed": False,
        "obsolete_6kb_folds_reused": False,
        "mandatory_task_native_baselines_complete": False,
        "missing_baselines": list(MISSING_BASELINES),
        "shortlist_created": False,
        "champion_claim": False,
        "sealed_evaluation": False,
        "sealed_assets_read": False,
        "confirmatory_p_values_reported": False,
        "fit_implementation_imported": False,
    }
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--task-spec-tree", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--fit-tree", type=Path, required=True)
    parser.add_argument("--fit-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evaluate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
