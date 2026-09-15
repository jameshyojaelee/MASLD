#!/usr/bin/env python3
"""Independently evaluate Enformer/Sei development OOF predictions by family."""

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

from masld_bench.artifacts import verify_frozen_tree
from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.freeze_gse281364_enformer_sei_head_taskspec import (
    CANDIDATES,
    CONTEXTS,
    MISSING_BASELINES,
    PREDICTION_FIELDS,
    SEEDS,
    file_sha256,
    load_config,
    validate_config,
)
from scripts.fit_gse281364_enformer_sei_heads import ROW_FIELDS


SCHEMA = "masld-bench-gse281364-enformer-sei-head-evaluation-v1"


class EnformerSeiEvaluationError(RuntimeError):
    """Raised when predictions, identities, or independent metrics differ."""


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise EnformerSeiEvaluationError(f"invalid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise EnformerSeiEvaluationError(f"JSON object required: {path}")
    return value


def _read_rows(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise EnformerSeiEvaluationError(f"gzip TSV schema differs: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise EnformerSeiEvaluationError(f"empty gzip TSV: {path}")
    return rows


def _write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise EnformerSeiEvaluationError(f"cannot write empty TSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _finite(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise EnformerSeiEvaluationError("prediction is not numeric") from error
    if not math.isfinite(result):
        raise EnformerSeiEvaluationError("prediction is not finite")
    return result


def _fisher_macro(values: Sequence[float | None]) -> float | None:
    if any(value is None or not math.isfinite(value) for value in values):
        return None
    clipped = np.clip(np.asarray(values, dtype=np.float64), -0.999999, 0.999999)
    return float(np.tanh(np.mean(np.arctanh(clipped))))


def _metrics(observed: np.ndarray, prediction: np.ndarray) -> dict[str, float | None]:
    residual = observed - prediction
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    result: dict[str, float | None] = {
        "spearman": None,
        "pearson": None,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "mae": float(np.mean(np.abs(residual))),
        "r2": None if denominator == 0 else float(1.0 - np.sum(residual**2) / denominator),
        "calibration_intercept": None,
        "calibration_slope": None,
    }
    if not np.all(prediction == prediction[0]) and not np.all(observed == observed[0]):
        result["spearman"] = float(spearmanr(observed, prediction).statistic)
        result["pearson"] = float(pearsonr(observed, prediction).statistic)
        slope, intercept = np.polyfit(prediction, observed, 1)
        result["calibration_intercept"] = float(intercept)
        result["calibration_slope"] = float(slope)
    return result


def _fast_spearman(observed: np.ndarray, prediction: np.ndarray) -> float | None:
    if observed.size < 3 or np.all(observed == observed[0]) or np.all(prediction == prediction[0]):
        return None
    left = rankdata(observed)
    right = rankdata(prediction)
    left -= left.mean()
    right -= right.mean()
    denominator = float(np.sqrt(np.sum(left**2) * np.sum(right**2)))
    return None if denominator == 0 else float(np.sum(left * right) / denominator)


def _display(value: float | None) -> str:
    return "not_applicable_constant_prediction" if value is None else format(value, ".17g")


def _bootstrap(
    observed: Mapping[str, np.ndarray],
    predictions: Mapping[tuple[str, str], Mapping[str, np.ndarray]],
    block_ids: np.ndarray,
    *,
    resamples: int,
    seed: int,
) -> dict[tuple[str, str], np.ndarray]:
    blocks = sorted(set(block_ids.tolist()))
    by_block = {block: np.flatnonzero(block_ids == block) for block in blocks}
    rng = np.random.default_rng(seed)
    output = {candidate: np.full(resamples, np.nan) for candidate in predictions}
    for iteration in range(resamples):
        sampled = rng.integers(0, len(blocks), size=len(blocks))
        indices = np.concatenate([by_block[blocks[index]] for index in sampled])
        for candidate, context_values in predictions.items():
            macro = _fisher_macro(
                [
                    _fast_spearman(observed[context][indices], context_values[context][indices])
                    for context in CONTEXTS
                ]
            )
            if macro is not None:
                output[candidate][iteration] = macro
    return output


def evaluate(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config.resolve(strict=True))
    validate_config(config)
    if arguments.output.exists():
        raise EnformerSeiEvaluationError("evaluation output already exists")
    task_tree = arguments.task_spec_tree.resolve(strict=True)
    if file_sha256(task_tree / "ARTIFACTS.json") != arguments.task_spec_sha256:
        raise EnformerSeiEvaluationError("TaskSpec ARTIFACTS identity differs")
    verify_frozen_tree(task_tree)
    fit_root = arguments.fit_root.resolve(strict=True)
    if file_sha256(fit_root / "ARTIFACTS.json") != arguments.fit_sha256:
        raise EnformerSeiEvaluationError("fit ARTIFACTS identity differs")
    verify_frozen_tree(fit_root)
    fit_receipt = _load_json(fit_root / "receipt.json")
    if (
        fit_receipt.get("status") != "pass_restricted_development_oof_fit"
        or fit_receipt.get("candidate_count") != len(CANDIDATES)
        or fit_receipt.get("prediction_rows") != len(CANDIDATES) * 10330
        or fit_receipt.get("task_spec_artifacts_sha256") != arguments.task_spec_sha256
        or fit_receipt.get("metrics_calculated")
        or fit_receipt.get("family_native_feature_spaces_combined")
        or fit_receipt.get("sealed_assets_read")
        or fit_receipt.get("champion_claim")
    ):
        raise EnformerSeiEvaluationError("fit receipt firewall differs")
    static = config["static_authority"]
    static_tree = (root / static["tree_path"]).resolve(strict=True)
    outcome = config["outcome_authority"]
    outcome_tree = (root / outcome["tree_path"]).resolve(strict=True)
    if (
        file_sha256(static_tree / "ARTIFACTS.json") != static["artifacts_sha256"]
        or file_sha256(outcome_tree / "ARTIFACTS.json") != outcome["artifacts_sha256"]
        or file_sha256(static_tree / static["row_member"]) != static["row_sha256"]
        or file_sha256(outcome_tree / outcome["member"]) != outcome["member_sha256"]
    ):
        raise EnformerSeiEvaluationError("evaluation authority differs")
    verify_frozen_tree(static_tree)
    verify_frozen_tree(outcome_tree)
    row_rows = _read_rows(static_tree / static["row_member"], ROW_FIELDS)
    if len(row_rows) != 10330:
        raise EnformerSeiEvaluationError("row denominator differs")
    row_lookup = {(int(row["seed"]), row["row_hash"]): row for row in row_rows}
    if len(row_lookup) != 10330:
        raise EnformerSeiEvaluationError("row identities differ")
    elements = sorted({row["element_id"] for row in row_rows})
    element_index = {element: index for index, element in enumerate(elements)}
    metadata: dict[str, tuple[str, str]] = {}
    for row in row_rows:
        pair = (row["long_range_block_id"], row["outer_fold"])
        if metadata.setdefault(row["element_id"], pair) != pair:
            raise EnformerSeiEvaluationError("element metadata differs")
    if len(elements) != 1033 or len({value[0] for value in metadata.values()}) != 239:
        raise EnformerSeiEvaluationError("element or block census differs")
    outcomes = load_outcomes(outcome_tree / outcome["member"], set(elements))
    observed = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements], dtype=np.float64)
        for context in CONTEXTS
    }
    block_ids = np.asarray([metadata[element][0] for element in elements])

    prediction_rows = _read_rows(fit_root / "oof_predictions.tsv.gz", PREDICTION_FIELDS)
    if len(prediction_rows) != len(CANDIDATES) * 10330:
        raise EnformerSeiEvaluationError("prediction denominator differs")
    values = {
        (model, head, seed, context): np.full(1033, np.nan)
        for model, head in CANDIDATES
        for seed in SEEDS
        for context in CONTEXTS
    }
    identities = set()
    authority_fields = PREDICTION_FIELDS[:11]
    for row in prediction_rows:
        candidate = (row["model_id"], row["head_id"])
        seed = int(row["seed"])
        identity = (*candidate, seed, row["row_hash"])
        authority = row_lookup.get((seed, row["row_hash"]))
        if candidate not in CANDIDATES or seed not in SEEDS or identity in identities or authority is None:
            raise EnformerSeiEvaluationError("prediction identity differs")
        identities.add(identity)
        if any(row[field] != authority[field] for field in authority_fields):
            raise EnformerSeiEvaluationError("prediction metadata differs")
        if row["experimental_replicates"] != "4" or row["biological_donors"] != "0" or row["outcome_role"] != "exposed_development_MPRA_only":
            raise EnformerSeiEvaluationError("prediction topology differs")
        index = element_index[row["element_id"]]
        slot = values[(*candidate, seed, row["assay_context_id"])]
        if math.isfinite(slot[index]):
            raise EnformerSeiEvaluationError("duplicate prediction")
        slot[index] = _finite(row["prediction"])
    if len(identities) != len(CANDIDATES) * 10330 or any(not np.isfinite(array).all() for array in values.values()):
        raise EnformerSeiEvaluationError("prediction coverage differs")

    candidate_metadata = {
        (row["model_id"], row["head_id"]): row for row in config["candidates"]
    }
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
                metric = _metrics(observed[context], values[(model, head, seed, context)])
                context_metrics.append(metric)
                per_seed_rows.append({
                    "family": candidate_metadata[candidate]["family"],
                    "feature_space": candidate_metadata[candidate]["feature_space"],
                    "model_id": model,
                    "head_id": head,
                    "seed": seed,
                    "assay_context_id": context,
                    "elements": 1033,
                    **{name: _display(value) for name, value in metric.items()},
                })
            primary = _fisher_macro([metric["spearman"] for metric in context_metrics])
            per_seed_primary[(*candidate, seed)] = primary
            per_seed_rows.append({
                "family": candidate_metadata[candidate]["family"],
                "feature_space": candidate_metadata[candidate]["feature_space"],
                "model_id": model,
                "head_id": head,
                "seed": seed,
                "assay_context_id": "macro",
                "elements": 2066,
                "spearman": _display(primary),
            })
        by_context = {
            context: np.mean(np.vstack([values[(model, head, seed, context)] for seed in SEEDS]), axis=0)
            for context in CONTEXTS
        }
        ensemble_predictions[candidate] = by_context
        context_metrics = []
        for context in CONTEXTS:
            metric = _metrics(observed[context], by_context[context])
            context_metrics.append(metric)
            ensemble_rows.append({
                "family": candidate_metadata[candidate]["family"],
                "feature_space": candidate_metadata[candidate]["feature_space"],
                "model_id": model,
                "head_id": head,
                "assay_context_id": context,
                "elements": 1033,
                **{name: _display(value) for name, value in metric.items()},
            })
        point_primary[candidate] = _fisher_macro([metric["spearman"] for metric in context_metrics])
        ensemble_rows.append({
            "family": candidate_metadata[candidate]["family"],
            "feature_space": candidate_metadata[candidate]["feature_space"],
            "model_id": model,
            "head_id": head,
            "assay_context_id": "macro",
            "elements": 2066,
            "spearman": _display(point_primary[candidate]),
        })

    controls = [candidate for candidate in CANDIDATES if candidate[0] == "available_simple_controls" and point_primary[candidate] is not None]
    if not controls:
        raise EnformerSeiEvaluationError("no nonconstant available control")
    reference = max(controls, key=lambda candidate: (float(point_primary[candidate]), candidate[1]))
    bootstrap = _bootstrap(
        observed,
        ensemble_predictions,
        block_ids,
        resamples=int(config["uncertainty"]["resamples"]),
        seed=int(config["uncertainty"]["seed"]),
    )
    reference_samples = bootstrap[reference]
    interval_rows = []
    family_rows = []
    for candidate in CANDIDATES:
        samples = bootstrap[candidate]
        valid = np.isfinite(samples)
        if point_primary[candidate] is None or int(valid.sum()) < 9500:
            low = high = standard_error = None
        else:
            low, high = (float(value) for value in np.quantile(samples[valid], (0.025, 0.975)))
            standard_error = float(np.std(samples[valid], ddof=1))
        gain = None if point_primary[candidate] is None else float(point_primary[candidate] - point_primary[reference])
        gain_samples = samples - reference_samples
        gain_valid = np.isfinite(gain_samples)
        if gain is None or int(gain_valid.sum()) < 9500:
            gain_low = gain_high = None
        else:
            gain_low, gain_high = (float(value) for value in np.quantile(gain_samples[gain_valid], (0.025, 0.975)))
        positive_seeds = sum(
            per_seed_primary[(*candidate, seed)] is not None
            and per_seed_primary[(*reference, seed)] is not None
            and per_seed_primary[(*candidate, seed)] > per_seed_primary[(*reference, seed)]
            for seed in SEEDS
        )
        common = {
            "family": candidate_metadata[candidate]["family"],
            "feature_space": candidate_metadata[candidate]["feature_space"],
            "model_id": candidate[0],
            "head_id": candidate[1],
            "primary_metric": "tanh_mean_Fisher_z_Spearman_across_contexts",
            "primary_score": _display(point_primary[candidate]),
            "primary_ci_low": _display(low),
            "primary_ci_high": _display(high),
            "primary_bootstrap_se": _display(standard_error),
            "available_control_reference": f"{reference[0]}/{reference[1]}",
            "gain_vs_available_control": _display(gain),
            "gain_ci_low": _display(gain_low),
            "gain_ci_high": _display(gain_high),
            "positive_gain_seeds": positive_seeds,
            "bootstrap_resamples": int(config["uncertainty"]["resamples"]),
            "bootstrap_unit": "long_range_block_id",
        }
        interval_rows.append(common)
        family_rows.append({
            **common,
            "restricted_comparator": str(candidate_metadata[candidate].get("restricted", False)).lower(),
            "family_native_feature_space_preserved": "true",
            "cross_family_rank": "not_reported",
            "mandatory_baselines_complete": "false",
            "shortlist": "false",
            "champion_claim": "false",
            "native_equivalence_claim": "false",
            "external_claim": "false",
            "confirmatory_p_value": "not_reported_development_screen",
        })
    arguments.output.mkdir(parents=True, mode=0o750)
    _write_tsv(arguments.output / "per_seed_metrics.tsv", per_seed_rows)
    _write_tsv(arguments.output / "ensemble_metrics.tsv", ensemble_rows)
    _write_tsv(arguments.output / "paired_block_bootstrap.tsv", interval_rows)
    _write_tsv(arguments.output / "family_native_summary.tsv", family_rows)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_independent_restricted_development_evaluation",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "candidate_count": len(CANDIDATES),
        "prediction_rows": len(CANDIDATES) * 10330,
        "elements": 1033,
        "long_range_blocks": 239,
        "contexts": list(CONTEXTS),
        "fixed_seeds": list(SEEDS),
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "experimental_replicates_used_as_independent_donors": False,
        "available_control_reference": {"model_id": reference[0], "head_id": reference[1], "primary_score": point_primary[reference]},
        "bootstrap_resamples": int(config["uncertainty"]["resamples"]),
        "bootstrap_unit": "long_range_block_id",
        "family_native_results_reported_separately": True,
        "models_ranked_as_interchangeable": False,
        "overall_leader_reported": False,
        "mandatory_task_native_baselines_complete": False,
        "missing_baselines": list(MISSING_BASELINES),
        "shortlist_created": False,
        "champion_claim": False,
        "native_equivalence_claim": False,
        "external_evaluation": False,
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
    parser.add_argument("--fit-root", type=Path, required=True)
    parser.add_argument("--fit-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evaluate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
