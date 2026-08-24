#!/usr/bin/env python
"""Independent cross-campaign evaluator for the Cobolt RNA-to-ATAC smoke screen.

The Cobolt campaign reruns the mandatory baselines but not scPair. This
evaluator scores that exact seven-model campaign, then independently joins the
already frozen scPair OOF predictions on the identical row universe. It never
opens ATAC during model fit or prediction and never changes either campaign.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


def _load_module(filename: str, module_name: str) -> Any:
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load evaluator module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = _load_module(
    "rna_atac_development.py",
    "masld_bench_standalone_rna_atac_development_for_cobolt",
)
scpair_evaluator = _load_module(
    "rna_atac_scpair_development.py",
    "masld_bench_standalone_rna_atac_scpair_for_cobolt",
)


COBOLT_CAMPAIGN_MODEL_IDS = (
    "assay_native_pseudobulk",
    "cobolt",
    "mean_track",
    "nearest_context",
    "shrunken_pseudobulk",
    "shuffled_context",
    "trans_only",
)
MODEL_IDS = (
    "assay_native_pseudobulk",
    "cobolt",
    "mean_track",
    "nearest_context",
    "scpair",
    "shrunken_pseudobulk",
    "shuffled_context",
    "trans_only",
)
TASK_NATIVE_BASELINES = (
    "assay_native_pseudobulk",
    "mean_track",
    "shrunken_pseudobulk",
)
CANDIDATE_MODEL = "cobolt"
FROZEN_TRANSLATOR_COMPARATOR = "scpair"
NEGATIVE_CONTROL_MODEL = "shuffled_context"
LINEAR_CONTEXT_MODEL = "trans_only"
COBOLT_UPSTREAM_REVISION = "cf5a448c6539025346a6393c215ba329cb3a0183"
COBOLT_UPSTREAM_SOURCE_SHA256 = (
    "a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311"
)
COBOLT_ENVIRONMENT_SHA256 = (
    "c2308fc3b7524081c7a49ff1d840213306bcc20dff003010da2b898e2bdaa461"
)
ENDPOINT_FIELDS = scpair_evaluator.ENDPOINT_FIELDS
PROFILE_FIELDS = scpair_evaluator.PROFILE_FIELDS
PAIRED_BLOCK_FIELDS = scpair_evaluator.PAIRED_BLOCK_FIELDS


class RNAATACCoboltEvaluationError(RuntimeError):
    """Raised when the cross-campaign Cobolt evaluation contract is violated."""


def relative_deviance_reduction(reference: float, candidate: float) -> float:
    reference_value = float(reference)
    candidate_value = float(candidate)
    if (
        not math.isfinite(reference_value)
        or not math.isfinite(candidate_value)
        or reference_value <= 0
        or candidate_value < 0
    ):
        raise RNAATACCoboltEvaluationError("deviances must be finite and valid")
    return (reference_value - candidate_value) / reference_value


def select_strongest_task_native_baseline(
    models: Mapping[str, Mapping[str, Any]],
) -> str:
    totals: dict[str, float] = {}
    for model_id in TASK_NATIVE_BASELINES:
        record = models.get(model_id)
        if not isinstance(record, Mapping):
            raise RNAATACCoboltEvaluationError(
                f"missing task-native baseline result: {model_id}"
            )
        try:
            value = float(record["total_multinomial_deviance"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACCoboltEvaluationError(
                f"invalid task-native baseline deviance: {model_id}"
            ) from error
        if not math.isfinite(value) or value <= 0:
            raise RNAATACCoboltEvaluationError(
                f"invalid task-native baseline deviance: {model_id}"
            )
        totals[model_id] = value
    return min(TASK_NATIVE_BASELINES, key=lambda model_id: (totals[model_id], model_id))


def _artifact_record(path: Path, *, relative_to: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(relative_to).as_posix(),
        "sha256": base._sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _prediction_index(
    run_paths: Mapping[str, Mapping[int, Path]], model_id: str
) -> dict[str, dict[str, str]]:
    folds = run_paths.get(model_id)
    if not isinstance(folds, Mapping) or set(folds) != set(range(5)):
        raise RNAATACCoboltEvaluationError(f"{model_id} does not have five OOF folds")
    index: dict[str, dict[str, str]] = {}
    for fold, bundle_path in sorted(folds.items()):
        for row in base._load_bundle(bundle_path, model_id=model_id, fold=fold):
            row_hash = row["row_hash"]
            if row_hash in index:
                raise RNAATACCoboltEvaluationError(
                    f"duplicate OOF prediction row for {model_id}"
                )
            index[row_hash] = row
    return index


def _comparison_rows(
    outcomes: Sequence[Mapping[str, str]],
    *,
    candidate: Mapping[str, Mapping[str, str]],
    baseline: Mapping[str, Mapping[str, str]],
) -> list[dict[str, str]]:
    expected = {row["row_hash"] for row in outcomes}
    if len(expected) != len(outcomes) or set(candidate) != expected or set(baseline) != expected:
        raise RNAATACCoboltEvaluationError(
            "comparison prediction rows differ from evaluator outcomes"
        )
    return [
        {
            "row_hash": row["row_hash"],
            "donor_hash": row["donor_hash"],
            "block_hash": row["block_hash"],
            "stratum": row["stratum"],
            "observed": row["observed"],
            "candidate": candidate[row["row_hash"]]["predicted"],
            "baseline": baseline[row["row_hash"]]["predicted"],
        }
        for row in outcomes
    ]


def score_external_model(
    outcomes: Sequence[Mapping[str, str]],
    predictions: Mapping[str, Mapping[str, str]],
    *,
    model_id: str,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Re-derive base-evaluator metrics for an immutable external campaign."""

    import numpy as np
    from scipy.stats import spearmanr

    expected = {row["row_hash"] for row in outcomes}
    if len(expected) != len(outcomes) or set(predictions) != expected:
        raise RNAATACCoboltEvaluationError(
            f"{model_id} prediction universe differs from evaluator outcomes"
        )
    groups: dict[tuple[str, str], list[Mapping[str, str]]] = {}
    for row in outcomes:
        groups.setdefault((row["donor_hash"], row["stratum"]), []).append(row)

    block_deviances: list[float] = []
    auprcs: list[float] = []
    correlations: list[float] = []
    lineage_deviances: dict[str, list[float]] = {
        lineage: [] for lineage in base.LINEAGES
    }
    profile_rows: list[dict[str, str]] = []
    for (donor_hash, stratum), rows in sorted(groups.items()):
        observed = np.asarray([float(row["observed"]) for row in rows], dtype=np.float64)
        predicted = np.asarray(
            [float(predictions[row["row_hash"]]["predicted"]) for row in rows],
            dtype=np.float64,
        )
        if (
            observed.sum() <= 0
            or np.any(observed < 0)
            or np.any(predicted <= 0)
            or not np.all(np.isfinite(predicted))
        ):
            raise RNAATACCoboltEvaluationError(
                f"{model_id} has an invalid donor-lineage profile"
            )
        auprc = base._average_precision(observed > 0, predicted)
        correlation = float(spearmanr(observed, predicted).statistic)
        if not math.isfinite(correlation):
            correlation = 0.0
        auprcs.append(auprc)
        correlations.append(correlation)
        block_rows: dict[str, list[int]] = {}
        for index, row in enumerate(rows):
            block_rows.setdefault(row["block_hash"], []).append(index)
        for block_hash, positions_list in sorted(block_rows.items()):
            positions = np.asarray(positions_list, dtype=np.int64)
            block_observed = observed[positions]
            if block_observed.sum() <= 0:
                continue
            deviance = base._multinomial_deviance(
                block_observed, predicted[positions]
            )
            if not math.isfinite(deviance):
                raise RNAATACCoboltEvaluationError(
                    f"{model_id} block deviance is not finite"
                )
            block_deviances.append(deviance)
            lineage_deviances[stratum].append(deviance)
            profile_rows.append(
                {
                    "model_id": model_id,
                    "donor_hash": donor_hash,
                    "block_hash": block_hash,
                    "stratum": stratum,
                    "deviance": format(deviance, ".17g"),
                    "peak_auprc": format(auprc, ".17g"),
                    "profile_spearman": format(correlation, ".17g"),
                }
            )
    if not block_deviances or not auprcs or not correlations:
        raise RNAATACCoboltEvaluationError(f"{model_id} produced no scored profiles")
    summary = {
        "total_multinomial_deviance": sum(block_deviances),
        "mean_peak_auprc": sum(auprcs) / len(auprcs),
        "mean_profile_spearman": sum(correlations) / len(correlations),
        "lineage_mean_deviance": {
            lineage: sum(values) / len(values)
            for lineage, values in lineage_deviances.items()
            if values
        },
    }
    return summary, profile_rows


def paired_block_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    candidate_model: str,
    baseline_model: str,
) -> list[dict[str, str]]:
    required_models = {candidate_model, baseline_model}
    index: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in rows:
        model_id = row.get("model_id")
        if model_id not in required_models:
            continue
        key = (row["donor_hash"], row["block_hash"], row["stratum"])
        values = index.setdefault(key, {})
        if model_id in values:
            raise RNAATACCoboltEvaluationError(
                "profile metric appears more than once for a paired block"
            )
        try:
            value = float(row["deviance"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACCoboltEvaluationError("paired block deviance is invalid") from error
        if not math.isfinite(value) or value < 0:
            raise RNAATACCoboltEvaluationError("paired block deviance is invalid")
        values[str(model_id)] = value
    if not index or any(set(values) != required_models for values in index.values()):
        raise RNAATACCoboltEvaluationError(
            "candidate and baseline block universes differ"
        )
    return [
        {
            "donor_hash": donor_hash,
            "block_hash": block_hash,
            "stratum": stratum,
            "candidate_model": candidate_model,
            "baseline_model": baseline_model,
            "candidate_deviance": format(values[candidate_model], ".17g"),
            "baseline_deviance": format(values[baseline_model], ".17g"),
            "candidate_minus_baseline_deviance": format(
                values[candidate_model] - values[baseline_model], ".17g"
            ),
            "relative_deviance_reduction": format(
                relative_deviance_reduction(
                    values[baseline_model], values[candidate_model]
                ),
                ".17g",
            ),
        }
        for (donor_hash, block_hash, stratum), values in sorted(index.items())
    ]


def _validate_local_artifact(path: Path, record: Any, *, label: str) -> None:
    if not isinstance(record, Mapping) or set(record) != {
        "path",
        "sha256",
        "size_bytes",
    }:
        raise RNAATACCoboltEvaluationError(f"{label} artifact record differs")
    artifact = path.parent / str(record["path"])
    if (
        artifact.parent != path.parent
        or artifact.is_symlink()
        or not artifact.is_file()
        or artifact.stat().st_size != record["size_bytes"]
        or base._sha256_file(artifact) != record["sha256"]
    ):
        raise RNAATACCoboltEvaluationError(f"{label} artifact changed")


def _cobolt_fit_receipts(
    run_paths: Mapping[str, Mapping[int, Path]],
) -> list[dict[str, Any]]:
    result = []
    invariant_fields: dict[str, set[str]] = {
        "environment_lock_sha256": set(),
        "parameter_sha256": set(),
        "state_key_shape_dtype_sha256": set(),
    }
    for fold, bundle_path in sorted(run_paths[CANDIDATE_MODEL].items()):
        fitted = bundle_path.parent.parent / "002-fit" / "fitted_model.json"
        try:
            payload = json.loads(fitted.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise RNAATACCoboltEvaluationError(
                f"invalid Cobolt fitted-model receipt: {fitted}"
            ) from error
        exact = {
            "schema_version": "masld-bench-rna-atac-cobolt-model-v1",
            "model_id": CANDIDATE_MODEL,
            "held_out_fold": fold,
            "held_atac_used_for_fit": False,
            "query_rna_used_for_fit": False,
            "smoke_only": True,
            "champion_claim_allowed": False,
            "whole_module_pickle_saved": False,
            "dataset_adjustments_enabled": False,
            "query_fitted_latent_correction": False,
            "upstream_revision": COBOLT_UPSTREAM_REVISION,
            "upstream_model_source_sha256": COBOLT_UPSTREAM_SOURCE_SHA256,
            "environment_lock_sha256": COBOLT_ENVIRONMENT_SHA256,
            "training_objective": "cobolt_v1_0_1_three_elbo_combinations_raw_counts",
            "prediction_output": "rna_posterior_topic_proportions_times_atac_beta_softmax",
        }
        if any(payload.get(field) != expected for field, expected in exact.items()):
            raise RNAATACCoboltEvaluationError(
                f"Cobolt fitted-model binding differs for fold {fold}"
            )
        try:
            best_epoch = int(payload["best_epoch"])
            epochs_completed = int(payload["epochs_completed"])
            validation_nll = float(
                payload["best_validation_atac_fragment_nll_from_rna"]
            )
            training_donors = int(payload["inner_training_donor_count"])
            validation_donors = int(payload["inner_validation_donor_count"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACCoboltEvaluationError(
                f"Cobolt fit summary differs for fold {fold}"
            ) from error
        if (
            best_epoch < 0
            or epochs_completed < best_epoch
            or not math.isfinite(validation_nll)
            or validation_nll <= 0
            or training_donors <= 0
            or validation_donors <= 0
        ):
            raise RNAATACCoboltEvaluationError(
                f"Cobolt fit summary differs for fold {fold}"
            )
        _validate_local_artifact(fitted, payload.get("state_dict"), label="state_dict")
        _validate_local_artifact(
            fitted, payload.get("state_dict_manifest"), label="state_dict_manifest"
        )
        for field in invariant_fields:
            value = payload.get(field)
            if not isinstance(value, str) or len(value) != 64:
                raise RNAATACCoboltEvaluationError(
                    f"Cobolt {field} differs for fold {fold}"
                )
            invariant_fields[field].add(value)
        result.append(
            {
                "held_out_fold": fold,
                "best_epoch": best_epoch,
                "epochs_completed": epochs_completed,
                "best_validation_atac_fragment_nll_from_rna": validation_nll,
                "inner_training_donor_count": training_donors,
                "inner_validation_donor_count": validation_donors,
                "parameter_sha256": payload["parameter_sha256"],
                "state_dict_sha256": payload["state_dict"]["sha256"],
                "state_key_shape_dtype_sha256": payload[
                    "state_key_shape_dtype_sha256"
                ],
            }
        )
    if len(result) != 5 or any(len(values) != 1 for values in invariant_fields.values()):
        raise RNAATACCoboltEvaluationError("Cobolt five-fold fit invariants differ")
    return result


def evaluate(
    *,
    candidate: Path,
    execution_root: Path,
    scpair_candidate: Path,
    scpair_execution_root: Path,
    h5_path: Path,
    h5_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise RNAATACCoboltEvaluationError(
            f"evaluation output already exists: {output}"
        )
    output.mkdir(parents=True, exist_ok=False)

    base.MODEL_IDS = COBOLT_CAMPAIGN_MODEL_IDS
    raw_output = output / "seven_model_cobolt_campaign_evaluation"
    raw_result = base.evaluate(
        candidate=candidate,
        execution_root=execution_root,
        h5_path=h5_path,
        h5_sha256=h5_sha256,
        output=raw_output,
    )
    campaign_models = raw_result.get("models")
    if not isinstance(campaign_models, Mapping) or set(campaign_models) != set(
        COBOLT_CAMPAIGN_MODEL_IDS
    ):
        raise RNAATACCoboltEvaluationError(
            "base evaluator did not score the exact Cobolt campaign roster"
        )
    cobolt_run_paths = base._discover_runs(candidate, execution_root)
    base.MODEL_IDS = (FROZEN_TRANSLATOR_COMPARATOR,)
    scpair_run_paths = base._discover_runs(scpair_candidate, scpair_execution_root)
    base.MODEL_IDS = COBOLT_CAMPAIGN_MODEL_IDS

    endpoint_path = raw_output / "nearest_context_vs_strongest_baseline_endpoint.tsv"
    endpoint_fields, outcomes = base._read_tsv(endpoint_path)
    if endpoint_fields != ENDPOINT_FIELDS:
        raise RNAATACCoboltEvaluationError("base endpoint schema differs")
    prediction_indices = {
        model_id: _prediction_index(cobolt_run_paths, model_id)
        for model_id in (
            CANDIDATE_MODEL,
            NEGATIVE_CONTROL_MODEL,
            LINEAR_CONTEXT_MODEL,
            *TASK_NATIVE_BASELINES,
        )
    }
    prediction_indices[FROZEN_TRANSLATOR_COMPARATOR] = _prediction_index(
        scpair_run_paths, FROZEN_TRANSLATOR_COMPARATOR
    )
    scpair_summary, scpair_profile_rows = score_external_model(
        outcomes,
        prediction_indices[FROZEN_TRANSLATOR_COMPARATOR],
        model_id=FROZEN_TRANSLATOR_COMPARATOR,
    )
    models = {model_id: dict(record) for model_id, record in campaign_models.items()}
    models[FROZEN_TRANSLATOR_COMPARATOR] = scpair_summary
    if set(models) != set(MODEL_IDS):
        raise RNAATACCoboltEvaluationError("combined model roster differs")
    strongest = select_strongest_task_native_baseline(models)

    endpoint_artifacts: dict[str, Path] = {}
    for label, baseline_model in (
        ("task_native", strongest),
        ("shuffled_context", NEGATIVE_CONTROL_MODEL),
        ("trans_only", LINEAR_CONTEXT_MODEL),
        ("scpair", FROZEN_TRANSLATOR_COMPARATOR),
    ):
        path = output / f"cobolt_vs_{label}_endpoint.tsv"
        base._write_tsv(
            path,
            ENDPOINT_FIELDS,
            _comparison_rows(
                outcomes,
                candidate=prediction_indices[CANDIDATE_MODEL],
                baseline=prediction_indices[baseline_model],
            ),
        )
        endpoint_artifacts[label] = path

    raw_profile_path = raw_output / "profile_metrics.tsv"
    profile_fields, raw_profile_rows = base._read_tsv(raw_profile_path)
    if profile_fields != PROFILE_FIELDS:
        raise RNAATACCoboltEvaluationError("base profile metric schema differs")
    combined_profile_rows = sorted(
        [*raw_profile_rows, *scpair_profile_rows],
        key=lambda row: (
            row["model_id"],
            row["donor_hash"],
            row["block_hash"],
            row["stratum"],
        ),
    )
    profile_path = output / "eight_model_profile_metrics.tsv"
    base._write_tsv(profile_path, PROFILE_FIELDS, combined_profile_rows)

    paired_artifacts: dict[str, Path] = {}
    for label, baseline_model in (
        ("task_native", strongest),
        ("scpair", FROZEN_TRANSLATOR_COMPARATOR),
    ):
        path = output / f"cobolt_vs_{label}_donor_block_deviance.tsv"
        base._write_tsv(
            path,
            PAIRED_BLOCK_FIELDS,
            paired_block_rows(
                combined_profile_rows,
                candidate_model=CANDIDATE_MODEL,
                baseline_model=baseline_model,
            ),
        )
        paired_artifacts[label] = path

    candidate_deviance = float(models[CANDIDATE_MODEL]["total_multinomial_deviance"])
    comparator_deviances = {
        model_id: float(models[model_id]["total_multinomial_deviance"])
        for model_id in (
            strongest,
            NEGATIVE_CONTROL_MODEL,
            LINEAR_CONTEXT_MODEL,
            FROZEN_TRANSLATOR_COMPARATOR,
        )
    }
    lineage_gain = {
        lineage: relative_deviance_reduction(
            float(models[strongest]["lineage_mean_deviance"][lineage]),
            float(models[CANDIDATE_MODEL]["lineage_mean_deviance"][lineage]),
        )
        for lineage in base.LINEAGES
    }
    result = {
        **raw_result,
        "schema_version": "masld-bench-rna-atac-cobolt-development-evaluation-v1",
        "models": models,
        "evaluation_roster": list(MODEL_IDS),
        "cobolt_campaign_roster": list(COBOLT_CAMPAIGN_MODEL_IDS),
        "external_frozen_comparator_roster": [FROZEN_TRANSLATOR_COMPARATOR],
        "task_native_baseline_roster": list(TASK_NATIVE_BASELINES),
        "strongest_task_native_baseline": strongest,
        "cobolt_relative_deviance_reduction_vs_task_native": relative_deviance_reduction(
            comparator_deviances[strongest], candidate_deviance
        ),
        "cobolt_relative_deviance_reduction_vs_shuffled_context": relative_deviance_reduction(
            comparator_deviances[NEGATIVE_CONTROL_MODEL], candidate_deviance
        ),
        "cobolt_relative_deviance_reduction_vs_trans_only": relative_deviance_reduction(
            comparator_deviances[LINEAR_CONTEXT_MODEL], candidate_deviance
        ),
        "cobolt_relative_deviance_reduction_vs_scpair": relative_deviance_reduction(
            comparator_deviances[FROZEN_TRANSLATOR_COMPARATOR], candidate_deviance
        ),
        "cobolt_lineage_relative_deviance_reduction_vs_task_native": lineage_gain,
        "task_development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_major_lineages": 4,
            "maximum_any_lineage_worsening": 0.02,
            "relative_deviance_threshold_passed": relative_deviance_reduction(
                comparator_deviances[strongest], candidate_deviance
            )
            >= 0.05,
            "improved_lineage_count": sum(value > 0 for value in lineage_gain.values()),
            "no_lineage_worse_by_more_than_two_percent": all(
                value >= -0.02 for value in lineage_gain.values()
            ),
        },
        "fit_selections_by_outer_fold": scpair_evaluator._regularized_fit_selections(
            cobolt_run_paths
        ),
        "cobolt_fit_receipts_by_outer_fold": _cobolt_fit_receipts(cobolt_run_paths),
        "scpair_fit_receipts_by_outer_fold": scpair_evaluator._scpair_fit_receipts(
            scpair_run_paths
        ),
        "training_and_inference_firewall_all_folds_passed": True,
        "cross_campaign_join_contract": (
            "scPair is not retrained. Its immutable five-fold OOF predictions are "
            "joined only when their complete row-hash universe exactly matches the "
            "Cobolt evaluator outcome universe."
        ),
        "comparison_policy": (
            "Primary comparison uses minimum total multinomial deviance among "
            "the prospectively frozen mean-track, assay-native donor-lineage "
            "pseudobulk, and nested-CV shrunken pseudobulk baselines. The frozen "
            "scPair campaign is the direct-translator comparator; trans-only is "
            "the context-linear comparator and shuffled context is a negative control."
        ),
        "smoke_only": True,
        "champion_claim_allowed": False,
        "conditional_model_trigger_eligible": False,
        "conditional_model_trigger_ineligibility_reason": (
            "Single-study one-seed smoke fixture with source-wide consensus peaks."
        ),
        "campaign_artifact_sha256": base._sha256_file(candidate / "ARTIFACTS.json"),
        "scpair_campaign_artifact_sha256": base._sha256_file(
            scpair_candidate / "ARTIFACTS.json"
        ),
        "artifacts": {
            "base_cobolt_campaign_evaluation": _artifact_record(
                raw_output / "evaluation.json", relative_to=output
            ),
            "eight_model_profile_metrics": _artifact_record(
                profile_path, relative_to=output
            ),
            **{
                f"cobolt_vs_{label}_endpoint": _artifact_record(
                    path, relative_to=output
                )
                for label, path in endpoint_artifacts.items()
            },
            **{
                f"cobolt_vs_{label}_donor_block_deviance": _artifact_record(
                    path, relative_to=output
                )
                for label, path in paired_artifacts.items()
            },
        },
    }
    base._write_json(output / "evaluation.json", result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--execution-root", required=True, type=Path)
    parser.add_argument("--scpair-candidate", required=True, type=Path)
    parser.add_argument("--scpair-execution-root", required=True, type=Path)
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--h5-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    result = evaluate(
        candidate=arguments.candidate,
        execution_root=arguments.execution_root,
        scpair_candidate=arguments.scpair_candidate,
        scpair_execution_root=arguments.scpair_execution_root,
        h5_path=arguments.h5,
        h5_sha256=arguments.h5_sha256,
        output=arguments.output,
    )
    print(base._canonical_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
