#!/usr/bin/env python
"""Independent seven-model evaluator for the scPair RNA-to-ATAC smoke screen.

The base evaluator performs the evaluator-only ATAC outcome join. This module
adds scPair to the frozen regularized roster and evaluates it against the
prospectively defined strongest task-native baseline without changing either
of the completed earlier evaluator implementations.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


def _load_base_evaluator() -> Any:
    path = Path(__file__).with_name("rna_atac_development.py")
    spec = importlib.util.spec_from_file_location(
        "masld_bench_standalone_rna_atac_development_for_scpair", path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load base RNA-ATAC evaluator: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = _load_base_evaluator()


MODEL_IDS = (
    "assay_native_pseudobulk",
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
CANDIDATE_MODEL = "scpair"
NEGATIVE_CONTROL_MODEL = "shuffled_context"
LINEAR_CONTEXT_MODEL = "trans_only"
SCPAIR_UPSTREAM_REVISION = "52d1ca1e26c8127b2023f3d80cfd1c140a0ea80a"
SCPAIR_ENVIRONMENT_SHA256 = (
    "c2308fc3b7524081c7a49ff1d840213306bcc20dff003010da2b898e2bdaa461"
)
ENDPOINT_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "observed",
    "candidate",
    "baseline",
)
PROFILE_FIELDS = (
    "model_id",
    "donor_hash",
    "block_hash",
    "stratum",
    "deviance",
    "peak_auprc",
    "profile_spearman",
)
PAIRED_BLOCK_FIELDS = (
    "donor_hash",
    "block_hash",
    "stratum",
    "candidate_model",
    "baseline_model",
    "candidate_deviance",
    "baseline_deviance",
    "candidate_minus_baseline_deviance",
    "relative_deviance_reduction",
)


class RNAATACScPairEvaluationError(RuntimeError):
    """Raised when the seven-model evaluation requirements are not met."""


def relative_deviance_reduction(reference: float, candidate: float) -> float:
    """Return relative deviance reduction, where larger is better."""

    reference_value = float(reference)
    candidate_value = float(candidate)
    if (
        not math.isfinite(reference_value)
        or not math.isfinite(candidate_value)
        or reference_value <= 0
        or candidate_value < 0
    ):
        raise RNAATACScPairEvaluationError("deviances must be finite and valid")
    return (reference_value - candidate_value) / reference_value


def select_strongest_task_native_baseline(
    models: Mapping[str, Mapping[str, Any]],
) -> str:
    """Select only from the prospectively frozen task-native baseline roster."""

    totals: dict[str, float] = {}
    for model_id in TASK_NATIVE_BASELINES:
        record = models.get(model_id)
        if not isinstance(record, Mapping):
            raise RNAATACScPairEvaluationError(
                f"missing task-native baseline result: {model_id}"
            )
        try:
            value = float(record["total_multinomial_deviance"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACScPairEvaluationError(
                f"invalid task-native baseline deviance: {model_id}"
            ) from error
        if not math.isfinite(value) or value <= 0:
            raise RNAATACScPairEvaluationError(
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
    index: dict[str, dict[str, str]] = {}
    folds = run_paths.get(model_id)
    if not isinstance(folds, Mapping) or set(folds) != set(range(5)):
        raise RNAATACScPairEvaluationError(
            f"{model_id} does not have five OOF folds"
        )
    for fold, bundle_path in sorted(folds.items()):
        for row in base._load_bundle(bundle_path, model_id=model_id, fold=fold):
            row_hash = row["row_hash"]
            if row_hash in index:
                raise RNAATACScPairEvaluationError(
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
    if len(expected) != len(outcomes):
        raise RNAATACScPairEvaluationError("outcome row appears more than once")
    if set(candidate) != expected or set(baseline) != expected:
        raise RNAATACScPairEvaluationError(
            "comparison prediction rows differ from evaluator outcomes"
        )
    result = []
    for outcome in outcomes:
        row_hash = outcome["row_hash"]
        result.append(
            {
                "row_hash": row_hash,
                "donor_hash": outcome["donor_hash"],
                "block_hash": outcome["block_hash"],
                "stratum": outcome["stratum"],
                "observed": outcome["observed"],
                "candidate": candidate[row_hash]["predicted"],
                "baseline": baseline[row_hash]["predicted"],
            }
        )
    return result


def paired_block_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    candidate_model: str,
    baseline_model: str,
) -> list[dict[str, str]]:
    """Pair block-specific deviances on the exact donor/block/lineage universe."""

    required_models = {candidate_model, baseline_model}
    index: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in rows:
        model_id = row.get("model_id")
        if model_id not in required_models:
            continue
        key = (row["donor_hash"], row["block_hash"], row["stratum"])
        model_values = index.setdefault(key, {})
        if model_id in model_values:
            raise RNAATACScPairEvaluationError(
                "profile metric appears more than once for a paired block"
            )
        try:
            deviance = float(row["deviance"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACScPairEvaluationError(
                "paired block deviance is invalid"
            ) from error
        if not math.isfinite(deviance) or deviance < 0:
            raise RNAATACScPairEvaluationError(
                "paired block deviance is invalid"
            )
        model_values[str(model_id)] = deviance
    if not index or any(set(values) != required_models for values in index.values()):
        raise RNAATACScPairEvaluationError(
            "candidate and baseline block universes differ"
        )
    result = []
    for (donor_hash, block_hash, stratum), values in sorted(index.items()):
        candidate = values[candidate_model]
        baseline = values[baseline_model]
        result.append(
            {
                "donor_hash": donor_hash,
                "block_hash": block_hash,
                "stratum": stratum,
                "candidate_model": candidate_model,
                "baseline_model": baseline_model,
                "candidate_deviance": format(candidate, ".17g"),
                "baseline_deviance": format(baseline, ".17g"),
                "candidate_minus_baseline_deviance": format(
                    candidate - baseline, ".17g"
                ),
                "relative_deviance_reduction": format(
                    relative_deviance_reduction(baseline, candidate), ".17g"
                ),
            }
        )
    return result


def _regularized_fit_selections(
    run_paths: Mapping[str, Mapping[int, Path]],
) -> dict[str, list[dict[str, Any]]]:
    fields = {
        "shrunken_pseudobulk": ("selected_shrinkage_weight",),
        "shuffled_context": ("selected_pca_rank", "selected_ridge_lambda"),
        "trans_only": ("selected_pca_rank", "selected_ridge_lambda"),
    }
    result: dict[str, list[dict[str, Any]]] = {}
    for model_id, names in fields.items():
        records = []
        for fold, bundle_path in sorted(run_paths[model_id].items()):
            fitted = bundle_path.parent.parent / "002-fit" / "fitted_model.json"
            try:
                payload = json.loads(fitted.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error:
                raise RNAATACScPairEvaluationError(
                    f"invalid fitted-model receipt: {fitted}"
                ) from error
            if (
                payload.get("model_id") != model_id
                or payload.get("held_out_fold") != fold
                or payload.get("held_atac_used_for_fit") is not False
            ):
                raise RNAATACScPairEvaluationError(
                    f"fitted-model binding differs for {model_id} fold {fold}"
                )
            record = {"held_out_fold": fold}
            for name in names:
                if name not in payload:
                    raise RNAATACScPairEvaluationError(
                        f"fitted-model selection is missing {name}"
                    )
                record[name] = payload[name]
            records.append(record)
        result[model_id] = records
    return result


def _validate_local_artifact(path: Path, record: Any, *, label: str) -> None:
    if not isinstance(record, Mapping) or set(record) != {
        "path",
        "sha256",
        "size_bytes",
    }:
        raise RNAATACScPairEvaluationError(f"{label} artifact record differs")
    artifact = path.parent / str(record["path"])
    if (
        artifact.parent != path.parent
        or artifact.is_symlink()
        or not artifact.is_file()
        or artifact.stat().st_size != record["size_bytes"]
        or base._sha256_file(artifact) != record["sha256"]
    ):
        raise RNAATACScPairEvaluationError(f"{label} artifact changed")


def _scpair_fit_receipts(
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
            raise RNAATACScPairEvaluationError(
                f"invalid scPair fitted-model receipt: {fitted}"
            ) from error
        exact = {
            "schema_version": "masld-bench-rna-atac-scpair-model-v1",
            "model_id": CANDIDATE_MODEL,
            "held_out_fold": fold,
            "held_atac_used_for_fit": False,
            "query_rna_used_for_fit": False,
            "smoke_only": True,
            "champion_claim_allowed": False,
            "whole_module_pickle_saved": False,
            "upstream_revision": SCPAIR_UPSTREAM_REVISION,
            "environment_lock_sha256": SCPAIR_ENVIRONMENT_SHA256,
        }
        if any(payload.get(field) != expected for field, expected in exact.items()):
            raise RNAATACScPairEvaluationError(
                f"scPair fitted-model binding differs for fold {fold}"
            )
        try:
            best_epoch = int(payload["best_epoch"])
            epochs_completed = int(payload["epochs_completed"])
            best_validation_bce = float(payload["best_validation_bce"])
            training_donors = int(payload["inner_training_donor_count"])
            validation_donors = int(payload["inner_validation_donor_count"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACScPairEvaluationError(
                f"scPair fit summary differs for fold {fold}"
            ) from error
        if (
            best_epoch < 0
            or epochs_completed <= best_epoch
            or not math.isfinite(best_validation_bce)
            or best_validation_bce <= 0
            or training_donors <= 0
            or validation_donors <= 0
        ):
            raise RNAATACScPairEvaluationError(
                f"scPair fit summary differs for fold {fold}"
            )
        _validate_local_artifact(fitted, payload.get("state_dict"), label="state_dict")
        _validate_local_artifact(
            fitted, payload.get("state_dict_manifest"), label="state_dict_manifest"
        )
        for field in invariant_fields:
            value = payload.get(field)
            if not isinstance(value, str) or len(value) != 64:
                raise RNAATACScPairEvaluationError(
                    f"scPair {field} differs for fold {fold}"
                )
            invariant_fields[field].add(value)
        result.append(
            {
                "held_out_fold": fold,
                "best_epoch": best_epoch,
                "epochs_completed": epochs_completed,
                "best_validation_bce": best_validation_bce,
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
        raise RNAATACScPairEvaluationError(
            "scPair five-fold fit invariants differ"
        )
    return result


def evaluate(
    *,
    candidate: Path,
    execution_root: Path,
    h5_path: Path,
    h5_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise RNAATACScPairEvaluationError(
            f"evaluation output already exists: {output}"
        )
    output.mkdir(parents=True, exist_ok=False)

    base.MODEL_IDS = MODEL_IDS
    raw_output = output / "seven_model_base_evaluation"
    raw_result = base.evaluate(
        candidate=candidate,
        execution_root=execution_root,
        h5_path=h5_path,
        h5_sha256=h5_sha256,
        output=raw_output,
    )
    models = raw_result.get("models")
    if not isinstance(models, Mapping) or set(models) != set(MODEL_IDS):
        raise RNAATACScPairEvaluationError(
            "base evaluator did not score the exact seven-model roster"
        )
    strongest = select_strongest_task_native_baseline(models)
    run_paths = base._discover_runs(candidate, execution_root)

    endpoint_path = raw_output / "nearest_context_vs_strongest_baseline_endpoint.tsv"
    fields, outcomes = base._read_tsv(endpoint_path)
    if fields != ENDPOINT_FIELDS:
        raise RNAATACScPairEvaluationError("base endpoint schema differs")
    prediction_indices = {
        model_id: _prediction_index(run_paths, model_id)
        for model_id in (
            CANDIDATE_MODEL,
            NEGATIVE_CONTROL_MODEL,
            LINEAR_CONTEXT_MODEL,
            strongest,
        )
    }

    task_native_path = output / "scpair_vs_strongest_task_native_endpoint.tsv"
    shuffled_path = output / "scpair_vs_shuffled_context_endpoint.tsv"
    linear_path = output / "scpair_vs_trans_only_endpoint.tsv"
    for path, baseline_model in (
        (task_native_path, strongest),
        (shuffled_path, NEGATIVE_CONTROL_MODEL),
        (linear_path, LINEAR_CONTEXT_MODEL),
    ):
        base._write_tsv(
            path,
            ENDPOINT_FIELDS,
            _comparison_rows(
                outcomes,
                candidate=prediction_indices[CANDIDATE_MODEL],
                baseline=prediction_indices[baseline_model],
            ),
        )

    profile_path = raw_output / "profile_metrics.tsv"
    profile_fields, profile_rows = base._read_tsv(profile_path)
    if profile_fields != PROFILE_FIELDS:
        raise RNAATACScPairEvaluationError("base profile metric schema differs")
    paired_block_path = output / "scpair_vs_task_native_donor_block_deviance.tsv"
    paired_rows = paired_block_rows(
        profile_rows,
        candidate_model=CANDIDATE_MODEL,
        baseline_model=strongest,
    )
    base._write_tsv(paired_block_path, PAIRED_BLOCK_FIELDS, paired_rows)

    candidate_deviance = float(models[CANDIDATE_MODEL]["total_multinomial_deviance"])
    baseline_deviance = float(models[strongest]["total_multinomial_deviance"])
    shuffled_deviance = float(
        models[NEGATIVE_CONTROL_MODEL]["total_multinomial_deviance"]
    )
    linear_deviance = float(models[LINEAR_CONTEXT_MODEL]["total_multinomial_deviance"])
    lineage_gain = {}
    for lineage in base.LINEAGES:
        reference = float(models[strongest]["lineage_mean_deviance"][lineage])
        candidate_value = float(
            models[CANDIDATE_MODEL]["lineage_mean_deviance"][lineage]
        )
        lineage_gain[lineage] = relative_deviance_reduction(
            reference, candidate_value
        )

    result = {
        **raw_result,
        "schema_version": "masld-bench-rna-atac-scpair-development-evaluation-v1",
        "evaluation_roster": list(MODEL_IDS),
        "task_native_baseline_roster": list(TASK_NATIVE_BASELINES),
        "strongest_task_native_baseline": strongest,
        "scpair_relative_deviance_reduction_vs_task_native": (
            relative_deviance_reduction(baseline_deviance, candidate_deviance)
        ),
        "scpair_relative_deviance_reduction_vs_shuffled_context": (
            relative_deviance_reduction(shuffled_deviance, candidate_deviance)
        ),
        "scpair_relative_deviance_reduction_vs_trans_only": (
            relative_deviance_reduction(linear_deviance, candidate_deviance)
        ),
        "scpair_lineage_relative_deviance_reduction_vs_task_native": lineage_gain,
        "task_development_gate": {
            "minimum_relative_deviance_reduction": 0.05,
            "minimum_improved_major_lineages": 4,
            "maximum_any_lineage_worsening": 0.02,
            "relative_deviance_threshold_passed": relative_deviance_reduction(
                baseline_deviance, candidate_deviance
            )
            >= 0.05,
            "improved_lineage_count": sum(value > 0 for value in lineage_gain.values()),
            "no_lineage_worse_by_more_than_two_percent": all(
                value >= -0.02 for value in lineage_gain.values()
            ),
        },
        "fit_selections_by_outer_fold": _regularized_fit_selections(run_paths),
        "scpair_fit_receipts_by_outer_fold": _scpair_fit_receipts(run_paths),
        "scpair_training_and_inference_firewall_all_folds_passed": True,
        "comparison_policy": (
            "Primary comparison uses minimum total multinomial deviance among "
            "the prospectively frozen mean-track, assay-native donor-lineage "
            "pseudobulk, and nested-CV shrunken pseudobulk baselines. Trans-only "
            "is the strongest context-linear comparator and shuffled context is "
            "a negative control."
        ),
        "smoke_only": True,
        "champion_claim_allowed": False,
        "conditional_model_trigger_eligible": False,
        "conditional_model_trigger_ineligibility_reason": (
            "Single-study one-seed smoke fixture with source-wide consensus peaks."
        ),
        "artifacts": {
            "base_evaluation": _artifact_record(
                raw_output / "evaluation.json", relative_to=output
            ),
            "profile_metrics": _artifact_record(profile_path, relative_to=output),
            "scpair_vs_task_native_endpoint": _artifact_record(
                task_native_path, relative_to=output
            ),
            "scpair_vs_shuffled_context_endpoint": _artifact_record(
                shuffled_path, relative_to=output
            ),
            "scpair_vs_trans_only_endpoint": _artifact_record(
                linear_path, relative_to=output
            ),
            "scpair_vs_task_native_donor_block_deviance": _artifact_record(
                paired_block_path, relative_to=output
            ),
        },
    }
    base._write_json(output / "evaluation.json", result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--execution-root", required=True, type=Path)
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--h5-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args(argv)
    result = evaluate(
        candidate=arguments.candidate,
        execution_root=arguments.execution_root,
        h5_path=arguments.h5,
        h5_sha256=arguments.h5_sha256,
        output=arguments.output,
    )
    print(base._canonical_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
