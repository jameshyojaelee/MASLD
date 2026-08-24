#!/usr/bin/env python
"""Independent six-model evaluator for RNA-conditioned ATAC smoke runs.

The base evaluator performs the evaluator-only ATAC outcome join. This module
adds the regularized baselines and prospectively named comparisons without
changing the completed three-model evaluation implementation.
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
        "masld_bench_standalone_rna_atac_development_for_regularized", path
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
    "shrunken_pseudobulk",
    "shuffled_context",
    "trans_only",
)
TASK_NATIVE_BASELINES = (
    "assay_native_pseudobulk",
    "mean_track",
    "shrunken_pseudobulk",
)
CANDIDATE_MODEL = "trans_only"
NEGATIVE_CONTROL_MODEL = "shuffled_context"
ENDPOINT_FIELDS = (
    "row_hash",
    "donor_hash",
    "block_hash",
    "stratum",
    "observed",
    "candidate",
    "baseline",
)


class RNAATACRegularizedEvaluationError(RuntimeError):
    """Raised when the six-model evaluation contract is violated."""


def relative_deviance_reduction(reference: float, candidate: float) -> float:
    """Return the relative deviance reduction, where larger is better."""

    reference_value = float(reference)
    candidate_value = float(candidate)
    if (
        not math.isfinite(reference_value)
        or not math.isfinite(candidate_value)
        or reference_value <= 0
        or candidate_value < 0
    ):
        raise RNAATACRegularizedEvaluationError("deviances must be finite and valid")
    return (reference_value - candidate_value) / reference_value


def select_strongest_task_native_baseline(
    models: Mapping[str, Mapping[str, Any]],
) -> str:
    """Select the prespecified task-native baseline with minimum deviance."""

    totals: dict[str, float] = {}
    for model_id in TASK_NATIVE_BASELINES:
        record = models.get(model_id)
        if not isinstance(record, Mapping):
            raise RNAATACRegularizedEvaluationError(
                f"missing task-native baseline result: {model_id}"
            )
        try:
            value = float(record["total_multinomial_deviance"])
        except (KeyError, TypeError, ValueError) as error:
            raise RNAATACRegularizedEvaluationError(
                f"invalid task-native baseline deviance: {model_id}"
            ) from error
        if not math.isfinite(value) or value <= 0:
            raise RNAATACRegularizedEvaluationError(
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
        raise RNAATACRegularizedEvaluationError(
            f"{model_id} does not have five OOF folds"
        )
    for fold, bundle_path in sorted(folds.items()):
        for row in base._load_bundle(bundle_path, model_id=model_id, fold=fold):
            row_hash = row["row_hash"]
            if row_hash in index:
                raise RNAATACRegularizedEvaluationError(
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
        raise RNAATACRegularizedEvaluationError("outcome row appears more than once")
    if set(candidate) != expected or set(baseline) != expected:
        raise RNAATACRegularizedEvaluationError(
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


def _fit_selections(
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
                raise RNAATACRegularizedEvaluationError(
                    f"invalid fitted-model receipt: {fitted}"
                ) from error
            if (
                payload.get("model_id") != model_id
                or payload.get("held_out_fold") != fold
                or payload.get("held_atac_used_for_fit") is not False
            ):
                raise RNAATACRegularizedEvaluationError(
                    f"fitted-model binding differs for {model_id} fold {fold}"
                )
            record = {"held_out_fold": fold}
            for name in names:
                if name not in payload:
                    raise RNAATACRegularizedEvaluationError(
                        f"fitted-model selection is missing {name}"
                    )
                record[name] = payload[name]
            records.append(record)
        result[model_id] = records
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
        raise RNAATACRegularizedEvaluationError(
            f"evaluation output already exists: {output}"
        )
    output.mkdir(parents=True, exist_ok=False)

    base.MODEL_IDS = MODEL_IDS
    raw_output = output / "six_model_base_evaluation"
    raw_result = base.evaluate(
        candidate=candidate,
        execution_root=execution_root,
        h5_path=h5_path,
        h5_sha256=h5_sha256,
        output=raw_output,
    )
    models = raw_result.get("models")
    if not isinstance(models, Mapping) or set(models) != set(MODEL_IDS):
        raise RNAATACRegularizedEvaluationError(
            "base evaluator did not score the exact six-model roster"
        )
    strongest = select_strongest_task_native_baseline(models)
    run_paths = base._discover_runs(candidate, execution_root)

    endpoint_path = raw_output / "nearest_context_vs_strongest_baseline_endpoint.tsv"
    fields, outcomes = base._read_tsv(endpoint_path)
    if fields != ENDPOINT_FIELDS:
        raise RNAATACRegularizedEvaluationError("base endpoint schema differs")
    prediction_indices = {
        model_id: _prediction_index(run_paths, model_id)
        for model_id in (CANDIDATE_MODEL, NEGATIVE_CONTROL_MODEL, strongest)
    }

    task_native_path = output / "trans_only_vs_strongest_task_native_endpoint.tsv"
    shuffled_path = output / "trans_only_vs_shuffled_context_endpoint.tsv"
    base._write_tsv(
        task_native_path,
        ENDPOINT_FIELDS,
        _comparison_rows(
            outcomes,
            candidate=prediction_indices[CANDIDATE_MODEL],
            baseline=prediction_indices[strongest],
        ),
    )
    base._write_tsv(
        shuffled_path,
        ENDPOINT_FIELDS,
        _comparison_rows(
            outcomes,
            candidate=prediction_indices[CANDIDATE_MODEL],
            baseline=prediction_indices[NEGATIVE_CONTROL_MODEL],
        ),
    )

    candidate_deviance = float(models[CANDIDATE_MODEL]["total_multinomial_deviance"])
    baseline_deviance = float(models[strongest]["total_multinomial_deviance"])
    shuffled_deviance = float(
        models[NEGATIVE_CONTROL_MODEL]["total_multinomial_deviance"]
    )
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
        "schema_version": "masld-bench-rna-atac-regularized-development-evaluation-v1",
        "evaluation_roster": list(MODEL_IDS),
        "task_native_baseline_roster": list(TASK_NATIVE_BASELINES),
        "strongest_task_native_baseline": strongest,
        "trans_only_relative_deviance_reduction_vs_task_native": (
            relative_deviance_reduction(baseline_deviance, candidate_deviance)
        ),
        "trans_only_relative_deviance_reduction_vs_shuffled_context": (
            relative_deviance_reduction(shuffled_deviance, candidate_deviance)
        ),
        "trans_only_lineage_relative_deviance_reduction_vs_task_native": (
            lineage_gain
        ),
        "context_signal_direction_consistent": candidate_deviance < shuffled_deviance,
        "fit_selections_by_outer_fold": _fit_selections(run_paths),
        "comparison_policy": (
            "Minimum total multinomial deviance among the prospectively frozen "
            "mean-track, assay-native donor-lineage pseudobulk, and nested-CV "
            "shrunken pseudobulk baselines; shuffled context is a negative control."
        ),
        "smoke_only": True,
        "champion_claim_allowed": False,
        "artifacts": {
            "base_evaluation": _artifact_record(
                raw_output / "evaluation.json", relative_to=output
            ),
            "profile_metrics": _artifact_record(
                raw_output / "profile_metrics.tsv", relative_to=output
            ),
            "trans_only_vs_task_native_endpoint": _artifact_record(
                task_native_path, relative_to=output
            ),
            "trans_only_vs_shuffled_context_endpoint": _artifact_record(
                shuffled_path, relative_to=output
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
