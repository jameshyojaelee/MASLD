#!/usr/bin/env python3
"""Independent evaluator for frozen GSE274114 within-platform predictions."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import tomllib

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.contracts import PredictionBundle
from masld_bench.hashing import sha256_file


REGISTERED_METRICS = (
    "participant_macro_f1",
    "participant_auprc",
    "participant_auroc",
    "participant_brier",
    "participant_log_loss",
)
MODEL_KINDS = (
    "training_class_prior",
    "gene_rank_nearest_centroid",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
STOCHASTIC_SEEDS = (1103, 2207, 3301, 4409, 5501)
BOOTSTRAP_SEED = 274_114
TASKS: dict[str, dict[str, Any]] = {
    "gse274114_hiseq_healthy_vs_hbv": {
        "groups": {"CTRL": 0, "ENEG": 1},
        "group_counts": {"CTRL": 9, "ENEG": 11},
        "participants": 20,
        "model_prefix": "hiseq_",
        "split_id": "gse274114_hiseq_ctrl_eneg_stratified_participant_outer_v1",
        "task_spec_sha256": "2a05bf09ac7fccd4b52cac787f089a956260bf4e150f29a53f0e4843c247d113",
        "interpretation": "development-only within-HiSeq healthy-versus-HBV etiology utility",
        "prohibited_claim": "not a MASLD case-control, pooled, cross-platform, external-transfer, diagnostic, or champion result",
    },
    "gse274114_novaseq_mash_vs_mash_hbv": {
        "groups": {"NASH": 0, "ENEG_NASH": 1},
        "group_counts": {"NASH": 10, "ENEG_NASH": 9},
        "participants": 19,
        "model_prefix": "novaseq_",
        "split_id": "gse274114_novaseq_nash_eneg_nash_stratified_participant_outer_v1",
        "task_spec_sha256": "6150de64a148d5b1fb8f07bdac783f70931924211a4486302cbfadc5ce471ce5",
        "interpretation": "development-only within-NovaSeq MASH-versus-MASH-plus-HBV comorbidity utility",
        "prohibited_claim": "not a MASH-versus-non-MASH, pooled, cross-platform, external-transfer, diagnostic, or champion result",
    },
}


class GSE274114IndependentEvaluatorError(RuntimeError):
    """Raised when prediction or evaluator requirement differs."""


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise GSE274114IndependentEvaluatorError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), [dict(row) for row in reader]


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def _expected_seeds(model_kind: str) -> tuple[int, ...]:
    if model_kind in {"hvg_pca_elastic_net", "hvg_pca_linear_svm"}:
        return STOCHASTIC_SEEDS
    if model_kind in {"training_class_prior", "gene_rank_nearest_centroid"}:
        return (0,)
    raise GSE274114IndependentEvaluatorError("model kind is not registered")


def _validate_task_spec(path: Path, task_id: str) -> None:
    task = TASKS[task_id]
    if sha256_file(path) != task["task_spec_sha256"]:
        raise GSE274114IndependentEvaluatorError("TaskSpec identity differs")
    spec = tomllib.loads(path.read_text(encoding="utf-8"))
    if (
        spec.get("task_id") != task_id
        or tuple(spec.get("metrics", ())) != REGISTERED_METRICS
        or spec.get("unit_of_inference") != "participant"
        or spec.get("split_id") != task["split_id"]
        or spec.get("bootstrap_replicates") != 10_000
        or tuple(spec.get("resampling_units", ())) != ("participant",)
        or spec.get("uncertainty_method")
        != "source_label_stratified_participant_bootstrap"
        or set(spec.get("baseline_model_ids", ()))
        != {task["model_prefix"] + model for model in MODEL_KINDS}
        or spec.get("evaluator_parameters", {}).get("four_class_performance_allowed")
        is not False
        or spec.get("evaluator_parameters", {}).get("global_ood_accuracy_allowed")
        is not False
        or "champion claim is eligible" not in spec.get("claim_gate", "")
    ):
        raise GSE274114IndependentEvaluatorError("TaskSpec evaluation contract differs")


def _verify_model_prediction_set(
    *, task_root: Path, task_id: str, model_kind: str
) -> dict[str, Any]:
    task = TASKS[task_id]
    index_path = task_root / "evaluator_bundle_indices" / f"{model_kind}.tsv"
    fields, index_rows = _read_tsv(index_path)
    if fields != ("outer_fold", "seed", "bundle_root", "artifacts_sha256"):
        raise GSE274114IndependentEvaluatorError("bundle index schema differs")
    seeds = _expected_seeds(model_kind)
    expected_runs = {(fold, seed) for fold in range(5) for seed in seeds}
    observed_runs = {
        (int(row["outer_fold"]), int(row["seed"])) for row in index_rows
    }
    if observed_runs != expected_runs or len(index_rows) != len(expected_runs):
        raise GSE274114IndependentEvaluatorError("prediction set is incomplete")

    model_id = task["model_prefix"] + model_kind
    by_row: dict[str, dict[int, tuple[int, float]]] = {}
    verified_bundles = 0
    for row in sorted(
        index_rows, key=lambda value: (int(value["outer_fold"]), int(value["seed"]))
    ):
        fold = int(row["outer_fold"])
        seed = int(row["seed"])
        expected_root = (
            task_root
            / "prediction_bundles"
            / model_kind
            / f"fold-{fold}"
            / f"seed-{seed}"
        )
        bundle_root = Path(row["bundle_root"])
        if bundle_root != expected_root or not bundle_root.is_absolute():
            raise GSE274114IndependentEvaluatorError("bundle path identity differs")
        if sha256_file(bundle_root / "ARTIFACTS.json") != row["artifacts_sha256"]:
            raise GSE274114IndependentEvaluatorError("bundle manifest identity differs")
        manifest = verify_frozen_tree(bundle_root)
        if (
            manifest["metadata"].get("artifact_class")
            != "gse274114_frozen_prediction_bundle"
            or manifest["metadata"].get("task_id") != task_id
            or manifest["metadata"].get("model_id") != model_id
            or manifest["metadata"].get("outer_fold") != fold
            or manifest["metadata"].get("seed") != seed
            or manifest["metadata"].get("query_labels_opened") is not False
            or manifest["metadata"].get("scoring_run") is not False
        ):
            raise GSE274114IndependentEvaluatorError("bundle firewall differs")
        bundle = PredictionBundle.load_json(bundle_root / "prediction_bundle.json")
        bundle.validate_artifacts(bundle_root)
        if (
            bundle.task_id != task_id
            or bundle.model_id != model_id
            or bundle.dataset_ids != ("gse274114_mash_hbv",)
            or bundle.split_id != task["split_id"]
            or bundle.biological_unit != "participant"
            or bundle.row_id_field != "row_hash"
            or bundle.unit_id_field != "unit_hash"
            or bundle.format_version != "gse274114_binary_probability_v1"
            or bundle.metadata.get("outer_fold") != fold
            or bundle.metadata.get("seed") != seed
            or bundle.metadata.get("query_labels_opened") is not False
            or bundle.metadata.get("scoring_run") is not False
        ):
            raise GSE274114IndependentEvaluatorError("PredictionBundle contract differs")
        prediction_fields, prediction_rows = _read_tsv(
            bundle_root / bundle.standardized_table.path
        )
        if prediction_fields != ("row_hash", "unit_hash", "predicted") or any(
            "label" in field.lower() or "observed" in field.lower()
            for field in prediction_fields
        ):
            raise GSE274114IndependentEvaluatorError("prediction table exposes outcomes")
        if len(prediction_rows) != bundle.n_predictions:
            raise GSE274114IndependentEvaluatorError("bundle prediction count differs")
        for prediction in prediction_rows:
            row_id = prediction["row_hash"]
            if not row_id or row_id != prediction["unit_hash"]:
                raise GSE274114IndependentEvaluatorError("participant identity differs")
            probability = float(prediction["predicted"])
            if not math.isfinite(probability) or not 0.0 <= probability <= 1.0:
                raise GSE274114IndependentEvaluatorError("probability is invalid")
            if seed in by_row.setdefault(row_id, {}):
                raise GSE274114IndependentEvaluatorError(
                    "participant prediction is duplicated across folds"
                )
            by_row[row_id][seed] = (fold, probability)
        verified_bundles += 1
    if (
        len(by_row) != task["participants"]
        or any(set(seed_records) != set(seeds) for seed_records in by_row.values())
        or any(
            len({fold for fold, _ in seed_records.values()}) != 1
            for seed_records in by_row.values()
        )
    ):
        raise GSE274114IndependentEvaluatorError(
            "participant prediction inventory differs"
        )
    return {
        "task_id": task_id,
        "model_kind": model_kind,
        "model_id": model_id,
        "seeds": seeds,
        "verified_bundles": verified_bundles,
        "by_row": by_row,
    }


def verify_all_predictions_before_labels(
    *,
    campaign_root: Path,
    campaign_artifacts_sha256: str,
    task_spec_paths: Mapping[str, Path],
) -> dict[str, dict[str, dict[str, Any]]]:
    if sha256_file(campaign_root / "ARTIFACTS.json") != campaign_artifacts_sha256:
        raise GSE274114IndependentEvaluatorError("bundled campaign identity differs")
    parent_manifest = verify_frozen_tree(campaign_root)
    if (
        parent_manifest["metadata"].get("artifact_class")
        != "gse274114_bundled_baseline_execution_v2"
        or parent_manifest["metadata"].get("status")
        != "passed_prediction_only_bundled_campaign"
        or parent_manifest["metadata"].get("task_count") != 2
        or parent_manifest["metadata"].get("prediction_bundle_count") != 120
        or parent_manifest["metadata"].get("scoring_or_evaluation_run") is not False
        or parent_manifest["metadata"].get("prediction_outputs_contain_observed_labels")
        is not False
        or parent_manifest["metadata"].get("prior_v1_attempts_reused") is not False
    ):
        raise GSE274114IndependentEvaluatorError("bundled campaign firewall differs")
    receipt = json.loads(
        (campaign_root / "bundled_campaign_receipt.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("status") != "passed_prediction_only_bundled_campaign"
        or receipt.get("task_count") != 2
        or receipt.get("prediction_bundle_count") != 120
        or receipt.get("scoring_or_evaluation_run") is not False
        or {row.get("task_id") for row in receipt.get("task_campaigns", ())}
        != set(TASKS)
    ):
        raise GSE274114IndependentEvaluatorError("bundled campaign receipt differs")

    verified: dict[str, dict[str, dict[str, Any]]] = {}
    task_rows = {row["task_id"]: row for row in receipt["task_campaigns"]}
    for task_id in TASKS:
        _validate_task_spec(task_spec_paths[task_id], task_id)
        task_root = campaign_root / "campaigns" / task_id
        row = task_rows[task_id]
        if (
            Path(row["campaign_root"]) != task_root
            or sha256_file(task_root / "ARTIFACTS.json")
            != row["campaign_artifacts_sha256"]
        ):
            raise GSE274114IndependentEvaluatorError("task campaign identity differs")
        task_manifest = verify_frozen_tree(task_root)
        if (
            task_manifest["metadata"].get("artifact_class")
            != "gse274114_prediction_only_baseline_task_campaign"
            or task_manifest["metadata"].get("task_id") != task_id
            or task_manifest["metadata"].get("prediction_bundle_count") != 60
            or task_manifest["metadata"].get("scoring_or_evaluation_run") is not False
            or task_manifest["metadata"].get("sealed_or_champion_claim_eligible")
            is not False
        ):
            raise GSE274114IndependentEvaluatorError("task campaign firewall differs")
        verified[task_id] = {
            model_kind: _verify_model_prediction_set(
                task_root=task_root, task_id=task_id, model_kind=model_kind
            )
            for model_kind in MODEL_KINDS
        }
        row_sets = [set(model["by_row"]) for model in verified[task_id].values()]
        if any(row_set != row_sets[0] for row_set in row_sets[1:]):
            raise GSE274114IndependentEvaluatorError(
                "models do not share a task participant inventory"
            )
    task_row_sets = [
        set(next(iter(models.values()))["by_row"]) for models in verified.values()
    ]
    if task_row_sets[0] & task_row_sets[1]:
        raise GSE274114IndependentEvaluatorError("task participant inventories overlap")
    return verified


def _registered_metrics(truth: Any, probability: Any) -> dict[str, float]:
    import numpy as np
    from sklearn.metrics import (
        average_precision_score,
        brier_score_loss,
        f1_score,
        log_loss,
        roc_auc_score,
    )

    y = np.asarray(truth, dtype=np.int64)
    p = np.asarray(probability, dtype=np.float64)
    if (
        y.ndim != 1
        or p.shape != y.shape
        or set(y.tolist()) != {0, 1}
        or not np.all(np.isfinite(p))
        or np.any((p < 0) | (p > 1))
    ):
        raise GSE274114IndependentEvaluatorError("metric inputs differ")
    hard = (p >= 0.5).astype(np.int64)
    return {
        "participant_macro_f1": float(
            f1_score(y, hard, average="macro", zero_division=0)
        ),
        "participant_auprc": float(average_precision_score(y, p)),
        "participant_auroc": float(roc_auc_score(y, p)),
        "participant_brier": float(brier_score_loss(y, p)),
        "participant_log_loss": float(
            log_loss(y, np.clip(p, 1e-8, 1 - 1e-8), labels=[0, 1])
        ),
    }


def _participant_stratified_bootstrap(
    truth: Any, probability: Any, *, replicates: int
) -> dict[str, dict[str, float]]:
    import numpy as np

    if replicates != 10_000:
        raise GSE274114IndependentEvaluatorError(
            "bootstrap replicate count is not registered"
        )
    y = np.asarray(truth, dtype=np.int64)
    p = np.asarray(probability, dtype=np.float64)
    class_indices = tuple(np.flatnonzero(y == value) for value in (0, 1))
    if any(len(indices) == 0 for indices in class_indices):
        raise GSE274114IndependentEvaluatorError("bootstrap class is empty")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    samples = {metric: [] for metric in REGISTERED_METRICS}
    for _ in range(replicates):
        indices = np.concatenate(
            [
                rng.choice(class_rows, size=len(class_rows), replace=True)
                for class_rows in class_indices
            ]
        )
        current = _registered_metrics(y[indices], p[indices])
        for metric in REGISTERED_METRICS:
            samples[metric].append(current[metric])
    return {
        metric: {
            "lower_2_5": float(np.quantile(values, 0.025)),
            "upper_97_5": float(np.quantile(values, 0.975)),
        }
        for metric, values in samples.items()
    }


def run_evaluation(
    *,
    campaign_root: Path,
    campaign_artifacts_sha256: str,
    labels_path: Path,
    labels_sha256: str,
    folds_path: Path,
    folds_sha256: str,
    task_spec_paths: Mapping[str, Path],
    bootstrap_replicates: int,
    output: Path,
) -> dict[str, Any]:
    import numpy as np

    if output.exists() or bootstrap_replicates != 10_000:
        raise GSE274114IndependentEvaluatorError("evaluation target or policy differs")

    # This must finish for both tasks and all models before evaluator labels open.
    verified = verify_all_predictions_before_labels(
        campaign_root=campaign_root,
        campaign_artifacts_sha256=campaign_artifacts_sha256,
        task_spec_paths=task_spec_paths,
    )
    if sha256_file(labels_path) != labels_sha256 or sha256_file(folds_path) != folds_sha256:
        raise GSE274114IndependentEvaluatorError("evaluator authority identity differs")
    label_fields, label_rows = _read_tsv(labels_path)
    fold_fields, fold_rows = _read_tsv(folds_path)
    if label_fields[:2] != ("row_id", "source_group") or fold_fields[:3] != (
        "row_id",
        "outer_fold",
        "source_group",
    ):
        raise GSE274114IndependentEvaluatorError("evaluator authority schema differs")
    labels_by_row = {row["row_id"]: row["source_group"] for row in label_rows}
    fold_by_row = {
        row["row_id"]: (int(row["outer_fold"]), row["source_group"])
        for row in fold_rows
    }
    if len(labels_by_row) != len(label_rows) or len(fold_by_row) != len(fold_rows):
        raise GSE274114IndependentEvaluatorError("evaluator participant is duplicated")

    output.mkdir(parents=True)
    result_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    task_receipts: list[dict[str, Any]] = []
    for task_id, task in TASKS.items():
        expected_rows = {
            row_id
            for row_id, source_group in labels_by_row.items()
            if source_group in task["groups"]
        }
        observed_group_counts = {
            source_group: sum(
                labels_by_row[row_id] == source_group for row_id in expected_rows
            )
            for source_group in task["groups"]
        }
        if (
            len(expected_rows) != task["participants"]
            or observed_group_counts != task["group_counts"]
            or not expected_rows <= set(fold_by_row)
            or any(
                fold_by_row[row_id][1] != labels_by_row[row_id]
                or fold_by_row[row_id][0] not in range(5)
                for row_id in expected_rows
            )
        ):
            raise GSE274114IndependentEvaluatorError("task evaluator inventory differs")
        task_output = output / "task_results" / task_id
        task_output.mkdir(parents=True)
        for model_kind in MODEL_KINDS:
            model = verified[task_id][model_kind]
            if set(model["by_row"]) != expected_rows:
                raise GSE274114IndependentEvaluatorError(
                    "prediction/outcome inventory differs"
                )
            rows = sorted(expected_rows)
            truth = np.asarray(
                [task["groups"][labels_by_row[row_id]] for row_id in rows],
                dtype=np.int64,
            )
            probability = np.asarray(
                [
                    np.mean(
                        [model["by_row"][row_id][seed][1] for seed in model["seeds"]]
                    )
                    for row_id in rows
                ],
                dtype=np.float64,
            )
            if any(
                model["by_row"][row_id][seed][0] != fold_by_row[row_id][0]
                for row_id in rows
                for seed in model["seeds"]
            ):
                raise GSE274114IndependentEvaluatorError(
                    "prediction does not follow evaluator-held fold"
                )
            metrics = _registered_metrics(truth, probability)
            intervals = _participant_stratified_bootstrap(
                truth, probability, replicates=bootstrap_replicates
            )
            hard = (probability >= 0.5).astype(np.int64)
            _write_tsv(
                task_output / f"{model_kind}.joined_evaluator_rows.tsv",
                (
                    "row_hash",
                    "outer_fold",
                    "observed",
                    "predicted_probability",
                    "predicted_class",
                ),
                (
                    {
                        "row_hash": row_id,
                        "outer_fold": fold_by_row[row_id][0],
                        "observed": int(observed),
                        "predicted_probability": format(float(predicted), ".17g"),
                        "predicted_class": int(predicted_class),
                    }
                    for row_id, observed, predicted, predicted_class in zip(
                        rows, truth, probability, hard
                    )
                ),
            )
            for metric in REGISTERED_METRICS:
                interval_rows.append(
                    {
                        "task_id": task_id,
                        "platform_specific_interpretation": task["interpretation"],
                        "model_kind": model_kind,
                        "model_id": model["model_id"],
                        "metric": metric,
                        "estimate": format(metrics[metric], ".17g"),
                        "lower_2_5": format(intervals[metric]["lower_2_5"], ".17g"),
                        "upper_97_5": format(intervals[metric]["upper_97_5"], ".17g"),
                    }
                )
            result_rows.append(
                {
                    "task_id": task_id,
                    "model_kind": model_kind,
                    "model_id": model["model_id"],
                    "participants": len(rows),
                    "prediction_bundles": model["verified_bundles"],
                    "seeds_ensembled": list(model["seeds"]),
                    "metrics": metrics,
                    "intervals": intervals,
                }
            )
        task_receipts.append(
            {
                "task_id": task_id,
                "participants": task["participants"],
                "source_group_counts": task["group_counts"],
                "platform_specific_interpretation": task["interpretation"],
                "prohibited_claim": task["prohibited_claim"],
                "model_count": len(MODEL_KINDS),
                "cross_platform_comparison_run": False,
                "task_champion_selected": False,
            }
        )
    _write_tsv(
        output / "registered_metric_intervals.tsv",
        (
            "task_id",
            "platform_specific_interpretation",
            "model_kind",
            "model_id",
            "metric",
            "estimate",
            "lower_2_5",
            "upper_97_5",
        ),
        interval_rows,
    )
    receipt = {
        "schema_version": "masld-bench-gse274114-independent-evaluation-v1",
        "status": "development_only_platform_specific_evaluation_complete",
        "campaign_artifacts_sha256": campaign_artifacts_sha256,
        "prediction_bundles_verified_before_labels_opened": True,
        "prediction_bundle_count": 120,
        "tasks": task_receipts,
        "model_results": result_rows,
        "registered_metrics": list(REGISTERED_METRICS),
        "biological_unit": "participant",
        "bootstrap_method": "source_label_stratified_participant_bootstrap",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "multiplicity_method": "primary_effect_and_interval_without_confirmatory_p_value",
        "nominal_or_confirmatory_p_values_computed": False,
        "global_four_class_scoring_run": False,
        "global_ood_scoring_run": False,
        "pooled_cross_platform_scoring_run": False,
        "mash_vs_non_mash_scoring_run": False,
        "task_champion_selected": False,
        "sealed_or_champion_claim_eligible": False,
        "diagnostic_or_clinical_claim_eligible": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    freeze_tree(
        output,
        {"artifact_class": "gse274114_independent_development_evaluation", **receipt},
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--campaign-artifacts-sha256", required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--labels-sha256", required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--folds-sha256", required=True)
    parser.add_argument("--hiseq-task-spec", type=Path, required=True)
    parser.add_argument("--novaseq-task-spec", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    receipt = run_evaluation(
        campaign_root=arguments.campaign_root,
        campaign_artifacts_sha256=arguments.campaign_artifacts_sha256,
        labels_path=arguments.labels,
        labels_sha256=arguments.labels_sha256,
        folds_path=arguments.folds,
        folds_sha256=arguments.folds_sha256,
        task_spec_paths={
            "gse274114_hiseq_healthy_vs_hbv": arguments.hiseq_task_spec,
            "gse274114_novaseq_mash_vs_mash_hbv": arguments.novaseq_task_spec,
        },
        bootstrap_replicates=arguments.bootstrap_replicates,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
