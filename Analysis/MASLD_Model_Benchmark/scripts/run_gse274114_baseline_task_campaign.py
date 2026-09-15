#!/usr/bin/env python3
"""Bundle one GSE274114 within-platform task into prediction-only baselines."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from scripts.gse274114_baseline_firewall import (
    DETERMINISTIC_SEED,
    FIXED_STOCHASTIC_SEEDS,
    MODEL_KINDS,
    STOCHASTIC_MODELS,
    TASKS,
    fit_predict,
    prepare_fold_view,
    write_tsv,
)


class GSE274114BaselineCampaignError(RuntimeError):
    """Raised before a campaign can mix tasks, folds, seeds, or roles."""


def planned_runs(task_id: str) -> list[dict[str, Any]]:
    if task_id not in TASKS:
        raise GSE274114BaselineCampaignError("task is not admitted")
    runs: list[dict[str, Any]] = []
    for outer_fold in range(5):
        for model_kind in MODEL_KINDS:
            seeds = (
                FIXED_STOCHASTIC_SEEDS
                if model_kind in STOCHASTIC_MODELS
                else (DETERMINISTIC_SEED,)
            )
            for seed in seeds:
                runs.append(
                    {
                        "task_id": task_id,
                        "model_kind": model_kind,
                        "model_id": TASKS[task_id]["model_prefix"] + model_kind,
                        "outer_fold": outer_fold,
                        "seed": seed,
                    }
                )
    identities = {
        (run["model_id"], run["outer_fold"], run["seed"]) for run in runs
    }
    if (
        len(runs) != 60
        or len(identities) != 60
        or len(set(FIXED_STOCHASTIC_SEEDS)) != 5
        or DETERMINISTIC_SEED in FIXED_STOCHASTIC_SEEDS
    ):
        raise GSE274114BaselineCampaignError("run/seed census differs")
    return runs


def _fit_one(run: dict[str, Any]) -> dict[str, Any]:
    bundle = fit_predict(
        view_root=Path(run["view_root"]),
        expected_view_sha256=run["view_sha256"],
        model_kind=run["model_kind"],
        seed=run["seed"],
        output=Path(run["output"]),
    )
    return {
        **run,
        "bundle_id": bundle["bundle_id"],
        "run_id": bundle["run_id"],
        "bundle_artifacts_sha256": sha256_file(Path(run["output"]) / "ARTIFACTS.json"),
    }


def run_campaign(
    *,
    feature_root: Path,
    feature_artifacts_sha256: str,
    labels_path: Path,
    folds_path: Path,
    task_id: str,
    task_spec_path: Path,
    task_spec_sha256: str,
    publication_root: Path,
    workers: int,
    output: Path,
) -> dict[str, Any]:
    if (
        output.exists()
        or not publication_root.is_absolute()
        or publication_root == output
        or task_id not in TASKS
        or not 1 <= workers <= 16
    ):
        raise GSE274114BaselineCampaignError("invalid or existing campaign target")
    if (
        sha256_file(feature_root / "ARTIFACTS.json") != feature_artifacts_sha256
        or sha256_file(task_spec_path) != task_spec_sha256
    ):
        raise GSE274114BaselineCampaignError("feature or TaskSpec identity differs")
    feature_manifest = verify_frozen_tree(feature_root)
    if (
        feature_manifest["metadata"].get("artifact_class")
        != "gse274114_unlabeled_feature_matrix"
        or feature_manifest["metadata"].get("labels_opened") is not False
        or feature_manifest["metadata"].get("folds_opened") is not False
        or feature_manifest["metadata"].get("participant_metadata_opened") is not False
        or feature_manifest["metadata"].get("normalization_run") is not False
        or feature_manifest["metadata"].get("model_fit_or_scoring_run") is not False
    ):
        raise GSE274114BaselineCampaignError("feature authority firewall differs")
    task_text = task_spec_path.read_text(encoding="utf-8")
    required_prohibitions = (
        "Global four-class",
        "cross-platform OOD",
        "task-champion claims are prohibited",
        "no diagnostic",
    )
    if any(phrase not in task_text for phrase in required_prohibitions):
        raise GSE274114BaselineCampaignError("TaskSpec prohibited-claim text differs")

    output.mkdir(parents=True)
    views_root = output / "trainer_views"
    predictions_root = output / "prediction_bundles"
    views_root.mkdir()
    predictions_root.mkdir()
    view_records: dict[int, dict[str, str]] = {}
    for outer_fold in range(5):
        view = views_root / f"fold-{outer_fold}"
        prepare_fold_view(
            feature_root=feature_root,
            labels_path=labels_path,
            folds_path=folds_path,
            task_id=task_id,
            outer_fold=outer_fold,
            output=view,
        )
        view_records[outer_fold] = {
            "view_root": view.as_posix(),
            "view_sha256": sha256_file(view / "ARTIFACTS.json"),
        }

    runs = planned_runs(task_id)
    worker_runs: list[dict[str, Any]] = []
    for run in runs:
        output_path = (
            predictions_root
            / run["model_kind"]
            / f"fold-{run['outer_fold']}"
            / f"seed-{run['seed']}"
        )
        worker_runs.append(
            {
                **run,
                **view_records[run["outer_fold"]],
                "output": output_path.as_posix(),
            }
        )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        completed = list(executor.map(_fit_one, worker_runs, chunksize=1))
    completed.sort(
        key=lambda row: (row["model_kind"], row["outer_fold"], row["seed"])
    )
    if len(completed) != 60 or {
        (row["model_id"], row["outer_fold"], row["seed"]) for row in completed
    } != {
        (row["model_id"], row["outer_fold"], row["seed"]) for row in runs
    }:
        raise GSE274114BaselineCampaignError("completed prediction census differs")
    write_tsv(
        output / "prediction_inventory.tsv",
        (
            "task_id",
            "model_kind",
            "model_id",
            "outer_fold",
            "seed",
            "bundle_root",
            "artifacts_sha256",
            "bundle_id",
            "run_id",
        ),
        (
            {
                "task_id": row["task_id"],
                "model_kind": row["model_kind"],
                "model_id": row["model_id"],
                "outer_fold": row["outer_fold"],
                "seed": row["seed"],
                "bundle_root": (
                    publication_root
                    / Path(row["output"]).relative_to(output)
                ).as_posix(),
                "artifacts_sha256": row["bundle_artifacts_sha256"],
                "bundle_id": row["bundle_id"],
                "run_id": row["run_id"],
            }
            for row in completed
        ),
    )
    index_root = output / "evaluator_bundle_indices"
    index_root.mkdir()
    model_census: dict[str, int] = {}
    for model_kind in MODEL_KINDS:
        selected = [row for row in completed if row["model_kind"] == model_kind]
        expected = 25 if model_kind in STOCHASTIC_MODELS else 5
        if len(selected) != expected:
            raise GSE274114BaselineCampaignError("per-model bundle census differs")
        model_census[model_kind] = len(selected)
        write_tsv(
            index_root / f"{model_kind}.tsv",
            ("outer_fold", "seed", "bundle_root", "artifacts_sha256"),
            (
                {
                    "outer_fold": row["outer_fold"],
                    "seed": row["seed"],
                    "bundle_root": (
                        publication_root
                        / Path(row["output"]).relative_to(output)
                    ).as_posix(),
                    "artifacts_sha256": row["bundle_artifacts_sha256"],
                }
                for row in selected
            ),
        )
    receipt = {
        "schema_version": "masld-bench-gse274114-baseline-task-campaign-v1",
        "status": "passed_prediction_only_baseline_campaign",
        "task_id": task_id,
        "task_spec_sha256": task_spec_sha256,
        "feature_artifacts_sha256": feature_artifacts_sha256,
        "participants": TASKS[task_id]["participants"],
        "outer_folds": 5,
        "model_kinds": list(MODEL_KINDS),
        "fixed_stochastic_seeds": list(FIXED_STOCHASTIC_SEEDS),
        "deterministic_seed": DETERMINISTIC_SEED,
        "prediction_bundle_count": len(completed),
        "prediction_bundle_census": model_census,
        "trainer_views": 5,
        "training_labels_in_views": True,
        "held_fold_labels_in_views": False,
        "source_group_strings_in_views": False,
        "participant_covariates_in_views": False,
        "age_sex_fibrosis_nas_inputs_used": False,
        "cohort_wide_normalization_run": False,
        "preprocessing_fit_inside_training_folds_only": True,
        "prediction_outputs_contain_observed_labels": False,
        "scoring_or_evaluation_run": False,
        "global_four_class_task_run": False,
        "global_ood_task_run": False,
        "mash_vs_non_mash_task_run": False,
        "sealed_or_champion_claim_eligible": False,
        "independent_evaluator_required_after_prediction_freeze": True,
    }
    write_json_exclusive(output / "campaign_receipt.json", receipt)
    freeze_tree(
        output,
        {"artifact_class": "gse274114_prediction_only_baseline_task_campaign", **receipt},
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature-root", type=Path, required=True)
    parser.add_argument("--feature-artifacts-sha256", required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--folds", type=Path, required=True)
    parser.add_argument("--task-id", choices=sorted(TASKS), required=True)
    parser.add_argument("--task-spec", type=Path, required=True)
    parser.add_argument("--task-spec-sha256", required=True)
    parser.add_argument("--publication-root", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run_campaign(
        feature_root=args.feature_root,
        feature_artifacts_sha256=args.feature_artifacts_sha256,
        labels_path=args.labels,
        folds_path=args.folds,
        task_id=args.task_id,
        task_spec_path=args.task_spec,
        task_spec_sha256=args.task_spec_sha256,
        publication_root=args.publication_root,
        workers=args.workers,
        output=args.output,
    )
    print(json.dumps(receipt, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
