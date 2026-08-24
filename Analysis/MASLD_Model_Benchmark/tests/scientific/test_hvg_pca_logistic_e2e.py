#!/usr/bin/env python
"""Exercise the real cell baseline and evaluator on synthetic raw counts."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Mapping

import anndata
import numpy as np
import pandas as pd
from scipy import sparse


PACKAGE_ROOT = Path(__file__).resolve().parents[2]


def _load_standalone(name: str, relative_path: str) -> Any:
    path = PACKAGE_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load standalone scientific module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_adapter = _load_standalone(
    "masld_bench_standalone_hvg_pca_logistic",
    "src/masld_bench/adapters/hvg_pca_logistic.py",
)
_evaluator = _load_standalone(
    "masld_bench_standalone_cell_state_development",
    "src/masld_bench/evaluators/cell_state_development.py",
)
fit = _adapter.fit
fold_index = _adapter.fold_index
MODEL_IDS = _adapter.MODEL_IDS
predict = _adapter.predict
prepare = _adapter.prepare
evaluate = _evaluator.evaluate


CLASS_ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
DATA_ROLE = "dataset_view_data:resource_atlas_geneformer_smoke_1000_v1"


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(path: Path, role: str, media_type: str) -> dict[str, Any]:
    return {
        "path": path.resolve(strict=True).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
        "media_type": media_type,
        "role": role,
    }


def _donors_by_fold() -> dict[int, list[str]]:
    result = {fold: [] for fold in range(5)}
    candidate = 0
    while any(len(items) < 20 for items in result.values()):
        donor = f"synthetic-donor-{candidate:05d}"
        assigned = fold_index(donor, seed=20260821, outer_folds=5)
        if len(result[assigned]) < 20:
            result[assigned].append(donor)
        candidate += 1
    return result


def _write_synthetic_h5ad(path: Path) -> None:
    donors_by_fold = _donors_by_fold()
    donors: list[str] = []
    labels: list[str] = []
    cell_ids: list[str] = []
    row_indices: list[int] = []
    column_indices: list[int] = []
    values: list[int] = []
    for fold in range(5):
        for donor in donors_by_fold[fold]:
            for class_index, class_id in enumerate(CLASS_ROSTER):
                for replicate in range(2):
                    row_index = len(donors)
                    donors.append(donor)
                    labels.append(class_id)
                    cell_ids.append(
                        f"cell::{donor}::{class_id}::{replicate}"
                    )
                    row_indices.extend((row_index, row_index))
                    column_indices.extend((class_index, 10))
                    values.extend((50, 1))
    matrix = sparse.csr_matrix(
        (values, (row_indices, column_indices)),
        shape=(1000, 37533),
        dtype=np.int32,
    )
    obs = pd.DataFrame(
        {
            "donor_id": donors,
            "broad_label": labels,
            "n_counts": np.asarray(matrix.sum(axis=1)).ravel(),
        },
        index=pd.Index(cell_ids, name="cell_id"),
    )
    var = pd.DataFrame(
        {
            "ensembl_id": [f"ENSG{index + 1:011d}" for index in range(37533)]
        },
        index=pd.Index(
            [f"feature-{index + 1:05d}" for index in range(37533)],
            name="feature_id",
        ),
    )
    adata = anndata.AnnData(X=matrix, obs=obs, var=var)
    adata.write_h5ad(path, compression="gzip")


def _base_request(
    *,
    action: str,
    model_id: str,
    run_id: str,
    environment_artifact: Mapping[str, Any],
    data_artifact: Mapping[str, Any],
    prior_outputs: list[dict[str, str]],
) -> dict[str, Any]:
    run_spec = {
        "model_id": model_id,
        "task_id": "cell_state_mapping",
        "split_id": "donor_outer",
        "stage": "smoke",
        "adaptation_regime": "native_lane",
        "runtime_id": "cpu_baseline_smoke",
        "dataset_ids": ["resource_atlas_current"],
        "seed": 1103,
        "fold": 0,
        "hyperparameters": {
            "calibration_c": 1.0,
            "calibration_max_iter": 5000,
            "classifier_max_iter": 5000,
            "elastic_net_l1_ratio": 0.5,
            "hvg_mean_bins": 20,
            "join_namespace": "synthetic:cell_state_mapping:donor_outer:v1",
            "knn_neighbors": 15,
            "normalization_target_sum": 10000.0,
            "n_pca_components": 50,
            "n_top_hvg": 2000,
            "outer_folds": 5,
            "pca_svd_solver": "full",
            "regularization_c": 1.0,
            "svm_calibration_folds": 3,
            "training_weight_policy": (
                "donor_class_balanced_rescaled_to_n_cells"
            ),
        },
        "inputs": [dict(environment_artifact), dict(data_artifact)],
        "metadata": {
            "split_contract": {"seed": 20260821},
            "evaluator_parameters": {"class_roster": list(CLASS_ROSTER)},
        },
    }
    return {
        "schema_version": "masld-bench-adapter-request-v1",
        "run_id": run_id,
        "action": action,
        "run_spec": run_spec,
        "prior_action_outputs": prior_outputs,
        "fit_dataset_ids": ["resource_atlas_current"],
        "development_prediction_dataset_ids": ["resource_atlas_current"],
        "prediction_first_stress_dataset_ids": [],
        "sealed_prediction_dataset_ids": [],
        "dataset_routing": {},
        "action_dataset_ids": ["resource_atlas_current"],
        "full_input_inventory_sha256": "b" * 64,
        "action_input_inventory_sha256": "c" * 64,
        "input_matrix": "raw_unselected_counts",
        "row_id_source": "obs_names",
        "unit_id_source": "obs.donor_id",
        "label_source": "obs.broad_label",
        "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
    }


def _write_request(path: Path, request: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(request, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def run_test(*, environment_lock: Path, output_root: Path) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=False)
    h5ad_path = output_root / "synthetic_atlas_smoke.h5ad"
    _write_synthetic_h5ad(h5ad_path)
    environment_artifact = _artifact(
        environment_lock,
        "environment:cpu_baseline_smoke",
        "application/json",
    )
    data_artifact = _artifact(h5ad_path, DATA_ROLE, "application/x-hdf5")

    results: dict[str, Any] = {}
    for model_id in MODEL_IDS:
        model_root = output_root / model_id
        requests = model_root / "requests"
        actions = model_root / "adapter_actions"
        requests.mkdir(parents=True)
        actions.mkdir()
        run_id = sha256(f"synthetic::{model_id}".encode("utf-8")).hexdigest()

        prepare_output = actions / "001-prepare"
        prepare_output.mkdir()
        prepare_request = _base_request(
            action="prepare",
            model_id=model_id,
            run_id=run_id,
            environment_artifact=environment_artifact,
            data_artifact=data_artifact,
            prior_outputs=[],
        )
        prepare_request_path = requests / "001-prepare.json"
        _write_request(prepare_request_path, prepare_request)
        prepare(prepare_request_path, prepare_request, prepare_output)

        fit_output = actions / "002-fit"
        fit_output.mkdir()
        fit_request = _base_request(
            action="fit",
            model_id=model_id,
            run_id=run_id,
            environment_artifact=environment_artifact,
            data_artifact=data_artifact,
            prior_outputs=[
                {"action": "prepare", "output_path": "adapter_actions/001-prepare"}
            ],
        )
        fit_request_path = requests / "002-fit.json"
        _write_request(fit_request_path, fit_request)
        fit(fit_request_path, fit_request, fit_output)

        predict_output = actions / "003-predict"
        predict_output.mkdir()
        predict_request = _base_request(
            action="predict",
            model_id=model_id,
            run_id=run_id,
            environment_artifact=environment_artifact,
            data_artifact=data_artifact,
            prior_outputs=[
                {"action": "prepare", "output_path": "adapter_actions/001-prepare"},
                {"action": "fit", "output_path": "adapter_actions/002-fit"},
            ],
        )
        predict_request_path = requests / "003-predict.json"
        _write_request(predict_request_path, predict_request)
        predict(predict_request_path, predict_request, predict_output)

        evaluation = evaluate(
            bundle_path=predict_output / "prediction_bundle.json",
            h5ad_path=h5ad_path,
            h5ad_sha256=data_artifact["sha256"],
            environment_lock=environment_lock,
            environment_sha256=environment_artifact["sha256"],
            class_roster=CLASS_ROSTER,
            expected_model_id=model_id,
            fold=0,
            outer_folds=5,
            split_seed=20260821,
            join_namespace="synthetic:cell_state_mapping:donor_outer:v1",
            output_root=output_root / "evaluations",
        )
        result = json.loads(
            (evaluation / "evaluation.json").read_text(encoding="utf-8")
        )
        metrics = result["metrics"]
        if result["independent_unit_count"] != 20 or result["row_count"] != 200:
            raise AssertionError("synthetic held-donor inventory differs")
        if metrics["donor_class_balanced_macro_f1"] < 0.99:
            raise AssertionError(
                f"synthetic {model_id} fit is unexpectedly weak: {metrics}"
            )
        if result["champion_eligible"] is not False or result["smoke_only"] is not True:
            raise AssertionError("synthetic smoke result crossed a promotion boundary")
        results[model_id] = {
            "evaluation": evaluation.as_posix(),
            "metrics": metrics,
        }
    summary = {
        "schema_version": "masld-bench-synthetic-cell-baselines-e2e-v2",
        "models": results,
        "held_out_donors": 20,
        "held_out_cells": 200,
        "synthetic_only": True,
    }
    summary_path = output_root / "synthetic_e2e_summary.json"
    summary_path.write_text(
        json.dumps(summary, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment-lock", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    arguments = parser.parse_args()
    summary = run_test(
        environment_lock=arguments.environment_lock.resolve(strict=True),
        output_root=arguments.output_root.resolve(),
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
