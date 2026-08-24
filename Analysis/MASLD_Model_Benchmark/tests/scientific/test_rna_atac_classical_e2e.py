#!/usr/bin/env python
"""Exercise all three RNA-to-ATAC baselines on one real GSE296875 fold."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Mapping


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_classical.py"
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_classical", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone scientific adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)


DATA_ROLE = "dataset_view_data:gse296875_rna_atac_smoke_1000_v1"


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(path: Path, *, role: str, media_type: str) -> dict[str, Any]:
    return {
        "path": path.resolve(strict=True).as_posix(),
        "sha256": _sha256_file(path),
        "size_bytes": path.stat().st_size,
        "role": role,
        "media_type": media_type,
    }


def _request(
    *,
    action: str,
    model_id: str,
    run_id: str,
    environment: Mapping[str, Any],
    data: Mapping[str, Any],
    prior: list[dict[str, str]],
) -> dict[str, Any]:
    inputs = [dict(environment)]
    if action == "prepare":
        inputs.append(dict(data))
    return {
        "schema_version": "masld-bench-adapter-request-v1",
        "run_id": run_id,
        "action": action,
        "run_spec": {
            "model_id": model_id,
            "task_id": "rna_conditioned_atac",
            "split_id": "donor_outer",
            "stage": "smoke",
            "adaptation_regime": "native_lane",
            "runtime_id": "cpu_baseline_smoke",
            "dataset_ids": ["gse296875"],
            "seed": 1103,
            "fold": 0,
            "inputs": inputs,
            "hyperparameters": {
                "join_namespace": "gse296875:rna_conditioned_atac:donor_outer:smoke:v1",
                "normalization_target_sum": 10000.0,
                "n_context_hvg": 2000,
                "n_smoke_peaks": 2000,
                "outer_folds": 5,
                "profile_pseudocount": 1e-8,
                "split_seed": 20260821,
            },
            "metadata": {
                "evaluator_parameters": {"strata": list(adapter.LINEAGES)}
            },
        },
        "prior_action_outputs": prior,
        "fit_dataset_ids": ["gse296875"],
        "development_prediction_dataset_ids": [],
        "prediction_first_stress_dataset_ids": [],
        "sealed_prediction_dataset_ids": [],
        "action_dataset_ids": ["gse296875"],
        "dataset_routing": {},
        "dataset_view_id": "gse296875_rna_atac_smoke_1000_v1",
        "withheld_input_roles_by_action": {
            "fit": [DATA_ROLE],
            "predict": [DATA_ROLE],
        },
    }


def _write_request(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def run_test(
    *, environment_lock: Path, data_path: Path, output_root: Path
) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=False)
    environment = _artifact(
        environment_lock,
        role="environment:cpu_baseline_smoke",
        media_type="application/json",
    )
    data = _artifact(data_path, role=DATA_ROLE, media_type="application/x-hdf5")
    results = {}
    for model_id in adapter.MODEL_IDS:
        model_root = output_root / model_id
        requests = model_root / "requests"
        actions = model_root / "adapter_actions"
        requests.mkdir(parents=True)
        actions.mkdir()
        run_id = sha256(f"rna-atac-e2e::{model_id}".encode()).hexdigest()

        prepare_path = requests / "001-prepare.json"
        prepare_payload = _request(
            action="prepare",
            model_id=model_id,
            run_id=run_id,
            environment=environment,
            data=data,
            prior=[],
        )
        _write_request(prepare_path, prepare_payload)
        prepare_output = actions / "001-prepare"
        prepare_output.mkdir()
        adapter.prepare(prepare_path, prepare_payload, prepare_output)
        exported_names = {path.name for path in prepare_output.iterdir()}
        if any("query_atac" in name or "held_atac" in name for name in exported_names):
            raise RuntimeError("prepare exported held-donor ATAC")
        axes = json.loads((prepare_output / "axes.json").read_text())
        if axes["held_atac_exported"] is not False:
            raise RuntimeError("prepare held-ATAC receipt differs")

        fit_path = requests / "002-fit.json"
        fit_payload = _request(
            action="fit",
            model_id=model_id,
            run_id=run_id,
            environment=environment,
            data=data,
            prior=[{"action": "prepare", "output_path": "adapter_actions/001-prepare"}],
        )
        if any(item["role"] == DATA_ROLE for item in fit_payload["run_spec"]["inputs"]):
            raise RuntimeError("fit request exposes held ATAC")
        _write_request(fit_path, fit_payload)
        fit_output = actions / "002-fit"
        fit_output.mkdir()
        adapter.fit(fit_path, fit_payload, fit_output)

        predict_path = requests / "003-predict.json"
        predict_payload = _request(
            action="predict",
            model_id=model_id,
            run_id=run_id,
            environment=environment,
            data=data,
            prior=[
                {"action": "prepare", "output_path": "adapter_actions/001-prepare"},
                {"action": "fit", "output_path": "adapter_actions/002-fit"},
            ],
        )
        if any(item["role"] == DATA_ROLE for item in predict_payload["run_spec"]["inputs"]):
            raise RuntimeError("predict request exposes held ATAC")
        _write_request(predict_path, predict_payload)
        predict_output = actions / "003-predict"
        predict_output.mkdir()
        adapter.predict(predict_path, predict_payload, predict_output)
        bundle = json.loads((predict_output / "prediction_bundle.json").read_text())
        if (
            bundle["model_id"] != model_id
            or bundle["n_predictions"] <= 0
            or bundle["metadata"]["held_atac_input_exposed"] is not False
        ):
            raise RuntimeError("prediction bundle contract differs")
        results[model_id] = {
            "n_predictions": bundle["n_predictions"],
            "query_donors": axes["query_donors"],
            "selected_peaks": axes["selected_peak_count"],
        }
    result_path = output_root / "result.json"
    result_path.write_text(
        json.dumps(results, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environment-lock", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(
        json.dumps(
            run_test(
                environment_lock=arguments.environment_lock,
                data_path=arguments.data,
                output_root=arguments.output,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
