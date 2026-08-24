#!/usr/bin/env python
"""Exercise the scPair direct tensor adapter on one real GSE296875 fold."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any, Mapping


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
ADAPTER_PATH = (
    PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_scpair.py"
)
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_scpair_e2e", ADAPTER_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load standalone scientific adapter: {ADAPTER_PATH}")
adapter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = adapter
SPEC.loader.exec_module(adapter)


DATA_ROLE = "dataset_view_data:gse296875_rna_atac_smoke_1000_v1"
RUNTIME_ID = "gpu_rna_atac_torch_smoke"


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
            "model_id": "scpair",
            "task_id": "rna_conditioned_atac",
            "split_id": "donor_outer",
            "stage": "smoke",
            "adaptation_regime": "native_lane",
            "runtime_id": RUNTIME_ID,
            "dataset_ids": ["gse296875"],
            "seed": 1103,
            "fold": 0,
            "inputs": inputs,
            "hyperparameters": {
                "batch_size": 130,
                "dropout_rate": 0.1,
                "early_stopping_patience": 2,
                "hidden_layers": [900, 40],
                "inference_batch_size": 1,
                "join_namespace": "gse296875:rna_conditioned_atac:donor_outer:smoke:v1",
                "learning_rate": 0.001,
                "max_epochs": 3,
                "n_hvg": 2000,
                "n_smoke_peaks": 2000,
                "outer_folds": 5,
                "profile_pseudocount": 1e-8,
                "split_seed": 20260821,
                "validation_fraction": 0.2,
                "weight_decay": 1e-9,
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
        role=f"environment:{RUNTIME_ID}",
        media_type="application/json",
    )
    data = _artifact(data_path, role=DATA_ROLE, media_type="application/x-hdf5")
    requests = output_root / "requests"
    actions = output_root / "adapter_actions"
    requests.mkdir()
    actions.mkdir()
    run_id = sha256(b"rna-atac-scpair-e2e").hexdigest()

    prepare_path = requests / "001-prepare.json"
    prepare_payload = _request(
        action="prepare", run_id=run_id, environment=environment, data=data, prior=[]
    )
    _write_request(prepare_path, prepare_payload)
    prepare_output = actions / "001-prepare"
    prepare_output.mkdir()
    adapter.prepare(prepare_path, prepare_payload, prepare_output)
    names = {path.name for path in prepare_output.iterdir()}
    if any("query_atac" in name or "held_atac" in name for name in names):
        raise RuntimeError("prepare exported held-query ATAC")
    axes = json.loads((prepare_output / "axes.json").read_text())
    if (
        axes["held_atac_exported"] is not False
        or axes["query_atac_placeholder_created"] is not False
        or axes["pairing"] != "same_nucleus_exact_row_identity"
    ):
        raise RuntimeError("prepared pairing or missingness contract differs")

    fit_path = requests / "002-fit.json"
    fit_payload = _request(
        action="fit",
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
    fitted = json.loads((fit_output / "fitted_model.json").read_text())
    if (
        fitted["held_atac_used_for_fit"] is not False
        or fitted["query_rna_used_for_fit"] is not False
        or fitted["whole_module_pickle_saved"] is not False
    ):
        raise RuntimeError("fitted-model firewall differs")

    predict_path = requests / "003-predict.json"
    predict_payload = _request(
        action="predict",
        run_id=run_id,
        environment=environment,
        data=data,
        prior=[
            {"action": "prepare", "output_path": "adapter_actions/001-prepare"},
            {"action": "fit", "output_path": "adapter_actions/002-fit"},
        ],
    )
    if any(
        item["role"] == DATA_ROLE for item in predict_payload["run_spec"]["inputs"]
    ):
        raise RuntimeError("predict request exposes held ATAC")
    _write_request(predict_path, predict_payload)
    predict_output = actions / "003-predict"
    predict_output.mkdir()
    adapter.predict(predict_path, predict_payload, predict_output)
    bundle = json.loads((predict_output / "prediction_bundle.json").read_text())
    metadata = bundle["metadata"]
    if (
        bundle["model_id"] != "scpair"
        or bundle["n_predictions"] <= 0
        or metadata["held_atac_input_exposed"] is not False
        or metadata["query_atac_tensor_created"] is not False
        or metadata["profile_fixture_passed"] is not True
        or metadata["output_px_distinct_from_output_result"] is not True
    ):
        raise RuntimeError("prediction bundle contract differs")

    profile_sums: dict[tuple[str, str], float] = {}
    with (predict_output / "predictions.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["donor_hash"], row["stratum"])
            profile_sums[key] = profile_sums.get(key, 0.0) + float(row["predicted"])
    if not profile_sums or any(abs(value - 1.0) > 1e-10 for value in profile_sums.values()):
        raise RuntimeError("predicted profiles do not sum to one")
    result = {
        "n_predictions": bundle["n_predictions"],
        "query_donors": axes["query_donors"],
        "query_cells": axes["query_cell_count"],
        "selected_genes": axes["selected_gene_count"],
        "selected_peaks": axes["selected_peak_count"],
        "best_epoch": fitted["best_epoch"],
        "profile_count": len(profile_sums),
        "profile_fixture_passed": True,
        "held_atac_exposed": False,
    }
    (output_root / "result.json").write_text(
        json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return result


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
