#!/usr/bin/env python
"""Exercise the Cobolt direct decoder on one real GSE296875 outer fold."""

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
ADAPTER_PATH = PACKAGE_ROOT / "src" / "masld_bench" / "adapters" / "rna_atac_cobolt.py"
SPEC = importlib.util.spec_from_file_location(
    "masld_bench_standalone_rna_atac_cobolt_e2e", ADAPTER_PATH
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
            "model_id": "cobolt",
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
                "alpha": 5.0,
                "annealing_epochs": 3,
                "batch_size": 128,
                "early_stopping_patience": 2,
                "hidden_dims": [128, 64],
                "inference_batch_size": 1,
                "intercept_adjustment": False,
                "join_namespace": "gse296875:rna_conditioned_atac:donor_outer:smoke:v1",
                "learning_rate": 0.005,
                "max_epochs": 3,
                "n_hvg": 2000,
                "n_latent": 10,
                "n_smoke_peaks": 2000,
                "outer_folds": 5,
                "slope_adjustment": False,
                "split_seed": 20260821,
                "validation_fraction": 0.2,
            },
            "metadata": {"evaluator_parameters": {"strata": list(adapter.LINEAGES)}},
        },
        "prior_action_outputs": prior,
        "fit_dataset_ids": ["gse296875"],
        "development_prediction_dataset_ids": [],
        "prediction_first_stress_dataset_ids": [],
        "sealed_prediction_dataset_ids": [],
        "action_dataset_ids": ["gse296875"],
        "dataset_routing": {},
        "dataset_view_id": "gse296875_rna_atac_smoke_1000_v1",
        "withheld_input_roles_by_action": {"fit": [DATA_ROLE], "predict": [DATA_ROLE]},
    }


def _write_request(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def _equation_parity(upstream_source: Path) -> None:
    import torch

    spec = importlib.util.spec_from_file_location(
        "cobolt_v1_0_1_exact_model_for_adapter_parity", upstream_source
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load exact Cobolt source")
    upstream = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = upstream
    spec.loader.exec_module(upstream)
    torch.manual_seed(1103)
    exact = upstream.CoboltModel(
        in_channels=[11, 13],
        latent_dim=10,
        n_dataset=[1, 1],
        hidden_dims=[128, 64],
        alpha=5.0,
        intercept_adj=False,
        slope_adj=False,
        log=True,
    ).cuda().eval()
    project = adapter._make_model(
        input_dims=[11, 13], latent_dim=10, hidden_dims=[128, 64], alpha=5.0
    )
    project.load_state_dict(exact.state_dict(), strict=True)
    project.cuda().eval()
    query = torch.arange(33, dtype=torch.float32, device="cuda").reshape(3, 11) % 5
    with torch.no_grad():
        exact_mu, _ = exact.get_posterior([query, None], [True, False])
        exact_profile = torch.softmax(
            torch.softmax(exact_mu, dim=1) @ exact.beta[1], dim=1
        )
        project_profile = project.decoded_profile(query)
    if not torch.allclose(exact_profile, project_profile, rtol=0.0, atol=1e-7):
        raise RuntimeError("project Cobolt decoder differs from exact v1.0.1 source")


def run_test(
    *, environment_lock: Path, data_path: Path, upstream_source: Path, output_root: Path
) -> dict[str, Any]:
    _equation_parity(upstream_source)
    output_root.mkdir(parents=True, exist_ok=False)
    environment = _artifact(
        environment_lock, role=f"environment:{RUNTIME_ID}", media_type="application/json"
    )
    data = _artifact(data_path, role=DATA_ROLE, media_type="application/x-hdf5")
    requests = output_root / "requests"
    actions = output_root / "adapter_actions"
    requests.mkdir()
    actions.mkdir()
    run_id = sha256(b"rna-atac-cobolt-e2e").hexdigest()

    prepare_path = requests / "001-prepare.json"
    prepare_payload = _request(
        action="prepare", run_id=run_id, environment=environment, data=data, prior=[]
    )
    _write_request(prepare_path, prepare_payload)
    prepare_output = actions / "001-prepare"
    prepare_output.mkdir()
    adapter.prepare(prepare_path, prepare_payload, prepare_output)
    axes = json.loads((prepare_output / "axes.json").read_text())
    if (
        axes["held_atac_exported"] is not False
        or axes["query_atac_placeholder_created"] is not False
        or axes["training_atac_scale"] != "raw_nonnegative_integer_fragment_counts"
    ):
        raise RuntimeError("prepared Cobolt firewall or count scale differs")

    fit_path = requests / "002-fit.json"
    fit_payload = _request(
        action="fit",
        run_id=run_id,
        environment=environment,
        data=data,
        prior=[{"action": "prepare", "output_path": "adapter_actions/001-prepare"}],
    )
    _write_request(fit_path, fit_payload)
    fit_output = actions / "002-fit"
    fit_output.mkdir()
    adapter.fit(fit_path, fit_payload, fit_output)
    fitted = json.loads((fit_output / "fitted_model.json").read_text())
    if (
        fitted["held_atac_used_for_fit"] is not False
        or fitted["query_rna_used_for_fit"] is not False
        or fitted["query_fitted_latent_correction"] is not False
        or fitted["dataset_adjustments_enabled"] is not False
        or fitted["whole_module_pickle_saved"] is not False
    ):
        raise RuntimeError("fitted Cobolt firewall differs")

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
    _write_request(predict_path, predict_payload)
    predict_output = actions / "003-predict"
    predict_output.mkdir()
    adapter.predict(predict_path, predict_payload, predict_output)
    bundle = json.loads((predict_output / "prediction_bundle.json").read_text())
    metadata = bundle["metadata"]
    if (
        bundle["model_id"] != "cobolt"
        or bundle["n_predictions"] <= 0
        or metadata["held_atac_input_exposed"] is not False
        or metadata["query_atac_tensor_created"] is not False
        or metadata["dataset_adjustments_enabled"] is not False
        or metadata["query_fitted_latent_correction"] is not False
        or metadata["profile_fixture_passed"] is not True
    ):
        raise RuntimeError("Cobolt prediction bundle contract differs")
    profile_sums: dict[tuple[str, str], float] = {}
    with (predict_output / "predictions.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["donor_hash"], row["stratum"])
            profile_sums[key] = profile_sums.get(key, 0.0) + float(row["predicted"])
    if not profile_sums or any(abs(value - 1.0) > 1e-10 for value in profile_sums.values()):
        raise RuntimeError("Cobolt predicted profiles do not sum to one")
    result = {
        "n_predictions": bundle["n_predictions"],
        "query_donors": axes["query_donors"],
        "query_cells": axes["query_cell_count"],
        "selected_genes": axes["selected_gene_count"],
        "selected_peaks": axes["selected_peak_count"],
        "best_epoch": fitted["best_epoch"],
        "profile_count": len(profile_sums),
        "equation_parity_passed": True,
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
    parser.add_argument("--upstream-source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    print(
        json.dumps(
            run_test(
                environment_lock=arguments.environment_lock,
                data_path=arguments.data,
                upstream_source=arguments.upstream_source,
                output_root=arguments.output,
            ),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
