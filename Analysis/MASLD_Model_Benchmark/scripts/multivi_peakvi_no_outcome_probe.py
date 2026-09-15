#!/usr/bin/env python3
"""Synthetic L40S inclusion probe for exact PeakVI and MultiVI modules."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Sequence


class MultiomeProbeError(ValueError):
    """Raised when the synthetic observed-multiome requirement differs."""


def validate_state_row_sums(
    row_sums: Sequence[float], states: Sequence[str]
) -> list[bool]:
    """Prove upstream row-sum masks equal explicit observed/missing states."""
    if len(row_sums) != len(states) or not row_sums:
        raise MultiomeProbeError("state and row-sum axes differ")
    explicit: list[bool] = []
    for value, state in zip(row_sums, states, strict=True):
        if state not in {"observed", "structurally_missing"}:
            raise MultiomeProbeError("state is not admitted to the model tensor")
        if value < 0:
            raise MultiomeProbeError("assay row sum is negative")
        if state == "observed" and value <= 0:
            raise MultiomeProbeError("observed all-zero assay must fail admission")
        if state == "structurally_missing" and value != 0:
            raise MultiomeProbeError("structurally missing assay contains values")
        explicit.append(state == "observed")
    upstream = [value > 0 for value in row_sums]
    if explicit != upstream:
        raise MultiomeProbeError("explicit and upstream modality masks differ")
    return explicit


def _tensor_sha256(value: Any) -> str:
    tensor = value.detach().contiguous().cpu()
    digest = sha256()
    digest.update(str(tensor.dtype).encode("ascii"))
    digest.update(json.dumps(list(tensor.shape)).encode("ascii"))
    digest.update(tensor.numpy().tobytes(order="C"))
    return digest.hexdigest()


def _state_dict_sha256(state: dict[str, Any]) -> str:
    digest = sha256()
    for name in sorted(state):
        digest.update(name.encode("utf-8"))
        digest.update(_tensor_sha256(state[name]).encode("ascii"))
    return digest.hexdigest()


def _device_contract(torch: Any) -> str:
    if torch.cuda.device_count() != 1:
        raise MultiomeProbeError("exactly one GPU is required")
    name = torch.cuda.get_device_name(0)
    if "L40S" not in name:
        raise MultiomeProbeError("exactly one L40S is required")
    return name


def _common_tensors(torch: Any, rna: Any, atac: Any) -> dict[str, Any]:
    from scvi import REGISTRY_KEYS

    n_rows = rna.shape[0]
    return {
        REGISTRY_KEYS.X_KEY: rna,
        REGISTRY_KEYS.ATAC_X_KEY: atac,
        REGISTRY_KEYS.BATCH_KEY: torch.zeros((n_rows, 1), dtype=torch.long, device=rna.device),
        REGISTRY_KEYS.INDICES_KEY: torch.arange(n_rows, dtype=torch.long, device=rna.device),
        REGISTRY_KEYS.LABELS_KEY: torch.zeros((n_rows, 1), dtype=torch.long, device=rna.device),
    }


def _make_multivi(torch: Any, n_obs: int) -> Any:
    from scvi.module import MULTIVAE

    return MULTIVAE(
        n_input_regions=32,
        n_input_genes=24,
        n_input_proteins=0,
        modality_weights="equal",
        modality_penalty="Jeffreys",
        n_batch=1,
        n_obs=n_obs,
        n_labels=1,
        gene_likelihood="poisson",
        gene_dispersion="gene",
        n_hidden=16,
        n_latent=8,
        n_layers_encoder=1,
        n_layers_decoder=1,
        dropout_rate=0.0,
        region_factors=True,
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    ).cuda()


def _make_peakvi() -> Any:
    from scvi.module import PEAKVAE

    return PEAKVAE(
        n_input_regions=32,
        n_batch=1,
        n_hidden=16,
        n_latent=8,
        n_layers_encoder=1,
        n_layers_decoder=1,
        dropout_rate=0.0,
        model_depth=True,
        region_factors=True,
        use_batch_norm="none",
        use_layer_norm="both",
        latent_distribution="normal",
    ).cuda()


def _train_one_step(torch: Any, model: Any, tensors: dict[str, Any]) -> float:
    model.train()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1.0e-3)
    optimizer.zero_grad(set_to_none=True)
    _inference, _generative, losses = model(tensors=tensors, compute_loss=True)
    if not torch.isfinite(losses.loss):
        raise MultiomeProbeError("training loss is non-finite")
    losses.loss.backward()
    if not any(
        parameter.grad is not None and torch.isfinite(parameter.grad).all()
        for parameter in model.parameters()
    ):
        raise MultiomeProbeError("finite gradients are absent")
    optimizer.step()
    return float(losses.loss.detach().cpu())


def _multivi_probe(torch: Any) -> dict[str, Any]:
    generator = torch.Generator(device="cuda").manual_seed(20260824)
    rna = torch.poisson(
        torch.full((16, 24), 1.5, device="cuda"), generator=generator
    )
    atac = torch.bernoulli(
        torch.full((16, 32), 0.20, device="cuda"), generator=generator
    )
    rna[:, 0] += 1
    atac[:, 0] = 1
    observed = ["observed"] * 16
    validate_state_row_sums(rna.sum(1).cpu().tolist(), observed)
    validate_state_row_sums(atac.sum(1).cpu().tolist(), observed)
    tensors = _common_tensors(torch, rna, atac)
    model = _make_multivi(torch, n_obs=16)
    loss = _train_one_step(torch, model, tensors)
    state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
    restored = _make_multivi(torch, n_obs=16)
    restored.load_state_dict(state, strict=True)
    restored.eval()

    query_rna = torch.poisson(
        torch.full((4, 24), 1.2, device="cuda"), generator=generator
    )
    query_rna[:, 0] += 1
    query_atac = torch.zeros((4, 32), device="cuda")
    validate_state_row_sums(query_rna.sum(1).cpu().tolist(), ["observed"] * 4)
    validate_state_row_sums(
        query_atac.sum(1).cpu().tolist(), ["structurally_missing"] * 4
    )
    query = _common_tensors(torch, query_rna, query_atac)
    with torch.inference_mode():
        first_inference, first_generative = restored(
            tensors=query,
            generative_kwargs={"use_z_mean": True},
            compute_loss=False,
        )
        second_inference, second_generative = restored(
            tensors=query,
            generative_kwargs={"use_z_mean": True},
            compute_loss=False,
        )
    first = first_generative["p"]
    second = second_generative["p"]
    if (
        first.shape != (4, 32)
        or not torch.isfinite(first).all()
        or torch.any(first < 0)
        or torch.any(first > 1)
        or not torch.equal(first, second)
        or not torch.equal(first_inference["qz_m"], second_inference["qz_m"])
    ):
        raise MultiomeProbeError("MultiVI deterministic RNA-only profile differs")

    distractor_rna = torch.poisson(
        torch.full((3, 24), 3.0, device="cuda"), generator=generator
    )
    distractor_rna[:, 0] += 1
    expanded_rna = torch.cat((query_rna, distractor_rna), dim=0)
    expanded_atac = torch.zeros((7, 32), device="cuda")
    validate_state_row_sums(expanded_rna.sum(1).cpu().tolist(), ["observed"] * 7)
    validate_state_row_sums(
        expanded_atac.sum(1).cpu().tolist(), ["structurally_missing"] * 7
    )
    with torch.inference_mode():
        expanded_inference, expanded_generative = restored(
            tensors=_common_tensors(torch, expanded_rna, expanded_atac),
            generative_kwargs={"use_z_mean": True},
            compute_loss=False,
        )
    if not torch.equal(first, expanded_generative["p"][:4]) or not torch.equal(
        first_inference["qz_m"], expanded_inference["qz_m"][:4]
    ):
        raise MultiomeProbeError("MultiVI query composition changes existing rows")
    return {
        "training_loss": loss,
        "state_dict_sha256": _state_dict_sha256(state),
        "state_key_count": len(state),
        "query_profile_shape": list(first.shape),
        "query_profile_sha256": _tensor_sha256(first),
        "rna_query_observed": True,
        "query_atac_state": "structurally_missing",
        "query_atac_nonzero_values": int(torch.count_nonzero(query_atac).cpu()),
        "explicit_state_bridge_matches_upstream_row_sum_masks": True,
        "observed_all_zero_rows_rejected": True,
        "query_composition_invariant": True,
        "repeat_forward_bit_identical": True,
        "forward_executed": True,
        "backward_executed": True,
        "state_dict_resume_executed": True,
        "rna_conditioned_atac_runtime_candidate": True,
    }


def _peakvi_probe(torch: Any) -> dict[str, Any]:
    from scvi import REGISTRY_KEYS

    generator = torch.Generator(device="cuda").manual_seed(20260825)
    atac = torch.bernoulli(
        torch.full((16, 32), 0.20, device="cuda"), generator=generator
    )
    atac[:, 0] = 1
    validate_state_row_sums(atac.sum(1).cpu().tolist(), ["observed"] * 16)
    tensors = {
        REGISTRY_KEYS.X_KEY: atac,
        REGISTRY_KEYS.BATCH_KEY: torch.zeros((16, 1), dtype=torch.long, device="cuda"),
    }
    model = _make_peakvi()
    loss = _train_one_step(torch, model, tensors)
    state = {name: value.detach().cpu().clone() for name, value in model.state_dict().items()}
    restored = _make_peakvi()
    restored.load_state_dict(state, strict=True)
    restored.eval()
    with torch.inference_mode():
        first_inference, first_generative = restored(
            tensors=tensors,
            generative_kwargs={"use_z_mean": True},
            compute_loss=False,
        )
        second_inference, second_generative = restored(
            tensors=tensors,
            generative_kwargs={"use_z_mean": True},
            compute_loss=False,
        )
    first = first_generative["px"]
    if (
        first.shape != (16, 32)
        or not torch.isfinite(first).all()
        or not torch.equal(first, second_generative["px"])
        or not torch.equal(first_inference["qz"].loc, second_inference["qz"].loc)
    ):
        raise MultiomeProbeError("PeakVI deterministic observed-ATAC output differs")
    return {
        "training_loss": loss,
        "state_dict_sha256": _state_dict_sha256(state),
        "state_key_count": len(state),
        "observed_atac_output_shape": list(first.shape),
        "observed_atac_output_sha256": _tensor_sha256(first),
        "observed_atac_required_at_query": True,
        "rna_conditioned_atac_eligible": False,
        "sealed_rna_only_eligible": False,
        "repeat_forward_bit_identical": True,
        "forward_executed": True,
        "backward_executed": True,
        "state_dict_resume_executed": True,
    }


def run(output: Path) -> dict[str, Any]:
    import scvi
    import torch

    if output.exists():
        raise MultiomeProbeError("output already exists")
    if scvi.__version__ != "1.5.0.post1":
        raise MultiomeProbeError("scvi-tools version differs")
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    gpu = _device_contract(torch)
    multivi = _multivi_probe(torch)
    peakvi = _peakvi_probe(torch)
    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-multivi-peakvi-no-outcome-probe-v1",
        "status": "pass",
        "python": __import__("platform").python_version(),
        "scvi_tools": scvi.__version__,
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "gpu": gpu,
        "peakvi": peakvi,
        "multivi": multivi,
        "synthetic_inputs": True,
        "project_data_read": False,
        "outcomes_read": False,
        "labels_read": False,
        "same_nucleus_project_fit_executed": False,
        "gse244832_false_pairing": False,
        "gse281367_false_pairing": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(arguments.output)


if __name__ == "__main__":
    main()
