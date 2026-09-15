#!/usr/bin/env python3
"""Full-window L40S step and resume probe for Regular Corgi FiLM-plus-head."""

from __future__ import annotations

import argparse
import csv
import gc
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import random
import sys
from typing import Any

import numpy as np
import torch
from torch.nn import functional as F

from scripts.corgi_predict_outcome_aligned_tiles import iter_fasta, one_hot
from scripts.corgi_probe_native_runtime import (
    _load_released_model_source,
    _load_state,
    _network_namespace_record,
)


CHECKPOINT_SHA256 = "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
EXPECTED_STATE_KEYS = 177
EXPECTED_STATE_NUMEL = 195_870_252
SEED = 20260825


class CorgiFilmProbeError(RuntimeError):
    """Raised when the full-window adaptation compatibility probe differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def tensor_digest(tensor: torch.Tensor) -> str:
    value = tensor.detach().cpu().contiguous().numpy()
    result = sha256()
    result.update(str(value.dtype).encode("ascii"))
    result.update(json.dumps(list(value.shape), separators=(",", ":")).encode("ascii"))
    result.update(value.tobytes(order="C"))
    return result.hexdigest()


def load_adapter(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("corgi_film_head_adapter", path)
    if spec is None or spec.loader is None:
        raise CorgiFilmProbeError("cannot construct adapter source spec")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_fixed_inputs(
    *, tile_contract: Path, tile_roster: Path, fixture: Path
) -> tuple[np.ndarray, np.ndarray, str, str]:
    tiles = read_tsv(tile_roster / "subset/tiles.tsv")
    chosen = min(tiles, key=lambda row: sha256(f"film-probe|{row['tile_id']}".encode()).hexdigest())
    sequence: str | None = None
    for tile_id, value in iter_fasta(tile_contract / "tiles/tiles.fa.gz"):
        if tile_id == chosen["tile_id"]:
            sequence = value
            break
    if sequence is None or len(sequence) != 524_288:
        raise CorgiFilmProbeError("fixed probe sequence is absent")
    if sha256(sequence.encode("ascii")).hexdigest() != chosen["sequence_sha256"]:
        raise CorgiFilmProbeError("fixed probe sequence checksum differs")
    records = read_tsv(fixture / "fixture/context_records.tsv")
    candidates = [
        row
        for row in records
        if row["valid_fold"] == "0"
        and row["context_arm"] == "actual_released_rank_masked"
    ]
    if not candidates:
        raise CorgiFilmProbeError("fixed actual context is absent")
    chosen_context = min(candidates, key=lambda row: row["donor_hash"])
    with np.load(fixture / "fixture/contexts.npz", allow_pickle=False) as archive:
        context = np.asarray(
            archive["context"][int(chosen_context["context_index"])], dtype=np.float32
        )
    if context.shape != (2_891,) or not np.isfinite(context).all():
        raise CorgiFilmProbeError("fixed context shape or values differ")
    return one_hot(sequence), context, chosen["tile_id"], chosen_context["donor_hash"]


def load_adapted_model(
    *, source: Path, checkpoint: Path, adapter: Any, device: torch.device
) -> torch.nn.Module:
    Corgi, _CorgiPlus, released_config = _load_released_model_source(source)
    config = dict(released_config)
    config["final_softplus"] = False
    state = _load_state(checkpoint)
    if len(state) != EXPECTED_STATE_KEYS or sum(
        int(value.numel()) for value in state.values()
    ) != EXPECTED_STATE_NUMEL:
        raise CorgiFilmProbeError("released checkpoint tensor schema differs")
    model = Corgi(config)
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise CorgiFilmProbeError("released checkpoint strict restoration differs")
    del state
    adapted = adapter.CorgiFilmAtacModel(model, require_exact=True)
    del model
    gc.collect()
    return adapted.to(device=device).eval()


def poisson_multinomial_loss(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    if prediction.shape != target.shape or prediction.ndim != 3 or prediction.shape[1] != 1:
        raise CorgiFilmProbeError("probe prediction/target shape differs")
    epsilon = 1.0e-5
    prediction = prediction.float() + epsilon
    target = target.float() + epsilon
    predicted_total = prediction.sum(dim=-1)
    target_total = target.sum(dim=-1)
    probability = prediction / predicted_total.unsqueeze(-1)
    profile = -(target * torch.log(probability)).sum(dim=-1) / prediction.shape[-1]
    count = F.poisson_nll_loss(
        predicted_total,
        target_total,
        log_input=False,
        eps=0.0,
        reduction="none",
    ) / prediction.shape[-1]
    loss = (profile + 0.25 * count).mean()
    if not torch.isfinite(loss):
        raise CorgiFilmProbeError("probe Poisson-multinomial loss is non-finite")
    return loss


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--adapter", type=Path, required=True)
    parser.add_argument("--tile-contract", type=Path, required=True)
    parser.add_argument("--tile-roster", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise CorgiFilmProbeError("refusing to overwrite L40S probe")
    if digest(arguments.checkpoint) != CHECKPOINT_SHA256:
        raise CorgiFilmProbeError("checkpoint checksum differs")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise CorgiFilmProbeError("exactly one CUDA device is required")
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.use_deterministic_algorithms(True)
    device = torch.device("cuda:0")
    adapter = load_adapter(arguments.adapter.resolve(strict=True))
    sequence_np, context_np, tile_id, donor_hash = load_fixed_inputs(
        tile_contract=arguments.tile_contract.resolve(strict=True),
        tile_roster=arguments.tile_roster.resolve(strict=True),
        fixture=arguments.fixture.resolve(strict=True),
    )
    sequence = torch.from_numpy(sequence_np).to(device=device, dtype=torch.bfloat16).unsqueeze(0)
    context = torch.from_numpy(context_np).to(device=device, dtype=torch.bfloat16).unsqueeze(0)

    Corgi, _CorgiPlus, released_config = _load_released_model_source(
        arguments.source.resolve(strict=True)
    )
    config = dict(released_config)
    config["final_softplus"] = False
    state = _load_state(arguments.checkpoint.resolve(strict=True))
    released = Corgi(config)
    incompatible = released.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise CorgiFilmProbeError("released checkpoint strict restoration differs")
    del state
    released = released.to(device=device).eval()
    torch.cuda.reset_peak_memory_stats(device)
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        released_atac_raw = released(sequence, context)[:, 1:2].float().cpu()
    released_forward_peak = int(torch.cuda.max_memory_allocated(device))
    released = released.to("cpu")
    adapted = adapter.CorgiFilmAtacModel(released, require_exact=True).to(device=device).eval()
    del released
    gc.collect()
    torch.cuda.empty_cache()
    manifest = adapter.trainable_parameter_manifest(adapted)
    if manifest["trainable_parameter_numel"] != adapter.EXPECTED_TRAINABLE_NUMEL:
        raise CorgiFilmProbeError("exact trainable parameter surface differs")
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        adapted_atac_raw = adapted.forward_raw(sequence, context).float().cpu()
    initialization_max_abs = float((released_atac_raw - adapted_atac_raw).abs().max().item())
    if initialization_max_abs > 0.05:
        raise CorgiFilmProbeError(
            f"copied ATAC head differs from released row: {initialization_max_abs}"
        )
    del released_atac_raw, adapted_atac_raw

    representative_names = (
        "conv_0.weight",
        "conv_5.conv_layer.weight",
        "final_conv.0.conv_layer.weight",
    )
    parameters = dict(adapted.named_parameters())
    frozen_before = {name: tensor_digest(parameters[name]) for name in representative_names}
    head_before = adapted.masld_atac_head.weight.detach().cpu().clone()
    optimizer = adapter.build_optimizer(adapted)
    target_axis = torch.linspace(0.0, 8.0, 6_144, device=device, dtype=torch.float32)
    target = (0.25 + torch.sin(target_axis).square()).reshape(1, 1, -1)
    optimizer.zero_grad(set_to_none=True)
    torch.cuda.reset_peak_memory_stats(device)
    with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        prediction = adapted(sequence, context)
        loss = poisson_multinomial_loss(prediction, target)
    loss_before = float(loss.detach().cpu().item())
    loss.backward()
    gradient_norms: dict[str, float] = {}
    for prefix in adapter.TRAINABLE_PREFIXES:
        gradients = [
            parameter.grad.detach().float().norm()
            for name, parameter in adapted.named_parameters()
            if name.startswith(prefix) and parameter.requires_grad and parameter.grad is not None
        ]
        if not gradients or any(not torch.isfinite(value) for value in gradients):
            raise CorgiFilmProbeError(f"finite gradients absent for {prefix}")
        gradient_norms[prefix] = float(torch.stack(gradients).norm().cpu().item())
    if any(
        parameter.grad is not None
        for parameter in adapted.parameters()
        if not parameter.requires_grad
    ):
        raise CorgiFilmProbeError("a frozen parameter received a gradient")
    torch.nn.utils.clip_grad_norm_(
        [parameter for parameter in adapted.parameters() if parameter.requires_grad], 1.0
    )
    optimizer.step()
    optimizer.zero_grad(set_to_none=True)
    training_peak = int(torch.cuda.max_memory_allocated(device))
    if torch.equal(head_before, adapted.masld_atac_head.weight.detach().cpu()):
        raise CorgiFilmProbeError("optimizer step did not change the ATAC head")
    frozen_after = {name: tensor_digest(parameters[name]) for name in representative_names}
    if frozen_after != frozen_before:
        raise CorgiFilmProbeError("representative frozen parameters changed")
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        after_step = adapted(sequence, context).float().cpu()
    delta = adapter.adaptation_delta(adapted)
    arguments.output.mkdir(mode=0o750)
    checkpoint_path = arguments.output / "synthetic_probe_delta.pt"
    torch.save(
        {
            "delta": delta,
            "optimizer_state_dict": optimizer.state_dict(),
            "probe_only_not_benchmark_model": True,
        },
        checkpoint_path,
    )
    del optimizer, adapted, parameters, prediction, loss
    gc.collect()
    torch.cuda.empty_cache()

    resumed = load_adapted_model(
        source=arguments.source.resolve(strict=True),
        checkpoint=arguments.checkpoint.resolve(strict=True),
        adapter=adapter,
        device=device,
    )
    loaded = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    if loaded.get("probe_only_not_benchmark_model") is not True:
        raise CorgiFilmProbeError("probe-only checkpoint boundary differs")
    adapter.load_adaptation_delta(resumed, loaded["delta"])
    resumed_optimizer = adapter.build_optimizer(resumed)
    resumed_optimizer.load_state_dict(loaded["optimizer_state_dict"])
    with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
        resumed_prediction = resumed(sequence, context).float().cpu()
    resume_max_abs = float((after_step - resumed_prediction).abs().max().item())
    if resume_max_abs > 1.0e-3:
        raise CorgiFilmProbeError(f"checkpoint resume prediction differs: {resume_max_abs}")

    receipt = {
        "schema_version": "masld-bench-corgi-film-head-l40s-probe-v1",
        "status": "pass_full_window_single_step_and_resume",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "adaptation_rung": "film_plus_head",
        "seed": SEED,
        "gpu_name": torch.cuda.get_device_name(device),
        "gpu_capability": list(torch.cuda.get_device_capability(device)),
        "torch": torch.__version__,
        "torch_cuda": torch.version.cuda,
        "network": _network_namespace_record(),
        "sequence_bp": 524_288,
        "output_shape": [1, 1, 6_144],
        "tile_id": tile_id,
        "donor_hash": donor_hash,
        "context_arm": "actual_released_rank_masked",
        "trainable_parameter_numel": manifest["trainable_parameter_numel"],
        "trainable_parameter_tensors": manifest["trainable_parameter_tensors"],
        "batchnorm_running_statistics_frozen": manifest[
            "batchnorm_running_statistics_frozen_by_eval_mode"
        ],
        "released_to_adapted_initialization_max_abs": initialization_max_abs,
        "released_to_adapted_initialization_tolerance": 0.05,
        "loss": "poisson_multinomial_synthetic_probe_target",
        "loss_before_step": loss_before,
        "gradient_norms": gradient_norms,
        "frozen_parameter_gradients_present": False,
        "representative_frozen_parameters_unchanged": True,
        "released_forward_peak_memory_bytes": released_forward_peak,
        "training_step_peak_memory_bytes": training_peak,
        "checkpoint_resume_prediction_max_abs": resume_max_abs,
        "checkpoint_resume_tolerance": 1.0e-3,
        "synthetic_probe_delta_sha256": digest(checkpoint_path),
        "synthetic_probe_delta_reusable_as_benchmark_model": False,
        "development_atac_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "benchmark_metrics_computed": False,
        "production_training_performed": False,
        "production_prediction_performed": False,
        "promotion_gate_evaluated": False,
        "native_numeric_parity_established": False,
        "open_champion_eligible": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
