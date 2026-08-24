#!/usr/bin/env python3
"""Run a tensor-only synthetic RNA/ATAC MIDAS forward/backward/resume probe."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable


class MIDASProbeError(RuntimeError):
    """Raised when the exact MIDAS training primitive differs."""


def tensors(value: Any) -> Iterable[Any]:
    import torch

    if torch.is_tensor(value):
        yield value
    elif isinstance(value, dict):
        for child in value.values():
            yield from tensors(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            yield from tensors(child)


def to_device(value: Any, device: Any) -> Any:
    import torch

    if torch.is_tensor(value):
        return value.to(device)
    if isinstance(value, dict):
        return {key: to_device(child, device) for key, child in value.items()}
    if isinstance(value, list):
        return [to_device(child, device) for child in value]
    if isinstance(value, tuple):
        return tuple(to_device(child, device) for child in value)
    return value


def state_sha256(state: dict[str, Any]) -> str:
    import numpy as np

    digest = sha256()
    for key in sorted(state):
        value = state[key].detach().cpu().contiguous().numpy()
        digest.update(key.encode())
        digest.update(str(value.dtype).encode())
        digest.update(np.asarray(value.shape, dtype=np.int64).tobytes())
        digest.update(value.tobytes())
    return digest.hexdigest()


def synthetic_mudata(seed: int) -> Any:
    import anndata as ad
    import mudata as mu
    import numpy as np
    import pandas as pd
    from scipy import sparse

    rng = np.random.default_rng(seed)
    cells = [f"synthetic_{index:03d}" for index in range(64)]
    batches = np.asarray(["batch_a"] * 32 + ["batch_b"] * 32)
    rna = rng.poisson(1.2, size=(64, 20)).astype(np.float32)
    atac = rng.binomial(1, 0.12, size=(64, 30)).astype(np.float32)
    rna[:, 0] += 1
    atac[:, 0] = 1
    obs = pd.DataFrame({"batch": batches}, index=cells)
    rna_adata = ad.AnnData(
        X=sparse.csr_matrix(rna),
        obs=obs.copy(),
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(20)]),
    )
    atac_adata = ad.AnnData(
        X=sparse.csr_matrix(atac),
        obs=obs.copy(),
        var=pd.DataFrame(index=[f"peak_{index}" for index in range(30)]),
    )
    return mu.MuData({"rna": rna_adata, "atac": atac_adata})


def small_config() -> dict[str, Any]:
    from scmidas.config import load_config

    config = dict(load_config())
    config.update(
        {
            "dim_c": 4,
            "dim_u": 2,
            "dims_shared_enc": [16, 8],
            "dims_shared_dec": [8, 16],
            "dims_before_enc_atac": [8, 4],
            "dims_after_dec_atac": [4, 8],
            "dims_enc_s": [8, 4],
            "dims_dec_s": [4, 8],
            "dims_dsc": [8, 4],
            "n_iter_disc": 1,
            "num_workers": 0,
            "pin_memory": False,
            "persistent_workers": False,
            "drop": 0.0,
        }
    )
    return config


def run(output: Path, seed: int) -> dict[str, Any]:
    import numpy as np
    import torch
    import scmidas
    from scmidas import MIDAS

    if output.exists():
        raise MIDASProbeError("output exists")
    output.mkdir(parents=True, mode=0o750)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    mdata = synthetic_mudata(seed)
    MIDAS.setup_mudata(
        mdata, batch_key="batch", dims_x={"rna": [20], "atac": [10, 20]}
    )
    config = small_config()
    model = MIDAS(
        mdata,
        configs=config,
        batch_size=16,
        n_save=1000,
        save_model_path=str(output / "unused_checkpoints"),
    )
    device = torch.device("cuda")
    model.to(device)
    batch = to_device(next(iter(model.train_dataloader())), device)
    model.net.train()
    optimizer = torch.optim.AdamW(model.net.parameters(), lr=1.0e-4)
    optimizer.zero_grad(set_to_none=True)
    torch.manual_seed(seed + 1)
    forward = model.net(batch)
    forward_tensors = [value for value in tensors(forward[0]) if value.is_floating_point()]
    if not forward_tensors:
        raise MIDASProbeError("MIDAS reconstruction forward is empty")
    loss = sum(value.square().mean() for value in forward_tensors)
    loss = loss + forward[2].square().mean()
    if not torch.isfinite(loss):
        raise MIDASProbeError("MIDAS synthetic loss is non-finite")
    loss.backward()
    gradient_parameters = sum(
        parameter.grad is not None and torch.isfinite(parameter.grad).all().item()
        for parameter in model.net.parameters()
    )
    if gradient_parameters == 0:
        raise MIDASProbeError("MIDAS backward produced no finite gradients")
    optimizer.step()

    model.net.eval()
    rna_masks = {}
    if isinstance(batch.get("e"), dict) and "rna" in batch["e"]:
        rna_masks["rna"] = batch["e"]["rna"]
    rna_only = {
        "x": {"rna": batch["x"]["rna"]},
        "s": batch["s"],
        "e": rna_masks,
    }
    torch.manual_seed(seed + 2)
    with torch.no_grad():
        translated = model.net(rna_only)[0]
    if "atac" not in translated:
        raise MIDASProbeError("RNA-only MIDAS forward has no ATAC decoder output")
    atac_shapes = [list(value.shape) for value in tensors(translated["atac"])]
    if not atac_shapes:
        raise MIDASProbeError("MIDAS ATAC decoder output is empty")

    state_path = output / "tensor_state.pt"
    torch.save(model.state_dict(), state_path)
    saved_sha = state_sha256(model.state_dict())
    resumed = MIDAS(
        mdata,
        configs=config,
        batch_size=16,
        n_save=1000,
        save_model_path=str(output / "unused_checkpoints_resumed"),
    ).to(device)
    loaded = torch.load(state_path, map_location=device, weights_only=True)
    resumed.load_state_dict(loaded, strict=True)
    resumed.net.eval()
    if state_sha256(resumed.state_dict()) != saved_sha:
        raise MIDASProbeError("MIDAS strict tensor-state resume differs")
    torch.manual_seed(seed + 2)
    with torch.no_grad():
        repeated = resumed.net(rna_only)[0]
    first = [value.detach().cpu().numpy() for value in tensors(translated)]
    second = [value.detach().cpu().numpy() for value in tensors(repeated)]
    if len(first) != len(second) or any(
        not np.array_equal(left, right) for left, right in zip(first, second, strict=True)
    ):
        raise MIDASProbeError("MIDAS repeated safe-state forward differs")

    receipt = {
        "schema_version": "masld-bench-midas-synthetic-training-probe-v1",
        "status": "pass",
        "model_id": "midas",
        "exact_distribution": f"scmidas {scmidas.__version__}",
        "synthetic_only": True,
        "project_data_read": False,
        "outcomes_read": False,
        "test_outcomes_read": False,
        "modalities": ["rna", "atac"],
        "atac_chromosome_chunks": [10, 20],
        "n_cells": 64,
        "forward_passed": True,
        "backward_passed": True,
        "finite_gradient_parameters": int(gradient_parameters),
        "loss": float(loss.detach().cpu()),
        "rna_only_atac_decoder_output_shapes": atac_shapes,
        "safe_weights_only_load": True,
        "strict_state_resume": True,
        "repeated_forward_bit_identical": True,
        "state_dict_sha256": saved_sha,
        "high_level_pickle_checkpoint_load_executed": False,
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
    parser.add_argument("--seed", type=int, default=7523)
    arguments = parser.parse_args()
    run(arguments.output, arguments.seed)


if __name__ == "__main__":
    main()
