#!/usr/bin/env python3
"""Compare restored gReLU Borzoi weights across two bounded implementations."""

from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import types
from typing import Any, Mapping

import numpy as np

from scripts.reconcile_borzoi_grelu_to_local_port_mapping import map_key


class BorzoiNumericParityError(RuntimeError):
    """Raised when the numeric fixture requirements or an input binding differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_tree(project_root: Path, binding: Mapping[str, Any]) -> Path:
    root = project_root / binding["path"]
    if (
        root.is_symlink()
        or not (root / "COMPLETE").is_file()
        or digest(root / "ARTIFACTS.json") != binding["artifacts_sha256"]
    ):
        raise BorzoiNumericParityError("bound artifact differs")
    return root


def validate_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("schema_version") != (
        "masld-bench-borzoi-grelu-local-numeric-parity-fixture-v1"
    ):
        raise BorzoiNumericParityError("numeric parity schema differs")
    fixture = contract["fixture"]
    history = contract["historical_implementation"]
    comparison = contract["comparison"]
    execution = contract["execution_contract"]
    claims = contract["claim_boundary"]
    expected_stages = [
        "conv_tower_output",
        "conv_skip_y0",
        "conv_skip_y1",
        "transformer_output",
        "unet_output",
        "embedding_output",
        "raw_human_head",
        "softplus_human_head",
    ]
    if (
        contract["model_id"] != "borzoi_grelu_converted_port"
        or fixture["replicate"] != 0
        or fixture["checkpoint_filename"] != "human_state_dict_rep0.h5"
        or fixture["input_length_bp"] != 1024
        or fixture["expected_transformer_bins"] != 8
        or fixture["expected_output_bins"] != 32
        or fixture["batch_size"] != 1
        or fixture["device"] != "cpu"
        or fixture["dtype"] != "float32"
        or fixture["threads"] != 8
        or history["source_revision_is_checkpoint_bound"] is not False
        or history["unpinned_enformer_dependency_reimplemented_only_for_imported_primitives"]
        is not True
        or history["flash_attention_allowed"] is not False
        or comparison["stages"] != expected_stages
        or comparison["absolute_tolerance"] != 1e-6
        or comparison["relative_tolerance"] != 1e-5
        or comparison["all_elements_must_be_close"] is not True
        or comparison["completion_records_pass_or_fail"] is not True
        or comparison["parity_failure_blocks_local_port_execution"] is not True
        or comparison["parity_success_does_not_establish_native_borzoi_parity"]
        is not True
        or execution["torch_load_weights_only"] is not True
        or execution["torch_load_unrestricted"] is not False
        or execution["checkpoint_export_allowed"] is not False
        or execution["prediction_export_allowed"] is not False
        or execution["biological_data_allowed"] is not False
        or execution["sealed_sources_allowed"] is not False
        or execution["global_frozen_census_mutation_allowed"] is not False
        or claims["native_borzoi_substitute"] is not False
        or claims["open_champion_eligible"] is not False
        or claims["external_champion_claim_allowed"] is not False
        or claims["universal_claim_allowed"] is not False
        or claims["historical_source_is_checkpoint_bound"] is not False
        or claims["biological_outcomes_accessed"] is not False
        or claims["sealed_sources_accessed"] is not False
    ):
        raise BorzoiNumericParityError("numeric parity contract opened or drifted")


def exponential_linspace_int(
    start: int, end: int, num: int, divisible_by: int = 1
) -> list[int]:
    """Minimal imported enformer-pytorch primitive used by historical gReLU."""

    def round_divisible(value: float) -> int:
        return int(round(value / divisible_by) * divisible_by)

    base = math.exp(math.log(end / start) / (num - 1))
    return [round_divisible(start * base**index) for index in range(num)]


def relative_shift(value):
    """Minimal imported enformer-pytorch relative-position primitive."""

    import torch

    padding = torch.zeros_like(value[..., :1])
    shifted = torch.cat((padding, value), dim=-1)
    _, heads, first, second = shifted.shape
    shifted = shifted.reshape(-1, heads, second, first)
    shifted = shifted[:, :, 1:, :]
    shifted = shifted.reshape(-1, heads, first, second - 1)
    return shifted[..., : ((second + 1) // 2)]


def install_dependency_shim() -> None:
    """Install only the four enformer primitives imported by frozen sources."""

    import torch
    from torch import nn
    import torch.nn.functional as functional

    package = types.ModuleType("enformer_pytorch")
    package.__path__ = []
    module = types.ModuleType("enformer_pytorch.modeling_enformer")

    class GELU(nn.Module):
        def forward(self, value):
            return torch.sigmoid(1.702 * value) * value

    class AttentionPool(nn.Module):
        def __init__(self, dim, pool_size=2):
            super().__init__()
            self.pool_size = pool_size
            self.to_attn_logits = nn.Conv2d(dim, dim, 1, bias=False)
            nn.init.dirac_(self.to_attn_logits.weight)
            with torch.no_grad():
                self.to_attn_logits.weight.mul_(2)

        def _pool(self, value):
            batch, channels, length = value.shape
            return value.reshape(batch, channels, length // self.pool_size, self.pool_size)

        def forward(self, value):
            batch, _, length = value.shape
            remainder = length % self.pool_size
            if remainder:
                pad = self.pool_size - remainder
                value = functional.pad(value, (0, pad), value=0)
                mask = torch.zeros((batch, 1, length), dtype=torch.bool, device=value.device)
                mask = functional.pad(mask, (0, pad), value=True)
            pooled = self._pool(value)
            logits = self.to_attn_logits(pooled)
            if remainder:
                logits = logits.masked_fill(
                    self._pool(mask), -torch.finfo(logits.dtype).max
                )
            return (pooled * logits.softmax(dim=-1)).sum(dim=-1)

    module.exponential_linspace_int = exponential_linspace_int
    module.relative_shift = relative_shift
    module.GELU = GELU
    module.AttentionPool = AttentionPool
    package.modeling_enformer = module
    sys.modules["enformer_pytorch"] = package
    sys.modules["enformer_pytorch.modeling_enformer"] = module


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise BorzoiNumericParityError(f"cannot load historical module: {name}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_historical_modules(source_root: Path):
    install_dependency_shim()
    grelu = types.ModuleType("grelu")
    grelu.__path__ = []
    model = types.ModuleType("grelu.model")
    model.__path__ = []
    trunks = types.ModuleType("grelu.model.trunks")
    trunks.__path__ = []
    sys.modules["grelu"] = grelu
    sys.modules["grelu.model"] = model
    sys.modules["grelu.model.trunks"] = trunks
    load_module("grelu.model.position", source_root / "position.py")
    load_module("grelu.model.layers", source_root / "layers.py")
    load_module("grelu.model.blocks", source_root / "blocks.py")
    borzoi = load_module("grelu.model.trunks.borzoi", source_root / "borzoi.py")
    heads = load_module("grelu.model.heads", source_root / "heads.py")
    return borzoi, heads


def build_historical_model(source_root: Path):
    import torch
    from torch import nn

    borzoi, heads = load_historical_modules(source_root)

    class HistoricalBorzoi(nn.Module):
        def __init__(self):
            super().__init__()
            self.embedding = borzoi.BorzoiTrunk(
                stem_channels=512,
                stem_kernel_size=15,
                init_channels=608,
                n_conv=7,
                kernel_size=5,
                channels=1536,
                n_transformers=8,
                key_len=64,
                value_len=192,
                pos_dropout=0.0,
                attn_dropout=0.0,
                n_heads=8,
                n_pos_features=32,
                crop_len=0,
                flash_attn=False,
            )
            self.head = heads.ConvHead(
                n_tasks=7611,
                in_channels=1920,
                norm=False,
                act_func=None,
                pool_func=None,
            )

        def forward(self, value):
            return self.head(self.embedding(value))

    torch.manual_seed(0)
    return HistoricalBorzoi()


def tensor_digest(value) -> str:
    array = value.detach().cpu().contiguous().numpy()
    return sha256(memoryview(array).cast("B")).hexdigest()


def compare_tensors(left, right, *, atol: float, rtol: float) -> dict[str, Any]:
    import torch

    if left.shape != right.shape or left.dtype != right.dtype:
        raise BorzoiNumericParityError("comparison tensor schema differs")
    left64 = left.detach().to(dtype=torch.float64, device="cpu")
    right64 = right.detach().to(dtype=torch.float64, device="cpu")
    delta = left64 - right64
    rmse = torch.sqrt(torch.mean(delta.square())).item()
    denominator = torch.sqrt(torch.mean(left64.square())).item()
    left_centered = left64.flatten() - left64.mean()
    right_centered = right64.flatten() - right64.mean()
    norm_product = torch.linalg.vector_norm(left_centered) * torch.linalg.vector_norm(
        right_centered
    )
    correlation = (
        float(torch.dot(left_centered, right_centered) / norm_product)
        if float(norm_product) > 0
        else float("nan")
    )
    close = torch.isclose(left64, right64, atol=atol, rtol=rtol, equal_nan=False)
    return {
        "shape": list(left.shape),
        "dtype": str(left.dtype),
        "left_sha256": tensor_digest(left),
        "right_sha256": tensor_digest(right),
        "bit_identical": bool(torch.equal(left64, right64)),
        "all_elements_close": bool(close.all()),
        "close_fraction": float(close.to(torch.float64).mean()),
        "max_absolute_difference": float(delta.abs().max()),
        "rmse": rmse,
        "relative_rmse": rmse / max(denominator, 1e-30),
        "pearson_correlation": correlation,
    }


def build_fixture(contract: Mapping[str, Any]):
    import torch

    fixture = contract["fixture"]
    generator = torch.Generator(device="cpu")
    generator.manual_seed(fixture["seed"])
    indices = torch.randint(
        0, 4, (fixture["input_length_bp"],), generator=generator, dtype=torch.long
    )
    one_hot = torch.nn.functional.one_hot(indices, num_classes=4).to(torch.float32)
    ambiguous = torch.randperm(
        fixture["input_length_bp"], generator=generator
    )[: fixture["ambiguous_position_count"]]
    one_hot[ambiguous] = 0
    return one_hot.transpose(0, 1).unsqueeze(0).contiguous(), ambiguous.sort().values


def run_probe(contract_path: Path, project_root: Path) -> dict[str, Any]:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    mapping_root = verify_tree(project_root, contract["semantic_mapping_artifact"])
    source_root = verify_tree(project_root, contract["source_artifact"])
    runtime = contract["runtime"]
    environment = Path(runtime["environment"])
    for key in ("python", "port_model_source"):
        path = environment / runtime[key]
        if path.is_symlink() or digest(path) != runtime[f"{key}_sha256"]:
            raise BorzoiNumericParityError(f"runtime binding differs: {key}")

    os.environ["HF_HUB_OFFLINE"] = "1"
    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    import torch
    from borzoi_pytorch import Borzoi

    fixture = contract["fixture"]
    torch.set_num_threads(fixture["threads"])
    torch.use_deterministic_algorithms(True)
    historical = build_historical_model(mapping_root / "historical_sources")
    local = Borzoi.from_hparams(
        return_center_bins_only=True,
        bins_to_return=fixture["expected_output_bins"],
        enable_mouse_head=False,
        flashed=False,
    )
    state = torch.load(
        source_root / "sources" / fixture["checkpoint_filename"],
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    historical.load_state_dict(state, strict=True, assign=True)
    target_state = local.state_dict()
    remapped = {map_key(key)[0]: value for key, value in state.items()}
    if set(remapped) != set(target_state):
        raise BorzoiNumericParityError("semantic mapping is not onto local port")
    local.load_state_dict(remapped, strict=True, assign=True)
    historical.eval()
    local.eval()
    input_tensor, ambiguous = build_fixture(contract)

    with torch.inference_mode():
        hist_conv, hist_y0, hist_y1 = historical.embedding.conv_tower(input_tensor)
        hist_transformer = historical.embedding.transformer_tower(hist_conv)
        hist_unet = historical.embedding.unet_tower(
            hist_transformer, [hist_y0, hist_y1]
        )
        hist_embedding = historical.embedding.act(
            historical.embedding.pointwise_conv(hist_unet)
        )
        hist_raw = historical.head(hist_embedding)

        local_stem = local.conv_dna(input_tensor)
        local_y1 = local.res_tower(local_stem)
        local_y0 = local.unet1(local_y1)
        local_conv = local._max_pool(local_y0)
        local_transformer = local.transformer(local_conv.permute(0, 2, 1)).permute(
            0, 2, 1
        )
        local_unet = local.upsampling_unet1(local_transformer) + local.horizontal_conv1(
            local_y0
        )
        local_unet = local.separable1(local_unet)
        local_unet = local.upsampling_unet0(local_unet) + local.horizontal_conv0(
            local_y1
        )
        local_unet = local.separable0(local_unet)
        local_cropped = local.crop(local_unet.permute(0, 2, 1)).permute(0, 2, 1)
        local_embedding = local.final_joined_convs(local_cropped)
        local_raw = local.human_head(local_embedding.float())

    if hist_conv.shape[-1] != fixture["expected_transformer_bins"]:
        raise BorzoiNumericParityError("transformer-bin cardinality differs")
    if hist_raw.shape[-1] != fixture["expected_output_bins"]:
        raise BorzoiNumericParityError("output-bin cardinality differs")
    atol = contract["comparison"]["absolute_tolerance"]
    rtol = contract["comparison"]["relative_tolerance"]
    pairs = {
        "conv_tower_output": (hist_conv, local_conv),
        "conv_skip_y0": (hist_y0, local_y0),
        "conv_skip_y1": (hist_y1, local_y1),
        "transformer_output": (hist_transformer, local_transformer),
        "unet_output": (hist_unet, local_unet),
        "embedding_output": (hist_embedding, local_embedding),
        "raw_human_head": (hist_raw, local_raw),
        "softplus_human_head": (
            torch.nn.functional.softplus(hist_raw),
            torch.nn.functional.softplus(local_raw),
        ),
    }
    comparisons = {
        name: compare_tensors(left, right, atol=atol, rtol=rtol)
        for name, (left, right) in pairs.items()
    }
    parity = all(row["all_elements_close"] for row in comparisons.values())
    terminal = (
        "bounded_numeric_parity_passed_but_native_and_checkpoint_bound_parity_"
        "remain_unproven"
        if parity
        else "bounded_numeric_parity_failed_and_local_port_execution_remains_blocked"
    )
    return {
        "schema_version": contract["schema_version"],
        "model_id": contract["model_id"],
        "torch_version": torch.__version__,
        "fixture": {
            **fixture,
            "input_sha256": tensor_digest(input_tensor),
            "ambiguous_positions_zero_based": ambiguous.tolist(),
        },
        "comparison": {
            "absolute_tolerance": atol,
            "relative_tolerance": rtol,
            "stages": comparisons,
            "all_stages_passed": parity,
        },
        "semantic_mapping_artifact_sha256": digest(
            mapping_root / "ARTIFACTS.json"
        ),
        "source_artifact_sha256": digest(source_root / "ARTIFACTS.json"),
        "weights_only_load": True,
        "unrestricted_torch_load_executed": False,
        "checkpoint_exported": False,
        "prediction_exported": False,
        "biological_data_accessed": False,
        "sealed_sources_accessed": False,
        "historical_source_is_checkpoint_bound": False,
        "native_borzoi_substitute": False,
        "open_champion_eligible": False,
        "universal_claim_allowed": False,
        "terminal_disposition": terminal,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise BorzoiNumericParityError("refusing to overwrite numeric parity output")
    receipt = run_probe(arguments.contract, arguments.project_root)
    arguments.output.mkdir(parents=True, exist_ok=False)
    (arguments.output / "numeric_parity_fixture.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
