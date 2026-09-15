#!/usr/bin/env python3
"""Run the included MPRALegNet checkpoint on an outcome-free allele fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
from safetensors.torch import load_file
import torch
from torch import nn


COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class MPRALegNetProbeError(ValueError):
    """Raised when checkpoint execution differs from its frozen requirements."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_module(path: Path):
    loader = SourceFileLoader("masld_frozen_mpralegnet_model", str(path))
    spec = importlib.util.spec_from_file_location(loader.name, path, loader=loader)
    if spec is None or spec.loader is None:
        raise MPRALegNetProbeError("cannot load frozen MPRALegNet source")
    module = importlib.util.module_from_spec(spec)
    sys.modules[loader.name] = module
    spec.loader.exec_module(module)
    return module


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                if name is not None:
                    records[name] = "".join(pieces)
                name, pieces = value[1:], []
            elif name is None:
                raise MPRALegNetProbeError("FASTA sequence precedes identifier")
            else:
                pieces.append(value.upper())
    if name is not None:
        records[name] = "".join(pieces)
    if len(records) != 12 or any(len(sequence) != 6000 for sequence in records.values()):
        raise MPRALegNetProbeError("fixture FASTA census differs")
    return records


def _one_hot(sequence: str) -> np.ndarray:
    output = np.empty((4, len(sequence)), dtype=np.float32)
    output.fill(0.0)
    for index, base in enumerate(sequence):
        if base in "AGCT":
            output["AGCT".index(base), index] = 1.0
        elif base == "N":
            output[:, index] = 0.25
        else:
            raise MPRALegNetProbeError(f"unsupported base: {base}")
    return output


def probe(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists():
        raise MPRALegNetProbeError("MPRALegNet probe output already exists")
    arguments.output.mkdir(mode=0o750)
    config = json.loads(arguments.config.read_text())
    if config.get("activation") != "SiLU" or config.get("in_ch") != 4:
        raise MPRALegNetProbeError("MPRALegNet config differs")
    module = _load_module(arguments.model_source)
    constructor = {
        key: value
        for key, value in config.items()
        if key not in {"model_type", "activation"}
    }
    model = module.LegNet(**constructor, activation=nn.SiLU)
    state = load_file(str(arguments.checkpoint), device="cpu")
    if len(state) != 134 or sum(tensor.numel() for tensor in state.values()) != 1_330_548:
        raise MPRALegNetProbeError("MPRALegNet checkpoint tensor census differs")
    prefix = "model."
    if any(not key.startswith(prefix) for key in state):
        raise MPRALegNetProbeError("MPRALegNet state-dict prefix differs")
    normalized = {key[len(prefix) :]: value for key, value in state.items()}
    incompatible = model.load_state_dict(normalized, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise MPRALegNetProbeError("strict MPRALegNet restore differs")
    del state, normalized

    manifest_rows: list[dict[str, str]] = []
    with arguments.manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["model_id"] == "hyenadna" and row["context_id"] == "common_6000":
                manifest_rows.append(row)
    if len(manifest_rows) != 3 or {int(row["genomic_fold"]) for row in manifest_rows} != {2, 3, 4}:
        raise MPRALegNetProbeError("MPRALegNet fixture manifest differs")
    records = _read_fasta(arguments.fasta)
    fixture_sequences: dict[str, dict[str, str]] = {}
    for row in manifest_rows:
        fixture_id = row["fixture_id"]
        reference = records[f"{fixture_id}|common_6000|REF"][2885:3115]
        alternative = records[f"{fixture_id}|common_6000|ALT"][2885:3115]
        differences = [
            index
            for index, values in enumerate(zip(reference, alternative, strict=True))
            if values[0] != values[1]
        ]
        if (
            len(reference) != 230
            or differences != [115]
            or reference[115] != row["ref"]
            or alternative[115] != row["alt"]
        ):
            raise MPRALegNetProbeError("230-bp allele construction differs")
        fixture_sequences[fixture_id] = {
            "REF": reference,
            "ALT": alternative,
            "REF_RC": reference.translate(COMPLEMENT)[::-1],
            "ALT_RC": alternative.translate(COMPLEMENT)[::-1],
        }

    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    device = torch.device("cuda")
    model.eval().to(device)
    if model.training or any(module.training for module in model.modules()):
        raise MPRALegNetProbeError("MPRALegNet BatchNorm eval-mode contract differs")
    ordered_names = [
        f"{fixture_id}|{allele}"
        for fixture_id in sorted(fixture_sequences)
        for allele in ("REF", "ALT", "REF_RC", "ALT_RC")
    ]
    inputs = torch.from_numpy(
        np.stack(
            [
                _one_hot(fixture_sequences[name.split("|")[0]][name.split("|")[1]])
                for name in ordered_names
            ]
        )
    ).to(device)
    running_before = {
        name: value.detach().clone()
        for name, value in model.named_buffers()
        if "running_" in name or "num_batches_tracked" in name
    }
    with torch.inference_mode():
        scores_a = model(inputs)
        scores_b = model(inputs)
    if scores_a.shape != (12,) or not torch.isfinite(scores_a).all():
        raise MPRALegNetProbeError("MPRALegNet output differs")
    deterministic_max_abs_diff = float(torch.max(torch.abs(scores_a - scores_b)).cpu())
    if deterministic_max_abs_diff != 0.0:
        raise MPRALegNetProbeError("MPRALegNet deterministic repeat differs")
    running_after = {
        name: value.detach().clone()
        for name, value in model.named_buffers()
        if "running_" in name or "num_batches_tracked" in name
    }
    if running_before.keys() != running_after.keys() or any(
        not torch.equal(running_before[name], running_after[name]) for name in running_before
    ):
        raise MPRALegNetProbeError("BatchNorm running state changed in eval mode")
    scores = scores_a.detach().cpu().numpy()
    by_name = dict(zip(ordered_names, scores, strict=True))
    allele_rows: list[dict[str, object]] = []
    for row in sorted(manifest_rows, key=lambda value: int(value["genomic_fold"])):
        fixture_id = row["fixture_id"]
        reference_score = 0.5 * (
            by_name[f"{fixture_id}|REF"] + by_name[f"{fixture_id}|REF_RC"]
        )
        alternative_score = 0.5 * (
            by_name[f"{fixture_id}|ALT"] + by_name[f"{fixture_id}|ALT_RC"]
        )
        allele_rows.append(
            {
                "fixture_id": fixture_id,
                "genomic_fold": int(row["genomic_fold"]),
                "forward_variant_index0": 115,
                "reverse_complement_variant_index0": 114,
                "reference_orientation_mean_score": float(reference_score),
                "alternative_orientation_mean_score": float(alternative_score),
                "alt_minus_ref_score": float(alternative_score - reference_score),
            }
        )

    for parameter in model.parameters():
        parameter.requires_grad_(False)
    backward_input = inputs[:1].detach().clone().requires_grad_(True)
    backward_score = model(backward_input).sum()
    backward_score.backward()
    if backward_input.grad is None or not torch.isfinite(backward_input.grad).all():
        raise MPRALegNetProbeError("MPRALegNet input-gradient check differs")
    input_gradient_max_abs = float(backward_input.grad.abs().max().cpu())

    scores_path = arguments.output / "native_scores.npz"
    np.savez_compressed(scores_path, names=np.asarray(ordered_names), scores=scores)
    allele_path = arguments.output / "allele_scores.json"
    allele_path.write_text(json.dumps(allele_rows, sort_keys=True, indent=2) + "\n")
    receipt = {
        "schema_version": "masld-bench-mpralegnet-native-probe-v1",
        "status": "pass",
        "checkpoint_sha256": _sha256_file(arguments.checkpoint),
        "checkpoint_provenance": "third_party_Hugging_Face_repackaging_author_equivalence_unresolved",
        "strict_checkpoint_restore": True,
        "parameter_and_buffer_elements": 1_330_548,
        "input_length_bp": 230,
        "channel_order": "A_G_C_T",
        "n_encoding": "0.25_each_channel",
        "forward_variant_index0": 115,
        "reverse_complement_variant_index0": 114,
        "strand_policy": "mean_forward_and_reverse_complement",
        "allele_effect_sign": "ALT_minus_REF",
        "output_shape": list(scores.shape),
        "deterministic_repeat_max_abs_diff": deterministic_max_abs_diff,
        "batchnorm_eval_mode": True,
        "batchnorm_running_state_unchanged": True,
        "input_gradient_backward_finite": True,
        "input_gradient_max_abs": input_gradient_max_abs,
        "scores_sha256": _sha256_file(scores_path),
        "allele_scores_sha256": _sha256_file(allele_path),
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0),
        "maximum_gpu_memory_bytes": torch.cuda.max_memory_allocated(),
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "standalone_champion_eligible": False,
    }
    (arguments.output / "probe_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-source", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    probe(parser.parse_args())


if __name__ == "__main__":
    main()
