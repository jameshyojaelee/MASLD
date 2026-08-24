#!/usr/bin/env python3
"""Run the exact converted Sei checkpoint on its outcome-free native fixture."""

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


class SeiProbeError(ValueError):
    """Raised when the native Sei fixture differs from its frozen contract."""


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(
        name, path, loader=SourceFileLoader(name, str(path))
    )
    if spec is None or spec.loader is None:
        raise SeiProbeError(f"cannot load source module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
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
                raise SeiProbeError("FASTA sequence precedes its identifier")
            else:
                pieces.append(value.upper())
    if name is not None:
        records[name] = "".join(pieces)
    if len(records) != 12 or any(len(value) != 4096 for value in records.values()):
        raise SeiProbeError("Sei FASTA census or sequence length differs")
    return records


def _one_hot(sequence: str) -> np.ndarray:
    output = np.zeros((4, len(sequence)), dtype=np.float32)
    for index, base in enumerate(sequence):
        if base in "ACGT":
            output["ACGT".index(base), index] = 1.0
        elif base != "N":
            raise SeiProbeError(f"unsupported FASTA base: {base}")
    return output


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def probe(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists():
        raise SeiProbeError("Sei probe output already exists")
    arguments.output.mkdir(mode=0o750)
    sei_module = _load_module("masld_frozen_sei_model", arguments.model_source)
    scoring_module = _load_module("masld_frozen_sei_scoring", arguments.scoring_source)
    model = sei_module.Sei(sequence_length=4096, n_genomic_features=21907)
    state = load_file(str(arguments.checkpoint), device="cpu")
    if len(state) != 38 or sum(tensor.numel() for tensor in state.values()) != 889_979_983:
        raise SeiProbeError("converted Sei tensor census differs")
    prefix = "module.model."
    if any(not key.startswith(prefix) for key in state):
        raise SeiProbeError("converted Sei state-dict prefix differs")
    normalized = {key[len(prefix) :]: value for key, value in state.items()}
    incompatible = model.load_state_dict(normalized, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise SeiProbeError("strict Sei checkpoint restore differs")
    del state, normalized

    manifest_rows: list[dict[str, str]] = []
    with arguments.manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["model_id"] == "sei":
                manifest_rows.append(row)
    if len(manifest_rows) != 3 or {int(row["genomic_fold"]) for row in manifest_rows} != {2, 3, 4}:
        raise SeiProbeError("Sei native fixture manifest differs")
    records = _read_fasta(arguments.fasta)
    expected_names = {
        f"{row['fixture_id']}|sei|{allele}{suffix}"
        for row in manifest_rows
        for allele in ("REF", "ALT")
        for suffix in ("", "_RC")
    }
    if set(records) != expected_names:
        raise SeiProbeError("Sei native FASTA identifiers differ")

    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    device = torch.device("cuda")
    model.eval().to(device)
    ordered_names = sorted(records)
    # The released module lazily creates its spline basis during the first
    # forward. Initialize that cached tensor under no_grad so a later input-
    # gradient check does not encounter an inference-mode tensor.
    warmup_input = torch.from_numpy(_one_hot(records[ordered_names[0]])).unsqueeze(0).to(device)
    with torch.no_grad():
        warmup_output = model(warmup_input)
    if warmup_output.shape != (1, 21907) or not torch.isfinite(warmup_output).all():
        raise SeiProbeError("Sei spline warmup output differs")
    outputs: list[np.ndarray] = []
    with torch.inference_mode():
        for name in ordered_names:
            value = torch.from_numpy(_one_hot(records[name])).unsqueeze(0).to(device)
            prediction = model(value)
            if prediction.shape != (1, 21907) or not torch.isfinite(prediction).all():
                raise SeiProbeError("Sei native profile output differs")
            outputs.append(prediction[0].cpu().numpy())
        repeat_input = torch.from_numpy(_one_hot(records[ordered_names[0]])).unsqueeze(0).to(device)
        repeat_a = model(repeat_input)
        repeat_b = model(repeat_input)
        deterministic_max_abs_diff = float(torch.max(torch.abs(repeat_a - repeat_b)).cpu())
    profiles = np.stack(outputs)
    if profiles.min() < 0 or profiles.max() > 1 or deterministic_max_abs_diff != 0.0:
        raise SeiProbeError("Sei probability or determinism contract differs")

    profile_by_name = dict(zip(ordered_names, profiles, strict=True))
    projection = np.load(arguments.projection, allow_pickle=False)
    histone_indices = np.load(arguments.histone_indices, allow_pickle=False)
    if projection.shape != (61, 21907) or histone_indices.shape != (10064,):
        raise SeiProbeError("Sei projection auxiliary array differs")
    variant_scores: list[np.ndarray] = []
    raw_projection_shapes: list[list[int]] = []
    for row in sorted(manifest_rows, key=lambda value: int(value["genomic_fold"])):
        prefix_name = f"{row['fixture_id']}|sei"
        reference = 0.5 * (
            profile_by_name[f"{prefix_name}|REF"]
            + profile_by_name[f"{prefix_name}|REF_RC"]
        )
        alternative = 0.5 * (
            profile_by_name[f"{prefix_name}|ALT"]
            + profile_by_name[f"{prefix_name}|ALT_RC"]
        )
        raw_projection_shapes.append(
            list(scoring_module.sc_projection(reference[None, :], projection).shape)
        )
        score = scoring_module.sc_hnorm_varianteffect(
            reference[None, :], alternative[None, :], projection, histone_indices
        )
        if score.shape != (1, 40) or not np.isfinite(score).all():
            raise SeiProbeError("Sei native variant projection differs")
        variant_scores.append(score[0])
    native_scores = np.stack(variant_scores)
    if raw_projection_shapes != [[1, 61], [1, 61], [1, 61]]:
        raise SeiProbeError("Sei raw projection shape differs")

    for parameter in model.parameters():
        parameter.requires_grad_(False)
    backward_input = torch.from_numpy(_one_hot(records[ordered_names[0]])).unsqueeze(0).to(device)
    backward_input.requires_grad_(True)
    backward_output = model(backward_input)
    backward_output[:, :16].sum().backward()
    if backward_input.grad is None or not torch.isfinite(backward_input.grad).all():
        raise SeiProbeError("Sei input-gradient backward check differs")
    input_gradient_max_abs = float(backward_input.grad.abs().max().cpu())

    profiles_path = arguments.output / "native_profiles.npz"
    scores_path = arguments.output / "native_variant_scores.npy"
    np.savez_compressed(profiles_path, names=np.asarray(ordered_names), profiles=profiles)
    np.save(scores_path, native_scores, allow_pickle=False)
    receipt = {
        "schema_version": "masld-bench-sei-native-probe-v1",
        "status": "pass",
        "checkpoint_sha256": _sha256_file(arguments.checkpoint),
        "parameter_count": 889_979_983,
        "strict_checkpoint_restore": True,
        "input_length_bp": 4096,
        "profile_shape": list(profiles.shape),
        "profile_min": float(profiles.min()),
        "profile_max": float(profiles.max()),
        "native_variant_score_shape": list(native_scores.shape),
        "native_variant_mapping": "first_40_of_61_projection_rows_histone_normalized_ALT_minus_REF_after_RC_mean",
        "deterministic_repeat_max_abs_diff": deterministic_max_abs_diff,
        "input_gradient_backward_finite": True,
        "input_gradient_max_abs": input_gradient_max_abs,
        "profiles_sha256": _sha256_file(profiles_path),
        "variant_scores_sha256": _sha256_file(scores_path),
        "torch": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0),
        "maximum_gpu_memory_bytes": torch.cuda.max_memory_allocated(),
        "model_terms": "academic_and_research_use_only_restricted_comparator",
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
    }
    (arguments.output / "probe_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-source", type=Path, required=True)
    parser.add_argument("--scoring-source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--histone-indices", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    probe(parser.parse_args())


if __name__ == "__main__":
    main()
