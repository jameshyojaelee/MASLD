#!/usr/bin/env python3
"""Run a bounded exact Caduceus 131-kb native-context fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import sys

import numpy as np
import torch
from safetensors.torch import load_file


TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}


class CaduceusLongProbeError(ValueError):
    """Raised when the exact long-context fixture differs."""


def _hash_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _hash_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise CaduceusLongProbeError("long FASTA header differs")
                records[name] = []
            elif name is None:
                raise CaduceusLongProbeError("long FASTA sequence precedes header")
            else:
                records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def _encode(sequence: str, device: torch.device) -> torch.Tensor:
    if len(sequence) != 131_072 or set(sequence) - set(TOKEN_IDS):
        raise CaduceusLongProbeError("long Caduceus sequence differs")
    return torch.tensor([TOKEN_IDS[base] for base in sequence], device=device)


def _feature(
    model: torch.nn.Module,
    sequence: str,
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> np.ndarray:
    ids = _encode(sequence, device).unsqueeze(0)
    with torch.inference_mode():
        hidden = model.caduceus(input_ids=ids, return_dict=True).last_hidden_state
        if tuple(hidden.shape) != (1, 131_072, 512):
            raise CaduceusLongProbeError("long Caduceus hidden shape differs")
        realigned = (hidden[..., :256] + torch.flip(hidden[..., 256:], dims=(-2, -1))) / 2
        value = realigned[0, pool_start0:pool_end0].float().mean(dim=0)
    return value.cpu().numpy()


def _masked_score(
    model: torch.nn.Module,
    reference: str,
    alternative: str,
    index0: int,
    device: torch.device,
) -> float:
    ref = _encode(reference, device)
    alt = _encode(alternative, device)
    if torch.nonzero(ref != alt, as_tuple=False).flatten().tolist() != [index0]:
        raise CaduceusLongProbeError("long allele difference differs")
    masked = ref.clone()
    masked[index0] = 3
    with torch.inference_mode():
        logits = model(input_ids=masked.unsqueeze(0), return_dict=True).logits
        probabilities = torch.log_softmax(logits[0, index0], dim=-1)
    return float((probabilities[alt[index0]] - probabilities[ref[index0]]).cpu())


def probe(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or any(
        path.is_symlink() for path in (arguments.code_root, arguments.checkpoint, arguments.fixture)
    ):
        raise CaduceusLongProbeError("long Caduceus request differs")
    arguments.output.mkdir(mode=0o750)
    random.seed(20260824)
    np.random.seed(20260824)
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    sys.path.insert(0, str(arguments.code_root.parent))
    from caduceus.configuration_caduceus import CaduceusConfig
    from caduceus.modeling_caduceus import CaduceusForMaskedLM

    config = CaduceusConfig.from_json_file(str(arguments.code_root / "config.json"))
    model = CaduceusForMaskedLM(config)
    state = load_file(str(arguments.checkpoint), device="cpu")
    aliases = {
        f"caduceus.backbone.layers.{layer}.mixer.submodule.mamba_rev.{projection}.weight"
        for layer in range(16)
        for projection in ("in_proj", "out_proj")
    } | {"lm_head.lm_head.weight"}
    incompatible = model.load_state_dict(state, strict=False)
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
        raise CaduceusLongProbeError("long Caduceus checkpoint restore differs")
    del state
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise CaduceusLongProbeError("long Caduceus probe requires CUDA")
    model.to(device).eval()
    with arguments.manifest.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "caduceus"
            and row["context_id"] == "caduceus_native_131072"
        ]
    if len(rows) != 3:
        raise CaduceusLongProbeError("long Caduceus manifest differs")
    row = sorted(rows, key=lambda value: int(value["genomic_fold"]))[0]
    records = _fasta(arguments.fasta)
    fixture_id = row["fixture_id"]
    expected = {
        "REF": row["reference_sequence_sha256"],
        "ALT": row["alternative_sequence_sha256"],
        "REF_RC": row["reverse_complement_reference_sha256"],
        "ALT_RC": row["reverse_complement_alternative_sha256"],
    }
    sequences: dict[str, str] = {}
    features: list[np.ndarray] = []
    names: list[str] = []
    for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
        name = f"{fixture_id}|caduceus_native_131072|{allele}"
        sequence = records.get(name)
        if sequence is None or _hash_text(sequence) != expected[allele]:
            raise CaduceusLongProbeError("long Caduceus FASTA identity differs")
        sequences[allele] = sequence
        names.append(name)
        features.append(
            _feature(
                model,
                sequence,
                int(row["pool_start0"]),
                int(row["pool_end0"]),
                device,
            )
        )
    forward = _masked_score(
        model,
        sequences["REF"],
        sequences["ALT"],
        int(row["forward_variant_index0"]),
        device,
    )
    reverse = _masked_score(
        model,
        sequences["REF_RC"],
        sequences["ALT_RC"],
        int(row["reverse_complement_variant_index0"]),
        device,
    )
    repeated = _masked_score(
        model,
        sequences["REF"],
        sequences["ALT"],
        int(row["forward_variant_index0"]),
        device,
    )
    repeat_difference = abs(repeated - forward)
    if repeat_difference > 1e-5:
        raise CaduceusLongProbeError("long Caduceus numeric repeat differs")
    feature_values = np.stack(features)
    np.savez_compressed(
        arguments.output / "native_131k_features.npz",
        names=np.asarray(names),
        features=feature_values,
        rc_averaged_allele_delta=(
            feature_values[1] - feature_values[0] + feature_values[3] - feature_values[2]
        )
        / 2,
    )
    receipt = {
        "schema_version": "masld-bench-caduceus-131k-probe-v1",
        "status": "pass",
        "checkpoint_sha256": _hash_file(arguments.checkpoint),
        "fixture_id": fixture_id,
        "genomic_fold": int(row["genomic_fold"]),
        "input_length_bp": 131_072,
        "profile_count": 4,
        "feature_width": 256,
        "forward_alt_minus_ref_masked_log_probability": forward,
        "reverse_alt_minus_ref_masked_log_probability": reverse,
        "rc_averaged_alt_minus_ref_masked_log_probability": (forward + reverse) / 2,
        "numeric_repeat_abs_difference": repeat_difference,
        "common_context_backward_probe_required_separately": True,
        "head_fit": False,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "runtime_network_allowed": False,
        "maximum_gpu_memory_bytes": torch.cuda.max_memory_allocated(),
        "gpu": torch.cuda.get_device_name(0),
        "torch": torch.__version__,
    }
    (arguments.output / "probe_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--code-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    probe(parser.parse_args())


if __name__ == "__main__":
    main()
