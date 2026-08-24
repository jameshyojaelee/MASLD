#!/usr/bin/env python3
"""Run one exact HyenaDNA one-megabase native-context fixture."""

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


class HyenaLongProbeError(ValueError):
    """Raised when the exact one-megabase fixture differs."""


def _hash_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def _hash_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                if not name or name in records:
                    raise HyenaLongProbeError("HyenaDNA native FASTA header differs")
                records[name] = []
            elif name is None:
                raise HyenaLongProbeError("HyenaDNA native sequence precedes header")
            else:
                records[name].append(line)
    return {key: "".join(value) for key, value in records.items()}


def _tokenize(sequence: str, device: torch.device) -> torch.Tensor:
    if len(sequence) != 1_000_000 or set(sequence) - set(TOKEN_IDS):
        raise HyenaLongProbeError("HyenaDNA native sequence differs")
    return torch.tensor(
        [TOKEN_IDS[base] for base in sequence] + [1], dtype=torch.long, device=device
    )


def _predict(
    model: torch.nn.Module,
    sequence: str,
    mutation_index0: int,
    pool_start0: int,
    pool_end0: int,
    device: torch.device,
) -> tuple[float, np.ndarray]:
    ids = _tokenize(sequence, device).unsqueeze(0)
    with torch.inference_mode():
        hidden = model.hyena(input_ids=ids, return_dict=True).last_hidden_state
        if tuple(hidden.shape[:2]) != (1, 1_000_001):
            raise HyenaLongProbeError("HyenaDNA native hidden-state shape differs")
        logits = model.lm_head(hidden).float()
        targets = torch.arange(mutation_index0, 1_000_000, device=device)
        log_probabilities = torch.log_softmax(logits[0, targets - 1], dim=-1)
        score = log_probabilities.gather(1, ids[0, targets].unsqueeze(1)).mean()
        pooled = hidden[0, pool_start0:pool_end0].float().mean(dim=0)
    return float(score.cpu()), pooled.cpu().numpy()


def probe(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or any(
        path.is_symlink() for path in (arguments.code_root, arguments.checkpoint, arguments.fixture)
    ):
        raise HyenaLongProbeError("HyenaDNA native request differs")
    arguments.output.mkdir(mode=0o750)
    random.seed(20260824)
    np.random.seed(20260824)
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    sys.path.insert(0, str(arguments.code_root.parent))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    config = HyenaConfig.from_json_file(str(arguments.code_root / "config.json"))
    if config.max_seq_len != 1_000_002:
        raise HyenaLongProbeError("HyenaDNA native config differs")
    model = HyenaDNAForCausalLM(config)
    state = load_file(str(arguments.checkpoint), device="cpu")
    aliases = {
        f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
        for layer in range(config.n_layer)
        for position in (3, 5)
    }
    incompatible = model.load_state_dict(state, strict=False)
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
        raise HyenaLongProbeError("HyenaDNA native checkpoint restore differs")
    del state
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise HyenaLongProbeError("HyenaDNA native probe requires CUDA")
    model.to(device).eval()
    with arguments.manifest.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["model_id"] == "hyenadna"
            and row["context_id"] == "hyenadna_native_1000000"
        ]
    if len(rows) != 3:
        raise HyenaLongProbeError("HyenaDNA native manifest differs")
    row = sorted(rows, key=lambda value: int(value["genomic_fold"]))[0]
    records = _read_fasta(arguments.fasta)
    fixture_id = row["fixture_id"]
    expected = {
        "REF": row["reference_sequence_sha256"],
        "ALT": row["alternative_sequence_sha256"],
        "REF_RC": row["reverse_complement_reference_sha256"],
        "ALT_RC": row["reverse_complement_alternative_sha256"],
    }
    sequences: dict[str, str] = {}
    predictions: dict[str, tuple[float, np.ndarray]] = {}
    for allele in ("REF", "ALT", "REF_RC", "ALT_RC"):
        name = f"{fixture_id}|hyenadna_native_1000000|{allele}"
        sequence = records.get(name)
        if sequence is None or _hash_text(sequence) != expected[allele]:
            raise HyenaLongProbeError("HyenaDNA native FASTA identity differs")
        sequences[allele] = sequence
        index0 = int(
            row[
                "reverse_complement_variant_index0"
                if allele.endswith("_RC")
                else "forward_variant_index0"
            ]
        )
        predictions[allele] = _predict(
            model,
            sequence,
            index0,
            int(row["pool_start0"]),
            int(row["pool_end0"]),
            device,
        )
    repeated = _predict(
        model,
        sequences["REF"],
        int(row["forward_variant_index0"]),
        int(row["pool_start0"]),
        int(row["pool_end0"]),
        device,
    )
    score_repeat = abs(repeated[0] - predictions["REF"][0])
    feature_repeat = float(np.max(np.abs(repeated[1] - predictions["REF"][1])))
    if score_repeat > 1e-6 or feature_repeat > 1e-5:
        raise HyenaLongProbeError("HyenaDNA native numeric repeat differs")
    forward_delta = predictions["ALT"][0] - predictions["REF"][0]
    reverse_delta = predictions["ALT_RC"][0] - predictions["REF_RC"][0]
    features = np.stack([predictions[allele][1] for allele in ("REF", "ALT", "REF_RC", "ALT_RC")])
    np.savez_compressed(
        arguments.output / "native_1m_features.npz",
        alleles=np.asarray(["REF", "ALT", "REF_RC", "ALT_RC"]),
        features=features,
        rc_averaged_allele_delta=(features[1] - features[0] + features[3] - features[2]) / 2,
    )
    receipt = {
        "schema_version": "masld-bench-hyenadna-1m-probe-v1",
        "status": "pass",
        "checkpoint_sha256": _hash_file(arguments.checkpoint),
        "fixture_id": fixture_id,
        "genomic_fold": int(row["genomic_fold"]),
        "input_length_bp": 1_000_000,
        "forward_count": 5,
        "feature_width": int(features.shape[1]),
        "forward_alt_minus_ref_mean_suffix_log_likelihood": forward_delta,
        "reverse_alt_minus_ref_mean_suffix_log_likelihood": reverse_delta,
        "rc_averaged_alt_minus_ref_mean_suffix_log_likelihood": (forward_delta + reverse_delta) / 2,
        "numeric_score_repeat_abs_difference": score_repeat,
        "numeric_feature_repeat_max_abs_difference": feature_repeat,
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
