#!/usr/bin/env python3
"""Extract HyenaDNA pooled embeddings for shortened GSE281364 input windows.

The published HyenaDNA MPRA result reads a 4,096 bp genomic window mean-pooled over the
central 1,536 bp, while the assayed reporter contains only the 107 bp oligo.  This script
re-extracts the same alleles under four input windows (4096, 1024, 512, 107 bp) and three
pooling readouts, so the downstream head can be refit under the frozen protocol and the
gain measured as a function of how much non-assayed flank the model is given.

Windows are sliced from the frozen 4,096 bp fixture sequences; no genome FASTA is read.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import sys

import numpy as np
import torch
from safetensors.torch import load_file

from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    NativeContractError,
    load_config,
    read_fixture,
)


TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
INPUT_LENGTH = 4_096
VARIANT_INDEX0 = 2_048
OLIGO_START0 = 1_986
OLIGO_LENGTH = 107
WINDOWS = (4_096, 1_024, 512, 107)
POOLINGS = ("proportional", "full", "oligo")
FORWARD_ALLELES = ("REF", "ALT")
REVERSE_ALLELES = ("REF_RC", "ALT_RC")


def forward_window(length: int) -> tuple[int, int]:
    """Forward-strand slice of the frozen 4,096 bp record for one input window."""
    if length == OLIGO_LENGTH:
        return OLIGO_START0, OLIGO_START0 + OLIGO_LENGTH
    if length % 2 or not 0 < length <= INPUT_LENGTH:
        raise NativeContractError("window length differs")
    return VARIANT_INDEX0 - length // 2, VARIANT_INDEX0 + length // 2


def reverse_window(length: int) -> tuple[int, int]:
    """Same physical interval expressed on the reverse-complement record."""
    start, end = forward_window(length)
    return INPUT_LENGTH - end, INPUT_LENGTH - start


def pool_span(length: int, pooling: str, *, reverse: bool) -> tuple[int, int]:
    """Pooled position span relative to the start of the window."""
    if pooling == "full":
        return 0, length
    if pooling == "proportional":
        return int(round(0.3125 * length)), int(round(0.6875 * length))
    if pooling != "oligo":
        raise NativeContractError("pooling differs")
    window_start, _ = reverse_window(length) if reverse else forward_window(length)
    oligo_start = (
        INPUT_LENGTH - (OLIGO_START0 + OLIGO_LENGTH) if reverse else OLIGO_START0
    )
    start = oligo_start - window_start
    if start < 0 or start + OLIGO_LENGTH > length:
        raise NativeContractError("oligo footprint escapes the window")
    return start, start + OLIGO_LENGTH


def load_hyenadna(root: Path, record: dict[str, object], device: torch.device):
    """Restore the frozen HyenaDNA checkpoint exactly as the frozen extraction did."""
    code = root / record["code_path"] / "code/hyenadna"
    sys.path.insert(0, str(code.parent))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    config = HyenaConfig.from_json_file(str(code / "config.json"))
    if config.d_model != 256 or config.max_seq_len != 1_000_002:
        raise NativeContractError("HyenaDNA model identity differs")
    model = HyenaDNAForCausalLM(config)
    checkpoint = root / record["weights_path"] / record["checkpoint_relative_path"]
    state = load_file(str(checkpoint), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    aliases = {
        f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
        for layer in range(config.n_layer)
        for position in (3, 5)
    }
    shared = all(
        model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[1]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[3]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[5]
        for layer in range(config.n_layer)
    )
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys or not shared:
        raise NativeContractError("HyenaDNA checkpoint restore differs")
    del state
    model.to(device).eval()
    return model, sorted(aliases)


def hidden_states(model, sequence: str, device: torch.device) -> torch.Tensor:
    ids = torch.tensor(
        [[TOKEN_IDS[base] for base in sequence] + [1]], dtype=torch.long, device=device
    )
    with torch.inference_mode():
        return model.hyena(input_ids=ids, return_dict=True).last_hidden_state[0]


def causality_probe(model, rows, fasta, device: torch.device, count: int) -> dict[str, object]:
    """ALT - REF hidden states must be identically zero before the variant if causal."""
    worst_before = 0.0
    smallest_after = float("inf")
    for row in rows[:count]:
        prefix = row["fasta_record_prefix"]
        reference = hidden_states(model, fasta[f"{prefix}|REF"], device)
        alternative = hidden_states(model, fasta[f"{prefix}|ALT"], device)
        difference = (alternative - reference).abs()
        worst_before = max(worst_before, float(difference[:VARIANT_INDEX0].max()))
        smallest_after = min(smallest_after, float(difference[VARIANT_INDEX0:].max()))
    return {
        "elements_probed": count,
        "max_abs_delta_strictly_before_variant": worst_before,
        "min_over_elements_of_max_abs_delta_at_or_after_variant": smallest_after,
        "causal_as_predicted": worst_before < 1.0e-5,
    }


def extract(*, root: Path, config_path: Path, output: Path, seed: int, probe: int) -> dict:
    if output.exists() or config_path.is_symlink():
        raise NativeContractError("ablation extraction request differs")
    config = load_config(config_path)
    fixture = root / config["fixture"]["path"]
    rows, fasta = read_fixture(fixture)
    record = config["models"]["hyenadna"]

    # Fail before touching the GPU if the window geometry does not hold on the data.
    complement = str.maketrans("ACGT", "TGCA")
    for row in rows:
        prefix = row["fasta_record_prefix"]
        forward = fasta[f"{prefix}|REF"]
        reverse = fasta[f"{prefix}|REF_RC"]
        oligo = forward[OLIGO_START0 : OLIGO_START0 + OLIGO_LENGTH]
        start, end = reverse_window(OLIGO_LENGTH)
        if reverse[start:end] != oligo.translate(complement)[::-1]:
            raise NativeContractError("reverse-complement oligo footprint differs")
        for length in WINDOWS:
            fstart, fend = forward_window(length)
            rstart, rend = reverse_window(length)
            if (
                fend - fstart != length
                or rend - rstart != length
                or reverse[rstart:rend] != forward[fstart:fend].translate(complement)[::-1]
            ):
                raise NativeContractError("window slice is not strand consistent")

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.use_deterministic_algorithms(True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda" or torch.cuda.device_count() != 1:
        raise NativeContractError("exactly one CUDA device is required")
    model, aliases = load_hyenadna(root, record, device)

    output.mkdir(parents=True)
    probe_receipt = causality_probe(model, rows, fasta, device, probe)

    spans = {
        (length, pooling, reverse): pool_span(length, pooling, reverse=reverse)
        for length in WINDOWS
        for pooling in POOLINGS
        for reverse in (False, True)
    }
    embeddings = {
        (length, pooling): np.empty((len(rows), len(ALLELES), 256), dtype=np.float32)
        for length in WINDOWS
        for pooling in POOLINGS
    }
    forwards = 0
    for row_index, row in enumerate(rows):
        prefix = row["fasta_record_prefix"]
        for allele_index, allele in enumerate(ALLELES):
            sequence = fasta[f"{prefix}|{allele}"]
            reverse = allele in REVERSE_ALLELES
            for length in WINDOWS:
                start, end = reverse_window(length) if reverse else forward_window(length)
                hidden = hidden_states(model, sequence[start:end], device)
                forwards += 1
                if hidden.shape[0] != length + 1:
                    raise NativeContractError("token count differs")
                for pooling in POOLINGS:
                    low, high = spans[(length, pooling, reverse)]
                    value = hidden[low:high].float().mean(dim=0).cpu().numpy()
                    if value.shape != (256,) or not np.isfinite(value).all():
                        raise NativeContractError("pooled embedding differs")
                    embeddings[(length, pooling)][row_index, allele_index] = value

    # Repeat one forward pass per window to record run-to-run determinism.
    repeat_diff = 0.0
    first = rows[0]["fasta_record_prefix"]
    for length in WINDOWS:
        start, end = forward_window(length)
        sequence = fasta[f"{first}|REF"][start:end]
        low, high = spans[(length, "full", False)]
        again = hidden_states(model, sequence, device)[low:high].float().mean(dim=0).cpu().numpy()
        repeat_diff = max(
            repeat_diff, float(np.max(np.abs(again - embeddings[(length, "full")][0, 0])))
        )
    if repeat_diff > 1.0e-6:
        raise NativeContractError("repeated pooled embedding differs")

    fixture_ids = np.asarray([row["fixture_id"] for row in rows])
    group_ids = np.asarray([row["outer_locus_sequence_group_id"] for row in rows])
    folds = np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8)
    arms = []
    for length in WINDOWS:
        for pooling in POOLINGS:
            arm = f"w{length}_{pooling}"
            arm_root = output / arm
            arm_root.mkdir()
            np.savez_compressed(
                arm_root / "allele_embeddings.npz",
                fixture_ids=fixture_ids,
                outer_locus_sequence_group_ids=group_ids,
                outer_folds=folds,
                allele_order=np.asarray(ALLELES),
                embeddings=embeddings[(length, pooling)],
            )
            arms.append(
                {
                    "arm_id": arm,
                    "input_length_bp": length,
                    "pooling": pooling,
                    "forward_window": list(forward_window(length)),
                    "reverse_complement_window": list(reverse_window(length)),
                    "forward_pool_span": list(spans[(length, pooling, False)]),
                    "reverse_complement_pool_span": list(spans[(length, pooling, True)]),
                    "pool_length_bp": spans[(length, pooling, False)][1]
                    - spans[(length, pooling, False)][0],
                }
            )

    receipt = {
        "schema_version": "masld-bench-gse281364-hyenadna-window-ablation-embeddings-v1",
        "status": "pass_outcome_blind_window_ablation_extraction",
        "dataset_id": "gse281364",
        "model_id": "hyenadna",
        "arms": arms,
        "elements": len(rows),
        "alleles": list(ALLELES),
        "sequence_forwards": forwards,
        "hidden_width": 256,
        "variant_index0": VARIANT_INDEX0,
        "oligo_start0": OLIGO_START0,
        "oligo_length_bp": OLIGO_LENGTH,
        "published_arm_id": "w4096_proportional",
        "causality_probe": probe_receipt,
        "repeat_max_abs_diff": repeat_diff,
        "fixture_artifacts_sha256": config["fixture"]["artifacts_sha256"],
        "checkpoint_sha256": record["checkpoint_sha256"],
        "license": record["license"],
        "safetensors_omitted_shared_parameter_aliases": aliases,
        "checkpoint_loaded_with_safetensors": True,
        "device": str(device),
        "torch": torch.__version__,
        "seed": seed,
        "runtime_network_allowed": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_labels_read": False,
        "downstream_head_fit": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--causality-probe-elements", type=int, default=8)
    arguments = parser.parse_args()
    extract(
        root=arguments.root.resolve(strict=True),
        config_path=arguments.config.resolve(strict=True),
        output=arguments.output,
        seed=arguments.seed,
        probe=arguments.causality_probe_elements,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
