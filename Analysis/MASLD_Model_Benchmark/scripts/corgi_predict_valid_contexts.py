#!/usr/bin/env python3
"""Predict fixed development sequences with valid-donor Corgi contexts and ablations."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np
import torch

from scripts.corgi_probe_native_runtime import _load_released_model_source, _load_state


CHECKPOINT_SHA256 = "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
STATE_KEYS = 177
STATE_NUMEL = 195_870_252
SEED = 20260824
CHANNEL_RC = np.asarray([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 12, 15, 14, 17, 16, 19, 18, 20, 21])
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")


class CorgiPredictionError(RuntimeError):
    """Raised when a checkpoint, context, or prediction contract differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    output: list[tuple[str, str]] = []
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            value = line.rstrip("\n")
            if value.startswith(">"):
                if name is not None:
                    output.append((name, "".join(pieces).upper()))
                name, pieces = value[1:], []
            elif name is None:
                raise CorgiPredictionError("FASTA sequence precedes header")
            else:
                pieces.append(value)
    if name is not None:
        output.append((name, "".join(pieces).upper()))
    if len(output) != 3 or any(len(sequence) != 524_288 for _, sequence in output):
        raise CorgiPredictionError("fixed Corgi sequence fixture differs")
    return output


def one_hot(sequence: str) -> np.ndarray:
    lookup = np.zeros((256, 4), dtype=np.float32)
    for index, base in enumerate("ACGT"):
        lookup[ord(base), index] = 1.0
    encoded = lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    if encoded.shape != (524_288, 4):
        raise CorgiPredictionError("one-hot shape differs")
    return encoded


def reverse_complement_one_hot(sequence: np.ndarray) -> np.ndarray:
    return sequence[::-1, ::-1].copy()


def restore_reverse_complement(prediction: np.ndarray) -> np.ndarray:
    if prediction.ndim != 4 or prediction.shape[2:] != (22, 6_144):
        raise CorgiPredictionError("reverse-complement prediction shape differs")
    return prediction[:, :, CHANNEL_RC, ::-1]


def build_records(
    units: list[dict[str, str]],
    mapped: dict[str, np.ndarray],
) -> tuple[list[dict[str, object]], np.ndarray]:
    valid = [index for index, row in enumerate(units) if row["outer_role"] == "valid"]
    if not valid:
        raise CorgiPredictionError("valid-role unit roster is empty")
    count_contexts = mapped["released_rank_masked_neutral"]
    tpm_contexts = mapped["length_adjusted_tpm_rank_masked_neutral"]
    means = mapped["training_lineage_mean_released_rank"]
    nearest = mapped["nearest_training_unit"]
    shuffled = mapped["shuffled_held_unit"]
    records: list[dict[str, object]] = []
    contexts: list[np.ndarray] = []
    for unit in valid:
        lineage = units[unit]["lineage_id"]
        lineage_index = LINEAGES.index(lineage)
        arms = (
            ("actual_released_rank_masked", count_contexts[unit], unit),
            ("actual_length_adjusted_tpm_rank_masked", tpm_contexts[unit], unit),
            ("training_lineage_mean_released_rank", means[lineage_index], -1),
            ("nearest_training_released_rank", count_contexts[int(nearest[unit])], int(nearest[unit])),
            ("shuffled_valid_released_rank", count_contexts[int(shuffled[unit])], int(shuffled[unit])),
        )
        for arm, context, source_unit in arms:
            records.append(
                {
                    "prediction_index": len(records),
                    "unit_index": unit,
                    "donor_id": units[unit]["donor_id"],
                    "lineage_id": lineage,
                    "outer_fold": units[unit]["outer_fold"],
                    "context_arm": arm,
                    "context_source_unit": source_unit,
                    "outcome_role": "valid",
                }
            )
            contexts.append(np.asarray(context, dtype=np.float32))
    output = np.vstack(contexts)
    if output.shape != (len(records), 2_891) or not np.isfinite(output).all():
        raise CorgiPredictionError("context record matrix differs")
    return records, output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--mapper", type=Path, required=True)
    parser.add_argument("--fixture-fasta", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if digest(args.checkpoint) != CHECKPOINT_SHA256:
        raise CorgiPredictionError("Corgi checkpoint checksum differs")
    if args.batch_size < 1 or args.batch_size > 8:
        raise CorgiPredictionError("batch size is outside the validated L40S ceiling")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise CorgiPredictionError("exactly one CUDA device is required")
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    units = read_tsv(args.mapper / "units_with_roles.tsv")
    with np.load(args.mapper / "mapped_contexts.npz", allow_pickle=False) as source:
        mapped = {key: source[key] for key in source.files}
    records, contexts = build_records(units, mapped)
    fasta = read_fasta(args.fixture_fasta)
    sequence_names = [name for name, _ in fasta]
    sequences = [one_hot(sequence) for _, sequence in fasta]
    Corgi, _CorgiPlus, released_config = _load_released_model_source(args.source)
    config = dict(released_config)
    if config["input_trans_regulators"] != 2_891 or config["output_channels"] != 22:
        raise CorgiPredictionError("released Corgi architecture differs")
    state = _load_state(args.checkpoint)
    if len(state) != STATE_KEYS or sum(int(value.numel()) for value in state.values()) != STATE_NUMEL:
        raise CorgiPredictionError("Corgi tensor schema differs")
    direct = "output_head.weight" in state
    sequential = "output_head.0.weight" in state
    if not direct or sequential:
        raise CorgiPredictionError("regular Corgi output head differs")
    config["final_softplus"] = False
    model = Corgi(config)
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise CorgiPredictionError("strict Corgi restore differs")
    del state
    device = torch.device("cuda:0")
    model = model.to(device=device).eval()
    forward = np.empty((len(records), len(sequences), 22, 6_144), dtype=np.float16)
    reverse_raw = np.empty_like(forward)
    repeat_max_abs = 0.0
    torch.cuda.reset_peak_memory_stats(device)
    for sequence_index, sequence in enumerate(sequences):
        encoded = torch.from_numpy(sequence).to(device=device, dtype=torch.bfloat16)
        encoded_rc = torch.from_numpy(reverse_complement_one_hot(sequence)).to(
            device=device, dtype=torch.bfloat16
        )
        for start in range(0, len(records), args.batch_size):
            end = min(start + args.batch_size, len(records))
            context = torch.from_numpy(contexts[start:end]).to(device=device, dtype=torch.bfloat16)
            sequence_batch = encoded.unsqueeze(0).expand(end - start, -1, -1)
            rc_batch = encoded_rc.unsqueeze(0).expand(end - start, -1, -1)
            with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                prediction = model(sequence_batch, context)
                rc_prediction = model(rc_batch, context)
                if sequence_index == 0 and start == 0:
                    repeat = model(sequence_batch, context)
            if prediction.shape != (end - start, 22, 6_144) or not torch.isfinite(prediction).all():
                raise CorgiPredictionError("forward prediction shape or finiteness differs")
            if rc_prediction.shape != prediction.shape or not torch.isfinite(rc_prediction).all():
                raise CorgiPredictionError("reverse-complement prediction differs")
            if sequence_index == 0 and start == 0:
                repeat_max_abs = float((prediction.float() - repeat.float()).abs().max().item())
                if repeat_max_abs > 1.0e-3:
                    raise CorgiPredictionError("repeat prediction tolerance failed")
            forward[start:end, sequence_index] = prediction.float().cpu().numpy().astype(np.float16)
            reverse_raw[start:end, sequence_index] = rc_prediction.float().cpu().numpy().astype(np.float16)
    restored = restore_reverse_complement(reverse_raw)
    tta = ((forward.astype(np.float32) + restored.astype(np.float32)) / 2.0).astype(np.float16)
    if not np.isfinite(tta).all():
        raise CorgiPredictionError("TTA prediction is non-finite")
    args.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        args.output / "predictions.npz",
        forward=forward,
        reverse_complement_restored=restored,
        reverse_complement_average=tta,
    )
    with (args.output / "prediction_records.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    with (args.output / "sequence_records.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sequence_index", "sequence_id"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {"sequence_index": index, "sequence_id": value} for index, value in enumerate(sequence_names)
        )
    receipt = {
        "schema_version": "masld-bench-corgi-valid-context-prediction-v1",
        "status": "pass_outcome_free_prediction",
        "model_id": "corgi_regular",
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "records": len(records),
        "valid_units": len(records) // 5,
        "context_arms": 5,
        "sequence_fixtures": len(sequences),
        "output_shape": list(tta.shape),
        "output_channels": 22,
        "output_bins": 6_144,
        "strand_TTA": "forward_plus_channel_swapped_bin_reversed_RC_average",
        "channel_RC_map": CHANNEL_RC.tolist(),
        "repeat_max_abs": repeat_max_abs,
        "repeat_tolerance": 1.0e-3,
        "gpu_name": torch.cuda.get_device_name(device),
        "peak_memory_bytes": int(torch.cuda.max_memory_allocated(device)),
        "batch_size": args.batch_size,
        "held_RNA_used": True,
        "held_ATAC_or_other_outcomes_used": False,
        "test_or_sealed_features_or_labels_read": False,
        "model_fitted_or_adapted": False,
        "native_numeric_parity_established": False,
        "open_champion_eligible": False,
        "terminal_disposition": "restricted_development_comparator_pending_license_reference_and_native_parity",
        "mapper_artifacts_sha256": digest(args.mapper / "ARTIFACTS.json"),
        "predictions_sha256": digest(args.output / "predictions.npz"),
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
