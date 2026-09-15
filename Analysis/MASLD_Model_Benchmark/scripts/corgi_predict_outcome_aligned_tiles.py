#!/usr/bin/env python3
"""Run outcome-free Corgi prediction on fixed held-chromosome development tiles."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
from typing import Iterator

import numpy as np
import torch

from scripts.corgi_predict_valid_contexts import (
    CHECKPOINT_SHA256,
    CHANNEL_RC,
    LINEAGES,
    SEED,
    STATE_KEYS,
    STATE_NUMEL,
    CorgiPredictionError,
    digest,
    one_hot,
    reverse_complement_one_hot,
)
from scripts.corgi_probe_native_runtime import _load_released_model_source, _load_state


ATAC_CHANNEL = 1
BIN_BP = 64
OUTPUT_BINS = 6_144
CONTEXT_ARMS = (
    "actual_released_rank_masked",
    "actual_length_adjusted_tpm_rank_masked",
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
)
MAPPER_SCHEMA_VERSION = "masld-bench-corgi-masked-context-mapper-v1"
MAPPER_STATUS = "pass_outcome_free_mapper"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_valid_lane_contract(
    mapper: Path, mapper_outer_fold: int
) -> tuple[dict[str, object], int, int]:
    receipt = json.loads((mapper / "receipt.json").read_text(encoding="utf-8"))
    expected_valid_fold = (mapper_outer_fold + 1) % 5
    if (
        receipt.get("schema_version") != MAPPER_SCHEMA_VERSION
        or receipt.get("status") != MAPPER_STATUS
        or receipt.get("outer_fold") != mapper_outer_fold
        or receipt.get("valid_fold") != expected_valid_fold
        or receipt.get("seed") != SEED
        or receipt.get("held_ATAC_or_other_outcomes_used") is not False
        or receipt.get("test_or_sealed_outcomes_read") is not False
    ):
        raise CorgiPredictionError("mapper valid-lane receipt differs")
    # The current development lane uses synchronized donor/genomic outer folds.
    # Therefore the mapper's held validation fold is also the genomic validation fold.
    valid_donor_fold = int(receipt["valid_fold"])
    valid_genomic_fold = int(receipt["valid_fold"])
    return receipt, valid_donor_fold, valid_genomic_fold


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        for line in handle:
            value = line.rstrip("\n")
            if value.startswith(">"):
                if name is not None:
                    yield name, "".join(pieces).upper()
                name, pieces = value[1:], []
            elif name is None:
                raise CorgiPredictionError("tile FASTA sequence precedes header")
            else:
                pieces.append(value)
    if name is not None:
        yield name, "".join(pieces).upper()


def build_context_records(
    units: list[dict[str, str]],
    mapped: dict[str, np.ndarray],
    lineage: str,
    valid_donor_fold: int,
) -> tuple[list[dict[str, object]], np.ndarray]:
    if lineage not in LINEAGES:
        raise CorgiPredictionError("lineage is outside the Corgi context contract")
    valid = [
        index
        for index, row in enumerate(units)
        if row["outer_role"] == "valid" and row["lineage_id"] == lineage
    ]
    if not valid:
        raise CorgiPredictionError("lineage has no valid-role context units")
    if any(int(units[index]["outer_fold"]) != valid_donor_fold for index in valid):
        raise CorgiPredictionError("valid-role context donor fold differs")
    count_contexts = mapped["released_rank_masked_neutral"]
    tpm_contexts = mapped["length_adjusted_tpm_rank_masked_neutral"]
    means = mapped["training_lineage_mean_released_rank"]
    nearest = mapped["nearest_training_unit"]
    shuffled = mapped["shuffled_held_unit"]
    lineage_index = LINEAGES.index(lineage)
    records: list[dict[str, object]] = []
    contexts: list[np.ndarray] = []
    for unit in valid:
        arms = (
            (CONTEXT_ARMS[0], count_contexts[unit], unit),
            (CONTEXT_ARMS[1], tpm_contexts[unit], unit),
            (CONTEXT_ARMS[2], means[lineage_index], -1),
            (CONTEXT_ARMS[3], count_contexts[int(nearest[unit])], int(nearest[unit])),
            (CONTEXT_ARMS[4], count_contexts[int(shuffled[unit])], int(shuffled[unit])),
        )
        for arm, context, source_unit in arms:
            records.append(
                {
                    "prediction_index": len(records),
                    "unit_index": unit,
                    "donor_id": units[unit]["donor_id"],
                    "lineage_id": lineage,
                    "outer_fold": units[unit]["outer_fold"],
                    "donor_fold": valid_donor_fold,
                    "context_arm": arm,
                    "context_source_unit": source_unit,
                    "outcome_role": "valid",
                }
            )
            contexts.append(np.asarray(context, dtype=np.float32))
    matrix = np.vstack(contexts)
    if matrix.shape != (len(records), 2_891) or not np.isfinite(matrix).all():
        raise CorgiPredictionError("lineage context matrix differs")
    return records, matrix


def overlap_weights(local_start: int, local_end: int) -> tuple[np.ndarray, np.ndarray]:
    if not 0 <= local_start < local_end <= OUTPUT_BINS * BIN_BP:
        raise CorgiPredictionError("regional interval differs")
    first = local_start // BIN_BP
    past = (local_end + BIN_BP - 1) // BIN_BP
    bins = np.arange(first, past, dtype=np.int64)
    overlap = np.minimum((bins + 1) * BIN_BP, local_end) - np.maximum(
        bins * BIN_BP, local_start
    )
    if overlap.sum() != local_end - local_start or np.any(overlap <= 0):
        raise CorgiPredictionError("regional bin-overlap weights differ")
    return bins, overlap.astype(np.float32)


def regional_mean(track: np.ndarray, local_start: int, local_end: int) -> np.ndarray:
    bins, overlap = overlap_weights(local_start, local_end)
    if track.ndim != 2 or track.shape[1] != OUTPUT_BINS:
        raise CorgiPredictionError("regional track shape differs")
    return (track[:, bins] * overlap[None, :]).sum(axis=1) / float(local_end - local_start)


def regional_softplus_sum(track: np.ndarray, local_start: int, local_end: int) -> np.ndarray:
    bins, overlap = overlap_weights(local_start, local_end)
    if track.ndim != 2 or track.shape[1] != OUTPUT_BINS:
        raise CorgiPredictionError("regional track shape differs")
    return (np.logaddexp(0.0, track[:, bins]) * (overlap[None, :] / BIN_BP)).sum(axis=1)


def regional_orientation_mean_softplus_sum(
    forward_track: np.ndarray,
    reverse_track: np.ndarray,
    local_start: int,
    local_end: int,
) -> np.ndarray:
    if forward_track.shape != reverse_track.shape:
        raise CorgiPredictionError("orientation track shapes differ")
    return (
        regional_softplus_sum(forward_track, local_start, local_end)
        + regional_softplus_sum(reverse_track, local_start, local_end)
    ) / 2.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--mapper", type=Path, required=True)
    parser.add_argument("--tiles", type=Path, required=True)
    parser.add_argument("--tile-roster", type=Path)
    parser.add_argument("--mapper-outer-fold", type=int, choices=range(5), required=True)
    parser.add_argument("--lineage", choices=LINEAGES, required=True)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if digest(args.checkpoint) != CHECKPOINT_SHA256:
        raise CorgiPredictionError("Corgi checkpoint checksum differs")
    if int(CHANNEL_RC[ATAC_CHANNEL]) != ATAC_CHANNEL:
        raise CorgiPredictionError("ATAC reverse-complement channel is not self-mapping")
    if not 1 <= args.batch_size <= 8:
        raise CorgiPredictionError("batch size is outside the validated L40S ceiling")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise CorgiPredictionError("exactly one CUDA device is required")
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)

    mapper_receipt, valid_donor_fold, valid_genomic_fold = read_valid_lane_contract(
        args.mapper, args.mapper_outer_fold
    )
    units = read_tsv(args.mapper / "units_with_roles.tsv")
    with np.load(args.mapper / "mapped_contexts.npz", allow_pickle=False) as source:
        mapped = {key: source[key] for key in source.files}
    records, contexts = build_context_records(
        units, mapped, args.lineage, valid_donor_fold
    )
    roster = args.tile_roster or args.tiles
    tile_rows = {
        row["tile_id"]: row
        for row in read_tsv(roster / "tiles.tsv")
        if int(row["genomic_fold"]) == valid_genomic_fold
    }
    windows = [
        row
        for row in read_tsv(roster / "window_map.tsv")
        if int(row["genomic_fold"]) == valid_genomic_fold
    ]
    if not tile_rows or not windows:
        raise CorgiPredictionError("held-fold tile or window roster is empty")
    by_tile: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for index, row in enumerate(windows):
        if row["tile_id"] not in tile_rows:
            raise CorgiPredictionError("window refers to a different held-fold tile")
        by_tile.setdefault(row["tile_id"], []).append((index, row))

    Corgi, _CorgiPlus, released_config = _load_released_model_source(args.source)
    config = dict(released_config)
    if config["input_trans_regulators"] != 2_891 or config["output_channels"] != 22:
        raise CorgiPredictionError("released Corgi architecture differs")
    state = _load_state(args.checkpoint)
    if len(state) != STATE_KEYS or sum(int(value.numel()) for value in state.values()) != STATE_NUMEL:
        raise CorgiPredictionError("Corgi tensor schema differs")
    if "output_head.weight" not in state or "output_head.0.weight" in state:
        raise CorgiPredictionError("regular Corgi output head differs")
    config["final_softplus"] = False
    model = Corgi(config)
    incompatible = model.load_state_dict(state, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise CorgiPredictionError("strict Corgi restore differs")
    del state
    device = torch.device("cuda:0")
    model = model.to(device=device).eval()
    score_names = (
        "forward_raw_regional_mean",
        "reverse_complement_raw_regional_mean",
        "strand_tta_raw_regional_mean",
        "strand_tta_orientation_mean_softplus_regional_sum",
    )
    scores = {
        name: np.full((len(records), len(windows)), np.nan, dtype=np.float32)
        for name in score_names
    }
    seen_tiles: set[str] = set()
    repeat_max_abs = 0.0
    torch.cuda.reset_peak_memory_stats(device)
    for tile_id, sequence in iter_fasta(args.tiles / "tiles.fa.gz"):
        if tile_id not in tile_rows:
            continue
        if tile_id in seen_tiles or len(sequence) != 524_288:
            raise CorgiPredictionError("tile FASTA identity or length differs")
        digest_sequence = sha256(sequence.encode("ascii")).hexdigest()
        if digest_sequence != tile_rows[tile_id]["sequence_sha256"]:
            raise CorgiPredictionError("tile sequence checksum differs")
        seen_tiles.add(tile_id)
        encoded_np = one_hot(sequence)
        encoded = torch.from_numpy(encoded_np).to(device=device, dtype=torch.bfloat16)
        encoded_rc = torch.from_numpy(reverse_complement_one_hot(encoded_np)).to(
            device=device, dtype=torch.bfloat16
        )
        tile_windows = by_tile[tile_id]
        for start in range(0, len(records), args.batch_size):
            end = min(start + args.batch_size, len(records))
            context = torch.from_numpy(contexts[start:end]).to(device=device, dtype=torch.bfloat16)
            sequence_batch = encoded.unsqueeze(0).expand(end - start, -1, -1)
            rc_batch = encoded_rc.unsqueeze(0).expand(end - start, -1, -1)
            with torch.no_grad(), torch.autocast(device_type="cuda", dtype=torch.bfloat16):
                forward = model(sequence_batch, context)
                reverse = model(rc_batch, context)
                if len(seen_tiles) == 1 and start == 0:
                    repeat = model(sequence_batch, context)
            if forward.shape != (end - start, 22, OUTPUT_BINS) or reverse.shape != forward.shape:
                raise CorgiPredictionError("Corgi tile output shape differs")
            if not torch.isfinite(forward).all() or not torch.isfinite(reverse).all():
                raise CorgiPredictionError("Corgi tile output is non-finite")
            if len(seen_tiles) == 1 and start == 0:
                repeat_max_abs = float((forward.float() - repeat.float()).abs().max().item())
                if repeat_max_abs > 1.0e-3:
                    raise CorgiPredictionError("repeat prediction tolerance failed")
            forward_track = forward[:, ATAC_CHANNEL].float().cpu().numpy()
            reverse_track = reverse[:, ATAC_CHANNEL].float().cpu().numpy()[:, ::-1].copy()
            tta_track = (forward_track + reverse_track) / 2.0
            for window_index, window in tile_windows:
                local_start, local_end = int(window["local_start"]), int(window["local_end"])
                scores[score_names[0]][start:end, window_index] = regional_mean(
                    forward_track, local_start, local_end
                )
                scores[score_names[1]][start:end, window_index] = regional_mean(
                    reverse_track, local_start, local_end
                )
                scores[score_names[2]][start:end, window_index] = regional_mean(
                    tta_track, local_start, local_end
                )
                scores[score_names[3]][start:end, window_index] = (
                    regional_orientation_mean_softplus_sum(
                        forward_track,
                        reverse_track,
                        local_start,
                        local_end,
                    )
                )
    if seen_tiles != set(tile_rows) or any(not np.isfinite(value).all() for value in scores.values()):
        raise CorgiPredictionError("tile coverage or regional prediction completeness differs")

    args.output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(args.output / "regional_predictions.npz", **scores)
    with (args.output / "prediction_records.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    window_fields = tuple(windows[0])
    with (args.output / "window_records.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=window_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(windows)
    receipt = {
        "schema_version": "masld-bench-corgi-outcome-aligned-prediction-v3",
        "status": "pass_outcome_free_prediction",
        "model_id": "corgi_regular",
        "dataset_id": "gse296875",
        "split_id": (
            f"donor{args.mapper_outer_fold}_genomic{args.mapper_outer_fold}"
        ),
        "outcome_role": "valid",
        "mapper_outer_fold": args.mapper_outer_fold,
        "donor_test_fold": args.mapper_outer_fold,
        "genomic_test_fold": args.mapper_outer_fold,
        "valid_donor_fold": valid_donor_fold,
        "valid_genomic_fold": valid_genomic_fold,
        "lineage_id": args.lineage,
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "valid_donors": len(records) // len(CONTEXT_ARMS),
        "context_arms": list(CONTEXT_ARMS),
        "tiles": len(tile_rows),
        "scoreable_windows": len(windows),
        "regional_score_names": list(score_names),
        "atac_channel": ATAC_CHANNEL,
        "atac_reverse_complement_channel": int(CHANNEL_RC[ATAC_CHANNEL]),
        "strand_tta": "forward_plus_bin_reversed_reverse_complement_average",
        "count_like_tta_order": "orientation_wise_softplus_then_average",
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
        "mapper_artifacts_sha256": digest(args.mapper / "ARTIFACTS.json"),
        "mapper_receipt_sha256": digest(args.mapper / "receipt.json"),
        "mapper_receipt_schema_version": mapper_receipt["schema_version"],
        "tile_contract_artifacts_sha256": digest(args.tiles.parent / "ARTIFACTS.json"),
        "tile_roster_artifacts_sha256": digest(roster.parent / "ARTIFACTS.json"),
        "predictions_sha256": digest(args.output / "regional_predictions.npz"),
        "terminal_disposition": "restricted_development_comparator_pending_license_reference_and_native_parity",
    }
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
