#!/usr/bin/env python3
"""Run frozen Sei native variant scoring on outcome-blind GSE281364 inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path
import random
import sys

import numpy as np
from safetensors.torch import load_file
import torch


PARAMETERS = 889_979_983
OUTPUTS = 21_907
SEQUENCE_CLASSES = 40
SEED = 20260824


class SeiPredictionError(ValueError):
    """Raised when native Sei prediction differs from its frozen requirements."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(
        name, path, loader=SourceFileLoader(name, str(path))
    )
    if spec is None or spec.loader is None:
        raise SeiPredictionError(f"cannot load source module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    pieces: list[str] = []
    with gzip.open(path, "rt", encoding="ascii") as handle:
        for line in handle:
            value = line.strip()
            if value.startswith(">"):
                if name is not None:
                    if name in records:
                        raise SeiPredictionError("duplicate FASTA identifier")
                    records[name] = "".join(pieces)
                name, pieces = value[1:], []
            elif name is None:
                raise SeiPredictionError("FASTA sequence precedes identifier")
            else:
                pieces.append(value.upper())
    if name is not None:
        if name in records:
            raise SeiPredictionError("duplicate FASTA identifier")
        records[name] = "".join(pieces)
    if len(records) != 4_132 or any(len(sequence) != 4_096 for sequence in records.values()):
        raise SeiPredictionError("FASTA census or sequence length differs")
    return records


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if (
        len(rows) != 1_033
        or len({row["fixture_id"] for row in rows}) != len(rows)
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != len(rows)
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
    ):
        raise SeiPredictionError("fixture manifest differs")
    return rows


def one_hot_many(sequences: list[str]) -> np.ndarray:
    output = np.zeros((len(sequences), 4, 4_096), dtype=np.float32)
    for row, sequence in enumerate(sequences):
        for column, base in enumerate(sequence):
            if base not in "ACGT":
                raise SeiPredictionError("unsupported FASTA base")
            output[row, "ACGT".index(base), column] = 1.0
    return output


def scalarize(scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(scores, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != SEQUENCE_CLASSES or not np.isfinite(values).all():
        raise SeiPredictionError("native score matrix differs")
    indices = np.argmax(np.abs(values), axis=1)
    signed = values[np.arange(values.shape[0]), indices]
    maximum = np.abs(signed)
    mean = values.mean(axis=1)
    return indices.astype(np.int64), signed, maximum, mean


def predict(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or not 1 <= arguments.batch_variants <= 16:
        raise SeiPredictionError("output or batch contract differs")
    manifest = read_manifest(arguments.manifest)
    records = read_fasta(arguments.fasta)
    expected_names = {
        f"{row['fixture_id']}|{arm}"
        for row in manifest
        for arm in ("REF", "ALT", "REF_RC", "ALT_RC")
    }
    if set(records) != expected_names:
        raise SeiPredictionError("FASTA identifiers differ from manifest")
    model_module = load_module("masld_gse281364_sei_model", arguments.model_source)
    scoring_module = load_module("masld_gse281364_sei_scoring", arguments.scoring_source)
    model = model_module.Sei(sequence_length=4096, n_genomic_features=OUTPUTS)
    state = load_file(str(arguments.checkpoint), device="cpu")
    if len(state) != 38 or sum(tensor.numel() for tensor in state.values()) != PARAMETERS:
        raise SeiPredictionError("converted checkpoint tensor census differs")
    prefix = "module.model."
    if any(not key.startswith(prefix) for key in state):
        raise SeiPredictionError("converted checkpoint prefix differs")
    normalized = {key[len(prefix) :]: value for key, value in state.items()}
    incompatible = model.load_state_dict(normalized, strict=True)
    if incompatible.missing_keys or incompatible.unexpected_keys:
        raise SeiPredictionError("strict checkpoint restore differs")
    del state, normalized
    projection = np.load(arguments.projection, allow_pickle=False)
    histone_indices = np.load(arguments.histone_indices, allow_pickle=False)
    if projection.shape != (61, OUTPUTS) or histone_indices.shape != (10_064,):
        raise SeiPredictionError("native scoring auxiliaries differ")

    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.use_deterministic_algorithms(True)
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise SeiPredictionError("exactly one CUDA device is required")
    device = torch.device("cuda:0")
    model.eval().to(device)
    first = records[f"{manifest[0]['fixture_id']}|REF"]
    with torch.no_grad():
        warmup = model(torch.from_numpy(one_hot_many([first])).to(device))
    if warmup.shape != (1, OUTPUTS) or not torch.isfinite(warmup).all():
        raise SeiPredictionError("warmup output differs")
    all_scores: list[np.ndarray] = []
    repeat_max_abs = -1.0
    torch.cuda.reset_peak_memory_stats(device)
    for start in range(0, len(manifest), arguments.batch_variants):
        chunk = manifest[start : start + arguments.batch_variants]
        sequences = [
            records[f"{row['fixture_id']}|{arm}"]
            for row in chunk
            for arm in ("REF", "ALT", "REF_RC", "ALT_RC")
        ]
        tensor = torch.from_numpy(one_hot_many(sequences)).to(device)
        with torch.inference_mode():
            profiles = model(tensor)
            if start == 0:
                repeated = model(tensor)
                repeat_max_abs = float((profiles - repeated).abs().max().cpu())
        if profiles.shape != (4 * len(chunk), OUTPUTS) or not torch.isfinite(profiles).all():
            raise SeiPredictionError("native profile batch differs")
        values = profiles.float().cpu().numpy().reshape(len(chunk), 4, OUTPUTS)
        reference = 0.5 * (values[:, 0] + values[:, 2])
        alternative = 0.5 * (values[:, 1] + values[:, 3])
        score = scoring_module.sc_hnorm_varianteffect(
            reference, alternative, projection, histone_indices
        )
        if score.shape != (len(chunk), SEQUENCE_CLASSES) or not np.isfinite(score).all():
            raise SeiPredictionError("native variant score batch differs")
        all_scores.append(np.asarray(score, dtype=np.float32))
        print(json.dumps({"completed": start + len(chunk), "total": len(manifest)}), flush=True)
    scores = np.vstack(all_scores)
    if scores.shape != (len(manifest), SEQUENCE_CLASSES) or repeat_max_abs != 0.0:
        raise SeiPredictionError("complete score or determinism contract differs")
    class_index, signed_max, max_abs, mean_signed = scalarize(scores)
    arguments.output.mkdir(parents=True)
    score_path = arguments.output / "native_sequence_class_scores.npy"
    np.save(score_path, scores, allow_pickle=False)
    with (arguments.output / "predictions.tsv").open("x", encoding="utf-8", newline="") as handle:
        fields = (
            "fixture_id",
            "element_id",
            "outer_locus_sequence_group_id",
            "outer_fold",
            "max_abs_sequence_class_index0",
            "signed_max_abs_score",
            "max_abs_score",
            "mean_signed_score",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(manifest):
            writer.writerow(
                {
                    "fixture_id": row["fixture_id"],
                    "element_id": row["element_id"],
                    "outer_locus_sequence_group_id": row["outer_locus_sequence_group_id"],
                    "outer_fold": row["outer_fold"],
                    "max_abs_sequence_class_index0": int(class_index[index]),
                    "signed_max_abs_score": format(float(signed_max[index]), ".12g"),
                    "max_abs_score": format(float(max_abs[index]), ".12g"),
                    "mean_signed_score": format(float(mean_signed[index]), ".12g"),
                }
            )
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-sei-gse281364-prediction-v1",
        "status": "pass_outcome_blind_native_prediction",
        "dataset_id": "gse281364",
        "model_id": "sei",
        "elements": len(manifest),
        "outer_locus_sequence_groups": len(manifest),
        "outer_folds": 5,
        "native_sequence_classes": SEQUENCE_CLASSES,
        "primary_signed_score": "signed_value_of_lowest_index_maximum_absolute_native_sequence_class",
        "native_unsigned_score": "maximum_absolute_native_sequence_class",
        "secondary_signed_score": "mean_of_40_native_sequence_class_scores",
        "checkpoint_sha256": digest(arguments.checkpoint),
        "fixture_artifacts_sha256": digest(arguments.fixture_artifacts),
        "repeat_max_abs": repeat_max_abs,
        "batch_variants": arguments.batch_variants,
        "gpu": torch.cuda.get_device_name(device),
        "peak_gpu_memory_bytes": int(torch.cuda.max_memory_allocated(device)),
        "scores_sha256": digest(score_path),
        "predictions_sha256": digest(arguments.output / "predictions.tsv"),
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "model_fitted_or_adapted": False,
        "champion_eligible": False,
        "terms": "academic_and_research_use_only_restricted_comparator",
    }
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-source", type=Path, required=True)
    parser.add_argument("--scoring-source", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--histone-indices", type=Path, required=True)
    parser.add_argument("--fixture-artifacts", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--batch-variants", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    predict(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
