#!/usr/bin/env python3
"""Predict Enformer allele effects for frozen, outcome-blind MPRA loci."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import random
from typing import Any

import numpy as np

from scripts.alphagenome_sei_build_fixture import IndexedFasta
from scripts.enformer_crested_probe_checkpoint import (
    ATTENTION_MEMBER,
    ATTENTION_SHA256,
    KERAS_NAME,
    KERAS_SHA256,
    KERAS_SIZE,
    LAYERS_MEMBER,
    LAYERS_SHA256,
    LAYERS_SIZE,
    OUTER_SHA256,
    _extract_exact_member,
    _import_custom_layers,
    _one_hot,
)


WINDOW = 196_608
CENTER = WINDOW // 2
TRACKS = 5_313
OUTPUT_BINS = 896
SEED = 20260824
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
HEPG2_ACCESSIBILITY = (
    (27, "ENCFF136DBS", "DNASE:HepG2"),
    (91, "ENCFF205TKQ", "DNASE:HepG2"),
    (234, "ENCFF577SOF", "DNASE:HepG2"),
)
LIVER_ACCESSIBILITY = (
    (26, "ENCFF818FXA", "DNASE:hepatocyte"),
    (448, "ENCFF462ZLK", "DNASE:right lobe of liver female adult (53 years)"),
)


class EnformerPredictionError(ValueError):
    """Raised when restricted Enformer inference differs from its frozen requirements."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if (
        len(rows) != 1_033
        or len({row["fixture_id"] for row in rows}) != 1_033
        or len({row["element_id"] for row in rows}) != 1_033
        or len({row["outer_locus_sequence_group_id"] for row in rows}) != 1_033
        or {int(row["outer_fold"]) for row in rows} != set(range(5))
        or any(int(row["variant_index0"]) != CENTER for row in rows)
    ):
        raise EnformerPredictionError("fixture manifest differs")
    return rows


def read_track_contract(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    if len(rows) != TRACKS or [int(row["position"]) for row in rows] != list(range(TRACKS)):
        raise EnformerPredictionError("track position contract differs")
    for index, identifier, description in HEPG2_ACCESSIBILITY + LIVER_ACCESSIBILITY:
        row = rows[index]
        if (
            int(row["canonical_index"]) != index
            or row["identifier"] != identifier
            or row["description"] != description
            or row["assay_family"] != "accessibility"
        ):
            raise EnformerPredictionError("prespecified accessibility track differs")
    if sum(row["assay_family"] == "accessibility" for row in rows) != 684:
        raise EnformerPredictionError("accessibility track census differs")
    return rows


def allele_sequences(reference: IndexedFasta, row: dict[str, str]) -> tuple[str, str, str, str]:
    start, end = int(row["input_start0"]), int(row["input_end0"])
    sequence = reference.fetch(row["contig"], start, end)
    if (
        len(sequence) != WINDOW
        or start + CENTER != int(row["variant_pos0"])
        or sequence[CENTER] != row["ref"]
        or sha256(sequence.encode()).hexdigest() != row["reference_sequence_sha256"]
    ):
        raise EnformerPredictionError("reference sequence identity differs")
    alternative = sequence[:CENTER] + row["alt"] + sequence[CENTER + 1 :]
    ref_rc, alt_rc = reverse_complement(sequence), reverse_complement(alternative)
    observed = (
        sha256(alternative.encode()).hexdigest(),
        sha256(ref_rc.encode()).hexdigest(),
        sha256(alt_rc.encode()).hexdigest(),
    )
    expected = (
        row["alternative_sequence_sha256"],
        row["reference_reverse_complement_sha256"],
        row["alternative_reverse_complement_sha256"],
    )
    if observed != expected:
        raise EnformerPredictionError("alternative or reverse-complement identity differs")
    return sequence, alternative, ref_rc, alt_rc


def score_sums(sums: np.ndarray) -> np.ndarray:
    values = np.asarray(sums, dtype=np.float64)
    if values.shape[0] != 4 or values.shape[1] != TRACKS or not np.isfinite(values).all():
        raise EnformerPredictionError("four-arm sum matrix differs")
    if np.any(values < 0):
        raise EnformerPredictionError("nonnegative Enformer sum contract differs")
    ref, alt, ref_rc, alt_rc = values
    forward_sad = alt - ref
    forward_sar = np.log2(alt + 1.0) - np.log2(ref + 1.0)
    tta_ref, tta_alt = 0.5 * (ref + ref_rc), 0.5 * (alt + alt_rc)
    tta_sad = tta_alt - tta_ref
    tta_sar = np.log2(tta_alt + 1.0) - np.log2(tta_ref + 1.0)
    result = np.stack((forward_sad, forward_sar, tta_sad, tta_sar)).astype(np.float32)
    if result.shape != (4, TRACKS) or not np.isfinite(result).all():
        raise EnformerPredictionError("allele summary differs")
    return result


def _cache_contract(arguments: argparse.Namespace) -> dict[str, object]:
    predictor_source = Path(__file__).resolve(strict=True)
    probe_source = predictor_source.with_name("enformer_crested_probe_checkpoint.py")
    return {
        "schema_version": "masld-bench-enformer-gse281364-cache-v1",
        "cache_id": arguments.cache_id,
        "model_archive_sha256": OUTER_SHA256,
        "source_archive_sha256": digest(arguments.source_archive),
        "checkpoint_member_sha256": KERAS_SHA256,
        "predictor_source_sha256": digest(predictor_source),
        "probe_source_sha256": digest(probe_source),
        "fixture_artifacts_sha256": digest(arguments.fixture_artifacts),
        "manifest_sha256": digest(arguments.manifest),
        "fasta_artifacts_sha256": digest(arguments.fasta_artifacts),
        "track_crosswalk_artifacts_sha256": digest(arguments.track_crosswalk_artifacts),
        "track_crosswalk_sha256": digest(arguments.track_crosswalk),
        "chunk_size": arguments.chunk_size,
        "sequence_batch_size": 1,
        "arms": ["REF", "ALT", "REF_RC", "ALT_RC"],
        "features": ["forward_sad", "forward_sar", "rc_ensemble_sad", "rc_ensemble_sar"],
        "outcomes_read": False,
    }


def _prepare_cache(root: Path, contract: dict[str, object]) -> Path:
    if not root.is_dir() or root.is_symlink():
        raise EnformerPredictionError("cache root must be an existing regular directory")
    cache = root / str(contract["cache_id"])
    try:
        cache.mkdir(mode=0o750)
    except FileExistsError:
        if not cache.is_dir() or cache.is_symlink():
            raise EnformerPredictionError("cache path differs")
    contract_path = cache / "contract.json"
    encoded = json.dumps(contract, indent=2, sort_keys=True) + "\n"
    try:
        with contract_path.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
    except FileExistsError:
        if contract_path.read_text() != encoded:
            raise EnformerPredictionError("cache contract differs")
    return cache


def _load_chunk(path: Path, start: int, end: int) -> np.ndarray:
    receipt_path, array_path = path / "receipt.json", path / "features.npy"
    if not path.is_dir() or path.is_symlink() or not receipt_path.is_file() or not array_path.is_file():
        raise EnformerPredictionError("checkpoint chunk is incomplete")
    receipt = json.loads(receipt_path.read_text())
    if (
        receipt.get("start") != start
        or receipt.get("end") != end
        or receipt.get("features_sha256") != digest(array_path)
        or receipt.get("outcomes_read") is not False
        or (start == 0 and receipt.get("deterministic_repeat_max_abs") != 0.0)
    ):
        raise EnformerPredictionError("checkpoint chunk receipt differs")
    values = np.load(array_path, allow_pickle=False)
    if values.shape != (end - start, 4, TRACKS) or values.dtype != np.float32 or not np.isfinite(values).all():
        raise EnformerPredictionError("checkpoint chunk array differs")
    return values


def _write_chunk(
    cache: Path,
    start: int,
    end: int,
    values: np.ndarray,
    deterministic_repeat_max_abs: float | None = None,
) -> Path:
    final = cache / f"chunk-{start:04d}-{end:04d}"
    if final.exists():
        raise EnformerPredictionError("checkpoint chunk unexpectedly exists")
    staging = cache / f"chunk-{start:04d}-{end:04d}.staging-{os.getpid()}"
    staging.mkdir(mode=0o750)
    array_path = staging / "features.npy"
    with array_path.open("xb") as handle:
        np.save(handle, values, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    receipt = {
        "schema_version": "masld-bench-enformer-gse281364-cache-chunk-v1",
        "start": start,
        "end": end,
        "features_sha256": digest(array_path),
        "deterministic_repeat_max_abs": deterministic_repeat_max_abs,
        "outcomes_read": False,
    }
    if start == 0 and deterministic_repeat_max_abs != 0.0:
        raise EnformerPredictionError("first chunk lacks an exact repeat check")
    with (staging / "receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    staging.rename(final)
    return final


def _load_model(arguments: argparse.Namespace, scratch: Path) -> Any:
    import keras
    import tensorflow as tf

    if digest(arguments.model_archive) != OUTER_SHA256:
        raise EnformerPredictionError("outer model archive differs")
    model_path = scratch / KERAS_NAME
    attention_path = scratch / "crested_attention.py"
    layers_path = scratch / "crested_layers.py"
    _extract_exact_member(arguments.model_archive, KERAS_NAME, model_path, KERAS_SHA256, KERAS_SIZE)
    _extract_exact_member(arguments.source_archive, ATTENTION_MEMBER, attention_path, ATTENTION_SHA256)
    _extract_exact_member(arguments.source_archive, LAYERS_MEMBER, layers_path, LAYERS_SHA256, LAYERS_SIZE)
    random.seed(SEED)
    np.random.seed(SEED)
    tf.keras.utils.set_random_seed(SEED)
    model = keras.models.load_model(
        model_path,
        custom_objects=_import_custom_layers(attention_path, layers_path),
        compile=False,
        safe_mode=True,
    )
    model.trainable = False
    if list(model.input_shape) != [None, WINDOW, 4] or list(model.output_shape) != [None, OUTPUT_BINS, TRACKS]:
        raise EnformerPredictionError("model input or output shape differs")
    if len(tf.config.list_physical_devices("GPU")) != 1:
        raise EnformerPredictionError("exactly one GPU is required")
    return model


def _predict_sums(model: Any, sequences: tuple[str, str, str, str]) -> np.ndarray:
    import tensorflow as tf

    rows: list[np.ndarray] = []
    for sequence in sequences:
        tensor = tf.convert_to_tensor(_one_hot(sequence)[None, :, :], dtype=tf.float32)
        output = model(tensor, training=False)
        if tuple(output.shape) != (1, OUTPUT_BINS, TRACKS):
            raise EnformerPredictionError("model output shape differs")
        value = np.asarray(tf.reduce_sum(output[0], axis=0).numpy(), dtype=np.float32)
        if value.shape != (TRACKS,) or not np.isfinite(value).all() or np.any(value < 0):
            raise EnformerPredictionError("model track sum differs")
        rows.append(value)
    return np.stack(rows)


def predict(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or not 1 <= arguments.chunk_size <= 64:
        raise EnformerPredictionError("output or chunk-size contract differs")
    rows = read_manifest(arguments.manifest)
    tracks = read_track_contract(arguments.track_crosswalk)
    contract = _cache_contract(arguments)
    cache = _prepare_cache(arguments.checkpoint_root, contract)
    reference = IndexedFasta(arguments.fasta, arguments.fai)
    model = _load_model(arguments, arguments.scratch)
    completed: list[np.ndarray] = []
    repeat_max_abs = -1.0
    chunks_computed = 0
    for start in range(0, len(rows), arguments.chunk_size):
        end = min(start + arguments.chunk_size, len(rows))
        chunk_path = cache / f"chunk-{start:04d}-{end:04d}"
        if chunk_path.exists():
            values = _load_chunk(chunk_path, start, end)
            if start == 0:
                repeat_max_abs = 0.0
        else:
            values = np.empty((end - start, 4, TRACKS), dtype=np.float32)
            for offset, row in enumerate(rows[start:end]):
                sequences = allele_sequences(reference, row)
                sums = _predict_sums(model, sequences)
                if start == 0 and offset == 0:
                    repeat = _predict_sums(model, (sequences[0],))[0]
                    repeat_max_abs = float(np.max(np.abs(sums[0] - repeat)))
                values[offset] = score_sums(sums)
                print(json.dumps({"completed": start + offset + 1, "total": len(rows)}), flush=True)
            chunk_path = _write_chunk(
                cache,
                start,
                end,
                values,
                deterministic_repeat_max_abs=repeat_max_abs if start == 0 else None,
            )
            values = _load_chunk(chunk_path, start, end)
            chunks_computed += 1
        completed.append(values)
    features = np.concatenate(completed, axis=0)
    if features.shape != (len(rows), 4, TRACKS):
        raise EnformerPredictionError("complete feature matrix differs")

    arguments.output.mkdir(parents=True)
    names = ("forward_sad", "forward_sar", "rc_ensemble_sad", "rc_ensemble_sar")
    feature_hashes: dict[str, str] = {}
    for index, name in enumerate(names):
        path = arguments.output / f"{name}.npy"
        with path.open("xb") as handle:
            np.save(handle, features[:, index, :], allow_pickle=False)
        feature_hashes[name] = digest(path)
    access_indices = [i for i, row in enumerate(tracks) if row["assay_family"] == "accessibility"]
    hepg2_indices = [index for index, _, _ in HEPG2_ACCESSIBILITY]
    liver_indices = [index for index, _, _ in LIVER_ACCESSIBILITY]
    tta_sad, tta_sar = features[:, 2, :], features[:, 3, :]
    with (arguments.output / "predictions.tsv").open("x", encoding="utf-8", newline="") as handle:
        fields = (
            "fixture_id", "element_id", "outer_locus_sequence_group_id", "outer_fold",
            "hepg2_accessibility_sad", "hepg2_accessibility_sar",
            "liver_accessibility_sad", "liver_accessibility_sar",
            "all_accessibility_mean_sad", "all_accessibility_mean_sar",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows):
            writer.writerow({
                "fixture_id": row["fixture_id"],
                "element_id": row["element_id"],
                "outer_locus_sequence_group_id": row["outer_locus_sequence_group_id"],
                "outer_fold": row["outer_fold"],
                "hepg2_accessibility_sad": format(float(tta_sad[index, hepg2_indices].mean()), ".12g"),
                "hepg2_accessibility_sar": format(float(tta_sar[index, hepg2_indices].mean()), ".12g"),
                "liver_accessibility_sad": format(float(tta_sad[index, liver_indices].mean()), ".12g"),
                "liver_accessibility_sar": format(float(tta_sar[index, liver_indices].mean()), ".12g"),
                "all_accessibility_mean_sad": format(float(tta_sad[index, access_indices].mean()), ".12g"),
                "all_accessibility_mean_sar": format(float(tta_sar[index, access_indices].mean()), ".12g"),
            })
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-enformer-gse281364-prediction-v1",
        "status": "pass_outcome_blind_restricted_prediction",
        "dataset_id": "gse281364",
        "model_id": "enformer_crested_restricted_port",
        "registered_scientific_identity": "restricted_conversion_not_native_sonnet",
        "elements": len(rows),
        "outer_locus_sequence_groups": len(rows),
        "outer_folds": 5,
        "tracks": TRACKS,
        "accessibility_tracks": len(access_indices),
        "features": list(names),
        "allele_delta_sign": "ALT_minus_REF",
        "reverse_complement": "forward_and_reverse_complement_sum_ensemble_track_axis_identity",
        "track_scalarization": "outcome_blind_exact_position_identifier_bound_accessibility_sets",
        "checkpoint_member_sha256": KERAS_SHA256,
        "fixture_artifacts_sha256": digest(arguments.fixture_artifacts),
        "track_crosswalk_artifacts_sha256": digest(arguments.track_crosswalk_artifacts),
        "feature_sha256": feature_hashes,
        "predictions_sha256": digest(arguments.output / "predictions.tsv"),
        "cache_id": arguments.cache_id,
        "chunks": len(completed),
        "chunks_computed_this_attempt": chunks_computed,
        "deterministic_repeat_max_abs": repeat_max_abs,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "sealed_outcomes_read": False,
        "model_fitted_or_adapted": False,
        "native_sonnet_parity_established": False,
        "champion_eligible": False,
        "terms": "internal_academic_noncommercial_nontransferable_restricted_comparator",
    }
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-archive", type=Path, required=True)
    parser.add_argument("--source-archive", type=Path, required=True)
    parser.add_argument("--fixture-artifacts", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta-artifacts", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--track-crosswalk-artifacts", type=Path, required=True)
    parser.add_argument("--track-crosswalk", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--cache-id", required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--output", type=Path, required=True)
    predict(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
