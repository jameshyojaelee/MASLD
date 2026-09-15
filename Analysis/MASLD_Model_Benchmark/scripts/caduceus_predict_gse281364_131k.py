#!/usr/bin/env python3
"""Extract outcome-blind Caduceus 131-kb REF/ALT/RC embeddings."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import random
import sys
from typing import Any, Mapping

import numpy as np
import torch
from safetensors.torch import load_file

from scripts.alphagenome_sei_build_fixture import IndexedFasta


WINDOW = 131_072
CENTER = WINDOW // 2
POOL_START = 64_768
POOL_END = 66_304
WIDTH = 256
ALLELES = ("REF", "ALT", "REF_RC", "ALT_RC")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
TOKEN_IDS = {"A": 7, "C": 8, "G": 9, "T": 10, "N": 11}
SEED = 20260824


class CaduceusPredictionError(ValueError):
    """Raised when the frozen long-context prediction requirement differs."""


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
        or any(int(row["input_length_bp"]) != WINDOW for row in rows)
        or any(int(row["forward_variant_index0"]) != CENTER for row in rows)
        or any(int(row["reverse_complement_variant_index0"]) != CENTER - 1 for row in rows)
        or any(int(row["pool_start0"]) != POOL_START for row in rows)
        or any(int(row["pool_end0"]) != POOL_END for row in rows)
    ):
        raise CaduceusPredictionError("fixture manifest differs")
    return rows


def allele_sequences(
    reference: IndexedFasta, row: Mapping[str, str]
) -> tuple[str, str, str, str]:
    start, end = int(row["input_start0"]), int(row["input_end0"])
    sequence = reference.fetch(row["contig"], start, end)
    if (
        len(sequence) != WINDOW
        or start + CENTER != int(row["variant_pos0"])
        or sequence[CENTER] != row["ref"]
        or digest_text(sequence) != row["reference_sequence_sha256"]
    ):
        raise CaduceusPredictionError("reference sequence identity differs")
    alternative = sequence[:CENTER] + row["alt"] + sequence[CENTER + 1 :]
    ref_rc, alt_rc = reverse_complement(sequence), reverse_complement(alternative)
    observed = (
        digest_text(alternative),
        digest_text(ref_rc),
        digest_text(alt_rc),
    )
    expected = (
        row["alternative_sequence_sha256"],
        row["reverse_complement_reference_sha256"],
        row["reverse_complement_alternative_sha256"],
    )
    if (
        observed != expected
        or sum(left != right for left, right in zip(sequence, alternative)) != 1
        or sum(left != right for left, right in zip(ref_rc, alt_rc)) != 1
        or ref_rc[CENTER - 1] != reverse_complement(row["ref"])
        or alt_rc[CENTER - 1] != reverse_complement(row["alt"])
    ):
        raise CaduceusPredictionError("alternative or reverse-complement identity differs")
    return sequence, alternative, ref_rc, alt_rc


def digest_text(value: str) -> str:
    return sha256(value.encode("ascii")).hexdigest()


def encode_batch(sequences: tuple[str, ...], device: torch.device) -> torch.Tensor:
    lookup = np.full(256, -1, dtype=np.int16)
    for base, token in TOKEN_IDS.items():
        lookup[ord(base)] = token
    rows = []
    for sequence in sequences:
        encoded = lookup[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
        if encoded.shape != (WINDOW,) or np.any(encoded < 0):
            raise CaduceusPredictionError("sequence tokenization differs")
        rows.append(encoded)
    values = np.stack(rows).astype(np.int64, copy=False)
    return torch.from_numpy(values).to(device=device, non_blocking=False)


def pooled_features(
    model: torch.nn.Module, sequences: tuple[str, ...], device: torch.device
) -> np.ndarray:
    ids = encode_batch(sequences, device)
    with torch.inference_mode():
        hidden = model.caduceus(input_ids=ids, return_dict=True).last_hidden_state
        if tuple(hidden.shape) != (len(sequences), WINDOW, 2 * WIDTH):
            raise CaduceusPredictionError("Caduceus RCPS hidden-state shape differs")
        forward = hidden[..., :WIDTH]
        reverse = torch.flip(hidden[..., WIDTH:], dims=(-2, -1))
        pooled = ((forward + reverse) / 2.0)[:, POOL_START:POOL_END].float().mean(dim=1)
    values = pooled.cpu().numpy().astype(np.float32, copy=False)
    if values.shape != (len(sequences), WIDTH) or not np.isfinite(values).all():
        raise CaduceusPredictionError("pooled embedding differs")
    return values


def load_model(
    code_root: Path, checkpoint: Path, device: torch.device
) -> torch.nn.Module:
    sys.path.insert(0, str(code_root.parent))
    from caduceus.configuration_caduceus import CaduceusConfig
    from caduceus.modeling_caduceus import CaduceusForMaskedLM

    config = CaduceusConfig.from_json_file(str(code_root / "config.json"))
    if (
        config.model_type != "caduceus"
        or config.d_model != WIDTH
        or config.n_layer != 16
        or config.vocab_size != 16
        or not config.rcps
        or not config.bidirectional
        or not config.bidirectional_weight_tie
        or config.bidirectional_strategy != "add"
    ):
        raise CaduceusPredictionError("checkpoint configuration differs")
    model = CaduceusForMaskedLM(config)
    state = load_file(str(checkpoint), device="cpu")
    aliases = {
        f"caduceus.backbone.layers.{layer}.mixer.submodule.mamba_rev.{projection}.weight"
        for layer in range(16)
        for projection in ("in_proj", "out_proj")
    } | {"lm_head.lm_head.weight"}
    incompatible = model.load_state_dict(state, strict=False)
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
        raise CaduceusPredictionError("checkpoint restore differs")
    del state
    if not all(
        layer.mixer.submodule.mamba_rev.in_proj.weight
        is layer.mixer.submodule.mamba_fwd.in_proj.weight
        and layer.mixer.submodule.mamba_rev.out_proj.weight
        is layer.mixer.submodule.mamba_fwd.out_proj.weight
        for layer in model.caduceus.backbone.layers
    ) or model.lm_head.weight is not model.get_input_embeddings().weight:
        raise CaduceusPredictionError("tied parameter identity differs")
    return model.to(device).eval()


def cache_contract(arguments: argparse.Namespace) -> dict[str, object]:
    return {
        "schema_version": "masld-bench-caduceus-gse281364-cache-v1",
        "cache_id": arguments.cache_id,
        "predictor_source_sha256": digest(Path(__file__).resolve(strict=True)),
        "config_sha256": digest(arguments.config),
        "fixture_artifacts_sha256": digest(arguments.fixture_artifacts),
        "manifest_sha256": digest(arguments.manifest),
        "fasta_artifacts_sha256": digest(arguments.fasta_artifacts),
        "checkpoint_sha256": digest(arguments.checkpoint),
        "chunk_size": arguments.chunk_size,
        "sequence_batch_size": 4,
        "alleles": list(ALLELES),
        "outcomes_read": False,
    }


def prepare_cache(root: Path, contract: Mapping[str, object]) -> Path:
    if not root.is_dir() or root.is_symlink():
        raise CaduceusPredictionError("cache root differs")
    cache = root / str(contract["cache_id"])
    try:
        cache.mkdir(mode=0o750)
    except FileExistsError:
        if not cache.is_dir() or cache.is_symlink():
            raise CaduceusPredictionError("cache path differs")
    return cache


def _ensure_contract(path: Path, contract: Mapping[str, object]) -> None:
    encoded = json.dumps(contract, indent=2, sort_keys=True) + "\n"
    try:
        with path.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
    except FileExistsError:
        if path.read_text(encoding="utf-8") != encoded:
            raise CaduceusPredictionError("cache contract differs")


def load_chunk(path: Path, start: int, end: int) -> tuple[np.ndarray, dict[str, Any]]:
    receipt_path, array_path = path / "receipt.json", path / "embeddings.npy"
    if (
        not path.is_dir()
        or path.is_symlink()
        or not receipt_path.is_file()
        or not array_path.is_file()
    ):
        raise CaduceusPredictionError("checkpoint chunk is incomplete")
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        values = np.load(array_path, allow_pickle=False)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise CaduceusPredictionError("checkpoint chunk is unreadable") from error
    if (
        receipt.get("start") != start
        or receipt.get("end") != end
        or receipt.get("embeddings_sha256") != digest(array_path)
        or receipt.get("outcomes_read") is not False
        or values.shape != (end - start, 4, WIDTH)
        or values.dtype != np.float32
        or not np.isfinite(values).all()
        or (start == 0 and float(receipt.get("deterministic_repeat_max_abs", 1.0)) > 1e-5)
    ):
        raise CaduceusPredictionError("checkpoint chunk contract differs")
    return values, receipt


def write_chunk(
    cache: Path,
    start: int,
    end: int,
    values: np.ndarray,
    deterministic_repeat_max_abs: float | None,
) -> Path:
    final = cache / f"chunk-{start:04d}-{end:04d}"
    if final.exists():
        raise CaduceusPredictionError("checkpoint chunk unexpectedly exists")
    staging = cache / f"chunk-{start:04d}-{end:04d}.staging-{os.getpid()}"
    staging.mkdir(mode=0o750)
    array_path = staging / "embeddings.npy"
    with array_path.open("xb") as handle:
        np.save(handle, values, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    receipt = {
        "schema_version": "masld-bench-caduceus-gse281364-cache-chunk-v1",
        "start": start,
        "end": end,
        "embeddings_sha256": digest(array_path),
        "deterministic_repeat_max_abs": deterministic_repeat_max_abs,
        "max_external_rc_feature_abs_difference": float(
            max(
                np.max(np.abs(values[:, 0] - values[:, 2])),
                np.max(np.abs(values[:, 1] - values[:, 3])),
            )
        ),
        "outcomes_read": False,
    }
    if start == 0 and (deterministic_repeat_max_abs is None or deterministic_repeat_max_abs > 1e-5):
        raise CaduceusPredictionError("first chunk deterministic repeat differs")
    with (staging / "receipt.json").open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    staging.rename(final)
    return final


def predict(arguments: argparse.Namespace) -> dict[str, object]:
    if arguments.output.exists() or not 1 <= arguments.chunk_size <= 32:
        raise CaduceusPredictionError("prediction request differs")
    config = json.loads(arguments.config.read_text(encoding="utf-8"))
    if (
        config.get("schema_version")
        != "masld-bench-gse281364-caduceus-131k-contract-v1"
        or config.get("model_id") != "caduceus"
        or any(config.get("firewall", {}).values())
        or config.get("checkpoint", {}).get("license") != "Apache-2.0"
        or config.get("exposure", {}).get("checkpoint_state") != "target_label_unexposed"
        or config.get("exposure", {}).get("allowed_supervision")
        != "locus-cross-fitted_assay-native_MPRA_activity_only"
    ):
        raise CaduceusPredictionError("admission configuration differs")
    rows = read_manifest(arguments.manifest)
    contract = cache_contract(arguments)
    cache = prepare_cache(arguments.checkpoint_root, contract)
    _ensure_contract(cache / "contract.json", contract)
    reference = IndexedFasta(arguments.fasta, arguments.fai)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda" or torch.cuda.device_count() != 1:
        raise CaduceusPredictionError("exactly one CUDA device is required")
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed_all(SEED)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    model = load_model(arguments.code_root, arguments.checkpoint, device)
    arrays: list[np.ndarray] = []
    chunk_receipts: list[dict[str, Any]] = []
    for start in range(0, len(rows), arguments.chunk_size):
        end = min(len(rows), start + arguments.chunk_size)
        chunk_path = cache / f"chunk-{start:04d}-{end:04d}"
        if chunk_path.exists():
            values, receipt = load_chunk(chunk_path, start, end)
        else:
            locus_values: list[np.ndarray] = []
            first_sequences: tuple[str, ...] | None = None
            first_values: np.ndarray | None = None
            for row in rows[start:end]:
                sequences = allele_sequences(reference, row)
                values = pooled_features(model, sequences, device)
                locus_values.append(values)
                if start == 0 and first_sequences is None:
                    first_sequences, first_values = sequences, values.copy()
            chunk_values = np.stack(locus_values).astype(np.float32, copy=False)
            repeat = None
            if start == 0:
                assert first_sequences is not None and first_values is not None
                repeated = pooled_features(model, first_sequences, device)
                repeat = float(np.max(np.abs(repeated - first_values)))
            write_chunk(cache, start, end, chunk_values, repeat)
            values, receipt = load_chunk(chunk_path, start, end)
        arrays.append(values)
        chunk_receipts.append(receipt)
        print(json.dumps({"cache_chunk": [start, end], "status": "complete"}), flush=True)
    embeddings = np.concatenate(arrays, axis=0)
    if embeddings.shape != (1_033, 4, WIDTH) or not np.isfinite(embeddings).all():
        raise CaduceusPredictionError("complete embedding tensor differs")
    max_rc = max(float(row["max_external_rc_feature_abs_difference"]) for row in chunk_receipts)
    repeat = float(chunk_receipts[0]["deterministic_repeat_max_abs"])
    if max_rc > 1e-5 or repeat > 1e-5:
        raise CaduceusPredictionError("RCPS or deterministic tolerance differs")
    arguments.output.mkdir(parents=True)
    np.savez_compressed(
        arguments.output / "allele_embeddings.npz",
        fixture_ids=np.asarray([row["fixture_id"] for row in rows]),
        element_ids=np.asarray([row["element_id"] for row in rows]),
        outer_locus_sequence_group_ids=np.asarray(
            [row["outer_locus_sequence_group_id"] for row in rows]
        ),
        outer_folds=np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8),
        allele_order=np.asarray(ALLELES),
        embeddings=embeddings,
    )
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-caduceus-131k-embeddings-v1",
        "status": "pass_outcome_blind_embedding_extraction",
        "dataset_id": "gse281364",
        "model_id": "caduceus",
        "checkpoint_sha256": digest(arguments.checkpoint),
        "fixture_artifacts_sha256": digest(arguments.fixture_artifacts),
        "fasta_artifacts_sha256": digest(arguments.fasta_artifacts),
        "license": "Apache-2.0",
        "exposure_state": "target_label_unexposed",
        "restricted_comparator": False,
        "open_champion_eligible_after_task_and_external_evaluation_gates": True,
        "elements": len(rows),
        "outer_locus_sequence_groups": len(rows),
        "outer_folds": 5,
        "input_length_bp": WINDOW,
        "pool_start0": POOL_START,
        "pool_end0": POOL_END,
        "pool_length_bp": POOL_END - POOL_START,
        "hidden_width": WIDTH,
        "alleles": list(ALLELES),
        "sequence_forwards": len(rows) + 1,
        "sequence_batch_size": 4,
        "cache_chunks": len(chunk_receipts),
        "feature_contract": "RCPS_internal_channels_realigned_and_averaged_then_1536bp_mean_pool",
        "max_external_rc_feature_abs_difference": max_rc,
        "deterministic_repeat_max_abs_difference": repeat,
        "checkpoint_loaded_with_safetensors": True,
        "runtime_network_allowed": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_labels_read": False,
        "downstream_head_fit": False,
        "allowed_future_head": "locus-cross-fitted_assay-native_MPRA_activity_only",
        "prohibited_head": "cell-type_eQTL_or_ieQTL_effect",
        "device": str(device),
        "gpu": torch.cuda.get_device_name(0),
        "torch": torch.__version__,
        "maximum_gpu_memory_bytes": torch.cuda.max_memory_allocated(),
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--code-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--fixture-artifacts", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--fasta-artifacts", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--cache-id", required=True)
    parser.add_argument("--chunk-size", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    predict(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
