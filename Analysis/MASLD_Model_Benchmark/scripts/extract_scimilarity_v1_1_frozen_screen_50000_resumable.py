#!/usr/bin/env python3
"""Resumably extract raw SCimilarity embeddings from the frozen 50k fixture."""

from __future__ import annotations

import argparse
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
from typing import Any

from scripts.extract_scimilarity_v1_1_frozen_screen_50000 import (
    EMBEDDING_POLICY,
    ENCODER_SHA256,
    INPUT_DIMENSION,
    LATENT_DIMENSION,
    ROWS,
    SCimilarityExtractionError,
    _load_encoder,
    embed_batch,
    load_fixture,
    sha256_file,
)


STATE_SCHEMA = "masld-bench-scimilarity-v1-1-resumable-state-v1"
RECEIPT_SCHEMA = "masld-bench-scimilarity-v1-1-resumable-extraction-v1"
_STOP_REQUESTED = False


def _request_stop(_signum: int, _frame: object) -> None:
    global _STOP_REQUESTED
    _STOP_REQUESTED = True


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _state_contract(
    fixture_root: Path,
    fixture: dict[str, Any],
    *,
    rows: int,
    dimension: int,
    batch_size: int,
    chunk_size: int,
) -> dict[str, Any]:
    return {
        "schema_version": STATE_SCHEMA,
        "model_id": "scimilarity_v1_1",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "checkpoint_sha256": ENCODER_SHA256,
        "fixture_receipt_sha256": sha256_file(fixture_root / "fixture_receipt.json"),
        "fixture_row_order_sha256": fixture["receipt"]["row_order_sha256"],
        "rows": rows,
        "input_dimension": dimension,
        "embedding_width": LATENT_DIMENSION,
        "batch_size": batch_size,
        "chunk_size": chunk_size,
        "embedding_policy": EMBEDDING_POLICY,
        "aggregate_development_exposure": "encoder_seen",
        "released_reference_index_used": False,
        "downstream_head_fit": False,
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
    }


def _load_or_create_contract(path: Path, expected: dict[str, Any]) -> None:
    if path.exists():
        observed = json.loads(path.read_text(encoding="utf-8"))
        if observed != expected:
            raise SCimilarityExtractionError("resumable extraction state contract differs")
        return
    _atomic_json(path, expected)


def _chunk_paths(root: Path, start: int, stop: int) -> tuple[Path, Path]:
    stem = f"chunk-{start:08d}-{stop:08d}"
    return root / f"{stem}.npy", root / f"{stem}.json"


def _chunk_receipt(path: Path, values: Any, start: int, stop: int) -> dict[str, Any]:
    return {
        "schema_version": "masld-bench-scimilarity-v1-1-embedding-chunk-v1",
        "start": start,
        "stop": stop,
        "shape": [stop - start, LATENT_DIMENSION],
        "dtype": "float32",
        "sha256": sha256_file(path),
        "unit_l2": True,
    }


def _validate_chunk(data_path: Path, receipt_path: Path, start: int, stop: int):
    import numpy as np

    if receipt_path.exists() and not data_path.exists():
        raise SCimilarityExtractionError("resumable chunk receipt lacks data")
    if not data_path.exists():
        return None
    values = np.load(data_path, allow_pickle=False)
    if (
        values.shape != (stop - start, LATENT_DIMENSION)
        or values.dtype != np.dtype("float32")
        or not np.isfinite(values).all()
        or not np.allclose(np.linalg.norm(values, axis=1), 1.0, rtol=1e-5, atol=1e-6)
    ):
        raise SCimilarityExtractionError("resumable embedding chunk differs")
    expected = _chunk_receipt(data_path, values, start, stop)
    if receipt_path.exists():
        observed = json.loads(receipt_path.read_text(encoding="utf-8"))
        if observed != expected:
            raise SCimilarityExtractionError("resumable embedding chunk receipt differs")
    else:
        _atomic_json(receipt_path, expected)
    return values


def _write_chunk(
    data_path: Path,
    receipt_path: Path,
    values: Any,
    start: int,
    stop: int,
) -> None:
    import numpy as np

    temporary = data_path.with_name(f".{data_path.name}.{os.getpid()}.tmp")
    with temporary.open("xb") as handle:
        np.save(handle, values, allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, data_path)
    _atomic_json(receipt_path, _chunk_receipt(data_path, values, start, stop))


def extract_resumable(
    fixture_root: Path,
    checkpoint: Path,
    state_root: Path,
    output: Path,
    *,
    expected_rows: int = ROWS,
    expected_dimension: int = INPUT_DIMENSION,
    batch_size: int = 512,
    chunk_size: int = 4096,
    device: str = "cuda",
    stop_after_new_chunks: int | None = None,
) -> dict[str, Any]:
    import numpy as np
    import torch

    if (
        output.exists()
        or batch_size < 1
        or chunk_size < batch_size
        or chunk_size % batch_size
        or device not in {"cpu", "cuda"}
        or (stop_after_new_chunks is not None and stop_after_new_chunks < 1)
    ):
        raise SCimilarityExtractionError("resumable output/batch/device contract differs")
    if device == "cuda" and (
        not torch.cuda.is_available() or torch.cuda.device_count() != 1
    ):
        raise SCimilarityExtractionError("exactly one CUDA device is required")

    state_root.mkdir(parents=True, exist_ok=True)
    chunks_root = state_root / "chunks"
    chunks_root.mkdir(exist_ok=True)
    lock_path = state_root / "extraction.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        fixture = load_fixture(
            fixture_root,
            expected_rows=expected_rows,
            expected_dimension=expected_dimension,
        )
        contract = _state_contract(
            fixture_root,
            fixture,
            rows=expected_rows,
            dimension=expected_dimension,
            batch_size=batch_size,
            chunk_size=chunk_size,
        )
        contract_path = state_root / "contract.json"
        _load_or_create_contract(contract_path, contract)
        model, loading = _load_encoder(checkpoint, device=device)
        if device == "cuda":
            torch.cuda.reset_peak_memory_stats()

        chunks: list[tuple[int, int, Path, Path]] = []
        new_chunks = 0
        reused_chunks = 0
        for start in range(0, expected_rows, chunk_size):
            stop = min(start + chunk_size, expected_rows)
            data_path, receipt_path = _chunk_paths(chunks_root, start, stop)
            values = _validate_chunk(data_path, receipt_path, start, stop)
            if values is not None:
                reused_chunks += 1
                chunks.append((start, stop, data_path, receipt_path))
                continue
            values = np.empty((stop - start, LATENT_DIMENSION), dtype=np.float32)
            for batch_start in range(start, stop, batch_size):
                batch_stop = min(batch_start + batch_size, stop)
                values[batch_start - start : batch_stop - start] = embed_batch(
                    model,
                    fixture["matrix"][batch_start:batch_stop],
                    device=device,
                )
            _write_chunk(data_path, receipt_path, values, start, stop)
            chunks.append((start, stop, data_path, receipt_path))
            new_chunks += 1
            _atomic_json(
                state_root / "progress.json",
                {
                    "schema_version": "masld-bench-scimilarity-v1-1-resumable-progress-v1",
                    "chunks_complete": len(chunks),
                    "rows_complete": stop,
                    "last_completed_chunk": data_path.name,
                    "sealed_outcomes_read": False,
                },
            )
            if _STOP_REQUESTED or (
                stop_after_new_chunks is not None
                and new_chunks >= stop_after_new_chunks
            ):
                return {
                    "status": "paused_after_atomic_chunk",
                    "chunks_complete": len(chunks),
                    "rows_complete": stop,
                    "new_chunks": new_chunks,
                    "reused_chunks": reused_chunks,
                }

        output.mkdir(parents=True)
        embeddings = np.empty((expected_rows, LATENT_DIMENSION), dtype=np.float32)
        chunk_manifest = []
        for start, stop, data_path, receipt_path in chunks:
            embeddings[start:stop] = _validate_chunk(
                data_path, receipt_path, start, stop
            )
            chunk_manifest.append(json.loads(receipt_path.read_text(encoding="utf-8")))
        representation = output / "raw_embeddings.npz"
        np.savez_compressed(
            representation,
            embeddings=embeddings,
            outer_folds=fixture["outer_folds"],
            row_ids=np.asarray(fixture["rows"]),
        )
        np.save(output / "observed_gene_mask.npy", fixture["observed_mask"], allow_pickle=False)
        np.save(
            output / "study_gene_observability.npy",
            fixture["study_observability"],
            allow_pickle=False,
        )
        _atomic_json(output / "state_contract.json", contract)
        _atomic_json(
            output / "chunk_manifest.json",
            {
                "schema_version": "masld-bench-scimilarity-v1-1-chunk-manifest-v1",
                "chunks": chunk_manifest,
            },
        )
        with np.load(representation, allow_pickle=False) as observed:
            if set(observed.files) != {"embeddings", "outer_folds", "row_ids"}:
                raise SCimilarityExtractionError("raw embedding NPZ field contract differs")
            if (
                observed["embeddings"].shape != (expected_rows, LATENT_DIMENSION)
                or observed["embeddings"].dtype != np.dtype("float32")
                or observed["outer_folds"].shape != (expected_rows,)
                or observed["row_ids"].shape != (expected_rows,)
            ):
                raise SCimilarityExtractionError("raw embedding output differs")
        receipt = {
            "schema_version": RECEIPT_SCHEMA,
            "status": "pass_raw_outcome_blind_resumable_embeddings",
            "model_id": "scimilarity_v1_1",
            "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
            "rows": expected_rows,
            "embedding_shape": [expected_rows, LATENT_DIMENSION],
            "embedding_dtype": "float32",
            "representation_sha256": sha256_file(representation),
            "checkpoint_sha256": ENCODER_SHA256,
            "strict_weights_only_restore": True,
            "loading": loading,
            "batch_size": batch_size,
            "chunk_size": chunk_size,
            "chunks": len(chunks),
            "chunks_reused_this_attempt": reused_chunks,
            "chunks_created_this_attempt": new_chunks,
            "state_contract_sha256": sha256_file(output / "state_contract.json"),
            "aggregate_development_exposure": "encoder_seen",
            "development_exposure_by_study": fixture["receipt"]["exposure_by_study"],
            "observed_genes": int(fixture["observed_mask"].sum()),
            "study_gene_observability_order": fixture["receipt"][
                "study_gene_observability_order"
            ],
            "study_gene_observability_shape": list(fixture["study_observability"].shape),
            "study_gene_observability_is_global": fixture["receipt"][
                "study_gene_observability_is_global"
            ],
            "native_model_accepts_observability_mask": False,
            "native_fixed_vocabulary_exception": fixture["receipt"][
                "native_fixed_vocabulary_exception"
            ],
            "device": torch.cuda.get_device_name(0) if device == "cuda" else "cpu",
            "cuda": torch.version.cuda if device == "cuda" else None,
            "cuda_peak_allocated_bytes": (
                int(torch.cuda.max_memory_allocated()) if device == "cuda" else 0
            ),
            "cuda_peak_reserved_bytes": (
                int(torch.cuda.max_memory_reserved()) if device == "cuda" else 0
            ),
            "raw_unrefined_embeddings_only": True,
            "released_reference_index_used": False,
            "downstream_head_fit": False,
            "common_head_fit": False,
            "metrics_run": False,
            "evaluation_label_columns_read": [],
            "histology_read": False,
            "sealed_outcomes_read": False,
            "sealed_champion_eligible": False,
        }
        _atomic_json(output / "receipt.json", receipt)
        return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, default=ROWS)
    parser.add_argument("--expected-dimension", type=int, default=INPUT_DIMENSION)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--chunk-size", type=int, default=4096)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    arguments = parser.parse_args()
    signal.signal(signal.SIGUSR1, _request_stop)
    result = extract_resumable(
        arguments.fixture_root,
        arguments.checkpoint,
        arguments.state_root,
        arguments.output,
        expected_rows=arguments.expected_rows,
        expected_dimension=arguments.expected_dimension,
        batch_size=arguments.batch_size,
        chunk_size=arguments.chunk_size,
        device=arguments.device,
    )
    print(json.dumps(result, sort_keys=True))
    return 75 if result["status"] == "paused_after_atomic_chunk" else 0


if __name__ == "__main__":
    raise SystemExit(main())
