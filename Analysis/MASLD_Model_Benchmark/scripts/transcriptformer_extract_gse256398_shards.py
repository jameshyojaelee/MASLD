#!/usr/bin/env python3
"""Extract label-blind TranscriptFormer embeddings from GSE256398 shards."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.transcriptformer_extract_unlabeled_embeddings import (
    BATCH_SIZE,
    DETERMINISM_TOLERANCE,
    EMBEDDING_WIDTH,
    SEQUENCE_LENGTH,
    _batch,
    _embedding_only_forward,
    _load_model,
)


EXPECTED_ROWS = 172_997
EXPECTED_DONORS = 26
POLICIES = ("native", "common")
FIXTURE_FIELDS = {
    "assay_token_indices",
    "common_gene_counts",
    "common_gene_token_indices",
    "in_source_exact",
    "native_gene_counts",
    "native_gene_token_indices",
    "outer_folds",
    "row_ids",
}


class GSE256398TranscriptFormerEmbeddingError(RuntimeError):
    """Raised when sharded inference inputs or outputs differ."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pad_last_batch(batch: Any, row_count: int) -> Any:
    if row_count == BATCH_SIZE:
        return batch
    from transcriptformer.data.dataclasses import BatchData

    pad = BATCH_SIZE - row_count
    return BatchData(
        gene_counts=torch.cat(
            [batch.gene_counts, torch.zeros((pad, SEQUENCE_LENGTH), device="cuda")]
        ),
        gene_token_indices=torch.cat(
            [
                batch.gene_token_indices,
                torch.ones((pad, SEQUENCE_LENGTH), dtype=torch.int64, device="cuda"),
            ]
        ),
        aux_token_indices=torch.cat(
            [
                batch.aux_token_indices,
                torch.zeros((pad, 1), dtype=torch.int64, device="cuda"),
            ]
        ),
    )


def validate_shard_arrays(
    *,
    tokens: np.ndarray,
    counts: np.ndarray,
    assay: np.ndarray,
    folds: np.ndarray,
    row_ids: np.ndarray,
    source_mask: np.ndarray,
) -> None:
    rows = len(row_ids)
    if (
        rows < 1
        or tokens.shape != (rows, SEQUENCE_LENGTH)
        or counts.shape != tokens.shape
        or assay.shape != (rows, 1)
        or folds.shape != (rows,)
        or source_mask.shape != (rows,)
        or len(set(row_ids.astype(str).tolist())) != rows
        or int(tokens.min()) < 0
        or int(tokens.max()) >= 23_830
        or float(counts.min()) < 0.0
        or float(counts.max()) > 30.0
        or np.any(assay != 0)
        or source_mask.dtype != np.bool_
    ):
        raise GSE256398TranscriptFormerEmbeddingError("token shard content differs")


@torch.inference_mode()
def extract_arrays(
    model: Any,
    tokens: np.ndarray,
    counts: np.ndarray,
    assay: np.ndarray,
    *,
    check_repeat: bool,
) -> tuple[np.ndarray, float | None, float, float]:
    rows = tokens.shape[0]
    repeat_diff: float | None = None
    torch.cuda.reset_peak_memory_stats()
    if check_repeat:
        first = _batch(tokens, counts, assay, 0, min(BATCH_SIZE, rows))
        first = pad_last_batch(first, min(BATCH_SIZE, rows))
        with torch.autocast(device_type="cuda", dtype=torch.float16):
            repeat_a = _embedding_only_forward(model, first)
            repeat_b = _embedding_only_forward(model, first)
        repeat_diff = float(torch.max(torch.abs(repeat_a - repeat_b)).item())
        if not np.isfinite(repeat_diff) or repeat_diff > DETERMINISM_TOLERANCE:
            raise GSE256398TranscriptFormerEmbeddingError(
                f"deterministic repeat differs by {repeat_diff}"
            )

    embeddings = np.empty((rows, EMBEDDING_WIDTH), dtype=np.float32)
    started = time.monotonic()
    for start in range(0, rows, BATCH_SIZE):
        end = min(start + BATCH_SIZE, rows)
        current = _batch(tokens, counts, assay, start, end)
        current = pad_last_batch(current, end - start)
        with torch.autocast(device_type="cuda", dtype=torch.float16):
            values = _embedding_only_forward(model, current)
        embeddings[start:end] = values[: end - start].numpy()
    elapsed = time.monotonic() - started
    if not np.all(np.isfinite(embeddings)):
        raise GSE256398TranscriptFormerEmbeddingError("embeddings contain nonfinite values")
    peak_gib = torch.cuda.max_memory_allocated() / (1024**3)
    return embeddings, repeat_diff, peak_gib, elapsed


def run(
    *, checkpoint_root: Path, fixture_root: Path, fixture_sha256: str, output: Path
) -> dict[str, Any]:
    if output.exists():
        raise GSE256398TranscriptFormerEmbeddingError("embedding output already exists")
    verify_frozen_tree(fixture_root)
    if sha256_file(fixture_root / "ARTIFACTS.json") != fixture_sha256:
        raise GSE256398TranscriptFormerEmbeddingError("frozen fixture differs")
    fixture_receipt = json.loads(
        (fixture_root / "fixture_receipt.json").read_text(encoding="utf-8")
    )
    shards = fixture_receipt.get("shards", [])
    if (
        fixture_receipt.get("status") != "pass_outcome_blind_sharded_fixture"
        or fixture_receipt.get("rows") != EXPECTED_ROWS
        or fixture_receipt.get("donors") != EXPECTED_DONORS
        or fixture_receipt.get("phenotype_metadata_read") is not False
        or fixture_receipt.get("barcode_cell_labels_read") is not False
        or fixture_receipt.get("sealed_outcomes_read") is not False
        or len(shards) != EXPECTED_DONORS
    ):
        raise GSE256398TranscriptFormerEmbeddingError("fixture firewall differs")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise GSE256398TranscriptFormerEmbeddingError("exactly one CUDA device is required")
    device_name = torch.cuda.get_device_name(0)
    capability = list(torch.cuda.get_device_capability(0))
    if "L40S" not in device_name or capability != [8, 9]:
        raise GSE256398TranscriptFormerEmbeddingError("inference device is not L40S sm89")

    torch.manual_seed(20260824)
    np.random.seed(20260824)
    torch.set_float32_matmul_precision("high")
    model = _load_model(checkpoint_root).to("cuda")
    output.mkdir(parents=True, exist_ok=False)
    records: list[dict[str, Any]] = []
    policy_totals = {
        policy: {"elapsed_seconds": 0.0, "peak_gpu_memory_gib": 0.0, "rows": 0}
        for policy in POLICIES
    }
    repeat_by_policy: dict[str, float] = {}
    seen_row_ids: set[str] = set()
    for shard_index, shard in enumerate(shards):
        shard_path = fixture_root / shard["path"]
        if sha256_file(shard_path) != shard["sha256"]:
            raise GSE256398TranscriptFormerEmbeddingError("fixture shard hash differs")
        with np.load(shard_path, allow_pickle=False) as values:
            if set(values.files) != FIXTURE_FIELDS:
                raise GSE256398TranscriptFormerEmbeddingError("fixture shard fields differ")
            assay = values["assay_token_indices"]
            folds = values["outer_folds"]
            row_ids = values["row_ids"].astype(str)
            source_mask = values["in_source_exact"]
            if any(row_id in seen_row_ids for row_id in row_ids):
                raise GSE256398TranscriptFormerEmbeddingError("row ID crosses shards")
            seen_row_ids.update(row_ids.tolist())
            for policy in POLICIES:
                tokens = values[f"{policy}_gene_token_indices"]
                counts = values[f"{policy}_gene_counts"]
                validate_shard_arrays(
                    tokens=tokens,
                    counts=counts,
                    assay=assay,
                    folds=folds,
                    row_ids=row_ids,
                    source_mask=source_mask,
                )
                embeddings, repeat_diff, peak_gib, elapsed = extract_arrays(
                    model,
                    tokens,
                    counts,
                    assay,
                    check_repeat=shard_index == 0,
                )
                if repeat_diff is not None:
                    repeat_by_policy[policy] = repeat_diff
                destination = output / f"{shard['source_sample_id'].lower()}_{policy}_embeddings.npz"
                np.savez_compressed(
                    destination,
                    embeddings=embeddings,
                    in_source_exact=source_mask,
                    outer_folds=folds,
                    row_ids=row_ids,
                )
                records.append(
                    {
                        "source_sample_id": shard["source_sample_id"],
                        "policy": policy,
                        "rows": len(row_ids),
                        "source_exact_rows": int(source_mask.sum()),
                        "path": destination.name,
                        "sha256": sha256_file(destination),
                        "size_bytes": destination.stat().st_size,
                        "elapsed_seconds": elapsed,
                        "peak_gpu_memory_gib": peak_gib,
                    }
                )
                policy_totals[policy]["elapsed_seconds"] += elapsed
                policy_totals[policy]["peak_gpu_memory_gib"] = max(
                    policy_totals[policy]["peak_gpu_memory_gib"], peak_gib
                )
                policy_totals[policy]["rows"] += len(row_ids)
        print(
            json.dumps(
                {
                    "completed_donor": shard["source_sample_id"],
                    "completed_shards": shard_index + 1,
                    "total_shards": len(shards),
                },
                sort_keys=True,
            ),
            flush=True,
        )

    if len(seen_row_ids) != EXPECTED_ROWS or any(
        values["rows"] != EXPECTED_ROWS for values in policy_totals.values()
    ):
        raise GSE256398TranscriptFormerEmbeddingError("prediction row census differs")
    receipt = {
        "schema_version": "masld-bench-transcriptformer-gse256398-embeddings-v1",
        "status": "pass_outcome_blind_cross_cohort_embeddings",
        "dataset_id": "gse256398",
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "embedding_width": EMBEDDING_WIDTH,
        "policies": policy_totals,
        "repeat_max_abs_diff": repeat_by_policy,
        "prediction_shards": records,
        "fixture_artifacts_sha256": fixture_sha256,
        "device_name": device_name,
        "device_capability": capability,
        "torch_version": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "downstream_head_fit": False,
        "phenotype_metadata_read": False,
        "barcode_cell_labels_read": False,
        "sealed_outcomes_read": False,
        "supervised_accuracy_claim_allowed": False,
        "exposure_status": "unknown",
        "sealed_champion_eligible": False,
    }
    write_json_exclusive(output / "embedding_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "transcriptformer_gse256398_outcome_blind_embeddings",
            "dataset_id": "gse256398",
            "donors": EXPECTED_DONORS,
            "rows": EXPECTED_ROWS,
            "prediction_shards": len(records),
            "phenotype_metadata_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--fixture-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run(
        checkpoint_root=args.checkpoint_root,
        fixture_root=args.fixture_root,
        fixture_sha256=args.fixture_sha256,
        output=args.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
