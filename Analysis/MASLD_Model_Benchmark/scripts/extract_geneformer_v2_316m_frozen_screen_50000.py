#!/usr/bin/env python3
"""Extract raw V2-316M embeddings for the outcome-blind 50k fixture."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any

from scripts.preflight_geneformer_v2_316m_common import (
    EMBEDDING_POLICY,
    HIDDEN_SIZE,
    load_prior_compatible_model,
    raw_mean_pool,
)


ROWS = 50_000
MAX_TOKENS = 4096


class GeneformerExtractionError(RuntimeError):
    """Raised when the frozen fixture or raw embedding requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_fixture(root: Path, expected_rows: int = ROWS) -> dict[str, Any]:
    import numpy as np

    receipt = json.loads((root / "fixture_receipt.json").read_text(encoding="utf-8"))
    if (
        receipt.get("status") != "pass_outcome_blind_fixture"
        or receipt.get("rows") != expected_rows
        or receipt.get("embedding_policy") != EMBEDDING_POLICY
        or receipt.get("max_tokens") != MAX_TOKENS
        or receipt.get("evaluation_label_columns_read")
        or receipt.get("histology_columns_read")
        or receipt.get("sealed_outcomes_read") is not False
    ):
        raise GeneformerExtractionError("outcome-blind 50k fixture receipt differs")
    tokens = np.load(root / "token_ids.npy", mmap_mode="r", allow_pickle=False)
    offsets = np.load(root / "token_offsets.npy", mmap_mode="r", allow_pickle=False)
    rows = (root / "embedding_row_order.txt").read_text(encoding="utf-8").splitlines()
    if (
        tokens.dtype != np.dtype("int32")
        or tokens.ndim != 1
        or offsets.dtype != np.dtype("int64")
        or offsets.shape != (expected_rows + 1,)
        or int(offsets[0]) != 0
        or int(offsets[-1]) != len(tokens)
        or bool(np.any(np.diff(offsets) <= 0))
        or int(np.diff(offsets).max()) > MAX_TOKENS
        or len(rows) != expected_rows
        or len(set(rows)) != expected_rows
    ):
        raise GeneformerExtractionError("outcome-blind token/row arrays differ")
    observed_row_hash = sha256(
        "\n".join(rows).encode("utf-8") + b"\n"
    ).hexdigest()
    if observed_row_hash != receipt.get("row_order_sha256"):
        raise GeneformerExtractionError("fixture row order hash differs")
    return {"receipt": receipt, "tokens": tokens, "offsets": offsets, "rows": rows}


def extract(
    fixture_root: Path,
    model_directory: Path,
    output: Path,
    *,
    expected_rows: int = ROWS,
    batch_size: int = 32,
) -> dict[str, Any]:
    import numpy as np
    import torch

    if output.exists() or batch_size < 1 or batch_size > 32:
        raise GeneformerExtractionError("output/batch-size contract differs")
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise GeneformerExtractionError("exactly one CUDA device is required")
    if "L40S" not in torch.cuda.get_device_name(0):
        raise GeneformerExtractionError("50k extraction requires an L40S")
    fixture = load_fixture(fixture_root, expected_rows=expected_rows)
    output.mkdir(parents=True)
    temporary_embeddings = output / ".embeddings.npy.incomplete"
    embeddings = np.lib.format.open_memmap(
        temporary_embeddings,
        mode="w+",
        dtype="float32",
        shape=(expected_rows, HIDDEN_SIZE),
    )
    torch.manual_seed(20260824)
    torch.cuda.manual_seed_all(20260824)
    torch.use_deterministic_algorithms(True)
    torch.cuda.reset_peak_memory_stats()
    model, loading = load_prior_compatible_model(model_directory)
    model.eval().to("cuda")
    tokens = fixture["tokens"]
    offsets = fixture["offsets"]
    with torch.inference_mode():
        for start in range(0, expected_rows, batch_size):
            stop = min(start + batch_size, expected_rows)
            sequences = [
                tokens[int(offsets[index]) : int(offsets[index + 1])]
                for index in range(start, stop)
            ]
            width = max(len(item) for item in sequences)
            ids = torch.zeros((stop - start, width), dtype=torch.long, device="cuda")
            attention = torch.zeros_like(ids)
            for index, values in enumerate(sequences):
                length = len(values)
                ids[index, :length] = torch.as_tensor(
                    np.array(values, dtype=np.int64, copy=True),
                    dtype=torch.long,
                    device="cuda",
                )
                attention[index, :length] = 1
            hidden = model(input_ids=ids, attention_mask=attention).last_hidden_state
            pooled = raw_mean_pool(hidden, attention)
            embeddings[start:stop] = pooled.cpu().numpy().astype("float32", copy=False)
    embeddings.flush()
    del embeddings
    os.rename(temporary_embeddings, output / "embeddings.npy")
    (output / "embedding_row_order.txt").write_text(
        "\n".join(fixture["rows"]) + "\n", encoding="utf-8"
    )
    observed = np.load(output / "embeddings.npy", mmap_mode="r", allow_pickle=False)
    if observed.shape != (expected_rows, HIDDEN_SIZE) or not np.isfinite(observed).all():
        raise GeneformerExtractionError("raw embedding output differs")
    receipt = {
        "schema_version": "masld-bench-geneformer-v2-316m-frozen-screen-extraction-v1",
        "status": "pass_raw_outcome_blind_extraction",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "rows": expected_rows,
        "embedding_shape": list(observed.shape),
        "embedding_dtype": str(observed.dtype),
        "embedding_policy": EMBEDDING_POLICY,
        "embedding_sha256": sha256_file(output / "embeddings.npy"),
        "row_order_sha256": sha256_file(output / "embedding_row_order.txt"),
        "fixture_row_order_sha256": fixture["receipt"]["row_order_sha256"],
        "batch_size": batch_size,
        "device": torch.cuda.get_device_name(0),
        "cuda": torch.version.cuda,
        "cuda_peak_allocated_bytes": int(torch.cuda.max_memory_allocated()),
        "cuda_peak_reserved_bytes": int(torch.cuda.max_memory_reserved()),
        "auto_model_restore": loading,
        "raw_unrefined_embeddings_only": True,
        "downstream_head_fit": False,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    (output / "extraction_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--model-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-rows", type=int, default=ROWS)
    parser.add_argument("--batch-size", type=int, default=32)
    arguments = parser.parse_args()
    result = extract(
        arguments.fixture_root,
        arguments.model_directory,
        arguments.output,
        expected_rows=arguments.expected_rows,
        batch_size=arguments.batch_size,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "embedding_shape")}))


if __name__ == "__main__":
    main()
