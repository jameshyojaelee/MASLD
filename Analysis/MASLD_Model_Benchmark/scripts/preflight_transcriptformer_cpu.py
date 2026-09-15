#!/usr/bin/env python3
"""CPU-only inclusion preflight for pinned TF-Sapiens source, state, and fixtures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch

from scripts.transcriptformer_extract_unlabeled_embeddings import (
    BATCH_SIZE,
    CHECKPOINT_SHA256,
    EMBEDDING_WIDTH,
    POLICIES,
    ROWS,
    SEQUENCE_LENGTH,
    _load_model,
    _read_rows,
    sha256_file,
)


class TranscriptFormerCPUPreflightError(ValueError):
    """Raised when the CPU inclusion requirement differs."""


def preflight(
    checkpoint_root: Path,
    fixture_root: Path,
    output: Path,
    *,
    expected_rows: int = ROWS,
) -> dict:
    if output.exists():
        raise TranscriptFormerCPUPreflightError("output already exists")
    receipt = json.loads((fixture_root / "fixture_receipt.json").read_text())
    if (
        receipt.get("status") != "pass_outcome_blind_fixture"
        or receipt.get("evaluation_label_columns_read") != []
        or receipt.get("sealed_outcomes_read") is not False
        or receipt.get("rows") != expected_rows
    ):
        raise TranscriptFormerCPUPreflightError("fixture firewall differs")
    rows = _read_rows(fixture_root / "row_contract.tsv", expected_rows)
    token_records = {}
    for policy in POLICIES:
        path = fixture_root / f"{policy}_tokens.npz"
        with np.load(path, allow_pickle=False) as values:
            if (
                values["gene_token_indices"].shape
                != (expected_rows, SEQUENCE_LENGTH)
                or values["gene_counts"].shape
                != (expected_rows, SEQUENCE_LENGTH)
                or values["assay_token_indices"].shape != (expected_rows, 1)
                or values["row_ids"].astype(str).tolist()
                != [row["row_id"] for row in rows]
                or values["outer_folds"].tolist()
                != [int(row["outer_fold"]) for row in rows]
            ):
                raise TranscriptFormerCPUPreflightError(f"{policy} token fixture differs")
        token_records[policy] = sha256_file(path)

    model = _load_model(checkpoint_root)
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    if (
        parameter_count != 429870637
        or model.gene_vocab.vocab_size != 23830
        or model.model_config.seq_len != SEQUENCE_LENGTH
        or model.model_config.embed_dim != EMBEDDING_WIDTH
        or model.inference_config.batch_size != BATCH_SIZE
    ):
        raise TranscriptFormerCPUPreflightError("strict-loaded architecture differs")
    result = {
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "checkpoint_strict_load": True,
        "development_rows": expected_rows,
        "downstream_head_fit": False,
        "evaluation_label_columns_read": [],
        "exposure_status": "unknown",
        "fixture_token_sha256": token_records,
        "model_parameter_count": parameter_count,
        "schema_version": "masld-bench-transcriptformer-cpu-preflight-v1",
        "sealed_champion_eligible": False,
        "sealed_outcomes_read": False,
        "status": "pass_cpu_preflight",
        "torch_cuda_available": bool(torch.cuda.is_available()),
        "torch_version": torch.__version__,
        "vocabulary_rows": model.gene_vocab.vocab_size,
    }
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", required=True, type=Path)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-rows", type=int, default=ROWS)
    arguments = parser.parse_args()
    result = preflight(
        arguments.checkpoint_root,
        arguments.fixture_root,
        arguments.output,
        expected_rows=arguments.expected_rows,
    )
    print(json.dumps({key: result[key] for key in ("status", "model_parameter_count")}))


if __name__ == "__main__":
    main()
