#!/usr/bin/env python3
"""Extract outcome-blind TF-Sapiens cell embeddings from frozen token fixtures."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import time
from typing import Any

import numpy as np
import torch
from torch.nn.attention.flex_attention import and_masks, create_block_mask

from transcriptformer.data.dataclasses import (
    BatchData,
    DataConfig,
    InferenceConfig,
    LossConfig,
    ModelConfig,
)
from transcriptformer.model.layers import mean_embeddings
from transcriptformer.model.masks import pad_mask_factory
from transcriptformer.model.model import Transcriptformer
from transcriptformer.tokenizer.vocab import construct_gene_embeddings


CHECKPOINT_SHA256 = "eff027e7393f92a32cadcbb2e13a917e24b62b71341571bf40f38ff6ec6d3021"
CONFIG_SHA256 = "2413d512e8fdd550389bc574dc465357e943e23b4411cab6d23423af08095c2f"
GENE_VOCAB_SHA256 = "676640be87ba4ba7ea7af8d9dc470a8ae8be02768989478027dec5fe61f465e2"
ASSAY_VOCAB_SHA256 = "0c405e4ead45a4b8350d8e874f834273efdd3f7ad4669b52d3e3727fb4fe70af"
SPECIAL_TOKENS = ["unknown", "[PAD]", "[START]", "[END]", "[RD]", "[CELL]", "[MASK]"]
POLICIES = ("native", "common")
ROWS = 1000
SEQUENCE_LENGTH = 2047
EMBEDDING_WIDTH = 2048
BATCH_SIZE = 8
DETERMINISM_TOLERANCE = 1.0e-3


class TranscriptFormerEmbeddingError(ValueError):
    """Raised when frozen inference inputs or outputs differ."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _strip_target(value: dict[str, Any]) -> dict[str, Any]:
    return {key: item for key, item in value.items() if key != "_target_"}


def _read_rows(path: Path, expected_rows: int = ROWS) -> list[dict[str, str]]:
    if expected_rows < 1:
        raise TranscriptFormerEmbeddingError("expected row count is invalid")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["row_position", "row_id", "donor_id", "dataset_id", "outer_fold"]
        if reader.fieldnames != expected:
            raise TranscriptFormerEmbeddingError("row contract fields differ")
        rows = list(reader)
    if (
        len(rows) != expected_rows
        or len({row["row_id"] for row in rows}) != expected_rows
    ):
        raise TranscriptFormerEmbeddingError("row contract cardinality differs")
    if [int(row["row_position"]) for row in rows] != list(range(expected_rows)):
        raise TranscriptFormerEmbeddingError("row contract order differs")
    fold_by_donor: dict[str, int] = {}
    for row in rows:
        fold = int(row["outer_fold"])
        previous = fold_by_donor.setdefault(row["donor_id"], fold)
        if previous != fold:
            raise TranscriptFormerEmbeddingError("one donor crosses outer folds")
    return rows


def _load_model(checkpoint_root: Path) -> Transcriptformer:
    config_path = checkpoint_root / "config.json"
    checkpoint_path = checkpoint_root / "model_weights.pt"
    gene_vocab_path = checkpoint_root / "vocabs/homo_sapiens_gene.h5"
    assay_vocab_path = checkpoint_root / "vocabs/assay_vocab.json"
    expected = {
        config_path: CONFIG_SHA256,
        checkpoint_path: CHECKPOINT_SHA256,
        gene_vocab_path: GENE_VOCAB_SHA256,
        assay_vocab_path: ASSAY_VOCAB_SHA256,
    }
    for path, checksum in expected.items():
        if not path.is_file() or sha256_file(path) != checksum:
            raise TranscriptFormerEmbeddingError(f"checkpoint artifact differs: {path.name}")

    config = json.loads(config_path.read_text(encoding="utf-8"))["model"]
    data_values = _strip_target(config["data_config"])
    data_values.update(
        clip_counts=30,
        filter_to_vocabs=True,
        filter_outliers=0.0,
        sort_genes=False,
        randomize_genes=False,
        min_expressed_genes=0,
        esm2_mappings_path=(checkpoint_root / "vocabs").as_posix(),
        aux_vocab_path=(checkpoint_root / "vocabs").as_posix(),
    )
    inference = InferenceConfig(
        output_keys=["embeddings"],
        batch_size=BATCH_SIZE,
        obs_keys=[],
        data_files=None,
        load_checkpoint=checkpoint_path.as_posix(),
        output_path=None,
        device="cuda",
        precision="16-mixed",
        emb_type="cell",
    )
    genes, embedding_matrix = construct_gene_embeddings(
        [gene_vocab_path.as_posix()], SPECIAL_TOKENS
    )
    assay_vocab = json.loads(assay_vocab_path.read_text(encoding="utf-8"))
    if len(genes) != 23830 or len(assay_vocab) != 40 or assay_vocab.get("unknown") != 0:
        raise TranscriptFormerEmbeddingError("checkpoint vocabulary cardinality differs")
    model = Transcriptformer(
        data_config=DataConfig(**data_values),
        model_config=ModelConfig(**_strip_target(config["model_config"])),
        loss_config=LossConfig(**_strip_target(config["loss_config"])),
        inference_config=inference,
        gene_vocab_dict=genes,
        aux_vocab_dict={"assay": assay_vocab},
        emb_matrix=torch.as_tensor(embedding_matrix),
    )
    state = torch.load(
        checkpoint_path,
        weights_only=True,
        mmap=True,
        map_location="cpu",
    )
    model.load_state_dict(state, strict=True)
    model.eval()
    return model


@torch.inference_mode()
def _embedding_only_forward(model: Transcriptformer, batch: BatchData) -> torch.Tensor:
    gene_tokens = batch.gene_token_indices
    gene_counts = batch.gene_counts
    aux_tokens = batch.aux_token_indices
    right_shifted_tokens = torch.cat(
        [model.gene_vocab.start_idx * torch.ones_like(gene_tokens[:, :1]), gene_tokens[:, :-1]],
        dim=1,
    )
    right_shifted_counts = torch.cat(
        [torch.ones_like(gene_counts[:, :1]), gene_counts[:, :-1]], dim=1
    )
    right_shifted_embeddings, _ = model.get_gene_embeddings(batch)
    if aux_tokens is not None:
        aux_embeddings = model.get_aux_embeddings(aux_tokens)
        right_shifted_embeddings = torch.cat(
            [aux_embeddings, right_shifted_embeddings], dim=1
        )
        right_shifted_counts = torch.cat(
            [torch.ones_like(aux_tokens), right_shifted_counts], dim=1
        )
    pad_mask = model._pad_mask(right_shifted_tokens, aux_tokens, dtype="bool")
    mask_mod = and_masks(pad_mask_factory(pad_mask), model.causal_mask)
    block_mask = create_block_mask(
        mask_mod,
        right_shifted_embeddings.shape[0],
        H=None,
        Q_LEN=right_shifted_embeddings.shape[1],
        KV_LEN=right_shifted_embeddings.shape[1],
        device=right_shifted_embeddings.device,
        BLOCK_SIZE=model.model_config.block_len,
        _compile=model.model_config.compile_block_mask,
    )
    score_mod = model._score_mod_factory(
        torch.log1p(right_shifted_counts + model.model_config.log_counts_eps),
        emb_mode=True,
    )
    output = model.transformer_encoder(
        x=right_shifted_embeddings,
        score_mod=score_mod,
        block_mask=block_mask,
    )
    if aux_tokens is not None:
        output = output[:, aux_tokens.shape[1] :, :]
        pad_mask = pad_mask[:, aux_tokens.shape[1] :]
    return mean_embeddings(output, pad_mask).detach().float().cpu()


def _batch(tokens: np.ndarray, counts: np.ndarray, assay: np.ndarray, start: int, end: int) -> BatchData:
    return BatchData(
        gene_counts=torch.from_numpy(counts[start:end]).to("cuda", non_blocking=True),
        gene_token_indices=torch.from_numpy(tokens[start:end].astype(np.int64)).to(
            "cuda", non_blocking=True
        ),
        aux_token_indices=torch.from_numpy(assay[start:end]).to("cuda", non_blocking=True),
    )


def _extract_policy(
    model: Transcriptformer,
    fixture: Path,
    rows: list[dict[str, str]],
    expected_rows: int,
) -> tuple[np.ndarray, float, float, float]:
    with np.load(fixture, allow_pickle=False) as values:
        required = {
            "assay_token_indices",
            "gene_counts",
            "gene_token_indices",
            "outer_folds",
            "row_ids",
        }
        if set(values.files) != required:
            raise TranscriptFormerEmbeddingError("token fixture fields differ")
        tokens = values["gene_token_indices"]
        counts = values["gene_counts"]
        assay = values["assay_token_indices"]
        folds = values["outer_folds"]
        row_ids = values["row_ids"].astype(str)
    if (
        tokens.shape != (expected_rows, SEQUENCE_LENGTH)
        or counts.shape != tokens.shape
        or assay.shape != (expected_rows, 1)
        or row_ids.tolist() != [row["row_id"] for row in rows]
        or folds.tolist() != [int(row["outer_fold"]) for row in rows]
        or int(tokens.min()) < 0
        or int(tokens.max()) >= 23830
        or float(counts.min()) < 0.0
        or float(counts.max()) > 30.0
        or np.any(assay != 0)
    ):
        raise TranscriptFormerEmbeddingError("token fixture content differs")

    torch.cuda.reset_peak_memory_stats()
    first = _batch(tokens, counts, assay, 0, BATCH_SIZE)
    with torch.autocast(device_type="cuda", dtype=torch.float16):
        repeat_a = _embedding_only_forward(model, first)
        repeat_b = _embedding_only_forward(model, first)
    repeat_diff = float(torch.max(torch.abs(repeat_a - repeat_b)).item())
    if not np.isfinite(repeat_diff) or repeat_diff > DETERMINISM_TOLERANCE:
        raise TranscriptFormerEmbeddingError(
            f"deterministic repeat differs by {repeat_diff}"
        )

    embeddings = np.empty((expected_rows, EMBEDDING_WIDTH), dtype=np.float32)
    started = time.monotonic()
    for start in range(0, expected_rows, BATCH_SIZE):
        end = min(start + BATCH_SIZE, expected_rows)
        current = _batch(tokens, counts, assay, start, end)
        if end - start < BATCH_SIZE:
            pad = BATCH_SIZE - (end - start)
            current = BatchData(
                gene_counts=torch.cat(
                    [current.gene_counts, torch.zeros((pad, SEQUENCE_LENGTH), device="cuda")]
                ),
                gene_token_indices=torch.cat(
                    [
                        current.gene_token_indices,
                        torch.ones((pad, SEQUENCE_LENGTH), dtype=torch.int64, device="cuda"),
                    ]
                ),
                aux_token_indices=torch.cat(
                    [current.aux_token_indices, torch.zeros((pad, 1), dtype=torch.int64, device="cuda")]
                ),
            )
        with torch.autocast(device_type="cuda", dtype=torch.float16):
            output = _embedding_only_forward(model, current)
        embeddings[start:end] = output[: end - start].numpy()
        if start % (10 * BATCH_SIZE) == 0:
            print(f"embedded {start}/{expected_rows}", flush=True)
    elapsed = time.monotonic() - started
    if not np.all(np.isfinite(embeddings)):
        raise TranscriptFormerEmbeddingError("embeddings contain nonfinite values")
    peak_gib = torch.cuda.max_memory_allocated() / (1024**3)
    return embeddings, repeat_diff, peak_gib, elapsed


def run(
    checkpoint_root: Path,
    fixture_root: Path,
    output: Path,
    *,
    expected_rows: int = ROWS,
) -> dict[str, Any]:
    if output.exists():
        raise TranscriptFormerEmbeddingError("output already exists")
    fixture_receipt = json.loads(
        (fixture_root / "fixture_receipt.json").read_text(encoding="utf-8")
    )
    if (
        fixture_receipt.get("status") != "pass_outcome_blind_fixture"
        or fixture_receipt.get("rows") != expected_rows
        or fixture_receipt.get("evaluation_label_columns_read") != []
        or fixture_receipt.get("sealed_outcomes_read") is not False
    ):
        raise TranscriptFormerEmbeddingError("fixture firewall differs")
    rows = _read_rows(fixture_root / "row_contract.tsv", expected_rows)
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise TranscriptFormerEmbeddingError("exactly one CUDA device is required")
    device_name = torch.cuda.get_device_name(0)
    capability = list(torch.cuda.get_device_capability(0))
    if "L40S" not in device_name or capability != [8, 9]:
        raise TranscriptFormerEmbeddingError("inference device is not L40S sm89")

    torch.manual_seed(20260824)
    np.random.seed(20260824)
    torch.set_float32_matmul_precision("high")
    model = _load_model(checkpoint_root).to("cuda")
    output.mkdir(parents=True)
    summaries = {}
    artifact_records = {}
    for policy in POLICIES:
        embeddings, repeat_diff, peak_gib, elapsed = _extract_policy(
            model,
            fixture_root / f"{policy}_tokens.npz",
            rows,
            expected_rows,
        )
        destination = output / f"{policy}_embeddings.npz"
        np.savez_compressed(
            destination,
            embeddings=embeddings,
            outer_folds=np.asarray([int(row["outer_fold"]) for row in rows], dtype=np.int8),
            row_ids=np.asarray([row["row_id"] for row in rows]),
        )
        artifact_records[destination.name] = {
            "sha256": sha256_file(destination),
            "size_bytes": destination.stat().st_size,
        }
        summaries[policy] = {
            "elapsed_seconds": elapsed,
            "peak_gpu_memory_gib": peak_gib,
            "repeat_max_abs_diff": repeat_diff,
        }
    receipt = {
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "common_and_native_reported_separately": True,
        "device_capability": capability,
        "device_name": device_name,
        "downstream_head_fit": False,
        "embedding_policy": "official_mean_over_full_2047_axis_with_padding_mask_zeroed",
        "embedding_width": EMBEDDING_WIDTH,
        "evaluation_label_columns_read": [],
        "exposure_status": "unknown",
        "policies": summaries,
        "rows": expected_rows,
        "schema_version": "masld-bench-transcriptformer-unlabeled-embeddings-v1",
        "sealed_champion_eligible": False,
        "sealed_outcomes_read": False,
        "status": "pass_outcome_blind_embeddings",
        "torch_cuda_build": torch.version.cuda,
        "torch_version": torch.__version__,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    artifact_records["receipt.json"] = {
        "sha256": sha256_file(output / "receipt.json"),
        "size_bytes": (output / "receipt.json").stat().st_size,
    }
    (output / "ARTIFACTS.json").write_text(
        json.dumps(
            {
                "artifacts": artifact_records,
                "schema_version": "masld-bench-artifact-manifest-v1",
                "status": "pass",
            },
            sort_keys=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-root", required=True, type=Path)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-rows", type=int, default=ROWS)
    arguments = parser.parse_args()
    receipt = run(
        arguments.checkpoint_root,
        arguments.fixture_root,
        arguments.output,
        expected_rows=arguments.expected_rows,
    )
    print(json.dumps({key: receipt[key] for key in ("status", "rows", "device_name")}))


if __name__ == "__main__":
    main()
