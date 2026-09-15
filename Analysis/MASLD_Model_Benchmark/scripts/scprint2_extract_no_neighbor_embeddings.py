#!/usr/bin/env python3
"""Extract label-blind scPRINT-2 embeddings without neighbor coupling."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader

from scdataloader import Collator
from scdataloader.data import SimpleAnnDataset
from scprint2 import scPRINT2


CHECKPOINT_SHA256 = "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
CHECKPOINT_SIZE = 889_706_878
FIXTURE_ARTIFACTS_SHA256 = "edf71c212340ebfdcc266ef7c0ecb73449193b36a17147e53d3801804082c9f7"
RUNTIME_ARTIFACTS_SHA256 = "4474a0d069d09a9a805e7d6e33936af95eeab9fb72d6e3c131e7e7733350842d"
PREFLIGHT_ARTIFACTS_SHA256 = "9ed6d6db96bca280b2b2a99f8c3256be329978ad145e234aaac4815ef33295fa"
HUMAN = "NCBITaxon:9606"
ROWS = 1000
DONORS = 102
BATCH_SIZE = 4
SEED = 20260824
POLICIES = {
    "common": ("most expr", 3200),
    "random_checkpoint_max": ("random expr", 3200),
    "released_selection_default_no_neighbor": ("random expr", 2000),
}
DETERMINISM_TOLERANCE = 1.0e-4


class Scprint2EmbeddingError(RuntimeError):
    """Raised when a frozen input or inference requirement differs."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def verify_manifest(path: Path, expected: str) -> None:
    if not path.is_file() or digest(path) != expected:
        raise Scprint2EmbeddingError(f"registered artifact differs: {path}")


def load_model(checkpoint: Path) -> tuple[scPRINT2, dict[str, list[str]], list[str]]:
    if checkpoint.stat().st_size != CHECKPOINT_SIZE or digest(checkpoint) != CHECKPOINT_SHA256:
        raise Scprint2EmbeddingError("checkpoint identity differs")
    payload = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
    state = payload.get("state_dict")
    hyperparameters = payload.get("hyper_parameters")
    if not isinstance(state, dict) or not isinstance(hyperparameters, dict):
        raise Scprint2EmbeddingError("checkpoint state contract differs")
    constructor = dict(hyperparameters)
    if constructor.pop("_instantiator", None) is None:
        raise Scprint2EmbeddingError("Lightning instantiator metadata differs")
    constructor["precpt_gene_emb"] = None
    constructor["gene_pos_file"] = None
    model = scPRINT2(**constructor)
    model.on_load_checkpoint(payload)
    incompatibility = model.load_state_dict(state, strict=False)
    if incompatibility.missing_keys or incompatibility.unexpected_keys:
        raise Scprint2EmbeddingError("strict checkpoint restoration differs")
    genes = hyperparameters.get("genes")
    organisms = hyperparameters.get("organisms")
    if not isinstance(genes, dict) or not isinstance(organisms, list) or HUMAN not in genes:
        raise Scprint2EmbeddingError("checkpoint vocabulary contract differs")
    return model.eval(), genes, organisms


def make_collator(
    genes: dict[str, list[str]], organisms: list[str], how: str, max_len: int
) -> Collator:
    ordered_genes = [gene for organism in organisms for gene in genes[organism]]
    organism_column = [organism for organism in organisms for _ in genes[organism]]
    gene_frame = pd.DataFrame({"organism": organism_column}, index=ordered_genes)
    collator = Collator(
        organisms=organisms,
        valid_genes=ordered_genes,
        how=how,
        max_len=max_len,
        add_zero_genes=0,
        genedf=gene_frame,
    )
    # scdataloader 2.1.0 omits this assignment with an offline gene table.
    collator.organism_ids = set(organisms)
    return collator


def make_loader(
    dataset: SimpleAnnDataset,
    genes: dict[str, list[str]],
    organisms: list[str],
    how: str,
    max_len: int,
) -> DataLoader:
    return DataLoader(
        dataset,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=0,
        collate_fn=make_collator(genes, organisms, how, max_len),
    )


@torch.inference_mode()
def embed_batch(
    model: scPRINT2,
    batch: dict[str, torch.Tensor],
    source_depth: np.ndarray,
) -> torch.Tensor:
    if "knn_cells" in batch or "knn_cells_info" in batch:
        raise Scprint2EmbeddingError("neighbor tensors entered no-neighbor inference")
    gene_pos = batch["genes"].to("cuda", non_blocking=True)
    expression = batch["x"].to("cuda", non_blocking=True)
    aligned_depth = batch["depth"].numpy().astype(np.float64)
    if source_depth.shape != aligned_depth.shape or np.any(source_depth < aligned_depth):
        raise Scprint2EmbeddingError("prefilter sequencing-depth contract differs")
    depth = torch.from_numpy(source_depth.astype(np.float32)).to(
        "cuda", non_blocking=True
    )
    result = model._predict(
        gene_pos,
        expression,
        depth,
        knn_cells=None,
        knn_cells_info=None,
        pred_embedding=["all"],
        keep_output=False,
    )
    pieces = result["embs"]
    names = ["other", *model.classes]
    if list(pieces) != names:
        raise Scprint2EmbeddingError("official embedding component order differs")
    if model.compressor is None:
        output = torch.stack([pieces[name] for name in names], dim=1).mean(dim=1)
    else:
        output = torch.cat([pieces[name] for name in names], dim=1)
    return output.detach().float().cpu()


def extract_policy(
    model: scPRINT2,
    dataset: SimpleAnnDataset,
    genes: dict[str, list[str]],
    organisms: list[str],
    how: str,
    max_len: int,
    source_depths: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    first_a = next(iter(make_loader(dataset, genes, organisms, how, max_len)))
    with torch.autocast(device_type="cuda", dtype=torch.float16):
        repeat_a = embed_batch(model, first_a, source_depths[:BATCH_SIZE])
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    first_b = next(iter(make_loader(dataset, genes, organisms, how, max_len)))
    with torch.autocast(device_type="cuda", dtype=torch.float16):
        repeat_b = embed_batch(model, first_b, source_depths[:BATCH_SIZE])
    repeat_diff = float(torch.max(torch.abs(repeat_a - repeat_b)).item())
    if not np.isfinite(repeat_diff) or repeat_diff > DETERMINISM_TOLERANCE:
        raise Scprint2EmbeddingError(f"deterministic repeat differs by {repeat_diff}")

    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.reset_peak_memory_stats()
    outputs: list[np.ndarray] = []
    started = time.monotonic()
    offset = 0
    for index, batch in enumerate(
        make_loader(dataset, genes, organisms, how, max_len)
    ):
        batch_rows = int(batch["x"].shape[0])
        with torch.autocast(device_type="cuda", dtype=torch.float16):
            output = embed_batch(
                model, batch, source_depths[offset : offset + batch_rows]
            )
        outputs.append(output.numpy())
        offset += batch_rows
        if index % 25 == 0:
            print(f"{how}: {index * BATCH_SIZE}/{ROWS}", flush=True)
    embeddings = np.concatenate(outputs, axis=0)
    if offset != ROWS or embeddings.shape[0] != ROWS or not np.isfinite(embeddings).all():
        raise Scprint2EmbeddingError("embedding output differs")
    return embeddings, {
        "elapsed_seconds": time.monotonic() - started,
        "embedding_width": int(embeddings.shape[1]),
        "max_len": max_len,
        "peak_gpu_memory_gib": torch.cuda.max_memory_allocated() / (1024**3),
        "repeat_max_abs_diff": repeat_diff,
    }


def run(
    checkpoint: Path,
    fixture_root: Path,
    runtime_root: Path,
    preflight_root: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise Scprint2EmbeddingError("output already exists")
    verify_manifest(fixture_root / "ARTIFACTS.json", FIXTURE_ARTIFACTS_SHA256)
    verify_manifest(runtime_root / "ARTIFACTS.json", RUNTIME_ARTIFACTS_SHA256)
    verify_manifest(preflight_root / "ARTIFACTS.json", PREFLIGHT_ARTIFACTS_SHA256)
    receipt = json.loads((fixture_root / "fixture/fixture_receipt.json").read_text())
    if (
        receipt.get("status") != "pass"
        or receipt.get("rows") != ROWS
        or receipt.get("donors") != DONORS
        or receipt.get("evaluation_label_columns_read") != []
        or receipt.get("sealed_outcomes_read") is not False
        or receipt.get("neighbor_tensors_allowed") is not False
        or receipt.get("source_depth_preserved_for_inference") is not True
    ):
        raise Scprint2EmbeddingError("fixture firewall differs")
    adata = ad.read_h5ad(fixture_root / "fixture/aligned_counts.h5ad")
    if adata.n_obs != ROWS or adata.n_vars != 20_004:
        raise Scprint2EmbeddingError("fixture dimensions differ")
    row_ids = np.asarray(adata.obs["row_id"].astype(str).to_numpy(), dtype=np.str_)
    folds = adata.obs["outer_fold"].astype(int).to_numpy()
    source_depths = adata.obs["source_total_count"].to_numpy(dtype=np.float64)
    if len(set(row_ids)) != ROWS or set(folds) != set(range(5)):
        raise Scprint2EmbeddingError("fixture row or fold contract differs")
    aligned_depths = np.asarray(adata.X.sum(axis=1)).ravel().astype(np.float64)
    if (
        source_depths.shape != (ROWS,)
        or np.any(source_depths < aligned_depths)
        or int(np.sum(source_depths != aligned_depths))
        != receipt.get("cells_with_out_of_vocabulary_counts")
    ):
        raise Scprint2EmbeddingError("source sequencing-depth fixture differs")
    dataset = SimpleAnnDataset(
        adata,
        obs_to_output=["organism_ontology_term_id"],
        get_knn_cells=False,
    )
    if not torch.cuda.is_available() or torch.cuda.device_count() != 1:
        raise Scprint2EmbeddingError("exactly one CUDA device is required")
    device_name = torch.cuda.get_device_name(0)
    capability = list(torch.cuda.get_device_capability(0))
    if "L40S" not in device_name or capability != [8, 9]:
        raise Scprint2EmbeddingError("inference device is not L40S sm89")
    torch.set_float32_matmul_precision("high")
    model, genes, organisms = load_model(checkpoint)
    model = model.to("cuda")

    output.mkdir(parents=True)
    policy_records = {}
    artifacts = {}
    widths = set()
    for policy, (how, max_len) in POLICIES.items():
        embeddings, record = extract_policy(
            model,
            dataset,
            genes,
            organisms,
            how,
            max_len,
            source_depths,
        )
        widths.add(record["embedding_width"])
        destination = output / f"{policy}_embeddings.npz"
        np.savez_compressed(
            destination,
            embeddings=embeddings,
            outer_folds=folds.astype(np.int8),
            row_ids=row_ids,
        )
        artifacts[destination.name] = {
            "sha256": digest(destination),
            "size_bytes": destination.stat().st_size,
        }
        policy_records[policy] = {"collator_policy": how, **record}
    if len(widths) != 1:
        raise Scprint2EmbeddingError("native and common embedding widths differ")
    result = {
        "schema_version": "masld-bench-scprint2-no-neighbor-embeddings-v1",
        "status": "pass_outcome_blind_embeddings",
        "model_id": "scprint2_small_v2",
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "rows": ROWS,
        "donors": DONORS,
        "embedding_width": widths.pop(),
        "embedding_policy": "released output_cell_emb; released random-expr/max_len defaults with neighbor coupling disabled, checkpoint-max random-expr 3200, and deterministic checkpoint-max most-expr 3200 reported separately",
        "policies": policy_records,
        "source_prefilter_depth_used": True,
        "aligned_count_fraction": receipt["aligned_count_fraction"],
        "cells_with_out_of_vocabulary_counts": receipt[
            "cells_with_out_of_vocabulary_counts"
        ],
        "neighbor_tensors_used": False,
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
        "exposure_status": "unknown",
        "sealed_champion_eligible": False,
        "device_name": device_name,
        "device_capability": capability,
        "torch_version": torch.__version__,
        "torch_cuda_build": torch.version.cuda,
        "source_runtime": "scPRINT-2 1.0.3 exact uv.lock",
    }
    receipt_path = output / "receipt.json"
    receipt_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    artifacts[receipt_path.name] = {
        "sha256": digest(receipt_path),
        "size_bytes": receipt_path.stat().st_size,
    }
    (output / "ARTIFACTS.json").write_text(
        json.dumps(
            {
                "artifacts": artifacts,
                "schema_version": "masld-bench-artifact-manifest-v1",
                "status": "pass",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--fixture-root", required=True, type=Path)
    parser.add_argument("--runtime-root", required=True, type=Path)
    parser.add_argument("--preflight-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = run(
        args.checkpoint,
        args.fixture_root,
        args.runtime_root,
        args.preflight_root,
        args.output,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "embedding_width", "device_name")}))


if __name__ == "__main__":
    main()
