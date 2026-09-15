#!/usr/bin/env python3
"""Build a label-blind, no-neighbor scPRINT-2 RNA fixture."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import torch
from torch.utils.data import DataLoader

from scdataloader import Collator
from scdataloader.data import SimpleAnnDataset


CHECKPOINT_SHA256 = "2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
CHECKPOINT_SIZE = 889_706_878
INPUT_SHA256 = "d6fa563b11570f37d25d9a62dd4626e8f63f0624adea873bb102a550fedc8854"
ROWS_SHA256 = "4996693e326661df563e741730709c45efd09b53cdb1d7fcbc9f223d06067620"
HUMAN = "NCBITaxon:9606"
ROWS = 1000
DONORS = 102
MAX_LEN = 3200
SEED = 20260824
FORBIDDEN = ("label", "fibrosis", "nas", "mash", "masld", "outcome", "stage")


class Scprint2FixtureError(ValueError):
    """Raised when the frozen input or no-neighbor requirement differs."""


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected = ["row_position", "row_id", "donor_id", "dataset_id", "outer_fold"]
        if reader.fieldnames != expected:
            raise Scprint2FixtureError("row contract fields differ")
        rows = list(reader)
    if len(rows) != ROWS or [int(row["row_position"]) for row in rows] != list(range(ROWS)):
        raise Scprint2FixtureError("row contract cardinality or order differs")
    return rows


def align_counts(matrix: sparse.spmatrix, source_genes: list[str], target_genes: list[str]) -> sparse.csr_matrix:
    source_index = {gene: offset for offset, gene in enumerate(source_genes)}
    target_columns = []
    source_columns = []
    for target_column, gene in enumerate(target_genes):
        source_column = source_index.get(gene)
        if source_column is not None:
            target_columns.append(target_column)
            source_columns.append(source_column)
    subset = sparse.csr_matrix(matrix)[:, source_columns].tocoo()
    mapped_columns = np.asarray(target_columns, dtype=np.int64)[subset.col]
    return sparse.csr_matrix(
        (subset.data, (subset.row, mapped_columns)),
        shape=(matrix.shape[0], len(target_genes)),
    )


def collator_probe(
    adata: ad.AnnData,
    genes: dict[str, list[str]],
    organisms: list[str],
    how: str,
    max_len: int,
) -> dict[str, Any]:
    ordered_genes = [gene for organism in organisms for gene in genes[organism]]
    organism_column = [organism for organism in organisms for _ in genes[organism]]
    gene_frame = pd.DataFrame({"organism": organism_column}, index=ordered_genes)
    dataset = SimpleAnnDataset(
        adata[:8],
        obs_to_output=["organism_ontology_term_id"],
        get_knn_cells=False,
    )
    collator = Collator(
        organisms=organisms,
        valid_genes=ordered_genes,
        how=how,
        max_len=max_len,
        add_zero_genes=0,
        genedf=gene_frame,
    )
    # scdataloader 2.1.0 initializes this attribute only when it downloads its
    # own gene table.  We provide the checkpoint-bound table to stay offline,
    # so restore the same organism set explicitly before collation.
    collator.organism_ids = set(organisms)
    np.random.seed(SEED)
    batch = next(iter(DataLoader(dataset, batch_size=8, shuffle=False, num_workers=0, collate_fn=collator)))
    if "knn_cells" in batch or "knn_cells_info" in batch:
        raise Scprint2FixtureError("no-neighbor collator emitted neighbor tensors")
    if batch["x"].shape != (8, max_len) or batch["genes"].shape != (8, max_len):
        raise Scprint2FixtureError("collated fixture shape differs")
    if int(batch["genes"].min()) < 0 or int(batch["genes"].max()) >= len(ordered_genes):
        raise Scprint2FixtureError("collated gene indices leave the checkpoint vocabulary")
    collator_depth = batch["depth"].numpy().astype(np.float64)
    source_depth = adata.obs["source_total_count"].to_numpy(dtype=np.float64)[:8]
    if np.any(source_depth < collator_depth):
        raise Scprint2FixtureError("source depth is below vocabulary-aligned depth")
    return {
        "batch_shape": list(batch["x"].shape),
        "max_len": max_len,
        "gene_index_min": int(batch["genes"].min()),
        "gene_index_max": int(batch["genes"].max()),
        "neighbor_tensor_keys": [],
        "expression_min": float(batch["x"].min()),
        "expression_max": float(batch["x"].max()),
        "collator_aligned_depth_sum": float(collator_depth.sum()),
        "source_prefilter_depth_sum": float(source_depth.sum()),
        "source_depth_override_required": bool(np.any(source_depth != collator_depth)),
    }


def build(source: Path, rows_path: Path, checkpoint: Path, output: Path) -> dict[str, Any]:
    source = source.resolve(strict=True)
    rows_path = rows_path.resolve(strict=True)
    checkpoint = checkpoint.resolve(strict=True)
    if output.exists():
        raise Scprint2FixtureError("output already exists")
    if digest(source) != INPUT_SHA256 or digest(rows_path) != ROWS_SHA256:
        raise Scprint2FixtureError("registered input artifact differs")
    if checkpoint.stat().st_size != CHECKPOINT_SIZE or digest(checkpoint) != CHECKPOINT_SHA256:
        raise Scprint2FixtureError("checkpoint identity differs")

    rows = read_rows(rows_path)
    adata = ad.read_h5ad(source)
    expected_obs = {"row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
    if set(adata.obs.columns) != expected_obs:
        raise Scprint2FixtureError("unlabeled input observation schema differs")
    if set(adata.obs["assay"].astype(str)) != {"unknown"}:
        raise Scprint2FixtureError("unlabeled input assay placeholder differs")
    if any(any(token in column.lower() for token in FORBIDDEN) for column in adata.obs.columns):
        raise Scprint2FixtureError("input fixture contains an evaluation field")
    if adata.n_obs != ROWS or adata.obs["row_id"].astype(str).tolist() != [row["row_id"] for row in rows]:
        raise Scprint2FixtureError("input row identity differs")
    if adata.obs["donor_id"].astype(str).nunique() != DONORS:
        raise Scprint2FixtureError("input donor cardinality differs")
    if not sparse.issparse(adata.X):
        raise Scprint2FixtureError("input counts are not sparse")
    matrix = sparse.csr_matrix(adata.X)
    if matrix.data.size and (
        float(matrix.data.min()) < 0
        or not np.allclose(matrix.data, np.rint(matrix.data), rtol=0.0, atol=1e-6)
    ):
        raise Scprint2FixtureError("input is not raw nonnegative integer counts")
    source_genes = adata.var_names.astype(str).tolist()
    if len(set(source_genes)) != len(source_genes) or any(not gene.startswith("ENSG") for gene in source_genes):
        raise Scprint2FixtureError("input feature identifiers differ")

    payload = torch.load(checkpoint, map_location="cpu", weights_only=True, mmap=True)
    hyperparameters = payload.get("hyper_parameters")
    if not isinstance(hyperparameters, dict):
        raise Scprint2FixtureError("checkpoint hyperparameters differ")
    genes = hyperparameters.get("genes")
    organisms = hyperparameters.get("organisms")
    if not isinstance(genes, dict) or not isinstance(organisms, list) or HUMAN not in genes:
        raise Scprint2FixtureError("checkpoint vocabulary metadata differs")
    if set(genes) != set(organisms) or len(set(gene for values in genes.values() for gene in values)) != sum(map(len, genes.values())):
        raise Scprint2FixtureError("checkpoint organism vocabulary is not one-to-one")

    source_total_counts = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)
    human_genes = genes[HUMAN]
    aligned = align_counts(matrix, source_genes, human_genes)
    aligned_total_counts = np.asarray(aligned.sum(axis=1)).ravel().astype(np.float64)
    if np.any(source_total_counts < aligned_total_counts):
        raise Scprint2FixtureError("source count totals are below aligned totals")
    output_obs = adata.obs.copy()
    output_obs["organism_ontology_term_id"] = HUMAN
    output_obs["source_total_count"] = source_total_counts
    aligned_adata = ad.AnnData(
        X=aligned,
        obs=output_obs,
        var=pd.DataFrame(index=pd.Index(human_genes, name="ensembl_gene_id")),
    )
    if aligned_adata.n_obs != ROWS or aligned_adata.n_vars != 20_004:
        raise Scprint2FixtureError("aligned matrix dimensions differ")

    probes = {
        "common_checkpoint_max_most_expr": collator_probe(
            aligned_adata, genes, organisms, "most expr", MAX_LEN
        ),
        "random_checkpoint_max": collator_probe(
            aligned_adata, genes, organisms, "random expr", MAX_LEN
        ),
        "native_default_random_expr": collator_probe(
            aligned_adata, genes, organisms, "random expr", 2000
        ),
    }
    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    stage.mkdir(parents=True)
    aligned_adata.uns["scprint2_input_firewall"] = {
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
        "neighbor_tensors_allowed": False,
    }
    aligned_adata.write_h5ad(stage / "aligned_counts.h5ad", compression="gzip")
    receipt = {
        "schema_version": "masld-bench-scprint2-no-neighbor-fixture-v1",
        "status": "pass",
        "dataset_view_id": "resource_atlas_geneformer_smoke_1000_v1",
        "rows": ROWS,
        "donors": DONORS,
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "checkpoint_total_gene_vocabulary": sum(map(len, genes.values())),
        "checkpoint_human_gene_vocabulary": len(human_genes),
        "source_gene_count": len(source_genes),
        "source_human_checkpoint_overlap": len(set(source_genes) & set(human_genes)),
        "aligned_nonzero_counts": int(aligned.nnz),
        "source_total_count_sum": float(source_total_counts.sum()),
        "aligned_total_count_sum": float(aligned_total_counts.sum()),
        "aligned_count_fraction": float(
            aligned_total_counts.sum() / source_total_counts.sum()
        ),
        "cells_with_out_of_vocabulary_counts": int(
            np.sum(source_total_counts != aligned_total_counts)
        ),
        "source_depth_preserved_for_inference": True,
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
        "neighbor_tensors_allowed": False,
        "offline_collator_organism_ids_restored": True,
        "collator_probes": probes,
    }
    (stage / "fixture_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    artifacts = {
        name: {"sha256": digest(stage / name), "size_bytes": (stage / name).stat().st_size}
        for name in ("aligned_counts.h5ad", "fixture_receipt.json")
    }
    (stage / "ARTIFACTS.json").write_text(
        json.dumps({"artifacts": artifacts, "schema_version": "masld-bench-artifact-manifest-v1", "status": "pass"}, indent=2, sort_keys=True) + "\n"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = build(args.source, args.rows, args.checkpoint, args.output)
    print(json.dumps({key: result[key] for key in ("status", "rows", "donors", "source_human_checkpoint_overlap")}))


if __name__ == "__main__":
    main()
