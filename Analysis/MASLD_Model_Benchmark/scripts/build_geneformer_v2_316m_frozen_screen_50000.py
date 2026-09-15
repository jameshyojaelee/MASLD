#!/usr/bin/env python3
"""Build outcome-blind Geneformer tokens for the frozen 50k Atlas screen."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any


ROWS = 50_000
DONORS = 102
STUDIES = 7
MAX_TOKENS = 4096
HIDDEN_SIZE = 1152
DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
EMBEDDING_POLICY = "last_hidden_state_mean_over_non_special_nonpadding_tokens"
ALLOWED_OBS = {"row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
FORBIDDEN_NAMES = {
    "broad_label",
    "source_cell_type",
    "cell_type",
    "label",
    "fibrosis",
    "nas",
    "mash",
    "masld",
    "sex",
    "age",
    "stage",
}


class GeneformerFixtureError(ValueError):
    """Raised when an input or outcome-blind token requirement differs."""


def sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _read_split(path: Path) -> dict[str, tuple[str, str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["row_id", "donor_id", "dataset", "outer_fold"]:
            raise GeneformerFixtureError("study-held-out row split schema differs")
        result: dict[str, tuple[str, str, int]] = {}
        for record in reader:
            row_id = record["row_id"]
            if row_id in result:
                raise GeneformerFixtureError("study-held-out row identifiers collide")
            fold = int(record["outer_fold"])
            if fold not in range(5):
                raise GeneformerFixtureError("study-held-out outer fold differs")
            result[row_id] = (record["donor_id"], record["dataset"], fold)
    return result


def _load_dictionary(path: Path, *, value_type: type) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or not value:
        raise GeneformerFixtureError("Geneformer canonical dictionary differs")
    try:
        return {str(key): value_type(item) for key, item in value.items()}
    except (TypeError, ValueError) as error:
        raise GeneformerFixtureError("Geneformer dictionary value differs") from error


def build(
    source: Path,
    split_path: Path,
    median_path: Path,
    token_path: Path,
    config_path: Path,
    output: Path,
    *,
    expected_source_sha256: str,
    expected_split_sha256: str,
    expected_rows: int = ROWS,
    expected_donors: int = DONORS,
    expected_studies: int = STUDIES,
) -> dict[str, Any]:
    import anndata
    import numpy as np
    from scipy import sparse

    if output.exists() or any(value < 1 for value in (expected_rows, expected_donors, expected_studies)):
        raise GeneformerFixtureError("fixture output/cardinality contract differs")
    if sha256_file(source) != expected_source_sha256:
        raise GeneformerFixtureError("outcome-blind count source hash differs")
    if sha256_file(split_path) != expected_split_sha256:
        raise GeneformerFixtureError("study-held-out row split hash differs")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if (
        config.get("architectures") != ["BertForMaskedLM"]
        or config.get("hidden_size") != HIDDEN_SIZE
        or config.get("max_position_embeddings") != MAX_TOKENS
        or config.get("vocab_size") != 20_275
    ):
        raise GeneformerFixtureError("Geneformer V2-316M config differs")
    medians = _load_dictionary(median_path, value_type=float)
    tokens = _load_dictionary(token_path, value_type=int)
    if len(tokens) != int(config["vocab_size"]):
        raise GeneformerFixtureError("Geneformer token vocabulary cardinality differs")

    adata = anndata.read_h5ad(source)
    observed_obs = set(map(str, adata.obs.columns))
    if not observed_obs.issubset(ALLOWED_OBS) or observed_obs.intersection(FORBIDDEN_NAMES):
        raise GeneformerFixtureError("outcome-blind observation firewall differs")
    if not {"donor_id", "dataset_id"}.issubset(observed_obs):
        raise GeneformerFixtureError("outcome-blind donor/study join is absent")
    if adata.n_obs != expected_rows or len(set(map(str, adata.obs_names))) != expected_rows:
        raise GeneformerFixtureError("outcome-blind row identity differs")
    matrix = adata.X
    if not sparse.issparse(matrix):
        raise GeneformerFixtureError("outcome-blind counts are not sparse")
    matrix = sparse.csr_matrix(matrix, copy=True)
    matrix.sort_indices()
    if (
        matrix.nnz < 1
        or float(matrix.data.min()) < 0.0
        or not np.allclose(matrix.data, np.rint(matrix.data), rtol=0.0, atol=1e-6)
    ):
        raise GeneformerFixtureError("outcome-blind matrix is not raw nonnegative UMI counts")
    row_ids = list(map(str, adata.obs_names))
    donors = list(map(str, adata.obs["donor_id"]))
    datasets = list(map(str, adata.obs["dataset_id"]))
    if len(set(donors)) != expected_donors or len(set(datasets)) != expected_studies:
        raise GeneformerFixtureError("outcome-blind donor/study cardinality differs")
    gene_ids = list(map(str, adata.var_names))
    if len(set(gene_ids)) != len(gene_ids) or not all(gene.startswith("ENSG") for gene in gene_ids):
        raise GeneformerFixtureError("stable Ensembl feature identity differs")

    split = _read_split(split_path)
    if set(split) != set(row_ids):
        raise GeneformerFixtureError("count rows and study split rows differ")
    outer_folds: list[int] = []
    donor_folds: dict[str, int] = {}
    study_folds: dict[str, int] = {}
    for row_id, donor, dataset in zip(row_ids, donors, datasets, strict=True):
        split_donor, split_dataset, fold = split[row_id]
        if (donor, dataset) != (split_donor, split_dataset):
            raise GeneformerFixtureError("count and study-split metadata join differs")
        if donor_folds.setdefault(donor, fold) != fold:
            raise GeneformerFixtureError("one donor crosses outer folds")
        if study_folds.setdefault(dataset, fold) != fold:
            raise GeneformerFixtureError("one study crosses outer folds")
        outer_folds.append(fold)
    if set(outer_folds) != set(range(5)):
        raise GeneformerFixtureError("study-held-out split lacks an outer fold")

    usable_by_feature = np.asarray(
        [gene in medians and gene in tokens and medians.get(gene, 0.0) > 0.0 for gene in gene_ids],
        dtype=bool,
    )
    median_by_feature = np.asarray(
        [medians.get(gene, 1.0) for gene in gene_ids], dtype=np.float64
    )
    token_by_feature = np.asarray(
        [tokens.get(gene, -1) for gene in gene_ids], dtype=np.int32
    )
    encoded: list[np.ndarray] = []
    offsets = np.zeros(expected_rows + 1, dtype=np.int64)
    truncated = 0
    maximum_length = 0
    for row_index in range(expected_rows):
        row = matrix.getrow(row_index)
        keep = usable_by_feature[row.indices] & (row.data > 0)
        indices = row.indices[keep]
        values = np.asarray(row.data[keep], dtype=np.float64)
        total = float(values.sum())
        if total <= 0.0:
            raise GeneformerFixtureError("a row has no counts in the frozen Geneformer vocabulary")
        normalized = (values / total * 1e4) / median_by_feature[indices]
        order = np.argsort(-normalized, kind="stable")
        if len(order) > MAX_TOKENS:
            truncated += 1
            order = order[:MAX_TOKENS]
        row_tokens = token_by_feature[indices[order]]
        if np.any(row_tokens < 0):
            raise GeneformerFixtureError("usable feature lacks a Geneformer token")
        encoded.append(np.asarray(row_tokens, dtype=np.int32))
        maximum_length = max(maximum_length, len(row_tokens))
        offsets[row_index + 1] = offsets[row_index] + len(row_tokens)
    token_values = np.concatenate(encoded)

    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    if stage.exists():
        raise GeneformerFixtureError("fixture staging path already exists")
    stage.mkdir(parents=True)
    np.save(stage / "token_ids.npy", token_values, allow_pickle=False)
    np.save(stage / "token_offsets.npy", offsets, allow_pickle=False)
    (stage / "embedding_row_order.txt").write_text(
        "\n".join(row_ids) + "\n", encoding="utf-8"
    )
    with (stage / "row_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_position", "row_id", "donor_id", "dataset_id", "outer_fold"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for position, (row_id, donor, dataset, fold) in enumerate(
            zip(row_ids, donors, datasets, outer_folds, strict=True)
        ):
            writer.writerow(
                {
                    "row_position": position,
                    "row_id": row_id,
                    "donor_id": donor,
                    "dataset_id": dataset,
                    "outer_fold": fold,
                }
            )
    receipt = {
        "schema_version": "masld-bench-geneformer-v2-316m-frozen-screen-fixture-v1",
        "status": "pass_outcome_blind_fixture",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": expected_rows,
        "donors": len(set(donors)),
        "studies": len(set(datasets)),
        "outer_folds": 5,
        "source_obs_fields_read": sorted(observed_obs),
        "evaluation_label_columns_read": [],
        "histology_columns_read": [],
        "sealed_outcomes_read": False,
        "embedding_policy": EMBEDDING_POLICY,
        "max_tokens": MAX_TOKENS,
        "maximum_observed_tokens": maximum_length,
        "cells_truncated": truncated,
        "genes_total": len(gene_ids),
        "genes_usable": int(usable_by_feature.sum()),
        "genes_dropped": int(len(gene_ids) - usable_by_feature.sum()),
        "total_tokens": int(offsets[-1]),
        "row_order_sha256": sha256("\n".join(row_ids).encode("utf-8") + b"\n").hexdigest(),
        "source_sha256": expected_source_sha256,
        "split_sha256": expected_split_sha256,
    }
    (stage / "fixture_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--split", type=Path, required=True)
    parser.add_argument("--median", type=Path, required=True)
    parser.add_argument("--token-dictionary", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-source-sha256", required=True)
    parser.add_argument("--expected-split-sha256", required=True)
    arguments = parser.parse_args()
    result = build(
        arguments.source,
        arguments.split,
        arguments.median,
        arguments.token_dictionary,
        arguments.config,
        arguments.output,
        expected_source_sha256=arguments.expected_source_sha256,
        expected_split_sha256=arguments.expected_split_sha256,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "total_tokens")}))


if __name__ == "__main__":
    main()
