#!/usr/bin/env python3
"""Build label-blind TranscriptFormer token fixtures from a registered RNA view."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


SOURCE_SHA256 = "d59c85b4fb0a1352a427de3b2a587c30aac12891a8d753e10d4442fc719c985f"
VOCAB_SHA256 = "676640be87ba4ba7ea7af8d9dc470a8ae8be02768989478027dec5fe61f465e2"
SPECIAL_TOKENS = ("unknown", "[PAD]", "[START]", "[END]", "[RD]", "[CELL]", "[MASK]")
PAD_TOKEN_ID = 1
ASSAY_UNKNOWN_TOKEN_ID = 0
SEQUENCE_LENGTH = 2047
COUNT_CLIP = 30.0
OUTER_FOLDS = 5
FOLD_NAMESPACE = "transcriptformer-resource-atlas-smoke-v1"
ALLOWED_SOURCE_OBS = ("_index", "donor_id", "dataset")
FORBIDDEN_EVALUATION_COLUMNS = (
    "broad_label",
    "source_cell_type",
    "cell_type",
    "label",
    "fibrosis",
    "nas",
    "mash",
    "masld",
)


class TranscriptFormerFixtureError(ValueError):
    """Raised when an outcome-blind fixture requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _strings(values: Any) -> list[str]:
    result = []
    for value in values:
        if isinstance(value, bytes):
            result.append(value.decode("utf-8"))
        else:
            result.append(str(value))
    return result


def _read_h5_column(group: h5py.Group, key: str) -> list[str]:
    if key not in group:
        raise TranscriptFormerFixtureError(f"source lacks allowed observation field {key}")
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return _strings(node[:])
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = _strings(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise TranscriptFormerFixtureError(f"source categorical field {key} has missing/invalid codes")
        return [categories[int(code)] for code in codes]
    raise TranscriptFormerFixtureError(f"unsupported H5AD encoding for allowed field {key}")


def _read_source(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str], list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        if not {"X", "obs", "var"}.issubset(handle):
            raise TranscriptFormerFixtureError("source is not the registered H5AD shape")
        obs = handle["obs"]
        index_key = obs.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode("utf-8")
        row_ids = _read_h5_column(obs, str(index_key))
        donors = _read_h5_column(obs, "donor_id")
        datasets = _read_h5_column(obs, "dataset")
        ensembl_ids = _read_h5_column(handle["var"], "ensembl_id")

        x = handle["X"]
        if not isinstance(x, h5py.Group) or not {"data", "indices", "indptr"}.issubset(x):
            raise TranscriptFormerFixtureError("registered source X is not sparse CSR")
        shape = tuple(int(value) for value in x.attrs["shape"])
        matrix = sparse.csr_matrix(
            (
                np.asarray(x["data"][:]),
                np.asarray(x["indices"][:], dtype=np.int64),
                np.asarray(x["indptr"][:], dtype=np.int64),
            ),
            shape=shape,
        )
    return matrix, row_ids, donors, datasets, ensembl_ids


def _outer_fold(donor_id: str, namespace: str = FOLD_NAMESPACE) -> int:
    if not namespace:
        raise TranscriptFormerFixtureError("outer-fold namespace is empty")
    payload = f"{namespace}\0{donor_id}".encode("utf-8")
    return int.from_bytes(sha256(payload).digest()[:8], "big") % OUTER_FOLDS


def _load_vocab(path: Path) -> tuple[list[str], dict[str, int]]:
    with h5py.File(path, "r") as handle:
        if "keys" not in handle or "arrays" not in handle:
            raise TranscriptFormerFixtureError("TF-Sapiens gene embedding vocabulary differs")
        genes = _strings(handle["keys"][:])
        if len(genes) != len(handle["arrays"]):
            raise TranscriptFormerFixtureError("TF-Sapiens gene names and embeddings differ")
    if len(set(genes)) != len(genes) or any(not gene.startswith("ENSG") for gene in genes):
        raise TranscriptFormerFixtureError("TF-Sapiens gene vocabulary is not unique stable Ensembl IDs")
    gene_to_token = {
        gene: offset + len(SPECIAL_TOKENS) for offset, gene in enumerate(genes)
    }
    return genes, gene_to_token


def _tokenize(
    matrix: sparse.csr_matrix,
    ensembl_ids: list[str],
    gene_to_token: dict[str, int],
    policy: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if policy not in {"native", "common"}:
        raise TranscriptFormerFixtureError("unknown gene-order policy")
    tokens = np.full((matrix.shape[0], SEQUENCE_LENGTH), PAD_TOKEN_ID, dtype=np.int32)
    counts = np.zeros((matrix.shape[0], SEQUENCE_LENGTH), dtype=np.float32)
    overlap = np.asarray([gene in gene_to_token for gene in ensembl_ids], dtype=bool)
    token_by_feature = np.asarray(
        [gene_to_token.get(gene, -1) for gene in ensembl_ids], dtype=np.int32
    )
    expressed_in_vocab: list[int] = []
    truncated = 0
    for row_offset in range(matrix.shape[0]):
        row = matrix.getrow(row_offset)
        if row.nnz and (
            float(row.data.min()) < 0.0
            or not np.allclose(row.data, np.rint(row.data), rtol=0.0, atol=1e-6)
        ):
            raise TranscriptFormerFixtureError("source X is not raw nonnegative integer counts")
        retained = [
            (int(index), float(value))
            for index, value in zip(row.indices, row.data, strict=True)
            if value > 0 and overlap[int(index)]
        ]
        if policy == "native":
            retained.sort(key=lambda item: item[0])
        else:
            retained.sort(key=lambda item: (-item[1], ensembl_ids[item[0]]))
        expressed_in_vocab.append(len(retained))
        if len(retained) > SEQUENCE_LENGTH:
            truncated += 1
            retained = retained[:SEQUENCE_LENGTH]
        length = len(retained)
        if length:
            feature_indices = np.asarray([item[0] for item in retained], dtype=np.int64)
            raw_counts = np.asarray([item[1] for item in retained], dtype=np.float32)
            tokens[row_offset, :length] = token_by_feature[feature_indices]
            counts[row_offset, :length] = np.clip(raw_counts, 0.0, COUNT_CLIP)
    return tokens, counts, {
        "cells_truncated": truncated,
        "expressed_vocab_genes_max": max(expressed_in_vocab),
        "expressed_vocab_genes_median": float(np.median(expressed_in_vocab)),
        "source_features_in_vocab": int(overlap.sum()),
        "source_feature_total": len(ensembl_ids),
    }


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(
    source: Path,
    vocab: Path,
    output: Path,
    *,
    expected_source_sha256: str = SOURCE_SHA256,
    expected_rows: int = 1000,
    expected_donors: int = 102,
    dataset_view_id: str = "resource_atlas_geneformer_smoke_1000_v1",
    fold_namespace: str = FOLD_NAMESPACE,
) -> dict[str, Any]:
    if output.exists():
        raise TranscriptFormerFixtureError("output already exists")
    if expected_rows < 1 or expected_donors < 1 or not dataset_view_id:
        raise TranscriptFormerFixtureError("registered view cardinality differs")
    if sha256_file(source) != expected_source_sha256 or sha256_file(vocab) != VOCAB_SHA256:
        raise TranscriptFormerFixtureError("source or vocabulary hash differs")
    matrix, row_ids, donors, datasets, ensembl_ids = _read_source(source)
    if matrix.shape != (expected_rows, len(ensembl_ids)) or len(set(row_ids)) != expected_rows:
        raise TranscriptFormerFixtureError("registered view row/feature contract differs")
    if (
        not (len(row_ids) == len(donors) == len(datasets))
        or len(set(donors)) != expected_donors
    ):
        raise TranscriptFormerFixtureError("registered view donor join differs")
    if len(set(ensembl_ids)) != len(ensembl_ids):
        raise TranscriptFormerFixtureError("source Ensembl feature IDs are not unique")

    genes, gene_to_token = _load_vocab(vocab)
    if len(genes) + len(SPECIAL_TOKENS) != 23830:
        raise TranscriptFormerFixtureError("TF-Sapiens checkpoint vocabulary row count differs")
    folds = np.asarray(
        [_outer_fold(donor, namespace=fold_namespace) for donor in donors],
        dtype=np.int8,
    )
    fold_by_donor = {}
    for donor, fold in zip(donors, folds, strict=True):
        previous = fold_by_donor.setdefault(donor, int(fold))
        if previous != int(fold):
            raise TranscriptFormerFixtureError("one donor crosses outer folds")
    if set(fold_by_donor.values()) != set(range(OUTER_FOLDS)):
        raise TranscriptFormerFixtureError("outcome-blind donor folds do not cover five folds")

    native_tokens, native_counts, native_summary = _tokenize(
        matrix, ensembl_ids, gene_to_token, "native"
    )
    common_tokens, common_counts, common_summary = _tokenize(
        matrix, ensembl_ids, gene_to_token, "common"
    )
    assay_tokens = np.full((len(row_ids), 1), ASSAY_UNKNOWN_TOKEN_ID, dtype=np.int64)

    stage = output.with_name(f".{output.name}.{os.getpid()}.staging")
    if stage.exists():
        raise TranscriptFormerFixtureError("staging output already exists")
    stage.mkdir(parents=True)
    obs = pd.DataFrame(
        {
            "row_id": row_ids,
            "donor_id": donors,
            "dataset_id": datasets,
            "outer_fold": folds,
            "assay": ["unknown"] * len(row_ids),
        },
        index=pd.Index(row_ids, name="row_id_index"),
    )
    var = pd.DataFrame(index=pd.Index(ensembl_ids, name="ensembl_id"))
    unlabeled = ad.AnnData(X=matrix, obs=obs, var=var)
    unlabeled.uns["transcriptformer_input_firewall"] = {
        "allowed_source_obs_fields_read": list(ALLOWED_SOURCE_OBS),
        "evaluation_label_columns_read": [],
        "sealed_outcomes_read": False,
    }
    unlabeled.write_h5ad(stage / "unlabeled_counts.h5ad", compression="gzip")

    common_policy = (
        "per_cell_nonzero_vocabulary_genes_sorted_by_descending_raw_count_"
        "then_ascending_stable_ensembl_id;truncate_2047;clip_counts_after_ordering_at_30"
    )
    native_policy = (
        "official_inference_default_source_var_order_of_nonzero_vocabulary_genes;"
        "sort_genes_false;randomize_genes_false;truncate_2047;clip_counts_at_30"
    )
    for policy, token_values, count_values in (
        ("native", native_tokens, native_counts),
        ("common", common_tokens, common_counts),
    ):
        np.savez_compressed(
            stage / f"{policy}_tokens.npz",
            assay_token_indices=assay_tokens,
            gene_counts=count_values,
            gene_token_indices=token_values,
            outer_folds=folds,
            row_ids=np.asarray(row_ids),
        )
    row_records = [
        {
            "row_position": offset,
            "row_id": row_id,
            "donor_id": donor,
            "dataset_id": dataset,
            "outer_fold": int(fold),
        }
        for offset, (row_id, donor, dataset, fold) in enumerate(
            zip(row_ids, donors, datasets, folds, strict=True)
        )
    ]
    _write_tsv(
        stage / "row_contract.tsv",
        ["row_position", "row_id", "donor_id", "dataset_id", "outer_fold"],
        row_records,
    )
    _write_tsv(
        stage / "gene_vocab.tsv",
        ["token_id", "ensembl_id"],
        [
            {"token_id": token_id, "ensembl_id": gene}
            for gene, token_id in gene_to_token.items()
        ],
    )
    receipt = {
        "assay_policy": "unknown_token_only_because_registered_view_has_no_assay_field",
        "common_gene_order_policy": common_policy,
        "common_summary": common_summary,
        "count_clip": COUNT_CLIP,
        "dataset_view_id": dataset_view_id,
        "donors": len(fold_by_donor),
        "evaluation_label_columns_read": [],
        "forbidden_evaluation_column_names": list(FORBIDDEN_EVALUATION_COLUMNS),
        "native_gene_order_policy": native_policy,
        "native_summary": native_summary,
        "outer_fold_donor_counts": {
            str(fold): sum(value == fold for value in fold_by_donor.values())
            for fold in range(OUTER_FOLDS)
        },
        "outer_fold_namespace": fold_namespace,
        "rows": len(row_ids),
        "schema_version": "masld-bench-transcriptformer-unlabeled-fixture-v1",
        "sealed_outcomes_read": False,
        "sequence_length": SEQUENCE_LENGTH,
        "source_obs_fields_read": list(ALLOWED_SOURCE_OBS),
        "status": "pass_outcome_blind_fixture",
        "vocabulary_rows_in_checkpoint": len(genes) + len(SPECIAL_TOKENS),
        "vocabulary_special_token_count": len(SPECIAL_TOKENS),
    }
    (stage / "fixture_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    artifact_names = (
        "common_tokens.npz",
        "fixture_receipt.json",
        "gene_vocab.tsv",
        "native_tokens.npz",
        "row_contract.tsv",
        "unlabeled_counts.h5ad",
    )
    artifacts = {
        "artifacts": {
            name: {"sha256": sha256_file(stage / name), "size_bytes": (stage / name).stat().st_size}
            for name in artifact_names
        },
        "schema_version": "masld-bench-artifact-manifest-v1",
        "status": "pass",
    }
    (stage / "ARTIFACTS.json").write_text(
        json.dumps(artifacts, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    os.rename(stage, output)
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--vocab", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-source-sha256", default=SOURCE_SHA256)
    parser.add_argument("--expected-rows", type=int, default=1000)
    parser.add_argument("--expected-donors", type=int, default=102)
    parser.add_argument(
        "--dataset-view-id", default="resource_atlas_geneformer_smoke_1000_v1"
    )
    parser.add_argument("--fold-namespace", default=FOLD_NAMESPACE)
    arguments = parser.parse_args()
    result = build(
        arguments.source,
        arguments.vocab,
        arguments.output,
        expected_source_sha256=arguments.expected_source_sha256,
        expected_rows=arguments.expected_rows,
        expected_donors=arguments.expected_donors,
        dataset_view_id=arguments.dataset_view_id,
        fold_namespace=arguments.fold_namespace,
    )
    print(json.dumps({key: result[key] for key in ("status", "rows", "donors")}))


if __name__ == "__main__":
    main()
