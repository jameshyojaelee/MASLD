#!/usr/bin/env python3
"""Build a label-free TF-Sapiens fixture from the frozen GSE296875 RNA view."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file
from scripts.build_transcriptformer_unlabeled_fixture import (
    ASSAY_UNKNOWN_TOKEN_ID,
    VOCAB_SHA256,
    _load_vocab,
    _tokenize,
)


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
EXPECTED_ROWS = 7_500
EXPECTED_GENES = 36_601


class QueryFixtureError(ValueError):
    """Raised when the feature-only query fixture violates its separation."""


def _decode(values: object) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def build(args: argparse.Namespace) -> dict[str, object]:
    if args.output.exists():
        raise QueryFixtureError("refusing to overwrite TF-Sapiens query fixture")
    features = args.features.resolve(strict=True)
    if sha256_file(features / "ARTIFACTS.json") != args.expected_features_sha256:
        raise QueryFixtureError("feature ARTIFACTS SHA-256 differs")
    meta = verify_frozen_tree(features)["metadata"]
    if (
        meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_features"
        or meta.get("view_id") != VIEW_ID
        or meta.get("rows") != EXPECTED_ROWS
        or meta.get("labels_present") is not False
        or meta.get("donor_ids_present") is not False
        or meta.get("sealed_outcomes_accessed") is not False
    ):
        raise QueryFixtureError("feature-only query authority differs")
    vocab = args.vocab.resolve(strict=True)
    if sha256_file(vocab) != VOCAB_SHA256:
        raise QueryFixtureError("TF-Sapiens gene vocabulary differs")

    with h5py.File(features / "rna_query.h5", "r") as handle:
        if set(handle["obs"].keys()) != {
            "row_id", "rna_status", "rna_observed_mask", "atac_observed_mask"
        }:
            raise QueryFixtureError("query observation firewall differs")
        row_ids = _decode(handle["obs/row_id"][:])
        gene_ids = _decode(handle["rna/ensembl_id"][:])
        group = handle["rna/counts_csr"]
        shape = tuple(int(value) for value in group["shape"][:])
        matrix = sparse.csr_matrix(
            (
                np.asarray(group["data"][:]),
                np.asarray(group["indices"][:], dtype=np.int64),
                np.asarray(group["indptr"][:], dtype=np.int64),
            ),
            shape=shape,
        )
        rna_mask = np.asarray(handle["obs/rna_observed_mask"][:])
        atac_mask = np.asarray(handle["obs/atac_observed_mask"][:])
    if (
        matrix.shape != (EXPECTED_ROWS, EXPECTED_GENES)
        or len(set(row_ids)) != EXPECTED_ROWS
        or len(set(gene_ids)) != EXPECTED_GENES
        or np.any(rna_mask != 1)
        or np.any(atac_mask != 0)
        or np.any(matrix.data < 0)
        or np.any(matrix.data != np.floor(matrix.data))
    ):
        raise QueryFixtureError("query RNA matrix or missingness differs")

    genes, gene_to_token = _load_vocab(vocab)
    if len(genes) + 7 != 23_830:
        raise QueryFixtureError("TF-Sapiens vocabulary cardinality differs")
    native_tokens, native_counts, native_summary = _tokenize(
        matrix, gene_ids, gene_to_token, "native"
    )
    common_tokens, common_counts, common_summary = _tokenize(
        matrix, gene_ids, gene_to_token, "common"
    )
    assay = np.full((EXPECTED_ROWS, 1), ASSAY_UNKNOWN_TOKEN_ID, dtype=np.int64)
    folds = np.full(EXPECTED_ROWS, -1, dtype=np.int8)

    args.output.mkdir(mode=0o750)
    for policy, tokens, counts in (
        ("native", native_tokens, native_counts),
        ("common", common_tokens, common_counts),
    ):
        with (args.output / f"{policy}_tokens.npz").open("xb") as handle:
            np.savez_compressed(
                handle,
                assay_token_indices=assay,
                gene_counts=counts,
                gene_token_indices=tokens,
                outer_folds=folds,
                row_ids=np.asarray(row_ids),
            )
    with (args.output / "row_contract.tsv").open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("row_position", "row_id", "donor_id", "dataset_id", "outer_fold"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, row_id in enumerate(row_ids):
            writer.writerow(
                {
                    "row_position": index,
                    "row_id": row_id,
                    "donor_id": "withheld_evaluator",
                    "dataset_id": "gse296875",
                    "outer_fold": -1,
                }
            )
    receipt = {
        "schema_version": "masld-bench-transcriptformer-gse296875-query-fixture-v1",
        "status": "pass_outcome_blind_fixture",
        "view_id": VIEW_ID,
        "rows": EXPECTED_ROWS,
        "genes": EXPECTED_GENES,
        "row_ids_sha256": canonical_sha256(row_ids),
        "feature_artifacts_sha256": args.expected_features_sha256,
        "vocab_sha256": VOCAB_SHA256,
        "policies": {"native": native_summary, "common": common_summary},
        "evaluation_label_columns_read": [],
        "donor_id_state": "withheld_evaluator",
        "outer_fold_state": "not_applicable_query",
        "atac_state": "structurally_missing",
        "atac_values_read": False,
        "target_adaptation_performed": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(args.output / "fixture_receipt.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_transcriptformer_outcome_blind_query_fixture",
            "view_id": VIEW_ID,
            "rows": EXPECTED_ROWS,
            "labels_read": False,
            "donor_ids_read": False,
            "atac_read": False,
            "target_adaptation_performed": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expected-features-sha256", required=True)
    parser.add_argument("--vocab", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(build(parse_args()), sort_keys=True))
