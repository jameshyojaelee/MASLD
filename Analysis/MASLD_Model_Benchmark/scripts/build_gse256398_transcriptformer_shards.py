#!/usr/bin/env python3
"""Build sharded, outcome-blind TranscriptFormer inputs for GSE256398."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from scripts.build_gse256398_unlabeled_smoke import load_allowed_genes, read_10x_counts
from scripts.build_transcriptformer_unlabeled_fixture import (
    ASSAY_UNKNOWN_TOKEN_ID,
    SEQUENCE_LENGTH,
    SPECIAL_TOKENS,
    VOCAB_SHA256,
    _load_vocab,
    _tokenize,
)


EXPECTED_SOURCE_ARTIFACTS = (
    "f54ac799b22cb02912b82432929255e6400d871cc0fd99dbf9b4756f8f977a20"
)
EXPECTED_CROSSWALK_ARTIFACTS = (
    "54ca7af31f6ae223defdee4a7c1d7931f7cba4362f621fe9646fd08530c4c6da"
)
EXPECTED_ROWS = 172_997
EXPECTED_DONORS = 26
OUTER_FOLDS = 5
FOLD_NAMESPACE = "transcriptformer-gse256398-unlabeled-v1"
H5_PATTERN = re.compile(
    r"^(GSM809\d{4})_(S\d+)_CB_raw_feature_bc_matrix_filtered\.h5$"
)
MEMBERSHIP_FIELDS = (
    "row_id",
    "gsm",
    "source_sample_id",
    "barcode",
    "in_source_exact",
    "in_prospective_capped",
)


class GSE256398TranscriptFormerShardError(RuntimeError):
    """Raised when a source or sharded fixture requirement differs."""


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def donor_sort_key(value: str) -> tuple[int, str]:
    if not re.fullmatch(r"S\d+", value):
        raise GSE256398TranscriptFormerShardError(f"unexpected donor ID: {value}")
    return int(value[1:]), value


def outer_fold(donor: str) -> int:
    payload = f"{FOLD_NAMESPACE}\0{donor}".encode("utf-8")
    return int.from_bytes(sha256(payload).digest()[:8], "big") % OUTER_FOLDS


def load_union_membership(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MEMBERSHIP_FIELDS:
            raise GSE256398TranscriptFormerShardError("union membership schema differs")
        rows = list(reader)
    if len(rows) != EXPECTED_ROWS or len({row["row_id"] for row in rows}) != EXPECTED_ROWS:
        raise GSE256398TranscriptFormerShardError("union membership cardinality differs")
    if any(row["in_prospective_capped"] != "true" for row in rows):
        raise GSE256398TranscriptFormerShardError("union includes non-prospective rows")
    expected_order = sorted(
        rows,
        key=lambda row: (donor_sort_key(row["source_sample_id"]), row["barcode"]),
    )
    if rows != expected_order:
        raise GSE256398TranscriptFormerShardError("union membership is not canonical")
    return rows


def source_h5_by_donor(source: Path) -> dict[str, tuple[str, Path]]:
    result: dict[str, tuple[str, Path]] = {}
    for path in sorted((source / "h5").glob("*.h5")):
        match = H5_PATTERN.fullmatch(path.name)
        if match is None:
            raise GSE256398TranscriptFormerShardError("source H5 filename differs")
        gsm, donor = match.groups()
        if donor in result:
            raise GSE256398TranscriptFormerShardError("duplicate source donor H5")
        result[donor] = (gsm, path)
    return result


def write_row_contract(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "row_position",
        "row_id",
        "source_sample_id",
        "gsm",
        "barcode",
        "outer_fold",
        "in_source_exact",
        "in_prospective_capped",
        "shard_name",
        "shard_row_position",
    )
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build(
    *,
    source: Path,
    crosswalk: Path,
    membership_lock: Path,
    membership_lock_sha256: str,
    vocab: Path,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise GSE256398TranscriptFormerShardError("fixture output already exists")
    expected = (
        (source, EXPECTED_SOURCE_ARTIFACTS),
        (crosswalk, EXPECTED_CROSSWALK_ARTIFACTS),
        (membership_lock, membership_lock_sha256),
    )
    for root, digest in expected:
        verify_frozen_tree(root)
        if sha256_file(root / "ARTIFACTS.json") != digest:
            raise GSE256398TranscriptFormerShardError(f"frozen input differs: {root}")
    if sha256_file(vocab) != VOCAB_SHA256:
        raise GSE256398TranscriptFormerShardError("TranscriptFormer vocabulary differs")

    rows = load_union_membership(
        membership_lock / "prospective_union_membership.tsv"
    )
    by_donor: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_donor.setdefault(row["source_sample_id"], []).append(row)
    h5_by_donor = source_h5_by_donor(source)
    if (
        len(by_donor) != EXPECTED_DONORS
        or set(by_donor) != set(h5_by_donor)
        or len({outer_fold(donor) for donor in by_donor}) != OUTER_FOLDS
    ):
        raise GSE256398TranscriptFormerShardError("donor roster or folds differ")

    allowed_indices, var = load_allowed_genes(crosswalk / "gene_crosswalk.tsv")
    ensembl_ids = var.index.astype(str).tolist()
    genes, gene_to_token = _load_vocab(vocab)
    if len(genes) + len(SPECIAL_TOKENS) != 23_830:
        raise GSE256398TranscriptFormerShardError("checkpoint vocabulary cardinality differs")

    output.mkdir(parents=True, exist_ok=False)
    shard_root = output / "shards"
    shard_root.mkdir()
    shard_records: list[dict[str, Any]] = []
    row_contract: list[dict[str, Any]] = []
    global_position = 0
    for donor in sorted(by_donor, key=donor_sort_key):
        donor_rows = by_donor[donor]
        expected_gsm, h5_path = h5_by_donor[donor]
        counts, barcodes, _ = read_10x_counts(h5_path)
        barcode_to_index = {barcode: index for index, barcode in enumerate(barcodes)}
        if len(barcode_to_index) != len(barcodes):
            raise GSE256398TranscriptFormerShardError("source barcode axis is duplicated")
        if any(row["gsm"] != expected_gsm for row in donor_rows):
            raise GSE256398TranscriptFormerShardError("membership GSM differs")
        try:
            selected_indices = [barcode_to_index[row["barcode"]] for row in donor_rows]
        except KeyError as error:
            raise GSE256398TranscriptFormerShardError(
                f"membership barcode absent for {donor}"
            ) from error
        matrix = counts[selected_indices, :][:, allowed_indices].tocsr()
        if (
            matrix.shape != (len(donor_rows), 35_455)
            or matrix.nnz < len(donor_rows)
            or float(matrix.data.min()) <= 0.0
            or not np.allclose(matrix.data, np.rint(matrix.data), rtol=0.0, atol=1e-6)
        ):
            raise GSE256398TranscriptFormerShardError("raw-count shard differs")

        native_tokens, native_counts, native_summary = _tokenize(
            matrix, ensembl_ids, gene_to_token, "native"
        )
        common_tokens, common_counts, common_summary = _tokenize(
            matrix, ensembl_ids, gene_to_token, "common"
        )
        row_ids = np.asarray([row["row_id"] for row in donor_rows])
        assay = np.full((len(donor_rows), 1), ASSAY_UNKNOWN_TOKEN_ID, dtype=np.int64)
        fold = outer_fold(donor)
        folds = np.full(len(donor_rows), fold, dtype=np.int8)
        source_mask = np.asarray(
            [row["in_source_exact"] == "true" for row in donor_rows], dtype=bool
        )
        shard_name = f"{donor.lower()}_tokens.npz"
        shard_path = shard_root / shard_name
        np.savez_compressed(
            shard_path,
            assay_token_indices=assay,
            common_gene_counts=common_counts,
            common_gene_token_indices=common_tokens,
            in_source_exact=source_mask,
            native_gene_counts=native_counts,
            native_gene_token_indices=native_tokens,
            outer_folds=folds,
            row_ids=row_ids,
        )
        shard_records.append(
            {
                "source_sample_id": donor,
                "gsm": expected_gsm,
                "rows": len(donor_rows),
                "outer_fold": fold,
                "source_exact_rows": int(source_mask.sum()),
                "path": f"shards/{shard_name}",
                "sha256": sha256_file(shard_path),
                "size_bytes": shard_path.stat().st_size,
                "native_summary": native_summary,
                "common_summary": common_summary,
            }
        )
        for shard_position, row in enumerate(donor_rows):
            row_contract.append(
                {
                    "row_position": global_position,
                    "row_id": row["row_id"],
                    "source_sample_id": donor,
                    "gsm": row["gsm"],
                    "barcode": row["barcode"],
                    "outer_fold": fold,
                    "in_source_exact": row["in_source_exact"],
                    "in_prospective_capped": "true",
                    "shard_name": shard_name,
                    "shard_row_position": shard_position,
                }
            )
            global_position += 1
        print(
            json.dumps(
                {"donor": donor, "rows": len(donor_rows), "completed_rows": global_position},
                sort_keys=True,
            ),
            flush=True,
        )

    if global_position != EXPECTED_ROWS or len(row_contract) != EXPECTED_ROWS:
        raise GSE256398TranscriptFormerShardError("sharded row census differs")
    row_contract_path = output / "row_contract.tsv"
    write_row_contract(row_contract_path, row_contract)
    receipt = {
        "schema_version": "masld-bench-gse256398-transcriptformer-shards-v1",
        "status": "pass_outcome_blind_sharded_fixture",
        "dataset_id": "gse256398",
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "features_allowed_gencode49": 35_455,
        "gene_order_policies": ["native", "common"],
        "sequence_length": SEQUENCE_LENGTH,
        "raw_integer_counts": True,
        "normalization_fitted": False,
        "hvg_selection_fitted": False,
        "phenotype_metadata_read": False,
        "barcode_cell_labels_read": False,
        "sealed_outcomes_read": False,
        "supervised_accuracy_claim_allowed": False,
        "membership_lock_artifacts_sha256": membership_lock_sha256,
        "row_contract_sha256": sha256_file(row_contract_path),
        "shards": shard_records,
    }
    write_json_exclusive(output / "fixture_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse256398_outcome_blind_transcriptformer_sharded_fixture",
            "dataset_id": "gse256398",
            "donors": EXPECTED_DONORS,
            "rows": EXPECTED_ROWS,
            "phenotype_metadata_read": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--membership-lock", type=Path, required=True)
    parser.add_argument("--membership-lock-sha256", required=True)
    parser.add_argument("--vocab", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build(
        source=args.source,
        crosswalk=args.crosswalk,
        membership_lock=args.membership_lock,
        membership_lock_sha256=args.membership_lock_sha256,
        vocab=args.vocab,
        output=args.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
