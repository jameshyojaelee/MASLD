#!/usr/bin/env python3
"""Apply the identical fold-safe PCA256 common-head requirements to Caduceus."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    NativeContractError,
    apply_projection,
    fit_projection,
    head_features,
)


def array_digest(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    result = sha256()
    result.update(str(array.dtype).encode())
    result.update(str(array.shape).encode())
    result.update(array.tobytes())
    return result.hexdigest()


def project(*, raw_root: Path, output: Path) -> dict[str, object]:
    if output.exists() or raw_root.is_symlink():
        raise NativeContractError("Caduceus projection request differs")
    with np.load(raw_root / "allele_embeddings.npz", allow_pickle=False) as source:
        fixture_ids = source["fixture_ids"]
        group_ids = source["outer_locus_sequence_group_ids"]
        folds = source["outer_folds"].astype(np.int64)
        allele_order = tuple(source["allele_order"].tolist())
        embeddings = source["embeddings"]
    if (
        fixture_ids.shape != (1_033,)
        or group_ids.shape != (1_033,)
        or folds.shape != (1_033,)
        or allele_order != ALLELES
        or embeddings.shape != (1_033, 4, 256)
        or {int(value) for value in folds} != set(range(5))
        or not np.isfinite(embeddings).all()
    ):
        raise NativeContractError("Caduceus raw embedding artifact differs")
    output.mkdir(parents=True)
    model_root = output / "caduceus"
    model_root.mkdir()
    fold_receipts = []
    for held_out_fold in range(5):
        parameters = fit_projection(embeddings, folds, held_out_fold, width=256)
        projected = apply_projection(embeddings, parameters)
        features = head_features(projected)
        if features.shape != (1_033, 1_024):
            raise NativeContractError("shared Caduceus head feature width differs")
        fold_root = model_root / f"heldout_fold{held_out_fold}"
        fold_root.mkdir()
        np.savez_compressed(fold_root / "projection_parameters.npz", **parameters)
        np.savez_compressed(
            fold_root / "head_features.npz",
            fixture_ids=fixture_ids,
            outer_locus_sequence_group_ids=group_ids,
            outer_folds=folds,
            feature_block_order=np.asarray(
                ["REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF"]
            ),
            features=features,
        )
        fold_receipts.append(
            {
                "held_out_fold": held_out_fold,
                "training_elements": int(parameters["training_elements"]),
                "held_out_elements": int(np.sum(folds == held_out_fold)),
                "projection_width": int(projected.shape[2]),
                "head_input_width": int(features.shape[1]),
                "features_sha256": array_digest(features),
                "projection_components_sha256": array_digest(parameters["components"]),
            }
        )
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-caduceus-projected-head-features-v1",
        "status": "pass_outcome_blind_fold_safe_projection",
        "dataset_id": "gse281364",
        "models": {
            "caduceus": {
                "hidden_width": 256,
                "restricted_comparator": False,
                "license": "Apache-2.0",
                "folds": fold_receipts,
            }
        },
        "elements": 1_033,
        "outer_locus_sequence_groups": 1_033,
        "outer_folds": 5,
        "projection_width": 256,
        "head_input_width": 1_024,
        "head_feature_blocks": [
            "REF",
            "ALT",
            "ALT_minus_REF",
            "absolute_ALT_minus_REF",
        ],
        "linear_head": "Linear_1024_to_1",
        "two_layer_head": "Linear_1024_to_256_GELU_Linear_256_to_1",
        "projection_fit_on_held_out_fold": False,
        "reporter_counts_read": False,
        "reporter_outcomes_read": False,
        "sealed_labels_read": False,
        "downstream_head_fit": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    print(json.dumps(project(**vars(parser.parse_args())), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
