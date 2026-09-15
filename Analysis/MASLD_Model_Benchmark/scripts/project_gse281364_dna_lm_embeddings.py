#!/usr/bin/env python3
"""Apply fold-safe outcome-blind normalization/PCA to GSE281364 DNA-LM embeddings."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from scripts.gse281364_dna_lm_native_contract import (
    ALLELES,
    MODELS,
    NativeContractError,
    apply_projection,
    fit_projection,
    head_features,
    load_config,
)


def array_digest(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str(array.dtype).encode())
    digest.update(str(array.shape).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def project(*, config_path: Path, raw_root: Path, output: Path) -> dict[str, object]:
    if output.exists() or config_path.is_symlink() or raw_root.is_symlink():
        raise NativeContractError("projection request differs")
    config = load_config(config_path)
    output.mkdir(parents=True)
    model_receipts: dict[str, object] = {}
    expected_widths = {
        model: int(config["models"][model]["hidden_width"]) for model in MODELS
    }
    reference_ids: np.ndarray | None = None
    reference_groups: np.ndarray | None = None
    reference_folds: np.ndarray | None = None
    for model in MODELS:
        source = np.load(raw_root / model / "allele_embeddings.npz", allow_pickle=False)
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
            or embeddings.shape != (1_033, 4, expected_widths[model])
            or not np.isfinite(embeddings).all()
        ):
            raise NativeContractError(f"{model} raw embedding artifact differs")
        if reference_ids is None:
            reference_ids, reference_groups, reference_folds = (
                fixture_ids,
                group_ids,
                folds,
            )
        elif (
            not np.array_equal(fixture_ids, reference_ids)
            or not np.array_equal(group_ids, reference_groups)
            or not np.array_equal(folds, reference_folds)
        ):
            raise NativeContractError("model embedding row alignment differs")
        model_root = output / model
        model_root.mkdir()
        fold_receipts = []
        for held_out_fold in range(5):
            parameters = fit_projection(embeddings, folds, held_out_fold)
            projected = apply_projection(embeddings, parameters)
            features = head_features(projected)
            if features.shape != (1_033, 1_024):
                raise NativeContractError("shared head feature width differs")
            fold_root = model_root / f"heldout_fold{held_out_fold}"
            fold_root.mkdir()
            np.savez_compressed(
                fold_root / "projection_parameters.npz",
                **parameters,
            )
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
                    "projection_components_sha256": array_digest(
                        parameters["components"]
                    ),
                }
            )
        model_receipts[model] = {
            "hidden_width": expected_widths[model],
            "restricted_comparator": bool(
                config["models"][model]["restricted_comparator"]
            ),
            "license": config["models"][model]["license"],
            "folds": fold_receipts,
        }
    receipt: dict[str, object] = {
        "schema_version": "masld-bench-gse281364-dna-lm-projected-head-features-v1",
        "status": "pass_outcome_blind_fold_safe_projection",
        "dataset_id": "gse281364",
        "models": model_receipts,
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
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    project(config_path=arguments.config, raw_root=arguments.raw_root, output=arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
