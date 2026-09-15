#!/usr/bin/env python3
"""Refit only the frozen shared control and require exact prior OOF rows."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path

import numpy as np

from scripts.evaluate_gse281364_dna_lm_common_lane import load_outcomes
from scripts.fit_gse281364_dna_language_seeded_heads import (
    ROW_FIELDS,
    _load_alleles,
    _load_frozen_shared_control,
    _read_tsv,
    _tree,
    fit_seeded_dual_projection_fold,
)
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    CONTEXTS,
    PREDICTION_FIELDS,
    SEEDS,
    file_sha256,
    load_config,
)
from scripts.verify_gse281364_shared_control_identity import _selected_rows


class SharedControlRefitError(RuntimeError):
    """Raised when the actual shared-control refit is not reproducible."""


def preflight(root: Path, config_path: Path) -> dict[str, object]:
    root = root.resolve(strict=True)
    config = load_config(config_path.resolve(strict=True))
    row_tree = _tree(root, config["row_authority"], label="row_authority")
    outcome_tree = _tree(root, config["outcome_authority"], label="outcome_authority")
    allele_tree = _tree(root, config["allele_authority"], label="allele_authority")
    reference_tree = _tree(root, config["shared_control_reference"], label="shared_control_reference")
    row_path = row_tree / config["row_authority"]["member"]
    outcome_path = outcome_tree / config["outcome_authority"]["member"]
    allele_path = allele_tree / config["allele_authority"]["member"]
    reference_path = reference_tree / config["shared_control_reference"]["member"]
    if (
        file_sha256(row_path) != config["row_authority"]["member_sha256"]
        or file_sha256(outcome_path) != config["outcome_authority"]["member_sha256"]
        or file_sha256(allele_path) != config["allele_authority"]["member_sha256"]
        or file_sha256(reference_path) != config["shared_control_reference"]["member_sha256"]
    ):
        raise SharedControlRefitError("shared-control authority member differs")
    rows = _read_tsv(row_path, ROW_FIELDS)
    if len(rows) != 10330:
        raise SharedControlRefitError("row denominator differs")
    metadata: dict[str, tuple[str, str]] = {}
    for row in rows:
        value = (row["long_range_block_id"], row["outer_fold"])
        if metadata.setdefault(row["element_id"], value) != value:
            raise SharedControlRefitError("element metadata differs")
    elements = list(metadata)
    element_index = {element: index for index, element in enumerate(elements)}
    folds = np.asarray(
        [int(metadata[element][1].removeprefix("fold-")) for element in elements],
        dtype=np.int64,
    )
    blocks = np.asarray([metadata[element][0] for element in elements])
    outcomes = load_outcomes(outcome_path, set(elements))
    targets = {
        context: np.asarray([outcomes[(element, context)]["mean"] for element in elements])
        for context in CONTEXTS
    }
    features = _load_alleles(allele_path, elements)
    predictions: dict[tuple[int, str], np.ndarray] = {}
    alphas = tuple(float(value) for value in config["fit"]["alpha_grid"])
    for context in CONTEXTS:
        for seed in SEEDS:
            oof = np.empty(1033)
            for held_fold in range(5):
                test = folds == held_fold
                prediction, _, _ = fit_seeded_dual_projection_fold(
                    inner_features=features,
                    outer_features=features,
                    outcomes=targets[context],
                    folds=folds,
                    block_ids=blocks,
                    held_fold=held_fold,
                    seed=seed,
                    top_k_grid=config["fit"]["allele_top_k_grid"],
                    alphas=alphas,
                )
                oof[test] = prediction
            predictions[(seed, context)] = oof
    selected = []
    for row in rows:
        value = float(predictions[(int(row["seed"]), row["assay_context_id"])][element_index[row["element_id"]]])
        if not math.isfinite(value):
            raise SharedControlRefitError("non-finite shared-control prediction")
        record = {
            **row,
            "model_id": "available_simple_controls",
            "head_id": "allele_identity_ridge",
            "prediction": format(value, ".17g"),
            "experimental_replicates": "4",
            "biological_donors": "0",
            "outcome_role": "exposed_development_MPRA_only",
        }
        selected.append(("\t".join(record[field] for field in PREDICTION_FIELDS) + "\n").encode("utf-8"))
    candidate = b"".join(selected)
    reference, count = _selected_rows(
        reference_path,
        model_id="available_simple_controls",
        head_id="allele_identity_ridge",
    )
    if count != 10330 or candidate != reference:
        raise SharedControlRefitError("actual shared-control OOF rows are not bit-identical")
    reused_rows, _ = _load_frozen_shared_control(root, config, rows)
    reused = b"".join(
        ("\t".join(reused_rows[(int(row["seed"]), row["row_hash"])][field] for field in PREDICTION_FIELDS) + "\n").encode("utf-8")
        for row in rows
    )
    if reused != reference:
        raise SharedControlRefitError("shared-control copy path transforms frozen rows")
    receipt = {
        "status": "pass_actual_shared_control_refit_bit_identical",
        "rows": 10330,
        "canonical_selected_rows_sha256": sha256(reference).hexdigest(),
        "reference_fit_artifacts_sha256": config["shared_control_reference"]["artifacts_sha256"],
        "bootstrap_draw_contract": "exact_frozen_enformer_sei_offsets",
        "element_order": "frozen_row_universe_insertion_order",
        "model_features_accessed": False,
        "shared_control_copy_path_bit_identical": True,
        "shared_control_rows_copied_without_numerical_transformation": True,
        "sealed_assets_read": False,
    }
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    arguments = parser.parse_args()
    preflight(arguments.root, arguments.config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
