#!/usr/bin/env python3
"""Verify exact frozen shared-control row copying without a numerical refit."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from scripts.fit_gse281364_dna_language_seeded_heads import (
    ROW_FIELDS,
    _load_frozen_shared_control,
    _read_tsv,
    _tree,
)
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    PREDICTION_FIELDS,
    file_sha256,
    load_config,
)
from scripts.verify_gse281364_shared_control_identity import _selected_rows


class SharedControlCopyError(RuntimeError):
    """Raised when authoritative control rows are transformed or rekeyed."""


def preflight(root: Path, config_path: Path) -> dict[str, object]:
    root = root.resolve(strict=True)
    config = load_config(config_path.resolve(strict=True))
    row_tree = _tree(root, config["row_authority"], label="row_authority")
    reference_tree = _tree(root, config["shared_control_reference"], label="shared_control_reference")
    row_path = row_tree / config["row_authority"]["member"]
    reference_path = reference_tree / config["shared_control_reference"]["member"]
    if (
        file_sha256(row_path) != config["row_authority"]["member_sha256"]
        or file_sha256(reference_path) != config["shared_control_reference"]["member_sha256"]
    ):
        raise SharedControlCopyError("shared-control authority member differs")
    rows = _read_tsv(row_path, ROW_FIELDS)
    reused_rows, selections = _load_frozen_shared_control(root, config, rows)
    reused = b"".join(
        ("\t".join(reused_rows[(int(row["seed"]), row["row_hash"])][field] for field in PREDICTION_FIELDS) + "\n").encode("utf-8")
        for row in rows
    )
    reference, count = _selected_rows(
        reference_path,
        model_id="available_simple_controls",
        head_id="allele_identity_ridge",
    )
    if count != 10330 or len(reused_rows) != 10330 or len(selections) != 50 or reused != reference:
        raise SharedControlCopyError("shared-control copy path is not bit-identical")
    receipt = {
        "status": "pass_authoritative_shared_control_copy_preflight",
        "rows": 10330,
        "selection_rows": 50,
        "canonical_selected_rows_sha256": sha256(reference).hexdigest(),
        "reference_fit_artifacts_sha256": config["shared_control_reference"]["artifacts_sha256"],
        "reference_prediction_member_sha256": config["shared_control_reference"]["member_sha256"],
        "copied_without_numerical_transformation": True,
        "control_refit_performed": False,
        "counted_as_new_evidence": False,
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
