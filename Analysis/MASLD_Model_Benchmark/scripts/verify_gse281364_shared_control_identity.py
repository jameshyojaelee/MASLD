#!/usr/bin/env python3
"""Require exact shared allele-control OOF rows before family evaluation."""

from __future__ import annotations

import argparse
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree
from scripts.reconcile_gse281364_dna_language_seeded_features import (
    PREDICTION_FIELDS,
    file_sha256,
    load_config,
)


SCHEMA = "masld-bench-gse281364-shared-control-identity-v1"


class SharedControlIdentityError(RuntimeError):
    """Raised when the shared OOF control is not byte-identical by row."""


def _tree(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    relative = Path(str(binding.get("tree_path", "")))
    if not relative.parts or relative.is_absolute() or ".." in relative.parts:
        raise SharedControlIdentityError(f"unsafe {label} tree path")
    path = (root / relative).resolve(strict=True)
    path.relative_to(root)
    if file_sha256(path / "ARTIFACTS.json") != binding.get("artifacts_sha256"):
        raise SharedControlIdentityError(f"{label} ARTIFACTS identity differs")
    verify_frozen_tree(path)
    return path


def _selected_rows(path: Path, *, model_id: str, head_id: str) -> tuple[bytes, int]:
    selected: list[bytes] = []
    identities: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != PREDICTION_FIELDS:
            raise SharedControlIdentityError(f"prediction schema differs: {path}")
        for row in reader:
            if row["model_id"] != model_id or row["head_id"] != head_id:
                continue
            identity = (row["seed"], row["row_hash"])
            if identity in identities:
                raise SharedControlIdentityError("duplicate shared-control identity")
            identities.add(identity)
            selected.append(("\t".join(row[field] for field in PREDICTION_FIELDS) + "\n").encode("utf-8"))
    return b"".join(selected), len(identities)


def verify(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config.resolve(strict=True))
    binding = config["shared_control_reference"]
    if arguments.output.exists():
        raise SharedControlIdentityError("identity output exists")
    reference_tree = _tree(root, binding, label="shared_control_reference")
    reference_member = reference_tree / binding["member"]
    if file_sha256(reference_member) != binding["member_sha256"]:
        raise SharedControlIdentityError("reference prediction member differs")
    candidate_tree = arguments.candidate_fit_tree.resolve(strict=True)
    if file_sha256(candidate_tree / "ARTIFACTS.json") != arguments.candidate_fit_sha256:
        raise SharedControlIdentityError("candidate fit ARTIFACTS identity differs")
    verify_frozen_tree(candidate_tree)
    reference, reference_count = _selected_rows(
        reference_member,
        model_id=binding["model_id"],
        head_id=binding["head_id"],
    )
    candidate, candidate_count = _selected_rows(
        candidate_tree / "oof_predictions.tsv.gz",
        model_id=binding["model_id"],
        head_id=binding["head_id"],
    )
    expected = int(binding["expected_rows"])
    if reference_count != expected or candidate_count != expected:
        raise SharedControlIdentityError("shared-control row denominator differs")
    if reference != candidate:
        raise SharedControlIdentityError("shared allele-control OOF rows are not bit-identical")
    digest = sha256(reference).hexdigest()
    arguments.output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_bit_identical_shared_allele_control_oof_rows",
        "dataset_id": "gse281364",
        "model_id": binding["model_id"],
        "head_id": binding["head_id"],
        "rows": expected,
        "canonical_selected_rows_sha256": digest,
        "reference_fit_artifacts_sha256": binding["artifacts_sha256"],
        "candidate_fit_artifacts_sha256": arguments.candidate_fit_sha256,
        "bit_identical": True,
        "comparison_performed_before_family_evaluation": True,
        "sealed_assets_read": False,
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidate-fit-tree", type=Path, required=True)
    parser.add_argument("--candidate-fit-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    verify(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
