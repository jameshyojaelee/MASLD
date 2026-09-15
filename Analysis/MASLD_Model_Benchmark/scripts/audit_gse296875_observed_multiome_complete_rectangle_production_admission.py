#!/usr/bin/env python3
"""Authorize complete-rectangle aggregation without deserializing values."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, reject_symlink_components, write_json_exclusive
from masld_bench.observed_multiome_aggregate import preflight_bundles
from scripts.aggregate_gse296875_observed_multiome_factorized_complete_rectangle import validate_manifest


class ProductionAdmissionError(ValueError):
    """Raised when complete production aggregation is not eligible."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ProductionAdmissionError("JSON object required")
    return value


def validate(root: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    bundles = validate_manifest(root, manifest)
    surfaces = preflight_bundles(root, bundles)
    claim = manifest.get("claim_boundary")
    if claim != {"development_baseline_selection_only": True, "partial_ranking_forbidden": True, "external_claim_allowed": False, "champion_claim_allowed": False, "sealed_outcomes_read": False}:
        raise ProductionAdmissionError("claim boundary differs")
    return {
        "schema_version": "masld-bench-observed-multiome-complete-rectangle-production-admission-receipt-v1",
        "campaign_id": manifest["campaign_id"],
        "dataset_id": "gse296875",
        "seed_bundles_verified": len(bundles),
        "surface_seed_runs_verified": len(surfaces),
        "prediction_matrices_expected": len(surfaces) * 5,
        "all_bundle_trees_hash_verified": True,
        "all_surface_trees_hash_verified": True,
        "production_prediction_arrays_deserialized": False,
        "production_evaluator_arrays_deserialized": False,
        "development_metric_calculated": False,
        "partial_rectangle_used": False,
        "partial_ranking_authorized": False,
        "promotion_authorized": False,
        "sealed_outcomes_read": False,
        "production_aggregation_authorized": True,
        "next_gate": "one_complete_rectangle_aggregation_then_independent_metric_rederivation"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    manifest_path = args.manifest.resolve(strict=True)
    manifest_path.relative_to(root)
    output = reject_symlink_components(args.output, label="production aggregation admission output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(manifest_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [manifest_path, Path(__file__).resolve(strict=True), root / "scripts/aggregate_gse296875_observed_multiome_factorized_complete_rectangle.py", root / "src/masld_bench/observed_multiome_aggregate.py", root / "tests/unit/test_observed_multiome_complete_rectangle_production_admission.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_complete_rectangle_production_admission", "production_aggregation_authorized": True, "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
