#!/usr/bin/env python3
"""Aggregate the exact complete observed-multiome rectangle once included."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive
from masld_bench.observed_multiome_aggregate import LINEAGES, aggregate_records, load_complete_records, preflight_bundles


SCHEMA = "masld-bench-observed-multiome-complete-rectangle-production-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_factorized_complete_rectangle_20260825"


class ProductionAggregateError(ValueError):
    """Raised when the included production aggregation surface differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ProductionAggregateError("JSON object required")
    return value


def validate_manifest(root: Path, manifest: dict[str, Any]) -> list[dict[str, str]]:
    if manifest.get("schema_version") != SCHEMA or manifest.get("campaign_id") != CAMPAIGN_ID or manifest.get("dataset_id") != "gse296875" or manifest.get("stage") != "development":
        raise ProductionAggregateError("production aggregation identity differs")
    fixture = manifest.get("aggregate_fixture")
    if not isinstance(fixture, dict):
        raise ProductionAggregateError("aggregate fixture binding is missing")
    fixture_path = reject_symlink_components(root / str(fixture.get("path", "")), label="aggregate fixture").resolve(strict=True)
    fixture_path.relative_to(root)
    verify_frozen_tree(fixture_path)
    if _digest(fixture_path / "ARTIFACTS.json") != fixture.get("artifacts_sha256") or _json(fixture_path / "receipt.json").get("production_aggregation_authorized") is not False:
        raise ProductionAggregateError("aggregate fixture differs")
    sources = manifest.get("source_bindings")
    if not isinstance(sources, dict) or set(sources) != {"aggregate_module", "production_driver", "production_sbatch"}:
        raise ProductionAggregateError("production source roster differs")
    for label, record in sources.items():
        path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record.get("sha256"):
            raise ProductionAggregateError(f"{label} drifted")
    bundles = manifest.get("seed_bundles")
    if not isinstance(bundles, list) or len(bundles) != 5 or {int(row.get("seed", -1)) for row in bundles} != {20260824, 20260825, 20260826, 20260827, 20260828}:
        raise ProductionAggregateError("seed bundle manifest differs")
    if manifest.get("revision") != 2 or manifest.get("supersedes_failed_job_id") != 21100166 or manifest.get("lineage_identifiers") != list(LINEAGES):
        raise ProductionAggregateError("production revision or lineage authority differs")
    inference = manifest.get("inference")
    if inference != {"seed_ensemble": "arithmetic_mean_prediction_across_five_fixed_seeds", "bootstrap_replicates": 10000, "bootstrap_seed": 20260825, "inferential_units": ["donor", "genomic_fold"], "lineage_macro_within_donor_block": True, "seeds_are_biological_replicates": False}:
        raise ProductionAggregateError("production inference contract differs")
    return bundles


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    manifest_path = args.manifest.resolve(strict=True)
    manifest_path.relative_to(root)
    manifest = _json(manifest_path)
    bundles = validate_manifest(root, manifest)
    surfaces = preflight_bundles(root, bundles)
    records = load_complete_records(surfaces)
    receipt, arrays = aggregate_records(records)
    output = reject_symlink_components(args.output, label="complete rectangle aggregation output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    with (output / "unit_results.npz").open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    receipt["campaign_id"] = CAMPAIGN_ID
    receipt["bundle_artifacts_sha256"] = {str(row["seed"]): str(row["artifacts_sha256"]) for row in bundles}
    receipt["production_aggregation_authorized_by_manifest"] = True
    receipt["partial_ranking_authorized"] = False
    receipt["promotion_authorized"] = False
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [manifest_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_aggregate.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_complete_rectangle_aggregation", "campaign_id": CAMPAIGN_ID, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
