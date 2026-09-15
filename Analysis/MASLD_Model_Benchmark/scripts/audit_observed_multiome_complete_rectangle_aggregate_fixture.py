#!/usr/bin/env python3
"""Freeze the synthetic-only complete-rectangle aggregator fixture."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from masld_bench.artifacts import freeze_tree, reject_symlink_components, write_json_exclusive
from masld_bench.observed_multiome_aggregate import BOOTSTRAP_REPLICATES, BOOTSTRAP_SEED, expected_surface_keys


class AggregateFixtureError(ValueError):
    """Raised when the prospective aggregator fixture differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def validate(root: Path) -> dict[str, object]:
    module = root / "src/masld_bench/observed_multiome_aggregate.py"
    test = root / "tests/unit/test_observed_multiome_aggregate.py"
    sbatch = root / "slurm/test_observed_multiome_complete_rectangle_aggregate_cpu.sbatch"
    for path in (module, test, sbatch):
        path.resolve(strict=True).relative_to(root)
    module_text = module.read_text(encoding="utf-8")
    for token in ("preflight_bundles", "load_complete_records", "aggregate_records", "arithmetic_mean_prediction_across_five_fixed_seeds", "two_way_donor_block_bootstrap", "partial_rectangle_used"):
        if token not in module_text:
            raise AggregateFixtureError(f"aggregator token is missing: {token}")
    test_text = test.read_text(encoding="utf-8")
    for token in ("test_complete_rectangle_uses_ensemble_and_two_way_units", "test_partial_rectangle_is_rejected_before_values_are_aggregated", "test_seed_join_drift_is_rejected"):
        if token not in test_text:
            raise AggregateFixtureError(f"fixture test is missing: {token}")
    if len(expected_surface_keys()) != 125 or BOOTSTRAP_REPLICATES != 10_000 or BOOTSTRAP_SEED != 20260825:
        raise AggregateFixtureError("rectangle or inference constant differs")
    forbidden = ("21100040", "21100041", "21100042", "21100043", "21100044", "model-cpu-train-202608")
    if any(token in module_text or token in test_text for token in forbidden):
        raise AggregateFixtureError("production bundle reference entered synthetic fixture")
    return {
        "schema_version": "masld-bench-observed-multiome-complete-rectangle-aggregate-fixture-receipt-v1",
        "dataset_id": "gse296875",
        "synthetic_surface_seed_runs_tested": 125,
        "synthetic_biological_donor_identifiers": 39,
        "synthetic_lineages": 5,
        "synthetic_genomic_blocks": 5,
        "arithmetic_seed_ensemble_tested": True,
        "partial_rectangle_rejection_tested": True,
        "seed_join_rejection_tested": True,
        "two_way_bootstrap_tested": True,
        "bootstrap_replicates_production": BOOTSTRAP_REPLICATES,
        "bootstrap_seed_production": BOOTSTRAP_SEED,
        "production_bundle_paths_bound": False,
        "production_prediction_values_read": False,
        "production_evaluator_values_read": False,
        "development_metric_calculated": False,
        "partial_ranking_authorized": False,
        "production_aggregation_authorized": False,
        "sealed_outcomes_read": False,
        "next_gate": "bind_all_five_complete_seed_bundles_then_independent_production_aggregation_admission"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    output = reject_symlink_components(args.output, label="aggregate fixture output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "receipt.json", validate(root))
    sources = [Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_aggregate.py", root / "tests/unit/test_observed_multiome_aggregate.py", root / "slurm/test_observed_multiome_complete_rectangle_aggregate_cpu.sbatch"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_complete_rectangle_aggregate_fixture", "production_aggregation_authorized": False, "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
