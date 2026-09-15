#!/usr/bin/env python3
"""Hold back the GSE267145 paired-contrast prespecification before any contrast runs.

This is the first half of a two-job freeze.  It copies the prespecification
into a read-only execution tree and records the digests of every input the
analysis job will later bind to.  No outcome, no prediction and no metric is
read here, so the held-back criteria cannot have been chosen after seeing a
result.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file


PRESPEC_REVISION_ID = "model-check-1101-paired-contrast-v1"
BOOTSTRAP_INDICES_SHA256 = (
    "fb7ea4b28ef37ee759a70bb3664868f6affa5accc8dadd96edf202238536b2c5"
)
BOUND_INPUTS = {
    "scorer_source": (
        "scripts/score_gse267145_histology_production_v4.py",
        "dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43",
    ),
    "reference_evaluator": (
        "scripts/evaluate_gse267145_histology_predictions.py",
        "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920",
    ),
    "task_contract": (
        "config/evaluation/gse267145_histology_state_task.toml",
        "53a85a843e75b1c0fa088af86cb2f2593976d28039807dfa17d1ec571603778c",
    ),
    "scoring_contract": (
        "config/evaluation/gse267145_histology_production_v4_scoring.json",
        "532f66eab44561b4b84124c745b24a839f7d267c970fe5fb1fa9206513be7971",
    ),
    "promotion_gate": (
        "config/evaluation/gse267145_histology_state_promotion_gate.json",
        "c8fd739e7d9cbaca9bf64729a8748dd0add42a9c76a5f7d4088362caf5add494",
    ),
    "analysis_source": (
        "scripts/analyze_gse267145_histology_paired_contrasts.py",
        None,
    ),
}
BOUND_TREES = {
    "scores": (
        "executions/model-scoring-078-21092130/scores",
        "a55aed3d9ac4cad2211ed3b319c45f911d05851b8559408ed99bf7fa59562393",
    ),
    "outcomes": (
        "executions/model-data-061-21079623/activation",
        "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a",
    ),
    "folds": (
        "executions/model-data-064-21079902/fixture/folds",
        "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3",
    ),
}


class PrespecSealError(RuntimeError):
    """Raised when a bound input does not authenticate."""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--prespec", type=Path, required=True)
    parser.add_argument("--prespec-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    root = arguments.root.resolve(strict=True)
    if sha256_file(arguments.prespec) != arguments.prespec_sha256:
        raise PrespecSealError("prespecification SHA-256 differs")
    prespec: dict[str, Any] = json.loads(arguments.prespec.read_text(encoding="utf-8"))
    if (
        prespec.get("status") != "sealed_before_any_contrast_is_computed"
        or prespec["resampling"]["bootstrap_indices_sha256"] != BOOTSTRAP_INDICES_SHA256
        or prespec["estimator"]["p_values_calculated"] is not False
        or prespec["estimator"]["bh_adjustment_calculated"] is not False
        or prespec["claim_boundary"]["model_ranking_allowed"] is not False
        or prespec["claim_boundary"]["champion_claim_allowed"] is not False
        or prespec["claim_boundary"]["confirmatory_inference_allowed"] is not False
        or not prespec["declared_outcome_cells"]
        or prespec["deposits"]["summaries_only"] is not False
    ):
        raise PrespecSealError("prespecification content differs")

    bound: dict[str, Any] = {}
    for label, (relative, expected) in BOUND_INPUTS.items():
        target = root / relative
        if not target.is_file() or target.is_symlink():
            raise PrespecSealError(f"bound input is missing: {relative}")
        digest = sha256_file(target)
        if expected is not None and digest != expected:
            raise PrespecSealError(f"bound input SHA-256 differs: {relative}")
        bound[label] = {"path": relative, "sha256": digest}
    for label, (relative, expected) in BOUND_TREES.items():
        target = root / relative
        manifest_digest = sha256_file(target / "ARTIFACTS.json")
        if manifest_digest != expected:
            raise PrespecSealError(f"bound tree SHA-256 differs: {relative}")
        verify_frozen_tree(target)
        bound[label] = {"path": relative, "artifacts_sha256": manifest_digest}

    arguments.output.mkdir(parents=True)
    write_json_exclusive(
        arguments.output / "prespec.json", prespec, mode=0o440
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-paired-contrast-prespec-seal-v1",
        "prespec_revision_id": PRESPEC_REVISION_ID,
        "prespec_contract_sha256": arguments.prespec_sha256,
        "bound_inputs": bound,
        "bootstrap_indices_sha256": BOOTSTRAP_INDICES_SHA256,
        "declared_contrasts": prespec["declared_contrast_totals"]["declared"],
        "computable_contrasts": prespec["declared_contrast_totals"]["computable"],
        "not_applicable_contrasts": prespec["declared_contrast_totals"]["not_applicable"],
        "declared_outcome_cell_count": len(prespec["declared_outcome_cells"]),
        "outcomes_read": False,
        "predictions_read": False,
        "metrics_calculated": False,
        "contrasts_calculated": False,
        "p_values_calculated": False,
        "bh_adjustment_calculated": False,
        "model_ranking_performed": False,
        "status": "sealed_prespec_no_outcome_read",
    }
    write_json_exclusive(
        arguments.output / "prespec_seal_receipt.json", receipt, mode=0o440
    )
    freeze_tree(
        arguments.output,
        {
            "artifact_class": "gse267145_histology_paired_contrast_prespec",
            "prespec_contract_sha256": arguments.prespec_sha256,
            "declared_contrasts": prespec["declared_contrast_totals"]["declared"],
            "outcomes_read": False,
            "metrics_calculated": False,
            "contrasts_calculated": False,
            "status": "sealed_prespec_no_outcome_read",
        },
    )
    print(json.dumps(receipt, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
