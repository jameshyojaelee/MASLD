#!/usr/bin/env python3
"""Summarize frozen StabMap seed evaluations without treating seeds as replicates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import mean, pstdev
from typing import Any, Mapping, Sequence


class StabMapStabilityError(ValueError):
    """Raised when the frozen StabMap seed evaluation census differs."""


def parse_evaluations(values: Sequence[str]) -> dict[int, Path]:
    roots: dict[int, Path] = {}
    for value in values:
        seed_text, separator, root_text = value.partition("=")
        if not separator or not seed_text.isdigit():
            raise StabMapStabilityError("seed-evaluation syntax differs")
        seed = int(seed_text)
        if seed in roots:
            raise StabMapStabilityError("seed is duplicated")
        roots[seed] = Path(root_text)
    if len(roots) != 3:
        raise StabMapStabilityError("exactly three seed evaluations are required")
    return roots


def summarize(evaluations: Mapping[int, Path], output: Path) -> dict[str, Any]:
    if output.exists():
        raise StabMapStabilityError("output exists")
    records: dict[str, Any] = {}
    gains: list[float] = []
    candidate_mrr: list[float] = []
    baseline_mrr: list[float] = []
    for seed, root in sorted(evaluations.items()):
        result = json.loads((root / "evaluation/evaluation.json").read_text())
        if (
            result.get("status") != "pass"
            or result.get("model_roster") != ["stabmap", "linear_cca"]
            or result.get("task_id") != "same_nucleus_cross_modal_retrieval"
            or result.get("outer_unit") != "donor"
            or result.get("n_donors") != 39
            or result.get("n_same_nucleus_pairs") != 1000
            or result.get("hidden_pair_map_available_to_model") is not False
            or result.get("cells_used_as_biological_replicates") is not False
            or result.get("test_outcomes_read") is not False
        ):
            raise StabMapStabilityError(f"seed {seed} evaluation differs")
        gain = float(result["mrr_gain_over_linear_cca"])
        candidate = float(
            result["summaries"]["stabmap"]["donor_macro_mean_reciprocal_rank"]
        )
        baseline = float(
            result["summaries"]["linear_cca"]["donor_macro_mean_reciprocal_rank"]
        )
        if not all(math.isfinite(value) for value in (gain, candidate, baseline)):
            raise StabMapStabilityError("non-finite stability value")
        gains.append(gain)
        candidate_mrr.append(candidate)
        baseline_mrr.append(baseline)
        records[str(seed)] = {
            "stabmap_donor_macro_mrr": candidate,
            "linear_cca_donor_macro_mrr": baseline,
            "mrr_gain_over_linear_cca": gain,
            "gain_direction_positive": gain > 0,
        }
    summary = {
        "schema_version": "masld-bench-stabmap-seed-stability-v1",
        "status": "pass",
        "model_id": "stabmap",
        "task_id": "same_nucleus_cross_modal_retrieval",
        "dataset_id": "gse296875",
        "seed_count": 3,
        "seeds": records,
        "mean_stabmap_donor_macro_mrr": mean(candidate_mrr),
        "sd_stabmap_donor_macro_mrr": pstdev(candidate_mrr),
        "mean_linear_cca_donor_macro_mrr": mean(baseline_mrr),
        "mean_mrr_gain_over_linear_cca": mean(gains),
        "sd_mrr_gain_over_linear_cca": pstdev(gains),
        "positive_gain_seed_count": sum(value > 0 for value in gains),
        "gain_direction_in_at_least_four_of_five_seeds": False,
        "interpretation": "three-seed screening stability; seeds are not biological replicates",
        "hidden_pair_map_available_to_models": False,
        "cells_used_as_biological_replicates": False,
        "test_outcomes_read": False,
        "rna_conditioned_atac_prediction_claim": False,
        "sealed_inference_eligible": False,
        "champion_claim_allowed": False,
    }
    output.mkdir(parents=True, mode=0o750)
    (output / "stability.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-evaluation", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    summarize(parse_evaluations(arguments.seed_evaluation), arguments.output)


if __name__ == "__main__":
    main()
