#!/usr/bin/env python3
"""Score frozen scVI/scANVI predictions with the independent cell evaluator.

This adapter changes only the included model roster.  It delegates all joins,
metrics, donor balancing, calibration summaries, and 10,000-donor-bootstrap
calculations to the existing independent evaluator.  Within-lane selection is
descriptive and does not establish superiority over the classical baselines.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import score_cell_baselines_study_50000 as evaluator


MODEL_IDS = ("scvi_baseline", "scanvi_baseline")
SEEDS = (20260824, 20260825, 20260826)


def run(
    *,
    source: Path,
    split: Path,
    predictions: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    expected_prediction_artifacts_sha256: str,
) -> None:
    evaluator.MODEL_IDS = MODEL_IDS
    evaluator.SCREEN_SEEDS = SEEDS
    evaluator.run(
        source=source,
        split=split,
        predictions=predictions,
        output=output,
        expected_source_artifacts_sha256=expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=expected_split_artifacts_sha256,
        expected_prediction_artifacts_sha256=expected_prediction_artifacts_sha256,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--split", required=True, type=Path)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-source-artifacts-sha256", required=True)
    parser.add_argument("--expected-split-artifacts-sha256", required=True)
    parser.add_argument("--expected-prediction-artifacts-sha256", required=True)
    arguments = parser.parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        predictions=arguments.predictions,
        output=arguments.output,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
        expected_prediction_artifacts_sha256=arguments.expected_prediction_artifacts_sha256,
    )
    print(json.dumps({"output": str(arguments.output.resolve()), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
