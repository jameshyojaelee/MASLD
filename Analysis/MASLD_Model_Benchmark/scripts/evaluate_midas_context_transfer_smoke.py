#!/usr/bin/env python3
"""Evaluate frozen inductive MIDAS held-donor RNA-to-ATAC profiles."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys


BASE_PATH = Path(__file__).with_name("evaluate_multivi_smoke.py")
SPEC = importlib.util.spec_from_file_location("midas_context_evaluator_base", BASE_PATH)
assert SPEC is not None and SPEC.loader is not None
base = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = base
SPEC.loader.exec_module(base)
base.CANDIDATE_MODEL_ID = "midas_inductive"
base.PAIRING_TOPOLOGY = "same_nucleus_training"
base.MODEL_IDS = (
    "midas_inductive",
    "training_lineage_mean",
    "training_global_mean",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--source-h5", type=Path, required=True)
    parser.add_argument("--source-h5-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    base.evaluate(
        arguments.prediction_root,
        arguments.prepared,
        arguments.source_h5,
        arguments.source_h5_sha256,
        arguments.output,
    )


if __name__ == "__main__":
    main()
