#!/usr/bin/env python3
"""Emit only the prespecified metrics permitted before selection."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.control_evaluation import evaluate_control_metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--candidate-embedding", required=True)
    parser.add_argument("--setting-id", required=True)
    parser.add_argument("--ewc-lambda", required=True, type=float)
    parser.add_argument("--replay-fraction", required=True, type=float)
    parser.add_argument("--reference-execution-record", required=True)
    parser.add_argument("--candidate-execution-record", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_control_metrics(
        load_config(args.config), args.reference_embedding, args.candidate_embedding,
        args.setting_id, args.ewc_lambda, args.replay_fraction,
        args.reference_execution_record, args.candidate_execution_record, args.output,
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
