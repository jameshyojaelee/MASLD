#!/usr/bin/env python3
"""Measure, never assume, same-seed tolerance for two selected GPU fits."""

from __future__ import annotations

import argparse
import json

from masld_cl.audits import measure_gpu_tolerance
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--left-embedding", required=True)
    parser.add_argument("--right-embedding", required=True)
    parser.add_argument("--left-execution-record", required=True)
    parser.add_argument("--right-execution-record", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = measure_gpu_tolerance(
        load_config(args.config), args.selection_lock,
        args.left_embedding, args.right_embedding,
        args.left_execution_record, args.right_execution_record,
        args.seed, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
