#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.distillation_pilot import evaluate_distillation_grid


def _candidate(value: str) -> tuple[float, str, str]:
    fields = value.split("|", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError(
            "candidate must be ALPHA|EMBEDDING|EXECUTION_RECORD"
        )
    return float(fields[0]), fields[1], fields[2]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--reference-execution-record", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--architecture-execution-record", required=True)
    parser.add_argument("--candidate", action="append", type=_candidate, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_distillation_grid(
        load_config(args.config), args.policy,
        args.reference_embedding, args.reference_execution_record,
        args.harmony_embedding, args.architecture_embedding,
        args.architecture_execution_record, args.candidate, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
