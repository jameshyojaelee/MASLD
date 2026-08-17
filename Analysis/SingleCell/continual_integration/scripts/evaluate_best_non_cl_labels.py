#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.non_cl_label_benchmark import evaluate_best_non_cl_labels


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--candidate-embedding", required=True)
    parser.add_argument("--embedding", action="append", default=[])
    parser.add_argument("--execution-record", action="append", default=[])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_best_non_cl_labels(
        load_config(args.config), args.policy, args.selection_lock,
        args.candidate_embedding, args.embedding, args.execution_record, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
