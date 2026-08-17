#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.benchmark import evaluate_benchmarks
from masld_cl.config import load_config


def _mapping(values):
    result = {}
    for value in values:
        key, separator, path = value.partition("=")
        if not separator or not key or not path or key in result:
            raise argparse.ArgumentTypeError("bundle arguments must be unique MODEL_KIND=MANIFEST")
        result[key] = path
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--non-cl", action="append", required=True)
    parser.add_argument("--candidate", action="append", required=True)
    parser.add_argument("--architecture", action="append", required=True)
    parser.add_argument("--execution-record", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_benchmarks(
        load_config(args.config), args.selection_lock, args.reference_embedding,
        args.harmony_embedding, _mapping(args.non_cl),
        _mapping(args.candidate), _mapping(args.architecture),
        args.execution_record, args.output,
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
