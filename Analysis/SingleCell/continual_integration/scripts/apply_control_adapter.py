#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.control_adapter import apply_control_adapter


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--base-embedding", required=True)
    parser.add_argument("--base-execution-record", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--reference-execution-record", required=True)
    parser.add_argument("--global-reference-weight", required=True, type=float)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = apply_control_adapter(
        load_config(args.config), args.policy, args.base_embedding,
        args.base_execution_record, args.reference_embedding,
        args.reference_execution_record, args.output, args.global_reference_weight,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
