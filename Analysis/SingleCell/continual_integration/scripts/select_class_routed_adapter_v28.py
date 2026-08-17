#!/usr/bin/env python3
"""Evaluate controls and freeze the prospective V28 selection."""

from __future__ import annotations

import argparse
import json

from masld_cl.class_routed_adapter import evaluate_and_select_class_routed_adapter
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--grid", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_and_select_class_routed_adapter(
        load_config(args.config), args.policy, args.grid,
        args.reference_embedding, args.harmony_embedding,
        args.architecture_embedding, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
