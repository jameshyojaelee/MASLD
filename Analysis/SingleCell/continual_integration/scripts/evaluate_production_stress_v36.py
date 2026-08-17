#!/usr/bin/env python3
"""Evaluate V35 technical stress gates and preserve the leakage audit."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.production_stress import evaluate_production_stress


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_production_stress(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
