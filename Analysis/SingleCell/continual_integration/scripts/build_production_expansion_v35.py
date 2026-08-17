#!/usr/bin/env python3
"""Build the six-model, stage-blind V35 production expansion."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.production_expansion import build_production_expansion


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_production_expansion(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
