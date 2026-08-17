#!/usr/bin/env python3
"""Build the prospective V28 class-routed adapter grid."""

from __future__ import annotations

import argparse
import json

from masld_cl.class_routed_adapter import build_class_routed_adapter_grid
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_class_routed_adapter_grid(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
