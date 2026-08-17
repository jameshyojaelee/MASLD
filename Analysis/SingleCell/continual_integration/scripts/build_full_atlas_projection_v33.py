#!/usr/bin/env python3
"""Build the immutable V33 full-atlas projection."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.full_atlas_projection import build_full_atlas_projection


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_full_atlas_projection(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
