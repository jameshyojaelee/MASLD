#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.orthogonal_bridge_retarget import build_retargeted_bridge_grid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--parent-embedding", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    result = build_retargeted_bridge_grid(
        load_config(args.config), args.policy, args.parent_embedding, args.output_root
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
