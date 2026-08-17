#!/usr/bin/env python3
"""Run the locked V32 confirmation bundle."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.routed_v9_confirmation import confirm_routed_v9


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = confirm_routed_v9(load_config(args.config), args.policy, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
