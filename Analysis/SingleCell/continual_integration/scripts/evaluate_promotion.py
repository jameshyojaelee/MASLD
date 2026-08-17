#!/usr/bin/env python3
"""Apply all prespecified gates and write promotion_decision.json."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.promotion import write_promotion_decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--metrics-lock", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = write_promotion_decision(
        load_config(args.config), args.metrics, args.metrics_lock, args.selection_lock, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
