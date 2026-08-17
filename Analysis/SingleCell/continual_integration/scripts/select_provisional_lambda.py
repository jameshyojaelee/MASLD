#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.selection import write_provisional_lambda_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--control-metrics", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = write_provisional_lambda_lock(
        load_config(args.config), args.control_metrics, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
