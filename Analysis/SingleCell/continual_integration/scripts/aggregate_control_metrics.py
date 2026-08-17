#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.control_aggregation import aggregate_control_metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = aggregate_control_metrics(load_config(args.config), args.input, args.output)
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
