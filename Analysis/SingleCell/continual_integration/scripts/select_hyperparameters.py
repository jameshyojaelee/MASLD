#!/usr/bin/env python3
"""Freeze a setting using reference and query-control metrics only."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.firewall import validate_program_firewall
from masld_cl.selection import write_selection_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--control-metrics", required=True)
    parser.add_argument("--provisional-lambda-lock", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    validate_program_firewall(config)
    lock = write_selection_lock(
        config, args.control_metrics, args.provisional_lambda_lock, args.output
    )
    print(json.dumps(lock, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
