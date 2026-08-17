#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.run_plan import create_run_plan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--best-lambda", type=float)
    parser.add_argument("--selection-lock")
    args = parser.parse_args()
    result = create_run_plan(
        load_config(args.config), args.output, best_lambda=args.best_lambda,
        selection_lock=args.selection_lock,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
