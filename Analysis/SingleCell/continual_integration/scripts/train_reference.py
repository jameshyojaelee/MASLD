#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.training import train_reference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model-kind", required=True)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    result = train_reference(
        load_config(args.config), args.prepared, args.contract_lock, args.prepared_lock, args.output,
        args.model_kind, args.seed,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
