#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.training import train_de_novo


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model-kind", required=True)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--selection-lock")
    parser.add_argument("--query-dataset", action="append")
    args = parser.parse_args()
    result = train_de_novo(
        load_config(args.config), args.prepared, args.contract_lock, args.prepared_lock, args.output,
        args.model_kind, args.seed, args.production,
        args.selection_lock,
        args.query_dataset,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
