#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.training import train_sequence


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--refit-lock", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model-kind", required=True)
    parser.add_argument("--order-name", required=True)
    parser.add_argument("--ewc-lambda", required=True, type=float)
    parser.add_argument("--replay-fraction", required=True, type=float)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    result = train_sequence(
        load_config(args.config), args.prepared, args.contract_lock,
        args.prepared_lock, args.reference_model, args.output, args.model_kind,
        args.order_name, args.ewc_lambda, args.replay_fraction, args.seed,
        args.selection_lock, args.refit_lock,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
