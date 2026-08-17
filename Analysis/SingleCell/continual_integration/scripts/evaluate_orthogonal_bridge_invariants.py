#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.orthogonal_bridge_invariants import evaluate_order_invariant, evaluate_prediction_invariant


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--outcome-decision", required=True)
    parser.add_argument("--prediction-output", required=True)
    parser.add_argument("--order-output", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    prediction = evaluate_prediction_invariant(config, args.selection_lock, args.prediction_output)
    order = evaluate_order_invariant(config, args.selection_lock, args.outcome_decision, args.order_output)
    print(json.dumps({"prediction": prediction, "order": order}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
