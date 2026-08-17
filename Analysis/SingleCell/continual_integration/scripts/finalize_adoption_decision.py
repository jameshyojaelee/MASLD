#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.adoption_decision import write_adoption_decision
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--rescue-decision", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--outcome-decision", required=True)
    parser.add_argument("--seed-confirmation", required=True)
    parser.add_argument("--order-invariant", required=True)
    parser.add_argument("--held-study", required=True)
    parser.add_argument("--label-benchmark", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = write_adoption_decision(
        load_config(args.config), args.rescue_decision, args.selection_lock,
        args.outcome_decision, args.seed_confirmation, args.order_invariant,
        args.held_study, args.label_benchmark, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
