#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.external_healthy_continual import train_external_healthy_continual_pilot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--common-policy", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--replay-fraction", type=float, default=0.2)
    parser.add_argument("--ewc-lambda", type=float, default=100.0)
    args = parser.parse_args()
    result = train_external_healthy_continual_pilot(
        load_config(args.config), args.prepared, args.prepared_lock,
        args.common_policy, args.policy, args.reference_model, args.output,
        seed=args.seed, replay_fraction=args.replay_fraction,
        ewc_lambda=args.ewc_lambda,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
