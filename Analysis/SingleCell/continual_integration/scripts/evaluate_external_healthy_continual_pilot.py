#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.external_healthy_continual import evaluate_external_healthy_continual_pilot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--common", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--anchor-policy")
    args = parser.parse_args()
    result = evaluate_external_healthy_continual_pilot(
        load_config(args.config), args.policy, args.common,
        args.candidate, args.output, args.anchor_policy,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
