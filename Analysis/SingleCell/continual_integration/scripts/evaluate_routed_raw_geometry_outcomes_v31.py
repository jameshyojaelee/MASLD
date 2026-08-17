#!/usr/bin/env python3
"""Run the post-lock V31 combined identity and disease audit."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.replay_ewc_latent_outcomes import (
    evaluate_replay_ewc_latent_adapter_outcomes,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_replay_ewc_latent_adapter_outcomes(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
