#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.external_reference_common import recover_external26_reference_export


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--roster", required=True,
        choices=["common_strict7", "external_clean26"],
    )
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    result = recover_external26_reference_export(
        load_config(args.config), args.prepared, args.prepared_lock,
        args.policy, args.source, args.output, args.roster, args.seed,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
