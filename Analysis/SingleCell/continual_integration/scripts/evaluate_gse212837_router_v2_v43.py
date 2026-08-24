#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.router_v2_evaluation import evaluate_gse212837_router_v2_diagnostic


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--router-manifest", required=True)
    parser.add_argument("--v40-policy", required=True)
    parser.add_argument("--v40-mapping-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_gse212837_router_v2_diagnostic(
        load_config(args.config), args.policy, args.router_manifest,
        args.v40_policy, args.v40_mapping_manifest, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
