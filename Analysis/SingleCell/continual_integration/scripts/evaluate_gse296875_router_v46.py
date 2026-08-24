#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.gse296875_router import evaluate_gse296875_router


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--mapping-policy", required=True)
    parser.add_argument("--routing-manifest", required=True)
    parser.add_argument("--author-labels", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_gse296875_router(
        load_config(args.config), args.mapping_policy, args.routing_manifest,
        args.author_labels, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
