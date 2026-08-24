#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.gse296875_router import prepare_gse296875_query_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--mapping-policy", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = prepare_gse296875_query_counts(
        load_config(args.config), args.mapping_policy, args.metadata, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
