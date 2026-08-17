#!/usr/bin/env python3
"""Build the locked, condition-blind GSE212837 mapping."""

import argparse
import json

from masld_cl.config import load_config
from masld_cl.external_gse212837 import build_external_gse212837_mapping


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = build_external_gse212837_mapping(
        load_config(args.config), args.policy, args.output
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
