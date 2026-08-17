#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.orthogonal_bridge import build_orthogonal_bridge_grid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--raw-pca-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--architecture-execution-record", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--reference-execution-record", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    result = build_orthogonal_bridge_grid(
        load_config(args.config), args.policy, args.raw_pca_embedding,
        args.architecture_embedding, args.architecture_execution_record,
        args.reference_embedding, args.reference_execution_record, args.output_root,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
