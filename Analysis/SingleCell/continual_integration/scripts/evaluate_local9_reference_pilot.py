#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.local9_reference_pilot import evaluate_local9_reference_pilot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--current-reference-embedding", required=True)
    parser.add_argument("--current-reference-record", required=True)
    parser.add_argument("--local9-reference-embedding", required=True)
    parser.add_argument("--local9-reference-record", required=True)
    parser.add_argument("--raw-pca-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_local9_reference_pilot(
        load_config(args.config), args.policy,
        args.current_reference_embedding, args.current_reference_record,
        args.local9_reference_embedding, args.local9_reference_record,
        args.raw_pca_embedding, args.harmony_embedding, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
