#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.control_adapter_pilot import evaluate_orthogonal_bridge_pilot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--architecture-embedding", required=True)
    parser.add_argument("--candidate-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = evaluate_orthogonal_bridge_pilot(
        load_config(args.config), args.reference_embedding, args.harmony_embedding,
        args.architecture_embedding, args.candidate_embedding, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
