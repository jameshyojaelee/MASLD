#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.outcome_evaluation import evaluate_control_adapter_outcome_preservation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--candidate-embedding", required=True)
    parser.add_argument("--raw-pca-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_control_adapter_outcome_preservation(
        load_config(args.config), args.selection_lock, args.candidate_embedding,
        args.raw_pca_manifest, args.output,
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
