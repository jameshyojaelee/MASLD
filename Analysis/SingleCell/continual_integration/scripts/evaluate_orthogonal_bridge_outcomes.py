#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.orthogonal_bridge_outcome import evaluate_orthogonal_bridge_outcomes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--candidate-embedding", required=True)
    parser.add_argument("--raw-pca-embedding", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_orthogonal_bridge_outcomes(
        load_config(args.config), args.selection_lock, args.candidate_embedding,
        args.raw_pca_embedding, args.output,
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
