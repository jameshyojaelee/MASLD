#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.confirmation_evaluation import evaluate_confirmation_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--bundle", required=True)
    parser.add_argument("--gpu-tolerance", required=True)
    parser.add_argument("--descriptive-projection-audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = evaluate_confirmation_bundle(
        load_config(args.config), args.selection_lock, args.bundle,
        args.gpu_tolerance, args.descriptive_projection_audit, args.output,
    )
    print(json.dumps(rows, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
