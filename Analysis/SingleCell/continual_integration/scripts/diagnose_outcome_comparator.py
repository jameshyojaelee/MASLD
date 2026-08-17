#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.outcome_comparator_diagnostic import diagnose_outcome_comparator


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--harmony-manifest", required=True)
    parser.add_argument("--raw-pca-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = diagnose_outcome_comparator(
        load_config(args.config), args.selection_lock, args.harmony_manifest,
        args.raw_pca_manifest, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
