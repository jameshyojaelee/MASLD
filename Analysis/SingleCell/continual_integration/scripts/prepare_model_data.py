#!/usr/bin/env python3
"""Create the 4,000-reference-HVG count object on a compute node."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.data import prepare_hvg_anndata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--library-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = prepare_hvg_anndata(
        load_config(args.config), args.contract_lock, args.library_manifest, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

