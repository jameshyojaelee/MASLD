#!/usr/bin/env python3
"""Build immutable atlas manifests and the production contract lock."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_cl.config import load_config
from masld_cl.contracts import write_contract_outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--development-skip-raw-hash",
        action="store_true",
        help="Write a non-production lock without streaming raw/X; training rejects it.",
    )
    parser.add_argument(
        "--skip-cell-manifest",
        action="store_true",
        help="Omit the 1.23M-row cell manifest for a quick contract audit.",
    )
    args = parser.parse_args()
    config = load_config(args.config)
    lock = write_contract_outputs(
        config,
        Path(args.output),
        hash_counts=not args.development_skip_raw_hash,
        write_cells=not args.skip_cell_manifest,
    )
    print(json.dumps(lock, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

