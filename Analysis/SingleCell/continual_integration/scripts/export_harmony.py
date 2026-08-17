#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.incumbent import export_harmony


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = export_harmony(
        load_config(args.config), args.contract_lock, args.output
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
