#!/usr/bin/env python3
"""Verify an existing production contract lock before model work."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.contracts import verify_contract_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--lock", required=True)
    parser.add_argument("--full-raw-hash", action="store_true")
    args = parser.parse_args()
    lock = verify_contract_lock(load_config(args.config), args.lock, full_hash=args.full_raw_hash)
    print(json.dumps({"verified": True, "lock": lock}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

