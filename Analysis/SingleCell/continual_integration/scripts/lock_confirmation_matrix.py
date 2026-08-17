#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.confirmation_matrix import write_confirmation_matrix_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--execution-record", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = write_confirmation_matrix_lock(
        load_config(args.config), args.selection_lock,
        args.execution_record, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
