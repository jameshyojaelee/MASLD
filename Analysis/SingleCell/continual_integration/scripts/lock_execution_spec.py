#!/usr/bin/env python3
"""Attach a pre-run source/input lock to an existing immutable execution spec."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_cl.execution import write_execution_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    args = parser.parse_args()
    spec = Path(args.spec).resolve()
    output = spec.with_suffix(spec.suffix + ".source-lock.json")
    lock = write_execution_lock(
        spec, Path(__file__).resolve().parents[1], output
    )
    print(json.dumps(lock, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
