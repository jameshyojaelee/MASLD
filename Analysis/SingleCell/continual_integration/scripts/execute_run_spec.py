#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.execution import execute_run_spec


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--pipeline-root", required=True)
    parser.add_argument("--record-dir", required=True)
    args = parser.parse_args()
    print(json.dumps(
        execute_run_spec(args.spec, args.pipeline_root, args.record_dir),
        indent=2, sort_keys=True,
    ))


if __name__ == "__main__":
    main()
