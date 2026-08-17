#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from masld_cl.config import canonical_json_bytes, write_json_exclusive
from masld_cl.execution import ALLOWED_GPU_SCRIPTS, write_execution_lock


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resource-class", required=True, choices=["small", "large"])
    parser.add_argument("--script", required=True, choices=sorted(ALLOWED_GPU_SCRIPTS))
    parser.add_argument("--argument", action="append", default=[])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    spec = {
        "schema_version": "masld-cl-execution-spec-v2",
        "resource_class": args.resource_class,
        "script": args.script,
        "arguments": args.argument,
    }
    spec["spec_sha256"] = hashlib.sha256(canonical_json_bytes(spec)).hexdigest()
    write_json_exclusive(args.output, spec)
    pipeline_root = Path(__file__).resolve().parents[1]
    lock_path = Path(args.output).with_suffix(Path(args.output).suffix + ".source-lock.json")
    lock = write_execution_lock(args.output, pipeline_root, lock_path)
    print(json.dumps({"spec": spec, "execution_lock": lock}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
