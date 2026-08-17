#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.audits import audit_descriptive_projection
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--update-manifest", required=True)
    parser.add_argument("--embedding-manifest", required=True)
    parser.add_argument("--execution-record", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = audit_descriptive_projection(
        load_config(args.config), args.selection_lock,
        args.update_manifest, args.embedding_manifest,
        args.execution_record, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
