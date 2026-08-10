#!/usr/bin/env python3
"""Fail-closed integrity check for resumable Plan 11 pipeline stages."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from yakubovsky_common import (
    RELEASE_ID,
    STAGE_ORDER,
    default_paths,
    validate_stage_chain,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--stage", choices=STAGE_ORDER, required=True)
    args = parser.parse_args()
    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    seals = validate_stage_chain(output, args.stage)
    result: dict[str, object] = {
        "release_id": RELEASE_ID,
        "stage": args.stage,
        "status": "seal_rederived",
    }
    result["stage_seal_sha256"] = seals[args.stage]
    result["validated_stage_chain"] = list(seals)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
