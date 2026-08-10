#!/usr/bin/env python3
"""Validate the signed Plan 60 workstream closure matrix without assembling it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from release_common import CANDIDATE_ID, ReleaseContractError, validate_closure


def build_summary(project_root: Path, closure_path: Path) -> dict[str, object]:
    bundle = validate_closure(project_root, closure_path)
    return {
        "status": "REL00_CLOSURE_READY",
        "candidate_id": CANDIDATE_ID,
        "closure_sha256": bundle.closure_sha256,
        "n_workstreams": len(bundle.workstreams),
        "n_artifacts": len(bundle.artifacts),
        "workstreams": [
            {
                "workstream_id": workstream.row["workstream_id"],
                "terminal_state": workstream.row["terminal_state"],
                "gate_verdict": workstream.row["gate_verdict"],
                "manifest_sha256": workstream.manifest_sha256,
                "n_artifacts": len(workstream.artifacts),
            }
            for workstream in bundle.workstreams
        ],
        "assembly_performed": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--closure", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        summary = build_summary(args.project_root.resolve(), args.closure)
    except ReleaseContractError as error:
        raise SystemExit(f"REL00_BLOCKED: {error}") from error
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
