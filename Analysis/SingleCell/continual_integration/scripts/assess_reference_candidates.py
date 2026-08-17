#!/usr/bin/env python3
"""Write the fail-closed healthy-reference roster assessment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from masld_cl.config import write_json_exclusive
from masld_cl.reference_assessment import build_assessment, sha256_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--repository-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    assessment = build_assessment(args.policy, args.registry, args.repository_root)
    assessment["policy_sha256"] = sha256_path(args.policy)
    assessment["registry_sha256"] = sha256_path(args.registry)
    write_json_exclusive(output / "reference_assessment.json", assessment)
    (output / "REFERENCE_ASSESSMENT_COMPLETE").touch(exist_ok=False)
    print(json.dumps(assessment, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
