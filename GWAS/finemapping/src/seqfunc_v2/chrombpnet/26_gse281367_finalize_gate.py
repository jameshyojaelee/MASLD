#!/usr/bin/env python3
"""Finalize the comparison-only GSE281367 gate without saturation authority."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=Path, required=True)
    args = ap.parse_args()
    base = args.base.resolve()
    internal = []
    for fold in range(5):
        path = base / f"fold{fold}/internal_qc.json"
        with path.open() as handle:
            internal.append(json.load(handle))
    with (base / "currin_calibration/currin_external_calibration.json").open() as handle:
        external = json.load(handle)
    all_internal_pass = all(bool(x.get("internal_qc_pass")) for x in internal)
    external_pass = bool(external.get("currin_external_calibration_pass"))
    verdict = {
        "all_internal_qc_pass": all_internal_pass,
        "internal_qc": internal,
        "currin_external_calibration_pass": external_pass,
        "currin_external_calibration": external,
        "comparison_ready": bool(all_internal_pass and external_pass),
        "saturation_authority": False,
        "saturation_unlocked": False,
        "rule": "comparison_ready = all_internal_qc_pass AND currin_external_calibration_pass; saturation remains review-blocked",
        "status": "comparison_ready" if all_internal_pass and external_pass else "comparison_gate_failed",
    }
    dest = base / "gate_verdict.json"
    tmp = dest.with_suffix(".json.tmp")
    with tmp.open("w") as handle:
        json.dump(verdict, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, dest)
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
