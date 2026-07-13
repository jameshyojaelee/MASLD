#!/usr/bin/env python3
"""Combine independent internal and external gates without conflating them."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"


def load(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--internal", type=Path, required=True, help="fold0 held-out internal-QC JSON")
    ap.add_argument("--external", type=Path, required=True, help="Currin calibration JSON with pass boolean")
    ap.add_argument("--out", type=Path, default=OUT / "gate_verdict.json")
    args = ap.parse_args()
    internal = load(args.internal)
    external = load(args.external)
    metrics = {
        "count_pearson": float(internal["count_pearson"]),
        "count_spearman": float(internal["count_spearman"]),
        "median_profile_jsd": float(internal["median_profile_jsd"]),
    }
    internal_pass = (
        metrics["count_pearson"] >= 0.55
        and metrics["count_spearman"] >= 0.50
        and metrics["median_profile_jsd"] <= 0.50
    )
    external_pass = bool(external["currin_external_calibration_pass"])
    verdict = {
        "fold0_internal_qc_pass": internal_pass,
        "fold0_internal_metrics": metrics,
        "fold0_internal_thresholds": {"count_pearson_min": 0.55, "count_spearman_min": 0.50, "median_profile_jsd_max": 0.50},
        "currin_external_calibration_pass": external_pass,
        "currin_external_report": str(args.external.resolve()),
        "saturation_unlocked": internal_pass and external_pass,
        "rule": "saturation_unlocked = fold0_internal_qc_pass AND currin_external_calibration_pass",
        "status": "pass" if internal_pass and external_pass else "locked",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(args.out.suffix + ".tmp")
    with tmp.open("w") as handle:
        json.dump(verdict, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, args.out)
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

