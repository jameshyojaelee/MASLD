#!/usr/bin/env python3
"""Normalize package-native held-out fold metrics and apply frozen thresholds."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
OUT = ROOT / "GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--out-root", type=Path, default=OUT)
    ap.add_argument("--require-pass", action="store_true")
    args = ap.parse_args()
    source = args.out_root / f"fold{args.fold}/chrombpnet_full/evaluation/chrombpnet_metrics.json"
    with source.open() as handle:
        raw = json.load(handle)
    peak_counts = raw["counts_metrics"]["peaks"]
    peak_profile = raw["profile_metrics"]["peaks"]
    result = {
        "fold": args.fold,
        "source": str(source.resolve()),
        "metric_scope": "package-native held-out test chromosomes; peaks",
        "count_pearson": float(peak_counts["pearsonr"]),
        "count_spearman": float(peak_counts["spearmanr"]),
        "median_profile_jsd": float(peak_profile["median_jsd"]),
    }
    result["internal_qc_pass"] = (
        result["count_pearson"] >= 0.55
        and result["count_spearman"] >= 0.50
        and result["median_profile_jsd"] <= 0.50
    )
    result["thresholds"] = {"count_pearson_min": 0.55, "count_spearman_min": 0.50, "median_profile_jsd_max": 0.50}
    dest = args.out_root / f"fold{args.fold}/internal_qc.json"
    tmp = dest.with_suffix(".json.tmp")
    with tmp.open("w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, dest)
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.require_pass and not result["internal_qc_pass"]:
        raise SystemExit(f"fold{args.fold} failed frozen internal QC thresholds")


if __name__ == "__main__":
    main()
