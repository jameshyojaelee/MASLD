#!/usr/bin/env python
"""
345b_merge_per_donor_parquets.py

One-shot consolidator: read every per_donor_lr/*.parquet from Script 345 and
write a single gzipped TSV that R can fread() without needing the `arrow`
package (which is not in the rnaseq env).

Output: all_donor_lr_scores.tsv.gz
"""

from __future__ import annotations
import os
import sys
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
PER_DONOR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/per_donor_lr"
OUT_TSV   = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/all_donor_lr_scores.tsv.gz"

def main():
    parquets = sorted(PER_DONOR.glob("*_lr_scores.parquet"))
    if not parquets:
        print("[345b] no per-donor parquets found", file=sys.stderr)
        sys.exit(1)
    print(f"[345b] consolidating {len(parquets)} per-donor parquets")
    dfs = [pd.read_parquet(p) for p in parquets]
    merged = pd.concat(dfs, ignore_index=True)
    print(f"[345b] merged shape: {merged.shape}")
    merged.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")
    print(f"[345b] wrote {OUT_TSV}")

if __name__ == "__main__":
    main()
