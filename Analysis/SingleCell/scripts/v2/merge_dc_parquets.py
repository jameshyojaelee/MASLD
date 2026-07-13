#!/usr/bin/env python
"""Consolidate donor-collapsed per-donor LIANA parquets -> one tsv.gz.

Mirrors the inline merge in 07_chain_346_349_v2.R but for the _dc outputs, so
the run-level all_donor_lr_scores_v2.tsv.gz is left untouched.
"""
import os
from pathlib import Path
import pandas as pd

BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
STAGE = BASE / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
PER_DONOR = STAGE / "per_donor_lr_v2_dc"
OUT = STAGE / "all_donor_lr_scores_v2_dc.tsv.gz"

parquets = sorted(PER_DONOR.glob("*_lr_scores.parquet"))
if not parquets:
    raise SystemExit(f"No parquets in {PER_DONOR}")
print(f"[merge] {len(parquets)} donor parquets")
dfs = [pd.read_parquet(p) for p in parquets]
merged = pd.concat(dfs, ignore_index=True)
print(f"[merge] shape={merged.shape}; unique donors={merged['sample'].nunique()}")
merged.to_csv(OUT, sep="\t", index=False, compression="gzip")
print(f"[merge] wrote {OUT}")
