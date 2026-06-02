#!/usr/bin/env python
"""
06_per_donor_liana_v2.py

Phase 0.5 v2 equivalent of Script 345.

Per-donor LIANA rank_aggregate on the v2 atlas (all cell types, not just
hepatocytes). LIANA needs hepatocyte + non-hepatocyte cells in the same
AnnData object so source/target lineage pairs are tested. We read Agent 2's
full v2 atlas:
  Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad

A pure hepatocyte atlas (from Agent 3) is NOT sufficient for CCC.

Outputs:
  Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2/{sample}_lr_scores.parquet
  + a lineage_cell_counts_v2.tsv sidecar (donor x lineage) written once
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

# v2 inputs
H5AD_V2_FULL = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/scalesc_human_annotated_celltypist_v2.h5ad"
)
H5AD_V2_HEP = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/atlas/hepatocyte_atlas_v2_annotated.h5ad"
)

# v1 fallback (for donor roster only; never used for cells)
DONOR_META_V1 = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
)
DONOR_META_V2 = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv"
)

V2_STAGE_DIR = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
OUT_DIR = V2_STAGE_DIR / "per_donor_lr_v2"
COUNTS_TSV = V2_STAGE_DIR / "lineage_cell_counts_v2.tsv"

KEEP_LINEAGES = {
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial cells",
    "Cholangiocytes",
    "T cells",
    "B cells",
    "Resident NK",
    "Plasma cells",
}
MIN_CELLS_PER_LINEAGE = 30
MIN_TOTAL_CELLS = 200
EXPR_PROP = 0.10
N_PERMS = 100


def resolve_atlas() -> Path:
    """Pick the full v2 atlas; halt if absent."""
    if H5AD_V2_FULL.exists():
        return H5AD_V2_FULL
    sys.exit(
        f"FATAL: v2 full atlas not found at {H5AD_V2_FULL}. "
        "LIANA requires all cell types; halting."
    )


def get_donor_roster() -> pd.DataFrame:
    """Prefer v2 donor metadata; fall back to v1 if not yet built."""
    if DONOR_META_V2.exists():
        roster = pd.read_csv(DONOR_META_V2, sep="\t")
        print(f"[roster] using v2 donor metadata: {len(roster)} donors")
    else:
        roster = pd.read_csv(DONOR_META_V1, sep="\t")
        print(f"[roster] v2 missing; falling back to v1: {len(roster)} donors")
    return roster


def run_donor(adata_full: ad.AnnData, sample: str) -> dict:
    t0 = time.time()
    sub = adata_full[adata_full.obs["sample"] == sample].to_memory()
    if "cell_type" not in sub.obs.columns:
        return {
            "sample": sample,
            "status": "missing_cell_type_col",
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    sub = sub[sub.obs["cell_type"].isin(KEEP_LINEAGES)].copy()
    if sub.n_obs < MIN_TOTAL_CELLS:
        return {
            "sample": sample,
            "status": "skipped_total_cells",
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    counts = sub.obs["cell_type"].value_counts()
    keep_cts = counts[counts >= MIN_CELLS_PER_LINEAGE].index.tolist()
    if len(keep_cts) < 2:
        return {
            "sample": sample,
            "status": "skipped_lineages",
            "n_cells": int(sub.n_obs),
            "n_lineages": len(keep_cts),
            "elapsed_s": time.time() - t0,
        }
    sub = sub[sub.obs["cell_type"].isin(keep_cts)].copy()
    sub.obs["cell_type"] = sub.obs["cell_type"].astype(str)

    import liana as li

    try:
        li.mt.rank_aggregate(
            sub,
            groupby="cell_type",
            resource_name="consensus",
            expr_prop=EXPR_PROP,
            min_cells=5,
            n_perms=N_PERMS,
            use_raw=False,
            verbose=False,
            seed=42,
        )
    except Exception as e:
        return {
            "sample": sample,
            "status": f"liana_failed:{type(e).__name__}:{e}",
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }

    res = sub.uns.get("liana_res")
    if res is None or len(res) == 0:
        return {
            "sample": sample,
            "status": "empty_result",
            "n_cells": int(sub.n_obs),
            "elapsed_s": time.time() - t0,
        }
    res = res.copy()
    res["sample"] = sample
    ct_counts = counts.to_dict()
    res["n_source_cells"] = res["source"].map(ct_counts)
    res["n_target_cells"] = res["target"].map(ct_counts)
    keep_cols = [
        "sample",
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "magnitude_rank",
        "specificity_rank",
        "n_source_cells",
        "n_target_cells",
    ]
    res = res[[c for c in keep_cols if c in res.columns]]
    out_parquet = OUT_DIR / f"{sample}_lr_scores.parquet"
    res.to_parquet(out_parquet, index=False)
    return {
        "sample": sample,
        "status": "ok",
        "n_cells": int(sub.n_obs),
        "n_lineages": len(keep_cts),
        "n_lr": int(len(res)),
        "out": str(out_parquet),
        "elapsed_s": round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor")
    ap.add_argument("--donor-list")
    ap.add_argument("--array-task-id", type=int)
    ap.add_argument("--array-batch-size", type=int, default=8)
    ap.add_argument("--write-counts", action="store_true")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    V2_STAGE_DIR.mkdir(parents=True, exist_ok=True)

    h5ad_path = resolve_atlas()
    print(f"[atlas] {h5ad_path}")
    adata_full = ad.read_h5ad(h5ad_path, backed="r")
    print(
        f"[atlas] shape={adata_full.shape}; "
        f"obs cols include cell_type={'cell_type' in adata_full.obs.columns}"
    )

    roster = get_donor_roster()
    all_donors = roster["sample"].astype(str).tolist()

    if args.donor:
        donors = [args.donor]
    elif args.donor_list:
        donors = [L.strip() for L in open(args.donor_list) if L.strip()]
    elif args.array_task_id is not None:
        i = args.array_task_id - 1
        donors = all_donors[
            i * args.array_batch_size : (i + 1) * args.array_batch_size
        ]
    else:
        donors = all_donors

    if not donors:
        print("[main] no donors to process; exiting", file=sys.stderr)
        sys.exit(0)
    print(f"[main] processing {len(donors)} donors")

    if args.write_counts and "cell_type" in adata_full.obs.columns:
        obs = adata_full.obs[["sample", "cell_type"]].copy()
        obs = obs[obs["cell_type"].isin(KEEP_LINEAGES)]
        wide = (
            obs.groupby(["sample", "cell_type"], observed=True)
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )
        wide.columns = [
            f"n_{c.replace(' ', '_')}" if c != "sample" else c
            for c in wide.columns
        ]
        wide.to_csv(COUNTS_TSV, sep="\t", index=False)
        print(f"[counts] wrote {COUNTS_TSV} ({wide.shape})")

    log_rows = []
    for d in donors:
        try:
            row = run_donor(adata_full, d)
        except Exception as e:
            row = {"sample": d, "status": f"unhandled:{type(e).__name__}:{e}"}
        print(
            f"[donor] {d}: {row.get('status')} "
            f"({row.get('n_cells', '?')} cells, "
            f"{row.get('n_lr', '?')} LR, {row.get('elapsed_s', '?')} s)",
            flush=True,
        )
        log_rows.append(row)
        gc.collect()

    log = pd.DataFrame(log_rows)
    log_csv = OUT_DIR / f"log_{donors[0]}_{donors[-1]}.csv"
    log.to_csv(log_csv, index=False)
    print(f"[main] log written to {log_csv}")


if __name__ == "__main__":
    main()
