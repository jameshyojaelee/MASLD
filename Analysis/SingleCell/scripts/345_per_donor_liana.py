#!/usr/bin/env python
"""
345_per_donor_liana.py

Per-donor LIANA rank_aggregate on the 9-lineage subset of the MASLD scRNA
atlas. Runs one donor (or a batch of donors) per SLURM array task. Outputs
a tidy parquet per donor.

9 lineages (CellTypist labels):
  Hepatocytes, Macrophages, Fibroblasts, Endothelial cells, Cholangiocytes,
  T cells, B cells, Resident NK, Plasma cells

Per-donor gating:
  - drop lineages with < MIN_CELLS_PER_LINEAGE (default 30)
  - skip donor if < 2 surviving lineages
  - skip donor if total cells < MIN_TOTAL_CELLS (default 200)

Output columns (per parquet):
  sample, source, target, ligand_complex, receptor_complex,
  magnitude_rank, specificity_rank, n_source_cells, n_target_cells

Usage:
  python 345_per_donor_liana.py --donor GSM4041150
  python 345_per_donor_liana.py --donor-list /path/to/donor_batch.txt
  python 345_per_donor_liana.py --array-task-id 3 --array-batch-size 8

SLURM array launch (see run_stage_ccc.sbatch).
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

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

H5AD_PATH  = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv"
OUT_DIR    = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/per_donor_lr"
COUNTS_TSV = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/lineage_cell_counts.tsv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

KEEP_LINEAGES = {
    "Hepatocytes", "Macrophages", "Fibroblasts", "Endothelial cells",
    "Cholangiocytes", "T cells", "B cells", "Resident NK", "Plasma cells",
}
MIN_CELLS_PER_LINEAGE = 30
MIN_TOTAL_CELLS       = 200
EXPR_PROP             = 0.10
N_PERMS               = 100  # LIANA permutations; 100 fast, 1000 slow

def run_donor(adata_full: ad.AnnData, sample: str) -> dict:
    """Run LIANA for one donor; return result summary."""
    t0 = time.time()
    sub = adata_full[adata_full.obs["sample"] == sample].to_memory()
    sub = sub[sub.obs["cell_type"].isin(KEEP_LINEAGES)].copy()

    if sub.n_obs < MIN_TOTAL_CELLS:
        return {"sample": sample, "status": "skipped_total_cells",
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    # Cell-type counts; drop lineages with too few cells (set their cells aside)
    counts = sub.obs["cell_type"].value_counts()
    keep_cts = counts[counts >= MIN_CELLS_PER_LINEAGE].index.tolist()
    if len(keep_cts) < 2:
        return {"sample": sample, "status": "skipped_lineages",
                "n_cells": int(sub.n_obs),
                "n_lineages": len(keep_cts), "elapsed_s": time.time() - t0}
    sub = sub[sub.obs["cell_type"].isin(keep_cts)].copy()
    sub.obs["cell_type"] = sub.obs["cell_type"].astype(str)

    # Run LIANA on this donor. Use logged + raw counts; depends on what's in .X
    import liana as li
    try:
        li.mt.rank_aggregate(
            sub,
            groupby="cell_type",
            resource_name="consensus",
            expr_prop=EXPR_PROP,
            min_cells=5,
            n_perms=N_PERMS,
            use_raw=False,   # .X is normalized log1p (CellTypist-annotated h5ad)
            verbose=False,
            seed=42,
        )
    except Exception as e:
        return {"sample": sample, "status": f"liana_failed:{type(e).__name__}:{e}",
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    res = sub.uns.get("liana_res")
    if res is None or len(res) == 0:
        return {"sample": sample, "status": "empty_result",
                "n_cells": int(sub.n_obs), "elapsed_s": time.time() - t0}

    # Annotate with donor + cell counts per cell-type
    res = res.copy()
    res["sample"] = sample
    ct_counts = counts.to_dict()
    res["n_source_cells"] = res["source"].map(ct_counts)
    res["n_target_cells"] = res["target"].map(ct_counts)
    keep_cols = ["sample", "source", "target", "ligand_complex", "receptor_complex",
                 "magnitude_rank", "specificity_rank",
                 "n_source_cells", "n_target_cells"]
    res = res[[c for c in keep_cols if c in res.columns]]

    out_parquet = OUT_DIR / f"{sample}_lr_scores.parquet"
    res.to_parquet(out_parquet, index=False)

    return {"sample": sample, "status": "ok",
            "n_cells": int(sub.n_obs),
            "n_lineages": len(keep_cts),
            "n_lr": int(len(res)),
            "out": str(out_parquet),
            "elapsed_s": round(time.time() - t0, 1)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", help="single donor (sample) ID")
    ap.add_argument("--donor-list", help="file with one donor per line")
    ap.add_argument("--array-task-id", type=int,
                    help="SLURM array task id (1-based)")
    ap.add_argument("--array-batch-size", type=int, default=8,
                    help="donors per array task")
    ap.add_argument("--write-counts", action="store_true",
                    help="also emit lineage_cell_counts.tsv (donor x lineage)")
    args = ap.parse_args()

    # Resolve donor list
    if args.donor:
        donors = [args.donor]
    elif args.donor_list:
        donors = [L.strip() for L in open(args.donor_list) if L.strip()]
    elif args.array_task_id is not None:
        full = pd.read_csv(DONOR_META, sep="\t")["sample"].tolist()
        i = args.array_task_id - 1
        donors = full[i * args.array_batch_size : (i + 1) * args.array_batch_size]
    else:
        # default: all donors
        donors = pd.read_csv(DONOR_META, sep="\t")["sample"].tolist()

    if not donors:
        print("[main] no donors to process; exiting", file=sys.stderr)
        sys.exit(0)

    print(f"[main] processing {len(donors)} donors")
    print(f"[main] h5ad: {H5AD_PATH}")
    adata_full = ad.read_h5ad(H5AD_PATH, backed="r")
    print(f"[main] full atlas shape: {adata_full.shape}")

    if args.write_counts:
        # Compute donor x lineage cell counts and write sidecar tsv
        obs = adata_full.obs[["sample", "cell_type"]].copy()
        obs = obs[obs["cell_type"].isin(KEEP_LINEAGES)]
        wide = (obs.groupby(["sample", "cell_type"]).size()
                  .unstack(fill_value=0).reset_index())
        wide.columns = [f"n_{c.replace(' ', '_')}"
                        if c != "sample" else c for c in wide.columns]
        wide.to_csv(COUNTS_TSV, sep="\t", index=False)
        print(f"[counts] wrote {COUNTS_TSV} ({wide.shape})")

    log_rows = []
    for d in donors:
        try:
            row = run_donor(adata_full, d)
        except Exception as e:
            row = {"sample": d, "status": f"unhandled:{type(e).__name__}:{e}"}
        print(f"[donor] {d}: {row.get('status')} "
              f"({row.get('n_cells', '?')} cells, "
              f"{row.get('n_lr', '?')} LR pairs, "
              f"{row.get('elapsed_s', '?')} s)",
              flush=True)
        log_rows.append(row)
        gc.collect()

    log = pd.DataFrame(log_rows)
    log_csv = OUT_DIR / f"log_{donors[0]}_{donors[-1]}.csv"
    log.to_csv(log_csv, index=False)
    print(f"[main] log written to {log_csv}")

if __name__ == "__main__":
    main()
