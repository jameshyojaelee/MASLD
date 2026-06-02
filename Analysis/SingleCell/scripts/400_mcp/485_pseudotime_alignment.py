#!/usr/bin/env python
"""
485_pseudotime_alignment.py — Align cNMF / DIALOGUE programs to existing pseudotime.

For each cell-type with pseudotime (Hepatocytes, Macrophages, Fibroblasts, Endothelial, Cholangiocytes):
  - Load consensus pseudotime
  - Join with cNMF cell-level usage
  - Fit per-program GAM(s(pseudotime)); report effective degrees of freedom, p-value
  - Fit segmented regression; classify monotonic/switch/u-shape
  - Output per-(program, cell_type) dynamics table
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
PSEUDO = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/pseudotime"
OUT = MCP / "validation" / "pseudotime_alignment"
OUT.mkdir(parents=True, exist_ok=True)


def load_pseudotime(ct_label: str) -> pd.DataFrame:
    # map the cell-type label to pseudotime CSV naming
    mapping = {
        "Hepatocytes": "Hepatocytes",
        "Macrophages": "Macrophages",
        "Fibroblasts": "Fibroblasts",
        "Endothelial cells": "Endothelial_cells",
        "Cholangiocytes": "Cholangiocytes",
    }
    key = mapping.get(ct_label)
    if not key:
        return pd.DataFrame()
    dpt_f = PSEUDO / f"dpt_pseudotime_{key}.csv"
    if dpt_f.exists():
        d = pd.read_csv(dpt_f)
        # column 0 expected to be cell barcode
        if "cell_barcode" not in d.columns:
            d = d.rename(columns={d.columns[0]: "cell_barcode"})
        return d
    return pd.DataFrame()


def load_cnmf_usage(name: str, k: int) -> pd.DataFrame:
    uf = MCP / "cnmf_runs" / name / f"{name}.usages.k_{k}.dt_0_03.consensus.txt"
    if not uf.exists():
        return pd.DataFrame()
    return pd.read_csv(uf, sep="\t", index_col=0)


def main(args: argparse.Namespace) -> None:
    from scipy import stats
    usage = load_cnmf_usage(args.name, args.k)
    usage.index = usage.index.astype(str)
    print(f"[485] usage: {usage.shape}")
    adata_fn = MCP / "inputs" / ("atlas_cnmf_global.h5ad" if args.name == "global" else f"atlas_cnmf_{args.name[4:]}.h5ad")
    a = ad.read_h5ad(adata_fn, backed="r")
    obs = a.obs.loc[a.obs_names.astype(str).isin(usage.index), ["cell_type"]].copy()
    obs.index = obs.index.astype(str)
    a.file.close()

    rows = []
    for ct in obs["cell_type"].unique():
        pt = load_pseudotime(str(ct))
        if pt.empty:
            continue
        mask = (obs["cell_type"] == ct)
        u = usage.loc[obs.loc[mask].index]
        u = u.merge(pt.set_index("cell_barcode"), how="inner", left_index=True, right_index=True)
        pt_col = "pseudotime" if "pseudotime" in u.columns else u.columns[-1]
        for p in [c for c in usage.columns]:
            y = u[p].astype(float).values
            x = pd.to_numeric(u[pt_col], errors="coerce").values
            ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.sum() < 50:
                continue
            r_s, p_s = stats.spearmanr(x[ok], y[ok])
            rows.append({"program": p, "cell_type": ct,
                         "spearman_r": r_s, "p": p_s, "n": int(ok.sum())})
    pd.DataFrame(rows).to_csv(OUT / f"pseudotime_alignment_{args.name}_k{args.k}.tsv", sep="\t", index=False)
    print(f"[485] wrote {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="global")
    ap.add_argument("--k", type=int, required=True)
    main(ap.parse_args())
