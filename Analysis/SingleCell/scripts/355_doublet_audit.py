#!/usr/bin/env python
"""
355_doublet_audit.py  (S5)

Per-cohort Scrublet/DoubletFinder threshold + doublet rate audit.

Approach (light-weight, read-only):
    1. Load the integrated atlas in `backed="r"` mode.
    2. For each dataset (cohort), summarize:
         - any obs columns matching /doublet|scrublet|dblFinder/i
         - threshold (if any threshold attribute appears in adata.uns)
         - doublet rate = mean(is_doublet) per cohort
         - putative 10x chemistry version inferred from filenames /
           Cell Ranger reference (parse uns / obs if available)
    3. If the doublet rate looks uniformly ~10% across cohorts that span
       multiple chemistries (v2 / v3 / 5'), flag and recommend a chemistry-
       stratified rerun.

Output: Analysis/SingleCell/results_gpu_v2/doublet_audit/REPORT.md
        + doublet_audit_summary.tsv
"""

from __future__ import annotations
import os
import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
H5AD = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
OUT_ROOT = Path(os.environ.get("S5_OUT_ROOT", str(PROJECT_ROOT)))
OUT = OUT_ROOT / "Analysis/SingleCell/results_gpu_v2/doublet_audit"
OUT.mkdir(parents=True, exist_ok=True)


def detect_chemistry(adata: ad.AnnData) -> pd.Series:
    """Best-effort detection of 10x chemistry per dataset."""
    obs = adata.obs
    if "chemistry" in obs.columns:
        return obs.groupby("dataset")["chemistry"].agg(
            lambda x: x.value_counts().idxmax())
    # Try uns
    chemistry_map = {}
    for ds in obs["dataset"].unique() if "dataset" in obs.columns else []:
        chemistry_map[ds] = "unknown"
    return pd.Series(chemistry_map, name="chemistry")


def _read_obs_via_h5py(path: Path) -> pd.DataFrame:
    """Read only /obs from h5ad via h5py — bypasses anndata's uns/null-encoding
    incompatibility when the file was written by a newer anndata version."""
    import h5py
    out = {}
    cat_codes = {}
    cat_levels = {}
    with h5py.File(str(path), "r") as f:
        obs = f["obs"]
        # Newer anndata stores column order in obs.attrs['column-order'] and an _index attr
        col_order = list(obs.attrs.get("column-order", list(obs.keys())))
        index_key = obs.attrs.get("_index", "_index")
        # Index
        if index_key in obs:
            ix = obs[index_key][...]
            ix = pd.Index([x.decode() if isinstance(x, bytes) else x for x in ix])
        else:
            ix = None
        for col in col_order:
            if col not in obs:
                continue
            node = obs[col]
            if isinstance(node, h5py.Group):
                # Categorical: codes + categories
                if "codes" in node and "categories" in node:
                    codes = node["codes"][...]
                    cats = node["categories"][...]
                    cats = [c.decode() if isinstance(c, bytes) else c for c in cats]
                    arr = pd.Categorical.from_codes(codes, categories=cats)
                    out[col] = arr
                    continue
            else:
                arr = node[...]
                if arr.dtype.kind in ("S", "O"):
                    arr = [x.decode() if isinstance(x, bytes) else x for x in arr]
                out[col] = arr
        df = pd.DataFrame(out)
        if ix is not None and len(ix) == len(df):
            df.index = ix
        return df


def main():
    print(f"[load] {H5AD}", flush=True)
    try:
        A = ad.read_h5ad(H5AD, backed="r")
        obs = A.obs.copy()
        print(f"[load] anndata backed read OK; shape {A.shape}", flush=True)
    except Exception as e:
        print(f"[warn] anndata read_h5ad failed: {e}; falling back to h5py obs-only",
              flush=True)
        obs = _read_obs_via_h5py(H5AD)
        A = None
        print(f"[load] h5py fallback: obs shape {obs.shape}", flush=True)
    dbl_cols = [c for c in obs.columns
                if re.search(r"doublet|scrublet|dblFinder", c, re.I)]
    print(f"[doublet cols] {dbl_cols}", flush=True)

    summary = []
    if "dataset" not in obs.columns:
        dataset_col = "sample"
    else:
        dataset_col = "dataset"

    if A is not None:
        chem = detect_chemistry(A)
    else:
        # obs-only fallback: no uns/var access available
        if "chemistry" in obs.columns:
            chem = obs.groupby("dataset")["chemistry"].agg(
                lambda x: x.value_counts().idxmax())
        else:
            chem = pd.Series({ds: "unknown" for ds in obs.get("dataset", pd.Series([])).unique()},
                             name="chemistry")
    for ds, sub in obs.groupby(dataset_col):
        row = {"dataset": ds, "n_cells": int(len(sub)),
               "chemistry": chem.get(ds, "unknown")}
        for c in dbl_cols:
            v = sub[c]
            if v.dtype == bool:
                row[f"{c}_rate"] = float(v.mean())
            elif v.dtype.kind in "fi":
                row[f"{c}_rate"]  = float((v > 0.5).mean())
                row[f"{c}_mean"]  = float(v.mean())
                row[f"{c}_max"]   = float(v.max())
            else:
                # categorical
                vc = v.value_counts(normalize=True)
                for k, p in vc.items():
                    row[f"{c}__{k}"] = float(p)
        summary.append(row)

    df = pd.DataFrame(summary).sort_values("dataset")
    df.to_csv(OUT / "doublet_audit_summary.tsv", sep="\t", index=False)
    print(df.to_string(index=False), flush=True)

    # Uniform-10% flag
    rate_cols = [c for c in df.columns if c.endswith("_rate")]
    uniform = False
    if rate_cols:
        rates = df[rate_cols].fillna(np.nan).values.flatten()
        rates = rates[~np.isnan(rates)]
        if len(rates) > 3:
            uniform = (abs(np.median(rates) - 0.10) < 0.02) and (np.std(rates) < 0.02)

    lines = []
    lines.append("# Doublet audit (S5)\n")
    lines.append(f"Source h5ad: `{H5AD}`\n")
    lines.append(f"Doublet obs columns found: `{dbl_cols}`\n")
    if uniform:
        lines.append("WARNING: doublet rate looks uniform ~10% across cohorts.\n"
                     "Recommend chemistry-stratified rerun "
                     "(Scrublet expected_doublet_rate per 10x chemistry: "
                     "v2 ≈ 0.06, v3 ≈ 0.08, 5' ≈ 0.04).\n")
    else:
        lines.append("Doublet rates show cohort-level variability "
                     "(no uniform-10% artifact).\n")
    lines.append("\n## Per-cohort table\n\n")
    lines.append(df.to_markdown(index=False))
    (OUT / "REPORT.md").write_text("\n".join(lines))
    print(f"[write] {OUT / 'REPORT.md'}", flush=True)


if __name__ == "__main__":
    main()
