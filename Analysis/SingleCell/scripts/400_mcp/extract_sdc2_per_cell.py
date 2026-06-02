#!/usr/bin/env python
"""
extract_sdc2_per_cell.py
==============================================================================
Extracts SDC2 (and the other 4 reciprocal-feedback genes APP / JAML / CXADR /
ACTR2) expression for hepatocytes and macrophages along consensus pseudotime
from the per-cell-type h5ad subsets in
  Analysis/SingleCell/results_gpu_v2/pseudotime/{Hepatocytes,Macrophages}_subset.h5ad

Joins barcode -> sample -> preparation_method via the row-aligned atlas
metadata CSV `atlas_umap_for_fig2.csv.gz` + integrated_atlas.h5ad obs.
Filters to preparation_method ∈ {unsorted, nuclei}.

Outputs (one CSV per cell type):
  Analysis/SingleCell/results_gpu_v2/ccc/sdc2_per_cell_<celltype>.csv
    columns: barcode, sample, dataset, condition, consensus_pseudotime,
             preparation_method, SDC2, APP, JAML, CXADR, ACTR2
And a binned aggregate CSV (10 quantile bins of consensus_pseudotime per
cell type) suitable for the Fig 2C4 inset:
  Analysis/SingleCell/results_gpu_v2/ccc/sdc2_pseudotime_binned.csv

Usage:
  micromamba run -n spatial python extract_sdc2_per_cell.py
"""

from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PSEUDO_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/pseudotime"
ATLAS_H5 = BASE / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
ATLAS_CSV = BASE / "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz"
OUT_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/ccc"

GENES = ["SDC2", "APP", "JAML", "CXADR", "ACTR2"]
PREP_KEEP = {"unsorted", "nuclei"}


def build_sample_prep_map() -> pd.Series:
    """Sample -> majority preparation_method via row-aligned atlas csv + h5ad obs."""
    print("[map] reading atlas h5ad obs (sample only)…", flush=True)
    a_atlas = ad.read_h5ad(ATLAS_H5, backed="r")
    samples = a_atlas.obs["sample"].astype(str).values
    print(f"[map]   n_atlas_cells={len(samples):,}", flush=True)
    print("[map] reading atlas csv (preparation_method only)…", flush=True)
    prep = pd.read_csv(ATLAS_CSV, usecols=["preparation_method"])["preparation_method"].astype(str).values
    if len(prep) != len(samples):
        raise RuntimeError(f"row count mismatch atlas csv {len(prep)} vs h5ad {len(samples)}")
    df = pd.DataFrame({"sample": samples, "preparation_method": prep})
    # one method per sample – take mode
    sample_prep = df.groupby("sample")["preparation_method"].agg(lambda x: x.value_counts().idxmax())
    print(f"[map]   unique samples={sample_prep.size}; "
          f"prep counts:\n{sample_prep.value_counts()}", flush=True)
    return sample_prep


def extract_one(label: str, h5ad_path: Path, sample_prep: pd.Series) -> pd.DataFrame:
    print(f"\n[{label}] reading h5ad backed: {h5ad_path}", flush=True)
    a = ad.read_h5ad(h5ad_path, backed="r")
    print(f"[{label}]   shape={a.shape}", flush=True)
    var_names = pd.Index(a.var_names.astype(str))
    have = [g for g in GENES if g in var_names]
    miss = [g for g in GENES if g not in var_names]
    if miss:
        print(f"[{label}]   missing genes: {miss}", flush=True)
    if not have:
        print(f"[{label}]   no target genes — skipping.", flush=True)
        return pd.DataFrame()
    gene_idx = [var_names.get_loc(g) for g in have]
    # h5py backed mode requires sorted indices — sort then permute back
    sort_order = sorted(range(len(gene_idx)), key=lambda i: gene_idx[i])
    sorted_idx = [gene_idx[i] for i in sort_order]
    sorted_have = [have[i] for i in sort_order]
    print(f"[{label}]   slicing {len(have)} genes (cols, sorted)…", flush=True)
    X = a.X[:, sorted_idx]
    if sp.issparse(X):
        X = X.toarray()
    expr = pd.DataFrame(np.asarray(X), index=a.obs_names.astype(str),
                        columns=sorted_have)
    # restore original gene order
    expr = expr[have]
    obs = a.obs[["sample", "dataset", "condition", "consensus_pseudotime"]].copy()
    obs.index = obs.index.astype(str)
    out = obs.join(expr, how="inner")
    out["preparation_method"] = out["sample"].astype(str).map(sample_prep).fillna("unknown")
    n0 = len(out)
    out = out[out["preparation_method"].isin(PREP_KEEP)]
    print(f"[{label}]   filter prep ∈ {PREP_KEEP}: {n0} -> {len(out)}", flush=True)
    out = out.dropna(subset=["consensus_pseudotime"])
    print(f"[{label}]   non-NA pseudotime: {len(out)}", flush=True)
    out["cell_type"] = label
    out.index.name = "barcode"
    return out.reset_index()


def bin_aggregate(df: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    """Per-(cell_type, gene, bin) mean / sem / n along consensus pseudotime."""
    if df.empty:
        return df
    rows = []
    for ct, sub in df.groupby("cell_type"):
        # quantile bins on pseudotime (within cell type) so each bin has ≈ equal n
        sub = sub.copy()
        sub["pt_bin"] = pd.qcut(sub["consensus_pseudotime"],
                                q=n_bins, labels=False, duplicates="drop")
        for g in GENES:
            if g not in sub.columns:
                continue
            grp = sub.groupby("pt_bin")[g].agg(["mean", "std", "count"]).reset_index()
            grp["sem"] = grp["std"] / np.sqrt(grp["count"].clip(lower=1))
            grp["pt_mean"] = sub.groupby("pt_bin")["consensus_pseudotime"].mean().values
            grp["cell_type"] = ct
            grp["gene"] = g
            rows.append(grp[["cell_type", "gene", "pt_bin", "pt_mean", "mean", "sem", "count"]])
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample_prep = build_sample_prep_map()

    parts = []
    for label, fn in (("Hepatocytes", "Hepatocytes_subset.h5ad"),
                      ("Macrophages", "Macrophages_subset.h5ad")):
        df = extract_one(label, PSEUDO_DIR / fn, sample_prep)
        if not df.empty:
            out_path = OUT_DIR / f"sdc2_per_cell_{label}.csv.gz"
            df.to_csv(out_path, index=False, compression="gzip")
            print(f"[{label}]   wrote {out_path}  ({len(df):,} cells)", flush=True)
            parts.append(df)
    if not parts:
        print("[fatal] no output produced", flush=True)
        return 1
    full = pd.concat(parts, ignore_index=True)
    binned = bin_aggregate(full, n_bins=10)
    binned_path = OUT_DIR / "sdc2_pseudotime_binned.csv"
    binned.to_csv(binned_path, index=False)
    print(f"[done] wrote {binned_path}  ({len(binned)} rows)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
