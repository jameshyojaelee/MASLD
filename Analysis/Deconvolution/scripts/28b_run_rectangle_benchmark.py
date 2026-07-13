#!/usr/bin/env python
"""
28b_run_rectangle_benchmark.py  --  Pseudobulk benchmark, step 2/4 (Rectangle).

For each fold: build a Rectangle reference AnnData from OUT-of-fold donors' cells
(subset of the 100k atlas), deconvolve the IN-fold donors' pseudobulk (CPM) with
rectanglepy.rectangle(), reindex to the 16 canonical cell types, and collect.

Rectangle consumes CPM (per-sample library-normalized ~1e6); the atlas is sc-derived
so there is no gene-length bias -> CPM (not TPM) is the correct unit.

Output: $BENCH_OUT_DIR/rectangle_estimates.tsv  (donors x 16 cell types)

Runs in the `rectangle` env python.
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad
import rectanglepy as rp

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5AD = os.path.join(ROOT, "Analysis/Deconvolution/reference/reference_human_scalesc_100k.h5ad")
OUT_DIR = os.environ.get("BENCH_OUT_DIR",
                         os.path.join(ROOT, "Analysis/Deconvolution/rectangle_comparison/pseudobulk"))
# optimize_cutoffs=True is Rectangle's recommended default (best accuracy); toggle off for speed.
RECT_OPTIMIZE = os.environ.get("BENCH_RECT_OPTIMIZE", "1") == "1"
NCPU = int(os.environ.get("SLURM_CPUS_PER_TASK", os.environ.get("BENCH_NCPU", "8")))


def read_categorical(g):
    cats = g["categories"][:]
    codes = g["codes"][:]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats], dtype=object)
    out = np.empty(codes.shape[0], dtype=object)
    mask = codes >= 0
    out[mask] = cats[codes[mask]]
    out[~mask] = None
    return out


def read_obs_col(f, name):
    node = f["obs"][name]
    if isinstance(node, h5py.Group):
        return read_categorical(node)
    vals = node[:]
    return np.array([v.decode() if isinstance(v, bytes) else v for v in vals], dtype=object)


def load_atlas():
    """Build AnnData (X raw counts, obs cell_type+sample, var symbols) via h5py.
    Avoids anndata.read_h5ad, which trips on this file's /uns/log1p null encoding."""
    f = h5py.File(H5AD, "r")
    Xg = f["X"]
    shape = tuple(int(x) for x in Xg.attrs["shape"])
    X = sp.csr_matrix((Xg["data"][:].astype(np.float32),
                       Xg["indices"][:].astype(np.int64),
                       Xg["indptr"][:].astype(np.int64)), shape=shape)
    genes = np.array([g.decode() if isinstance(g, bytes) else g
                      for g in f["var"]["_index"][:]], dtype=object)
    cell_type = read_obs_col(f, "cell_type")
    donor = read_obs_col(f, "sample")
    f.close()
    obs = pd.DataFrame({"cell_type": pd.Categorical(cell_type),
                        "sample": pd.Categorical(donor)})
    obs.index = [f"cell{i}" for i in range(X.shape[0])]
    var = pd.DataFrame(index=[str(g) for g in genes])
    a = ad.AnnData(X=X, obs=obs, var=var)
    return a


def main():
    canonical = [l.strip() for l in
                 open(os.path.join(OUT_DIR, "canonical_celltypes.txt")) if l.strip()]
    folds = pd.read_csv(os.path.join(OUT_DIR, "donor_folds.tsv"), sep="\t", dtype={"donor": str})
    cpm = pd.read_csv(os.path.join(OUT_DIR, "pseudobulk_cpm.tsv"), sep="\t",
                      index_col=0)                       # genes x donors
    cpm.columns = cpm.columns.astype(str)

    print(f"[28b] loading atlas AnnData ...", flush=True)
    atlas = load_atlas()
    atlas_donor = atlas.obs["sample"].astype(str).values
    print(f"[28b] atlas {atlas.shape}; folds={sorted(folds.fold.unique())}; "
          f"optimize_cutoffs={RECT_OPTIMIZE}; n_cpus={NCPU}", flush=True)

    est_rows = []
    for fold in sorted(folds.fold.unique()):
        ref_donors = set(folds.loc[folds.fold != fold, "donor"])
        test_donors = folds.loc[folds.fold == fold, "donor"].tolist()
        ref_mask = np.isin(atlas_donor, list(ref_donors))
        ref_adata = atlas[ref_mask].copy()
        bulks = cpm[test_donors].T                       # samples x genes (Rectangle layout)
        print(f"[28b] fold {fold}: ref cells={ref_mask.sum()} "
              f"(donors={len(ref_donors)}, types={ref_adata.obs.cell_type.nunique()}); "
              f"test donors={len(test_donors)}", flush=True)

        est, _sig = rp.rectangle(ref_adata, bulks, cell_type_col="cell_type",
                                 optimize_cutoffs=RECT_OPTIMIZE, n_cpus=NCPU)
        # est: samples x (cell types + 'Unknown'); keep 16 canonical, drop Unknown
        unknown = float(est["Unknown"].mean()) if "Unknown" in est.columns else 0.0
        est = est.reindex(columns=canonical, fill_value=0.0)
        est.index = [str(x) for x in est.index]
        est_rows.append(est)
        print(f"[28b] fold {fold}: done (mean Unknown mass dropped = {unknown:.4f})", flush=True)

    out = pd.concat(est_rows, axis=0)
    out = out.reindex(index=folds.donor.tolist(), columns=canonical)
    out.index.name = "donor"
    out.to_csv(os.path.join(OUT_DIR, "rectangle_estimates.tsv"), sep="\t")
    print(f"[28b] wrote rectangle_estimates.tsv {out.shape}", flush=True)
    print("[28b] DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
