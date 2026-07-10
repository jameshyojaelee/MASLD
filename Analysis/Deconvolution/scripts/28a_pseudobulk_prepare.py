#!/usr/bin/env python
"""
28a_pseudobulk_prepare.py  --  Pseudobulk ground-truth benchmark, step 1/4.

Builds a leakage-free pseudobulk deconvolution benchmark from the working
single-cell atlas (reference_human_scalesc_100k.h5ad; .X = raw integer counts).

Per DONOR (obs 'sample') with >= N_MIN cells:
  * TRUTH proportions      = that donor's cell_type fractions over the SAME cells
  * PSEUDOBULK expression  = SUM of raw counts over that donor's cells (a gene vector)

Donors are split into N_FOLDS grouped folds (seeded). For each fold the reference is
built from OUT-of-fold donors and the IN-fold donor pseudobulks are deconvolved. Truth
and pseudobulk come from the same cells, so they are exactly consistent.

Writes to $BENCH_OUT_DIR (default rectangle_comparison/pseudobulk):
  pseudobulk_counts.tsv   genes x donors  (raw summed integer counts)   -> MuSiC / InstaPrism
  pseudobulk_cpm.tsv      genes x donors  (per-donor CPM, sum 1e6)      -> Rectangle
  truth_proportions.tsv   donors x 16 cell types (rows sum to 1)
  donor_folds.tsv         donor, fold, n_cells
  canonical_celltypes.txt 16 cell-type strings, one per line (shared order for R + py)

Runs in the `rectangle` env python.
"""
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py

# ----------------------------------------------------------------------------- config
ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5AD = os.path.join(ROOT, "Analysis/Deconvolution/reference/reference_human_scalesc_100k.h5ad")
OUT_DIR = os.environ.get("BENCH_OUT_DIR",
                         os.path.join(ROOT, "Analysis/Deconvolution/rectangle_comparison/pseudobulk"))
N_MIN = int(os.environ.get("BENCH_MIN_CELLS", "200"))
N_FOLDS = int(os.environ.get("BENCH_N_FOLDS", "5"))
MAX_DONORS = int(os.environ.get("BENCH_MAX_DONORS", "0"))   # 0 = use all eligible donors
SEED = int(os.environ.get("BENCH_SEED", "42"))

# 16 canonical cell-type strings (EXACT; shared order for all four scripts)
CANONICAL_CELLTYPES = [
    "Endothelial cells", "Hepatocytes", "Plasma cells", "T cells", "Cholangiocytes",
    "Fibroblasts", "Macrophages", "Circulating NK/NKT", "Resident NK",
    "Mono+mono derived cells", "Basophils", "B cells", "cDC1s", "cDC2s", "pDCs",
    "Neutrophils",
]


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
    if isinstance(node, h5py.Group):          # categorical
        return read_categorical(node)
    vals = node[:]
    return np.array([v.decode() if isinstance(v, bytes) else v for v in vals], dtype=object)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"[28a] reading {H5AD}", flush=True)
    f = h5py.File(H5AD, "r")

    # ---- X as CSR (cells x genes) -------------------------------------------------
    Xg = f["X"]
    shape = tuple(int(x) for x in Xg.attrs["shape"])
    n_cells, n_genes = shape
    X = sp.csr_matrix((Xg["data"][:].astype(np.float64),
                       Xg["indices"][:].astype(np.int64),
                       Xg["indptr"][:].astype(np.int64)), shape=shape)
    print(f"[28a] X: {n_cells} cells x {n_genes} genes, nnz={X.nnz}", flush=True)

    genes = np.array([g.decode() if isinstance(g, bytes) else g
                      for g in f["var"]["_index"][:]], dtype=object)
    cell_type = read_obs_col(f, "cell_type")
    donor = read_obs_col(f, "sample")
    f.close()

    # ---- donor filter (>= N_MIN cells) -------------------------------------------
    donor_ser = pd.Series(donor)
    counts_per_donor = donor_ser.value_counts()
    eligible = counts_per_donor[counts_per_donor >= N_MIN].index.tolist()
    eligible = sorted(eligible)  # deterministic order before shuffle
    if MAX_DONORS > 0:
        # dry-run subset: seeded random pick to preserve cell-type diversity
        rng = np.random.default_rng(SEED)
        eligible = sorted(rng.choice(eligible, size=min(MAX_DONORS, len(eligible)),
                                     replace=False).tolist())
    print(f"[28a] donors with >= {N_MIN} cells: {len(eligible)} "
          f"(MAX_DONORS={MAX_DONORS or 'all'})", flush=True)

    donor_to_idx = {d: i for i, d in enumerate(eligible)}
    keep_cell = np.array([d in donor_to_idx for d in donor])
    print(f"[28a] cells retained: {keep_cell.sum()} / {n_cells}", flush=True)

    Xk = X[keep_cell]
    donor_k = donor[keep_cell]
    ct_k = cell_type[keep_cell]

    # ---- pseudobulk: donors x genes = (donors x cells indicator) @ (cells x genes) -
    row = np.array([donor_to_idx[d] for d in donor_k], dtype=np.int64)
    col = np.arange(len(donor_k), dtype=np.int64)
    ind = sp.csr_matrix((np.ones(len(donor_k)), (row, col)),
                        shape=(len(eligible), Xk.shape[0]))
    pb = np.asarray((ind @ Xk).todense())                     # donors x genes, raw sums
    pb = np.rint(pb).astype(np.int64)
    pb_df = pd.DataFrame(pb.T, index=genes, columns=eligible)  # genes x donors
    pb_df.index.name = "SYMBOL"

    # ---- truth proportions: donors x 16 canonical types --------------------------
    truth = pd.crosstab(pd.Series(donor_k, name="donor"),
                        pd.Series(ct_k, name="cell_type"))
    truth = truth.reindex(index=eligible, columns=CANONICAL_CELLTYPES, fill_value=0)
    truth = truth.div(truth.sum(axis=1), axis=0)               # row-normalize to 1
    truth.index.name = "donor"

    # ---- CPM for Rectangle (per-donor library size -> 1e6) -----------------------
    lib = pb_df.sum(axis=0).replace(0, np.nan)                 # per donor total counts
    cpm_df = pb_df.div(lib, axis=1) * 1e6
    cpm_df = cpm_df.fillna(0.0)
    cpm_df.index.name = "SYMBOL"

    # ---- grouped folds (seeded shuffle of donors) --------------------------------
    rng = np.random.default_rng(SEED)
    shuffled = np.array(eligible, dtype=object)
    perm = rng.permutation(len(shuffled))
    shuffled = shuffled[perm]
    fold_of = {}
    for fold, chunk in enumerate(np.array_split(shuffled, N_FOLDS)):
        for d in chunk:
            fold_of[d] = fold
    folds_df = pd.DataFrame({
        "donor": eligible,
        "fold": [fold_of[d] for d in eligible],
        "n_cells": [int(counts_per_donor[d]) for d in eligible],
    })

    # ---- write -------------------------------------------------------------------
    pb_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_counts.tsv"), sep="\t")
    cpm_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_cpm.tsv"), sep="\t")
    truth.to_csv(os.path.join(OUT_DIR, "truth_proportions.tsv"), sep="\t")
    folds_df.to_csv(os.path.join(OUT_DIR, "donor_folds.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT_DIR, "canonical_celltypes.txt"), "w") as fh:
        fh.write("\n".join(CANONICAL_CELLTYPES) + "\n")

    print(f"[28a] wrote pseudobulk_counts.tsv  {pb_df.shape} (genes x donors)")
    print(f"[28a] wrote pseudobulk_cpm.tsv     {cpm_df.shape}")
    print(f"[28a] wrote truth_proportions.tsv  {truth.shape}")
    print(f"[28a] wrote donor_folds.tsv        folds={sorted(folds_df.fold.unique())} "
          f"sizes={folds_df.fold.value_counts().sort_index().tolist()}")
    print(f"[28a] output dir: {OUT_DIR}")
    print("[28a] DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
