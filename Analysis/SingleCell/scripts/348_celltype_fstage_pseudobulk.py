#!/usr/bin/env python3
"""Per-(donor, cell_type) pseudobulk counts for the F-stage cell-type
mobilization analysis (fig2B1 ground-truth replacement).

Inputs:
  - .atlas_patched.h5ad   (full integrated atlas, raw counts in raw/X)
  - donor_metadata_extended.tsv  (carries F_stage_augmented, F_stage_documented)

Outputs (under results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk/):
  - pseudobulk_<celltype_safe>.tsv.gz  (genes x donors, raw count sums)
  - sample_celltype_metadata.tsv       (sample, cell_type, n_cells, F_stage_*, dataset)
  - genes.tsv                          (gene symbols, in row order of pseudobulks)
"""

import os
import re
import gc
import sys
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
ATLAS = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime/.atlas_patched.h5ad")
DONOR_META = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv")
OUT_DIR = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk")
os.makedirs(OUT_DIR, exist_ok=True)

MIN_CELLS_PER_GROUP = 30   # drop (donor, cell_type) cells with fewer cells


def safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def read_categorical(f, key):
    obj = f[f"obs/{key}"]
    cats = obj["categories"][:]
    if cats.dtype.kind in ("S", "O"):
        cats = np.array([c.decode("utf-8") if isinstance(c, bytes) else str(c) for c in cats])
    return cats[obj["codes"][:]]


def read_var_index(f, key="var"):
    idx = f[f"{key}/_index"][:]
    if idx.dtype.kind in ("S", "O"):
        idx = np.array([s.decode("utf-8") if isinstance(s, bytes) else str(s) for s in idx])
    return idx


def main():
    print(f"[1] Reading donor metadata: {DONOR_META}")
    donor_meta = pd.read_csv(DONOR_META, sep="\t")
    # Backwards-compat: protocol-contamination flag may be absent if Phase 1
    # atlas refresh hasn't completed yet. Treat absence as "include all".
    if "exclude_stage_analysis" not in donor_meta.columns:
        donor_meta["exclude_stage_analysis"] = False
    if "F_stage_augmented_clean" not in donor_meta.columns:
        donor_meta["F_stage_augmented_clean"] = donor_meta["F_stage_augmented"]
    # Drop protocol-contaminated donors (GSE136103, Liver_Atlas).
    n_before = len(donor_meta)
    donor_meta = donor_meta[~donor_meta["exclude_stage_analysis"].fillna(False).astype(bool)].copy()
    print(f"    excluded {n_before - len(donor_meta)} donors flagged "
          f"exclude_stage_analysis; {len(donor_meta)} donors remain")
    print(f"    donors: {len(donor_meta)}; F_stage_augmented_clean coverage: "
          f"{donor_meta['F_stage_augmented_clean'].notna().sum()}")

    print(f"[2] Reading atlas obs + raw counts from {ATLAS}")
    with h5py.File(ATLAS, "r") as f:
        cell_type = read_categorical(f, "cell_type")
        sample    = read_categorical(f, "sample")
        dataset   = read_categorical(f, "dataset")
        gene_names = read_var_index(f, "raw/var")  # raw.var matches main var
        n_cells = len(sample)
        print(f"    n_cells: {n_cells:,}, n_genes: {len(gene_names):,}")

        print("[3] Loading raw.X (sparse CSR) — this is the heavy step")
        rg = f["raw/X"]
        data = rg["data"][:]
        indices = rg["indices"][:]
        indptr  = rg["indptr"][:]
        shape = tuple(rg.attrs["shape"])
        print(f"    raw.X nnz: {len(data):,}  shape: {shape}  dtype: {data.dtype}")

    X = sparse.csr_matrix((data, indices, indptr), shape=shape)
    del data, indices, indptr
    gc.collect()

    # ----------------------------------------------------------------------
    # Restrict to donors with F_stage_augmented_clean (drop NA, contaminated,
    # and otherwise filtered donors). _clean is NA for GSE136103 + Liver_Atlas.
    # ----------------------------------------------------------------------
    keep_donors = donor_meta.dropna(subset=["F_stage_augmented_clean"]).copy()
    keep_donors["F_stage_augmented_clean"] = keep_donors["F_stage_augmented_clean"].astype(int)
    keep_set = set(keep_donors["sample"].astype(str))
    print(f"[4] Donors with F_stage_augmented_clean: {len(keep_set)}")

    cell_donor_keep = np.isin(sample, list(keep_set))
    print(f"    cells from these donors: {cell_donor_keep.sum():,} "
          f"({100 * cell_donor_keep.mean():.1f}%)")

    # ----------------------------------------------------------------------
    # Pseudobulk per (sample, cell_type) — iterate over cell types
    # ----------------------------------------------------------------------
    cell_type = np.asarray(cell_type)
    sample = np.asarray(sample)
    dataset = np.asarray(dataset)

    ct_levels = np.unique(cell_type)
    print(f"[5] Cell types ({len(ct_levels)}): {list(ct_levels)}")

    # Master metadata accumulator
    meta_rows = []

    # For donor-row order in output matrices: stable, sorted by sample id
    donor_order = sorted(keep_set)
    donor_idx_lookup = {s: i for i, s in enumerate(donor_order)}

    pd.DataFrame({"gene": gene_names}).to_csv(
        os.path.join(OUT_DIR, "genes.tsv"), sep="\t", index=False)

    for ct in ct_levels:
        print(f"\n[CT] {ct}")
        ct_mask = (cell_type == ct) & cell_donor_keep
        n_ct_cells = int(ct_mask.sum())
        if n_ct_cells == 0:
            print("    no cells, skipping")
            continue

        ct_sample = sample[ct_mask]

        # n_cells per donor in this cell type
        cnt = pd.Series(ct_sample).value_counts()
        valid_donors = cnt[cnt >= MIN_CELLS_PER_GROUP].index.tolist()
        if not valid_donors:
            print(f"    no donor with >= {MIN_CELLS_PER_GROUP} cells, skipping")
            continue

        donor_to_row = {d: i for i, d in enumerate(valid_donors)}
        keep_cell = np.isin(ct_sample, valid_donors)
        ct_idx_global = np.where(ct_mask)[0][keep_cell]
        ct_sample_kept = ct_sample[keep_cell]
        n_kept = len(ct_idx_global)
        print(f"    cells (kept >= {MIN_CELLS_PER_GROUP}/donor): {n_kept:,} across {len(valid_donors)} donors")

        # Indicator matrix: (n_donors x n_kept_cells), CSR for fast multiply
        rows = np.fromiter((donor_to_row[s] for s in ct_sample_kept), dtype=np.int32, count=n_kept)
        cols = np.arange(n_kept, dtype=np.int32)
        data = np.ones(n_kept, dtype=np.float32)
        ind = sparse.csr_matrix((data, (rows, cols)), shape=(len(valid_donors), n_kept))

        # Pseudobulk: sum cells per donor → (n_donors x n_genes) dense int
        X_ct = X[ct_idx_global]
        pb = ind @ X_ct                        # sparse (n_donors x n_genes)
        pb_dense = np.asarray(pb.todense(), dtype=np.int32)
        print(f"    pseudobulk shape: {pb_dense.shape}, total counts: {pb_dense.sum():,}")

        # Write matrix
        df = pd.DataFrame(pb_dense.T, index=gene_names, columns=valid_donors)
        df.index.name = "gene"
        out_path = os.path.join(OUT_DIR, f"pseudobulk_{safe(ct)}.tsv.gz")
        df.to_csv(out_path, sep="\t", compression="gzip")
        print(f"    wrote {out_path} ({df.shape[0]} genes x {df.shape[1]} donors)")

        # Metadata rows
        for d in valid_donors:
            meta_rows.append({
                "sample": d,
                "cell_type": ct,
                "n_cells": int(cnt[d]),
            })

        del ind, X_ct, pb, pb_dense, df
        gc.collect()

    print("\n[6] Writing sample x celltype metadata")
    meta = pd.DataFrame(meta_rows)
    meta = meta.merge(donor_meta[["sample", "dataset", "F_stage_documented",
                                  "F_stage_inferred", "F_stage_augmented",
                                  "F_stage_augmented_clean",
                                  "exclude_stage_analysis",
                                  "F_stage_source", "disease_stage_coarse",
                                  "diagnosis_harmonized", "sex", "age", "nas_score"]],
                      on="sample", how="left")
    meta_path = os.path.join(OUT_DIR, "sample_celltype_metadata.tsv")
    meta.to_csv(meta_path, sep="\t", index=False)
    print(f"    wrote {meta_path} ({len(meta)} rows)")
    print("\nDone.")


if __name__ == "__main__":
    main()
