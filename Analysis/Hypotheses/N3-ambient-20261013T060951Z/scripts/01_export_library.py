#!/usr/bin/env python
"""N3 ambient check, step 1: export one library for decontX.

Reads the atlas raw UMI counts of one library (all cell types) and the
empty-droplet profile from that library's Cell Ranger raw matrix. Writes
node-local temporary files for 02_decontx_library.R.

Usage: python -I 01_export_library.py <REPO> <LIBRARY> <DATASET> <TMPDIR>
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse

REPO, LIB, DATASET, TMP = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
ATLAS = REPO / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
RAW = REPO / f"data/{DATASET}/cellranger/{LIB}/outs/raw_feature_bc_matrix.h5"
FILT = REPO / f"data/{DATASET}/cellranger/{LIB}/outs/filtered_feature_bc_matrix.h5"
EMPTY_MIN, EMPTY_MAX = 1, 100


def decode(arr):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])


def main():
    out = TMP / LIB
    out.mkdir(parents=True, exist_ok=True)
    with h5py.File(ATLAS, "r") as f:
        lib_cats = decode(f["obs/sample/categories"][:])
        codes = f["obs/sample/codes"][:]
        code = int(np.where(lib_cats == LIB)[0][0])
        rows = np.flatnonzero(codes == code)
        s, e = int(rows[0]), int(rows[-1]) + 1
        assert e - s == len(rows), "library rows are not contiguous in the atlas"
        indptr = f["raw/X/indptr"][s:e + 1].astype(np.int64)
        a, b = int(indptr[0]), int(indptr[-1])
        data = f["raw/X/data"][a:b]
        idx = f["raw/X/indices"][a:b]
        gene_ids = decode(f["raw/var/gene_ids"][:])
        ct_cats = decode(f["obs/cell_type/categories"][:])
        ct = ct_cats[f["obs/cell_type/codes"][s:e]]
        names = decode(f["obs/_index"][s:e])
    n_genes = len(gene_ids)
    assert np.all(data == np.round(data)), "atlas raw/X is not integer-valued"
    M = sparse.csr_matrix((data.astype(np.int32), idx, indptr - a), shape=(e - s, n_genes))
    M.sort_indices()
    # cells x genes CSR has the same arrays as genes x cells CSC
    M.indices.astype("<i4").tofile(out / "indices.bin")
    M.indptr.astype("<i4").tofile(out / "indptr.bin")
    M.data.astype("<i4").tofile(out / "data.bin")
    prefix = LIB + "_"
    assert all(n.startswith(prefix) for n in names)
    atlas_bc = np.array([n[len(prefix):] for n in names])
    with open(out / "cells.tsv", "w") as fh:
        fh.write("barcode\tcell_type\n")
        for bc, c in zip(atlas_bc, ct):
            fh.write(f"{bc}\t{c}\n")

    # ------------------------------------------------ empty-droplet profile
    with h5py.File(FILT, "r") as h:
        called = set(decode(h["matrix/barcodes"][:]))
    with h5py.File(RAW, "r") as h:
        r_data = h["matrix/data"][:].astype(np.int32)
        r_idx = h["matrix/indices"][:]
        r_ptr = h["matrix/indptr"][:].astype(np.int64)
        r_bc = decode(h["matrix/barcodes"][:])
        r_ids = decode(h["matrix/features/id"][:])
    pos = {g: i for i, g in enumerate(r_ids)}
    feat_to_atlas = np.full(len(r_ids), -1, dtype=np.int64)
    for j, g in enumerate(gene_ids):
        feat_to_atlas[pos[g]] = j
    cs = np.concatenate([[0], np.cumsum(r_data, dtype=np.int64)])
    tot = cs[r_ptr[1:]] - cs[r_ptr[:-1]]
    is_called = np.fromiter((b in called for b in r_bc), dtype=bool, count=len(r_bc))
    atlas_set = set(atlas_bc)
    is_atlas = np.fromiter((b in atlas_set for b in r_bc), dtype=bool, count=len(r_bc))
    empty = (tot >= EMPTY_MIN) & (tot <= EMPTY_MAX) & ~is_called & ~is_atlas
    nnz_per_bc = np.diff(r_ptr)
    mask = np.repeat(empty, nnz_per_bc)
    feat_sum = np.bincount(r_idx[mask], weights=r_data[mask], minlength=len(r_ids))
    keep = feat_to_atlas >= 0
    bg = np.zeros(n_genes)
    bg[feat_to_atlas[keep]] = feat_sum[keep]
    bg.astype("<f8").tofile(out / "background.bin")

    # fidelity: atlas cell totals against the raw matrix restricted to atlas genes
    cs2 = np.concatenate([[0], np.cumsum(r_data * keep[r_idx], dtype=np.int64)])
    tot_atlas_genes = cs2[r_ptr[1:]] - cs2[r_ptr[:-1]]
    bpos = {b: i for i, b in enumerate(r_bc)}
    found = np.array([bpos.get(b, -1) for b in atlas_bc])
    atlas_tot = np.asarray(M.sum(1)).ravel()
    ok = found >= 0
    same = int((tot_atlas_genes[found[ok]] == atlas_tot[ok]).sum())
    stats = {
        "library": LIB, "dataset": DATASET, "n_cells": int(e - s), "n_genes": int(n_genes),
        "nnz": int(M.nnz), "atlas_umis": float(atlas_tot.sum()),
        "n_raw_barcodes": int(len(r_bc)), "n_cellranger_called": int(is_called.sum()),
        "atlas_cells_in_raw_matrix": int(ok.sum()), "atlas_cells_total_equal_raw_matrix": same,
        "n_empty_droplets": int(empty.sum()), "empty_umis_all_features": float(tot[empty].sum()),
        "empty_umis_atlas_genes": float(bg.sum()), "raw_matrix_umis": float(tot.sum()),
        "called_cell_umis": float(tot[is_called].sum()),
        "empty_min_umi": EMPTY_MIN, "empty_max_umi": EMPTY_MAX,
    }
    json.dump(stats, open(out / "stats.json", "w"), indent=1)
    print(json.dumps(stats), flush=True)


if __name__ == "__main__":
    main()
