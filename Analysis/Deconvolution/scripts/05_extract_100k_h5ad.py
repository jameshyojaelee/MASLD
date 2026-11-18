#!/usr/bin/env python3
"""Stratified 100,000-cell extraction from NUCLEUS-FIX unique_counts.h5ad (exclude == False)."""
import os
import json
import numpy as np
import h5py
from scipy import sparse
from scipy.io import mmwrite

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
H5AD = os.path.join(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts/unique_counts.h5ad")
DONOR_META = os.path.join(ROOT, "Analysis/Hypotheses/NUCLEUS-FIX-20261107/unique_counts/metadata/donor_metadata.tsv")
OUT_DIR = os.path.join(ROOT, "Analysis/Deconvolution/candidates/loco_20261118/tmp_sce")
os.makedirs(OUT_DIR, exist_ok=True)

TARGET_CELLS = 100000
SMALL_STRATUM_KEEP_ALL = 200
SEED = 20261118

# Parse donor_metadata.tsv without pandas
keep_donors = set()
with open(DONOR_META, "r", encoding="utf-8") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    i_donor = hdr.index("donor")
    i_excl = hdr.index("exclude")
    i_ds = hdr.index("dataset")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if parts[i_excl].strip().lower() == "false":
            keep_donors.add(parts[i_donor])

def decode_cat(group):
    if "categories" in group and "codes" in group:
        cats = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in group["categories"][:]]
        codes = group["codes"][:]
        return np.array([cats[c] if c >= 0 else "" for c in codes], dtype=object)
    raw = group[:]
    return np.array([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in raw], dtype=object)

with h5py.File(H5AD, "r") as f:
    donor = decode_cat(f["obs"]["donor"])
    dataset = decode_cat(f["obs"]["dataset"])
    cell_type = decode_cat(f["obs"]["cell_type"])
    obs_idx = decode_cat(f["obs"]["_index"])
    var_idx = decode_cat(f["var"]["_index"])
    var_ens = decode_cat(f["var"]["gene_ids"])

    keep_mask = np.array([d in keep_donors for d in donor], dtype=bool)
    cand_indices = np.where(keep_mask)[0]
    n_avail = len(cand_indices)
    assert n_avail >= TARGET_CELLS, f"Available non-excluded cells {n_avail} < {TARGET_CELLS}"

    # Group non-excluded indices by (dataset, cell_type)
    strata = {}
    for idx in cand_indices:
        key = (str(dataset[idx]), str(cell_type[idx]))
        strata.setdefault(key, []).append(idx)

    # Preserve all cells in small strata (<= 200), allocate remainder proportionally
    small_keys = [k for k, v in strata.items() if len(v) <= SMALL_STRATUM_KEEP_ALL]
    large_keys = [k for k, v in strata.items() if len(v) > SMALL_STRATUM_KEEP_ALL]
    alloc = {k: len(strata[k]) for k in small_keys}
    rem_target = TARGET_CELLS - sum(alloc.values())
    assert rem_target > 0
    large_total = sum(len(strata[k]) for k in large_keys)

    raw_quota = {k: rem_target * len(strata[k]) / large_total for k in large_keys}
    for k in large_keys:
        alloc[k] = min(len(strata[k]), max(1, int(np.floor(raw_quota[k]))))
    diff = TARGET_CELLS - sum(alloc.values())
    if diff > 0:
        order = sorted(large_keys, key=lambda k: raw_quota[k] - alloc[k], reverse=True)
        for k in order:
            if diff == 0:
                break
            if alloc[k] < len(strata[k]):
                alloc[k] += 1
                diff -= 1
    elif diff < 0:
        order = sorted(large_keys, key=lambda k: alloc[k] - raw_quota[k], reverse=True)
        for k in order:
            if diff == 0:
                break
            if alloc[k] > 1:
                alloc[k] -= 1
                diff += 1

    assert sum(alloc.values()) == TARGET_CELLS

    rng = np.random.default_rng(SEED)
    chosen = []
    for k in sorted(strata.keys()):
        arr = np.array(strata[k], dtype=np.int64)
        n_take = alloc[k]
        if n_take == len(arr):
            chosen.append(arr)
        else:
            sel = rng.choice(arr, size=n_take, replace=False)
            chosen.append(np.sort(sel))
    chosen_idx = np.sort(np.concatenate(chosen))
    assert len(chosen_idx) == TARGET_CELLS

    # Stream CSR rows cleanly from HDF5
    indptr = f["X"]["indptr"][:]
    indices_ds = f["X"]["indices"]
    data_ds = f["X"]["data"]
    n_genes = len(var_idx)

    new_indptr = np.zeros(TARGET_CELLS + 1, dtype=np.int64)
    col_chunks = []
    val_chunks = []
    pos = 0
    for i, r in enumerate(chosen_idx):
        s0, s1 = int(indptr[r]), int(indptr[r + 1])
        if s1 > s0:
            col_chunks.append(indices_ds[s0:s1])
            val_chunks.append(data_ds[s0:s1])
        pos += (s1 - s0)
        new_indptr[i + 1] = pos

    all_cols = np.concatenate(col_chunks)
    all_vals = np.concatenate(val_chunks)
    csr = sparse.csr_matrix((all_vals, all_cols, new_indptr), shape=(TARGET_CELLS, n_genes))
    csc_genes_by_cells = csr.T.tocsc()

mmwrite(os.path.join(OUT_DIR, "counts.mtx"), csc_genes_by_cells)

with open(os.path.join(OUT_DIR, "obs.tsv"), "w", encoding="utf-8") as fh:
    fh.write("cell_id\tdataset\tdonor\tcell_type\n")
    for r in chosen_idx:
        fh.write(f"{obs_idx[r]}\t{dataset[r]}\t{donor[r]}\t{cell_type[r]}\n")

with open(os.path.join(OUT_DIR, "var.tsv"), "w", encoding="utf-8") as fh:
    fh.write("symbol\tensembl_base\n")
    for s, e in zip(var_idx, var_ens):
        fh.write(f"{s}\t{e}\n")

summary = {
    "n_non_excluded_donors": len(keep_donors),
    "n_non_excluded_cells": int(n_avail),
    "n_sampled_cells": int(TARGET_CELLS),
    "n_genes": int(n_genes),
    "nnz": int(csc_genes_by_cells.nnz),
    "datasets": {str(ds): int(np.sum(dataset[chosen_idx] == ds)) for ds in np.unique(dataset[chosen_idx])},
    "cell_types": {str(ct): int(np.sum(cell_type[chosen_idx] == ct)) for ct in np.unique(cell_type[chosen_idx])},
}
with open(os.path.join(OUT_DIR, "extract_summary.json"), "w", encoding="utf-8") as fh:
    json.dump(summary, fh, indent=2, sort_keys=True)
print("EXTRACTED_100K_OK", json.dumps(summary["cell_types"], sort_keys=True))
