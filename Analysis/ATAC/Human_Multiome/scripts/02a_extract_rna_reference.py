#!/usr/bin/env python3
"""
02a_extract_rna_reference.py — Extract GSE244832 snRNA-seq reference subset
from the integrated atlas using h5py for maximum anndata version compatibility.

Reads the atlas h5ad via h5py (bypassing anndata IO which has version-specific
encoding issues), extracts obs metadata and expression matrix for GSE244832
donors, and saves as a clean h5ad readable by any anndata version.

Usage:
  python 02a_extract_rna_reference.py \
      --atlas ../../Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad \
      --donor_pairing ../../data/GSE244832/metadata/donor_pairing.csv \
      --output results/label_transfer/rna_reference_GSE244832.h5ad
"""

import argparse
import gc
import logging
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)


def read_categorical_h5py(group):
    """Read an anndata categorical column from h5py, handling missing 'ordered' attr."""
    codes = group["codes"][:]
    categories = group["categories"][:]
    if isinstance(categories[0], bytes):
        categories = [c.decode("utf-8") for c in categories]
    ordered = False
    if "ordered" in group.attrs:
        ordered = bool(group.attrs["ordered"])
    return pd.Categorical.from_codes(codes, categories=categories, ordered=ordered)


def read_string_array_h5py(dataset):
    """Read a string array from h5py."""
    data = dataset[:]
    if isinstance(data[0], bytes):
        return [s.decode("utf-8") for s in data]
    return list(data)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--atlas", required=True)
    p.add_argument("--donor_pairing", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    pairing = pd.read_csv(args.donor_pairing)
    rna_gsms = set(pairing["rna_gsm"].dropna().unique())
    log.info("Donor pairing: %d donors, %d GSMs", len(pairing), len(rna_gsms))

    # --- Read atlas via h5py ---
    log.info("Opening atlas: %s", args.atlas)
    f = h5py.File(args.atlas, "r")
    # Read obs_names
    obs_names = read_string_array_h5py(f["obs"]["_index"])
    n_cells = len(obs_names)
    log.info("Atlas: %d cells", n_cells)

    # Read var_names
    var_names = read_string_array_h5py(f["var"]["_index"])
    n_genes = len(var_names)
    log.info("Atlas: %d genes", n_genes)

    # Read obs columns we need
    obs_df = pd.DataFrame(index=obs_names)
    obs_group = f["obs"]

    for col_name in ["sample", "dataset", "cell_type", "cell_type_raw",
                     "cell_type_conf", "condition", "leiden"]:
        if col_name not in obs_group:
            log.warning("  obs column '%s' not found, skipping", col_name)
            continue
        col = obs_group[col_name]
        if isinstance(col, h5py.Group):
            # Categorical
            obs_df[col_name] = read_categorical_h5py(col)
            log.info("  %s: categorical (%d categories)", col_name,
                     len(obs_df[col_name].cat.categories))
        elif isinstance(col, h5py.Dataset):
            data = col[:]
            if data.dtype.kind in ("S", "U", "O"):
                obs_df[col_name] = read_string_array_h5py(col)
            else:
                obs_df[col_name] = data
            log.info("  %s: %s", col_name, obs_df[col_name].dtype)

    # --- Filter to GSE244832 ---
    gse_mask = None
    if "dataset" in obs_df.columns:
        gse_mask = obs_df["dataset"].astype(str).str.contains("GSE244832", case=False, na=False)
        log.info("Matched %d cells via dataset column", gse_mask.sum())

    if gse_mask is None or gse_mask.sum() == 0:
        if "sample" in obs_df.columns:
            gse_mask = obs_df["sample"].astype(str).isin(rna_gsms)
            log.info("Matched %d cells via GSM accessions", gse_mask.sum())

    if gse_mask is None or gse_mask.sum() == 0:
        log.error("No cells found for GSE244832!")
        f.close()
        sys.exit(1)

    indices = np.where(gse_mask.values)[0]
    n_ref = len(indices)
    log.info("Extracting %d reference cells...", n_ref)

    # --- Read expression matrix (X) for selected cells ---
    log.info("Reading expression matrix...")
    x_group = f["X"]

    if "encoding-type" in x_group.attrs:
        enc = x_group.attrs["encoding-type"]
        if isinstance(enc, bytes):
            enc = enc.decode()
    else:
        enc = "unknown"
    log.info("  X encoding: %s", enc)

    if enc in ("csr_matrix", "csc_matrix"):
        data_arr = x_group["data"][:]
        ind_arr = x_group["indices"][:]
        indptr_arr = x_group["indptr"][:]
        shape = tuple(x_group.attrs["shape"])

        if enc == "csr_matrix":
            # CSR: efficient row slicing — extract only the rows we need
            # without constructing the full matrix
            new_indptrs = [0]
            new_data = []
            new_indices = []
            for row_idx in indices:
                start, end = indptr_arr[row_idx], indptr_arr[row_idx + 1]
                new_data.append(data_arr[start:end])
                new_indices.append(ind_arr[start:end])
                new_indptrs.append(new_indptrs[-1] + (end - start))
            X_ref = sp.csr_matrix(
                (np.concatenate(new_data) if new_data else np.array([], dtype=data_arr.dtype),
                 np.concatenate(new_indices) if new_indices else np.array([], dtype=ind_arr.dtype),
                 np.array(new_indptrs, dtype=np.int64)),
                shape=(len(indices), shape[1]),
            )
            del new_data, new_indices, new_indptrs
        else:
            X_full = sp.csc_matrix((data_arr, ind_arr, indptr_arr), shape=shape)
            X_ref = X_full[indices].tocsr()
            del X_full
        del data_arr, ind_arr, indptr_arr
        gc.collect()
    elif isinstance(x_group, h5py.Dataset):
        # Dense
        X_ref = x_group[indices]
    else:
        log.error("Unknown X format: %s", type(x_group))
        f.close()
        sys.exit(1)

    log.info("  X shape: %s, nnz: %s", X_ref.shape,
             X_ref.nnz if sp.issparse(X_ref) else "dense")

    # --- Read obsm (PCA, UMAP) ---
    obsm = {}
    if "obsm" in f:
        for key in f["obsm"]:
            data = f["obsm"][key]
            if isinstance(data, h5py.Dataset):
                obsm[key] = data[:][indices]
                log.info("  obsm/%s: %s", key, obsm[key].shape)
            elif isinstance(data, h5py.Group):
                # Sparse obsm
                log.info("  obsm/%s: skipping (sparse group)", key)

    f.close()
    log.info("Atlas file closed.")

    # --- Build reference AnnData ---
    import anndata as ad

    obs_ref = obs_df.iloc[indices].copy()
    obs_ref.index = [obs_names[i] for i in indices]

    var_df = pd.DataFrame(index=var_names)

    adata_ref = ad.AnnData(
        X=X_ref,
        obs=obs_ref,
        var=var_df,
    )

    for key, val in obsm.items():
        adata_ref.obsm[key] = val

    log.info("Reference AnnData: %d cells x %d genes", adata_ref.n_obs, adata_ref.n_vars)
    log.info("Cell types:\n%s", adata_ref.obs["cell_type"].value_counts().to_string())

    # --- Save ---
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    log.info("Saving to %s...", args.output)
    adata_ref.write_h5ad(args.output)
    log.info("Done! File size: %.1f GB", os.path.getsize(args.output) / 1e9)


if __name__ == "__main__":
    main()
