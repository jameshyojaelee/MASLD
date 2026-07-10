#!/usr/bin/env python
# 25_build_rectangle_reference.py
# ---------------------------------------------------------------------------
# Build the two Rectangle (Wagner et al. 2024) single-cell references:
#   1. COARSE  : 16 canonical major cell types.
#   2. FINE    : multiscale — non-hepatocytes keep their major cell type;
#                hepatocytes are replaced by their 5-class meta-subtype
#                ("Hep_" + {Healthy, Disease-Neutral, Neutral,
#                 Disease-Progressor, Disease-Associated}).
#
# Both references keep .X = RAW INTEGER COUNTS and var_names = gene SYMBOLS,
# matching the shared Rectangle contract.
#
# Input:
#   - Analysis/Deconvolution/reference/reference_human_scalesc_100k.h5ad
#       .X = raw counts, obs['cell_type'] = 16 types, obs index = {sample}_{barcode}
#   - Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/
#       hepatocyte_subtype_metadata.csv  (index = {sample}_{barcode};
#                                          col 'hepatocyte_subtype' = 0..38)
#       meta_subtype_mapping.csv          (col 'subtype' 0..38 -> 'meta_subtype')
#
# Output:
#   - Analysis/Deconvolution/reference/rectangle_ref_coarse.h5ad
#   - Analysis/Deconvolution/reference/rectangle_ref_fine.h5ad
#
# Run with the rectangle env python:
#   $ROOT/RNA-seq/.mamba/rectangle/bin/python 25_build_rectangle_reference.py
# ---------------------------------------------------------------------------

import os
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import issparse

try:
    from anndata.experimental import read_elem
except ImportError:  # older/newer layout
    from anndata._io.specs import read_elem


def load_adata_core(path):
    """Load only X/obs/var from an h5ad, skipping uns/obsm/layers.

    The reference was written by a newer anndata that encodes a None scalar
    (uns/log1p/base) with an IOSpec ('null') that anndata 0.10.8 cannot read.
    Rectangle only needs the count matrix, cell-type labels, and gene symbols,
    so we read those three elements directly and rebuild a clean AnnData.
    """
    with h5py.File(path, "r") as f:
        X = read_elem(f["X"])
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
    return ad.AnnData(X=X, obs=obs, var=var)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

REF_IN = os.path.join(
    BASE, "Analysis/Deconvolution/reference/reference_human_scalesc_100k.h5ad"
)
HEP_META = os.path.join(
    BASE,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/"
    "hepatocyte_subtype_metadata.csv",
)
META_MAP = os.path.join(
    BASE,
    "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/"
    "meta_subtype_mapping.csv",
)
REF_DIR = os.path.join(BASE, "Analysis/Deconvolution/reference")
OUT_COARSE = os.path.join(REF_DIR, "rectangle_ref_coarse.h5ad")
OUT_FINE = os.path.join(REF_DIR, "rectangle_ref_fine.h5ad")

CANONICAL_CELLTYPES = {
    "Endothelial cells", "Hepatocytes", "Plasma cells", "T cells",
    "Cholangiocytes", "Fibroblasts", "Macrophages", "Circulating NK/NKT",
    "Resident NK", "Mono+mono derived cells", "Basophils", "B cells",
    "cDC1s", "cDC2s", "pDCs", "Neutrophils",
}

for p in (REF_IN, HEP_META, META_MAP):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing required input: {p}")

for p in (OUT_COARSE, OUT_FINE):
    if os.path.exists(p):
        sys.exit(f"ERROR: refusing to overwrite existing output: {p}")

print("=== 25: Build Rectangle references (coarse + fine) ===", flush=True)

# --- 1. Load reference ------------------------------------------------------
print(f"\nLoading reference (X/obs/var only): {REF_IN}", flush=True)
adata = load_adata_core(REF_IN)
print(f"  {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
print(f"  obs columns: {list(adata.obs.columns)}", flush=True)
print(f"  obs index examples: {list(adata.obs_names[:3])}", flush=True)

if "cell_type" not in adata.obs.columns:
    sys.exit("ERROR: obs has no 'cell_type' column")

# --- 2. Verify .X are raw integer counts ------------------------------------
X = adata.X
if issparse(X):
    data = X.data
else:
    data = np.asarray(X).ravel()
n_check = min(data.size, 5_000_000)
sample_vals = data[:n_check]
max_frac = float(np.max(np.abs(sample_vals - np.round(sample_vals)))) if sample_vals.size else 0.0
print(f"\n.X integer check: max |val - round(val)| over {n_check} nnz = {max_frac:.3g}",
      flush=True)
assert max_frac < 1e-6, ".X is not integer-valued raw counts"
assert float(np.min(data)) >= 0.0 if data.size else True, ".X has negative values"
print("  OK: .X are non-negative integer raw counts", flush=True)

# --- 3. COARSE reference ----------------------------------------------------
adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
ct_counts = adata.obs["cell_type"].value_counts()
print("\n[COARSE] obs['cell_type'] value_counts:", flush=True)
print(ct_counts.to_string(), flush=True)

observed = set(adata.obs["cell_type"].unique())
unexpected = observed - CANONICAL_CELLTYPES
if unexpected:
    print(f"  WARNING: cell types outside the 16 canonical set: {sorted(unexpected)}",
          flush=True)
missing_ct = CANONICAL_CELLTYPES - observed
if missing_ct:
    print(f"  NOTE: canonical types absent from this reference: {sorted(missing_ct)}",
          flush=True)

adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

adata.write_h5ad(OUT_COARSE)
print(f"\n[COARSE] wrote -> {OUT_COARSE}", flush=True)

# --- 4. FINE reference (hepatocyte multiscale) ------------------------------
print("\n[FINE] Building hepatocyte-substate labels...", flush=True)

hep_meta = pd.read_csv(HEP_META, index_col=0)
print(f"  hepatocyte_subtype_metadata: {hep_meta.shape[0]} cells, "
      f"cols={list(hep_meta.columns)}", flush=True)
if "hepatocyte_subtype" not in hep_meta.columns:
    sys.exit("ERROR: hepatocyte metadata missing 'hepatocyte_subtype' column")

meta_map = pd.read_csv(META_MAP)
print(f"  meta_subtype_mapping: cols={list(meta_map.columns)}", flush=True)
if not {"subtype", "meta_subtype"}.issubset(meta_map.columns):
    sys.exit("ERROR: mapping missing 'subtype'/'meta_subtype' columns")

# subtype (int 0..38) -> meta_subtype (str)
sub2meta = dict(zip(meta_map["subtype"].astype(int), meta_map["meta_subtype"].astype(str)))
print(f"  meta_subtype classes: {sorted(set(sub2meta.values()))}", flush=True)

# per-cell-barcode -> meta_subtype (via hepatocyte_subtype -> mapping)
hep_sub = hep_meta["hepatocyte_subtype"]
# drop any rows with missing subtype before int-casting
hep_sub = hep_sub.dropna().astype(int)
barcode2meta = hep_sub.map(sub2meta)
barcode2meta = barcode2meta.dropna()
print(f"  hepatocyte barcodes with a meta_subtype: {barcode2meta.shape[0]}", flush=True)

is_hep = adata.obs["cell_type"].astype(str) == "Hepatocytes"
n_hep_ref = int(is_hep.sum())
print(f"  reference hepatocytes: {n_hep_ref}", flush=True)

# Default fine label = coarse cell_type
cell_type_fine = adata.obs["cell_type"].astype(str).copy()

# For reference hepatocytes, look up meta_subtype by obs index
hep_idx = adata.obs_names[is_hep.values]
matched_meta = barcode2meta.reindex(hep_idx)  # NaN where unmatched
n_matched = int(matched_meta.notna().sum())
coverage = 100.0 * n_matched / n_hep_ref if n_hep_ref else 0.0
print(f"  matched to a meta_subtype: {n_matched} / {n_hep_ref} "
      f"({coverage:.1f}% coverage)", flush=True)

# Assign "Hep_<meta_subtype>" for matched hepatocytes; unmatched keep "Hepatocytes"
fine_hep_labels = matched_meta.apply(
    lambda m: f"Hep_{m}" if pd.notna(m) else "Hepatocytes"
)
cell_type_fine.loc[hep_idx] = fine_hep_labels.values

adata.obs["cell_type_fine"] = pd.Categorical(cell_type_fine.values)

print("\n[FINE] obs['cell_type_fine'] value_counts:", flush=True)
print(adata.obs["cell_type_fine"].value_counts().to_string(), flush=True)

adata.write_h5ad(OUT_FINE)
print(f"\n[FINE] wrote -> {OUT_FINE}", flush=True)

# --- 5. Summary -------------------------------------------------------------
print("\n=== SUMMARY ===", flush=True)
print(f"  coarse: {OUT_COARSE}  ({adata.n_obs} cells x {adata.n_vars} genes, "
      f"{adata.obs['cell_type'].nunique()} cell types)", flush=True)
print(f"  fine  : {OUT_FINE}  ({adata.obs['cell_type_fine'].nunique()} fine labels)",
      flush=True)
print(f"  hepatocyte meta-subtype coverage: {n_matched}/{n_hep_ref} ({coverage:.1f}%)",
      flush=True)
print("=== Done ===", flush=True)
