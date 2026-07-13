#!/usr/bin/env python
"""04c_fix_annotation.py -- repair cell_type on the fast-merge h5ad.

02's step5 fell back to leiden labels because
    ct_assignments.reindex(adata.obs_names)
raised "cannot reindex on an axis with duplicate labels" (10x barcodes collide
across the 12 donors). The gene-activity SCORING itself succeeded (all markers
found in the merge log) -- only the label transfer failed. Since gene_mat and the
processed AnnData are the SAME cells in the SAME order (both derive from one
AnnDataSet built by step2 in discover_donors order), the reindex is unnecessary:
recompute the gene matrix and assign cell_type POSITIONALLY, after verifying the
per-cell donor_id sequence is identical. Marker set + per-cell idxmax scoring are
copied VERBATIM from 02.step5 so the annotation method matches GSE244832 exactly.

Env: snapatac2.  Does NOT modify the shared 02 script.
"""
import os, sys, importlib.util, logging, gc
import numpy as np, pandas as pd, scipy.sparse as sp
import snapatac2 as snap
import anndata as ad

ROOT = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
S = os.path.join(ROOT, "Analysis/ATAC/Human_Multiome/scripts/02_snapatac2_processing.py")
spec = importlib.util.spec_from_file_location("snap02", S)
mod = importlib.util.module_from_spec(spec); sys.modules["snap02"] = mod
spec.loader.exec_module(mod)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("04c")

OUT       = os.path.join(ROOT, "Analysis/ATAC/Human_External/snapatac2_fast")
INPUT_DIR = os.path.join(ROOT, "Analysis/ATAC/Human_External/frag_links")
SAMPLE_CSV= os.path.join(ROOT, "data/GSE281367/metadata/donor_pairing.csv")
H5AD      = os.path.join(OUT, "snapatac2_processed.h5ad")

# copied verbatim from 02.step5_cell_type_annotation
LIVER_MARKERS = {
    "Hepatocyte":    ["ALB", "APOB", "CYP3A4", "CYP2E1", "HNF4A", "PCK1"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9"],
    "Stellate_Cell": ["ACTA2", "COL1A1", "PDGFRB", "DCN", "LUM"],
    "Endothelial":   ["PECAM1", "CDH5", "VWF", "ERG"],
    "LSEC":          ["CLEC4G", "CLEC4M", "FCN2", "FCN3", "STAB2"],
    "Kupffer_Cell":  ["CD68", "MARCO", "TIMD4", "CLEC4F"],
    "Macrophage":    ["CD14", "CD163", "CSF1R", "ITGAM"],
    "NK_T_Cell":     ["CD3E", "CD3D", "NKG7", "GNLY", "GZMB"],
    "B_Cell":        ["CD79A", "MS4A1", "CD19"],
    "Plasma_Cell":   ["JCHAIN", "MZB1", "SDC1"],
}

# rebuild the AnnDataSet in the SAME donor order the merge used (deterministic)
df = mod.load_sample_metadata(SAMPLE_CSV)
donors_info = mod.discover_donors(INPUT_DIR, df, logger)
processed = []
for d, cond, frag in donors_info:
    p = os.path.join(OUT, "per_donor", f"{d}.h5ad")
    processed.append((d, snap.read(p)))
logger.info("Reopened %d per-donor h5ads", len(processed))
dataset = mod.step2_create_dataset(processed, OUT, logger)
ds_donor = np.asarray(dataset.obs["donor_id"]).astype(str)

logger.info("Recomputing gene activity matrix (make_gene_matrix)...")
gene_mat = snap.pp.make_gene_matrix(dataset, gene_anno=snap.genome.hg38)
if hasattr(gene_mat, "to_memory"):
    gene_mat = gene_mat.to_memory()
logger.info("Gene activity matrix: %d cells x %d genes", gene_mat.n_obs, gene_mat.n_vars)

X = gene_mat.X
row_sums = np.maximum(np.asarray(X.sum(axis=1)).flatten(), 1)
gene_names = list(gene_mat.var_names)
scores = pd.DataFrame(index=np.arange(gene_mat.n_obs))          # integer index -> positional
for ct, markers in LIVER_MARKERS.items():
    present = [m for m in markers if m in gene_names]
    if not present:
        logger.warning("  No markers found for %s", ct); continue
    idx = [gene_names.index(m) for m in present]
    s = np.asarray(X[:, idx].sum(axis=1)).flatten()
    scores[ct] = s / row_sums
    logger.info("  %s: %d/%d markers", ct, len(present), len(markers))
ct_assign = scores.idxmax(axis=1).to_numpy()                    # per-cell, positional

# load processed h5ad, verify identical cell order, assign positionally
A = ad.read_h5ad(H5AD)
assert A.n_obs == len(ct_assign), f"cell mismatch {A.n_obs} vs {len(ct_assign)}"
a_donor = np.asarray(A.obs["donor_id"]).astype(str)
if not np.array_equal(a_donor, ds_donor):
    raise SystemExit("FATAL: donor_id order differs between processed h5ad and rebuilt "
                     "dataset -> positional assign UNSAFE. Aborting without write.")
logger.info("donor_id sequence verified identical (n=%d) -> positional assign safe", A.n_obs)

A.obs["cell_type"] = pd.Categorical(ct_assign)
A.write(H5AD)
logger.info("REPAIRED cell_type on %s\n%s", H5AD, A.obs["cell_type"].value_counts().to_string())

# also report the 4 DA-relevant types x condition
meta = pd.read_csv(SAMPLE_CSV)
cmap = dict(zip(meta.donor_id.astype(str), meta.condition.astype(str)))
A.obs["cond"] = a_donor  # temp
A.obs["cond"] = [cmap.get(d, "NA") for d in a_donor]
keep = ["Hepatocyte", "Stellate_Cell", "Macrophage", "Cholangiocyte"]
sub = A.obs[A.obs.cell_type.isin(keep)]
logger.info("DA cell types x condition:\n%s", pd.crosstab(sub.cell_type, sub.cond).to_string())
