#!/usr/bin/env python3
"""
41_govaere2026_cosmx_load.py — Load NanoString CosMx SMI 1000-plex
spatial transcriptomics for the Govaere 2026 MASLD/MASH deposit (GSE312698).

Inputs (data/external/govaere2026_natgenetics/cosmx/):
  GSM935170{4,5,6,7}_Leuven_{1-4}_exprMat_file.csv.gz        - cell × gene counts
  GSM935170{4,5,6,7}_Leuven_{1-4}_metadata_file.csv.gz       - per-cell features (PanCK/CD68/CD45/DAPI, Area, FOV, slide_ID, ...)
  GSM935170{4,5,6,7}_Leuven_{1-4}_fov_positions_file.csv.gz  - FOV grid coords (Slide / FOV / X_mm / Y_mm)

Per-Leuven cell counts (this deposit):
  Leuven_1 (GSM9351704):  405,940 cells   left + right explant livers (end-stage)
  Leuven_2 (GSM9351705):   70,886 cells   F3, MASL, Normal, F3
  Leuven_3 (GSM9351706):   35,537 cells   F1, F4
  Leuven_4 (GSM9351707):   40,313 cells   F3, F3, MASL

Output:
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad
  Analysis/Spatial/results/govaere2026/cosmx_cell_metadata.tsv

Conventions:
  - .X = raw integer counts (sparse CSR)
  - .var = TARGET genes only (Negative*, SystemControl*, FalseCode* probes excluded)
  - .obs columns: sample_id, cell_id, fov, slide_id, area_um2,
                  centerX_global_px, centerY_global_px, nCount_RNA, nFeature_RNA,
                  PanCK_intensity, CD68_intensity, CD45_intensity, DAPI_intensity,
                  sample_disease (slide-level composition string),
                  leiden, cell_type
  - .obsm["spatial"] = (centerX_global_px, centerY_global_px) microns
  - .uns["cosmx_panel"] = {"target_n", "neg_n", "sys_n", "false_n"}

Cell-type annotation:
  Pre-paper CosMx annotations (KC / MetMac / TransMac / preMac / Monocyte /
  cDC1 / cDC2 / migDC / Hepatocyte / Cholangiocyte / Endothelial / Mesenchymal /
  Lymphocyte) are NOT shipped in the GEO metadata CSVs. We run minimal QC +
  leiden clustering and assign each cluster to a canonical cell type via the
  marker-mean argmax of the paper's signatures.

SLURM:
  --partition=cpu --cpus-per-task=8 --mem=64G --time=8:00:00
"""

from __future__ import annotations

import argparse
import gc
import pathlib
import re
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
COSMX_DIR = PROJECT_ROOT / "data/external/govaere2026_natgenetics/cosmx"
OUT_PREPROC = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed"
OUT_GOVAERE = PROJECT_ROOT / "Analysis/Spatial/results/govaere2026"

# ── Sample registry ─────────────────────────────────────────────────────────
# Disease composition strings sourced from GEO Sample_description on
# 2026-05-21 for accessions GSM9351704–07 (Series GSE312698):
#   "n=1 normal biopsy, n=2 MASL biopsies, n=1 F1 biopsy, n=4 F3 biopsies,
#    n=1 F4 biopsy, n=2 end-stage explants."
SAMPLES = [
    {
        "sample_id": "Leuven_1",
        "gsm": "GSM9351704",
        "expr": "GSM9351704_Leuven_1_exprMat_file.csv.gz",
        "meta": "GSM9351704_Leuven_1_metadata_file.csv.gz",
        "fov":  "GSM9351704_Leuven_1_fov_positions_file.csv.gz",
        "disease_composition": "Explant_x2",   # 2 end-stage explants on same slide
        "expected_n_tissues": 2,
        "n_expected_cells": 405940,
    },
    {
        "sample_id": "Leuven_2",
        "gsm": "GSM9351705",
        "expr": "GSM9351705_Leuven_2_exprMat_file.csv.gz",
        "meta": "GSM9351705_Leuven_2_metadata_file.csv.gz",
        "fov":  "GSM9351705_Leuven_2_fov_positions_file.csv.gz",
        "disease_composition": "F3_MASL_Normal_F3",
        "expected_n_tissues": 4,
        "n_expected_cells": 70886,
    },
    {
        "sample_id": "Leuven_3",
        "gsm": "GSM9351706",
        "expr": "GSM9351706_Leuven_3_exprMat_file.csv.gz",
        "meta": "GSM9351706_Leuven_3_metadata_file.csv.gz",
        "fov":  "GSM9351706_Leuven_3_fov_positions_file.csv.gz",
        "disease_composition": "F1_F4",
        "expected_n_tissues": 2,
        "n_expected_cells": 35537,
    },
    {
        "sample_id": "Leuven_4",
        "gsm": "GSM9351707",
        "expr": "GSM9351707_Leuven_4_exprMat_file.csv.gz",
        "meta": "GSM9351707_Leuven_4_metadata_file.csv.gz",
        "fov":  "GSM9351707_Leuven_4_fov_positions_file.csv.gz",
        "disease_composition": "F3_F3_MASL",
        "expected_n_tissues": 3,
        "n_expected_cells": 40313,
    },
]

# Probe types to strip from the panel (CosMx Universal 1000-plex controls).
_NEG_RE   = re.compile(r"^Negative", re.IGNORECASE)
_SYS_RE   = re.compile(r"^SystemControl", re.IGNORECASE)
_FALSE_RE = re.compile(r"^FalseCode", re.IGNORECASE)

# Canonical cell-type markers (Govaere 2026 Fig 3f labels). Each entry is a
# tuple (label, marker_genes). Multiple markers per type → mean Z across
# cells, then per-cluster mean Z; cluster label = argmax type.
CELLTYPE_MARKERS = {
    "KC":            ["MARCO", "CD5L", "VSIG4", "CLEC4F", "TIMD4"],
    "MetMac":        ["GPNMB", "HS3ST2", "LPL", "FABP5", "TREM2", "CD9"],
    "TransMac":      ["CXCL10", "CXCL9", "CXCL11", "IL1B", "TNF"],
    "preMac":        ["PCNX2", "RUNX2", "PLTP", "CCL18"],
    "Monocyte":      ["VCAN", "FCN1", "S100A8", "S100A9", "CD14"],
    "cDC1":          ["CLEC9A", "XCR1", "IRF8"],
    "cDC2":          ["CLEC10A", "FCER1A", "CD1C"],
    "migDC":         ["LAMP3", "CCR7", "FSCN1"],
    "Hepatocyte":    ["ALB", "APOA1", "TTR", "TF", "SERPINA1"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9"],
    "Endothelial":   ["PECAM1", "VWF", "CDH5", "STAB1", "STAB2"],
    "Mesenchymal":   ["ACTA2", "COL1A1", "COL1A2", "COL3A1", "PDGFRA", "PDGFRB"],
    "Lymphocyte":    ["PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD79A", "NKG7"],
}

# Marker genes that MUST appear in the panel (used for the validation step).
REQUIRED_PANEL_GENES = [
    "GPNMB", "LPL", "FABP5", "HS3ST2",   # MetMac
    "MARCO", "CD5L",                      # KC
    "ALB",                                # Hepatocyte
    "KRT19",                              # Cholangiocyte
    "PECAM1",                             # Endothelial
    "ACTA2", "COL1A1",                    # Mesenchymal
    "PTPRC", "CD3D",                      # Lymphocyte
    "VCAN", "FCN1",                       # Monocyte
]


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _split_panel(gene_names: list[str]) -> dict:
    """Partition CosMx panel into target genes vs Neg/SystemControl/FalseCode probes."""
    is_neg   = np.array([bool(_NEG_RE.match(g))   for g in gene_names])
    is_sys   = np.array([bool(_SYS_RE.match(g))   for g in gene_names])
    is_false = np.array([bool(_FALSE_RE.match(g)) for g in gene_names])
    is_target = ~(is_neg | is_sys | is_false)
    return {
        "target_idx": np.where(is_target)[0],
        "neg_idx":    np.where(is_neg)[0],
        "sys_idx":    np.where(is_sys)[0],
        "false_idx":  np.where(is_false)[0],
        "n_target": int(is_target.sum()),
        "n_neg":    int(is_neg.sum()),
        "n_sys":    int(is_sys.sum()),
        "n_false":  int(is_false.sum()),
    }


def load_one_sample(samp: dict) -> ad.AnnData:
    """Read one Leuven slide's expression matrix + metadata + FOV grid into AnnData."""
    sid = samp["sample_id"]
    _log(f"  [{sid}] reading expression matrix")
    expr = pd.read_csv(COSMX_DIR / samp["expr"])
    n_cells = expr.shape[0]
    _log(f"    expr shape: {expr.shape}  (cells × cols)")

    # First two columns are (fov, cell_ID) per CosMx layout — gene columns follow.
    leading_cols = ["fov", "cell_ID"]
    assert list(expr.columns[:2]) == leading_cols, (
        f"Unexpected leading columns: {list(expr.columns[:2])}")
    gene_cols = list(expr.columns[2:])
    panel = _split_panel(gene_cols)
    _log(f"    panel: target={panel['n_target']}  Neg={panel['n_neg']}  "
         f"SystemControl={panel['n_sys']}  FalseCode={panel['n_false']}")

    # Strip controls. Convert to int32 (counts) sparse CSR.
    target_idx = panel["target_idx"]
    target_genes = [gene_cols[i] for i in target_idx]
    counts_dense = expr.iloc[:, 2:].values  # full dense matrix
    counts_dense = counts_dense[:, target_idx]
    X = sparse.csr_matrix(counts_dense.astype(np.int32))
    del counts_dense
    gc.collect()
    _log(f"    target counts matrix: {X.shape}, nnz={X.nnz}, "
         f"density={X.nnz/(X.shape[0]*X.shape[1]):.4f}")

    # Cell identifiers (per-FOV unique within a slide). Make globally unique by prefixing sample.
    fov_arr = expr["fov"].astype(int).values
    cellid_arr = expr["cell_ID"].astype(int).values
    cell_index = [f"{sid}_FOV{fov_arr[i]}_C{cellid_arr[i]}" for i in range(n_cells)]
    del expr
    gc.collect()

    # Metadata file
    _log(f"  [{sid}] reading per-cell metadata")
    meta = pd.read_csv(COSMX_DIR / samp["meta"])
    assert meta.shape[0] == n_cells, (
        f"Cell count mismatch in {sid}: meta={meta.shape[0]} vs expr={n_cells}")
    # Build a lookup from (fov, cell_ID) → row in meta, then realign to expr order
    meta_key = list(zip(meta["fov"].astype(int), meta["cell_ID"].astype(int)))
    expr_key = list(zip(fov_arr.tolist(), cellid_arr.tolist()))
    if meta_key == expr_key:
        meta_aligned = meta.reset_index(drop=True)
    else:
        _log(f"    [{sid}] metadata row order differs from expr — re-aligning")
        idx_lookup = {k: i for i, k in enumerate(meta_key)}
        row_idx = np.array([idx_lookup[k] for k in expr_key], dtype=np.int64)
        meta_aligned = meta.iloc[row_idx].reset_index(drop=True)
    del meta

    # FOV positions for spatial QC + plotting (joined to obs by fov column).
    fov_pos = pd.read_csv(COSMX_DIR / samp["fov"])
    fov_pos = fov_pos.set_index("FOV")[["X_mm", "Y_mm", "Slide"]]
    obs_fov_xy = fov_pos.reindex(fov_arr)

    # Assemble obs.
    # IMPORTANT: use `.values` for every pandas Series so pandas does NOT try to
    # align by index — `meta_aligned` has integer RangeIndex while `cell_index`
    # is string `f"{sid}_FOV{fov}_C{cell_id}"`; mismatched indices silently
    # produce all-NaN columns (centerX_global_px regression caught 2026-05-21).
    obs = pd.DataFrame({
        "sample_id":         sid,
        "gsm":               samp["gsm"],
        "cell_id":           cellid_arr.astype(np.int64),
        "fov":               fov_arr.astype(np.int32),
        "slide_id":          meta_aligned["slide_ID"].astype("int32").values,
        "area_um2":          meta_aligned["Area.um2"].astype("float32").values,
        "centerX_global_px": meta_aligned["CenterX_global_px"].astype("float32").values,
        "centerY_global_px": meta_aligned["CenterY_global_px"].astype("float32").values,
        "nCount_RNA":        meta_aligned["nCount_RNA"].astype("int32").values,
        "nFeature_RNA":      meta_aligned["nFeature_RNA"].astype("int32").values,
        "PanCK_intensity":   meta_aligned["Mean.PanCK"].astype("float32").values,
        "CD68_intensity":    meta_aligned["Mean.CD68"].astype("float32").values,
        "CD45_intensity":    meta_aligned["Mean.CD45"].astype("float32").values,
        "DAPI_intensity":    meta_aligned["Mean.DAPI"].astype("float32").values,
        "fov_X_mm":          obs_fov_xy["X_mm"].astype("float32").values,
        "fov_Y_mm":          obs_fov_xy["Y_mm"].astype("float32").values,
        "sample_disease":    samp["disease_composition"],
    }, index=pd.Index(cell_index, name="obs_name"))

    # Var
    var = pd.DataFrame(index=pd.Index(target_genes, name="gene_symbol"))
    var["target_panel"] = True

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["sample_id"] = sid
    adata.uns["cosmx_panel"] = {
        "target_n": panel["n_target"],
        "neg_n":    panel["n_neg"],
        "sys_n":    panel["n_sys"],
        "false_n":  panel["n_false"],
    }
    # Spatial coords (in pixels — 180 nm / px per CosMx SMI specs)
    adata.obsm["spatial"] = obs[["centerX_global_px", "centerY_global_px"]].to_numpy(dtype=np.float32)
    _log(f"  [{sid}] AnnData built: {adata.shape}")
    return adata


def concat_samples(adatas: list[ad.AnnData]) -> ad.AnnData:
    """Inner-join concatenation across the four Leuven slides."""
    _log("  Concatenating samples (inner-join on var)")
    merged = ad.concat(adatas, join="inner", merge="same", uns_merge="unique",
                       label="sample_id_concat", index_unique=None)
    # ad.concat creates a string suffix when index_unique is set; here we already
    # pre-prefixed cells so just reuse existing obs_names (must be globally unique).
    merged.obs_names_make_unique()
    _log(f"  Merged shape: {merged.shape}")
    return merged


def minimal_qc(adata: ad.AnnData,
               min_counts: int = 20,
               min_genes: int = 10,
               min_cells_per_gene: int = 5) -> ad.AnnData:
    """Soft QC for CosMx (per-cell totals are much lower than 10x scRNA)."""
    _log(f"  QC: min_counts={min_counts}, min_genes={min_genes}, "
         f"min_cells_per_gene={min_cells_per_gene}")
    n_before = adata.n_obs
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    mask = (adata.obs["total_counts"] >= min_counts) & (adata.obs["n_genes_by_counts"] >= min_genes)
    adata = adata[mask].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    _log(f"  QC kept {adata.n_obs}/{n_before} cells ({100*adata.n_obs/n_before:.1f}%), "
         f"{adata.n_vars} genes")
    return adata


def cluster_and_label(adata: ad.AnnData, n_pcs: int = 30,
                      n_neighbors: int = 15, resolution: float = 0.5,
                      hvg_n: int = 2000) -> ad.AnnData:
    """Leiden clustering + marker-based cell-type labelling."""
    _log("  Normalising (target=1e4) + log1p")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_target = adata.n_vars
    n_use = min(hvg_n, n_target)
    _log(f"  Selecting up to {n_use} HVGs")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_use, flavor="seurat")
    _log(f"  Scaling + PCA")
    sc.pp.scale(adata, max_value=10, zero_center=True)
    n_pcs_use = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_pcs_use)
    _log(f"  Neighbors + Leiden (resolution={resolution})")
    sc.pp.neighbors(adata, n_pcs=n_pcs_use, n_neighbors=n_neighbors)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    _log(f"  Found {n_clusters} leiden clusters")

    # Score each cell on each marker set (use scanpy score_genes on log-normalized data)
    _log("  Scoring canonical cell-type markers")
    available = set(adata.var_names)
    for label, markers in CELLTYPE_MARKERS.items():
        present = [m for m in markers if m in available]
        if not present:
            _log(f"    WARNING: no markers present for {label}: {markers}")
            adata.obs[f"score_{label}"] = 0.0
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{label}",
                          random_state=0, ctrl_size=50, use_raw=False)
    # Cluster-level mean score
    score_cols = [f"score_{l}" for l in CELLTYPE_MARKERS]
    cluster_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    # argmax per cluster
    cluster_label = cluster_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
    _log("  Cluster → cell-type assignment:")
    for c, lbl in cluster_label.items():
        nc = (adata.obs["leiden"] == c).sum()
        _log(f"    cluster {c:>3} → {lbl:>14}   ({nc:,} cells)")
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype("category")
    return adata


def validate(adata: ad.AnnData, samples: list[dict]) -> dict:
    """Run sanity checks per task spec; returns a dict of summary stats."""
    issues = []
    summary = {}

    sample_ids = list(adata.obs["sample_id"].unique())
    summary["n_samples"] = len(sample_ids)
    summary["sample_ids"] = sorted(sample_ids)
    expected_sids = sorted([s["sample_id"] for s in samples])
    if summary["sample_ids"] != expected_sids:
        issues.append(f"sample_id mismatch: got {summary['sample_ids']}, expected {expected_sids}")

    summary["sample_id_counts"] = adata.obs["sample_id"].value_counts().to_dict()
    summary["total_cells"] = int(adata.n_obs)
    summary["total_genes"] = int(adata.n_vars)

    if "cell_type" in adata.obs:
        ct_counts = adata.obs["cell_type"].value_counts()
        summary["cell_type_counts"] = ct_counts.to_dict()
        summary["n_cell_types"] = int(ct_counts.shape[0])
        if ct_counts.shape[0] < 8:
            issues.append(f"only {ct_counts.shape[0]} cell types (expected >= 8)")

    present = [g for g in REQUIRED_PANEL_GENES if g in adata.var_names]
    missing = [g for g in REQUIRED_PANEL_GENES if g not in adata.var_names]
    summary["required_markers_present"] = present
    summary["required_markers_missing"] = missing
    if missing:
        issues.append(f"missing required markers: {missing}")

    summary["issues"] = issues
    return summary


def main():
    parser = argparse.ArgumentParser(description="Build CosMx AnnData for GSE312698 (Govaere 2026)")
    parser.add_argument("--samples", type=str, default=None,
                        help="comma-separated subset (Leuven_1,...,Leuven_4); default all")
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--hvg-n", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-counts", type=int, default=20)
    parser.add_argument("--min-genes", type=int, default=10)
    parser.add_argument("--min-cells-per-gene", type=int, default=5)
    parser.add_argument("--output", type=str,
                        default=str(OUT_PREPROC / "cosmx_govaere2026.h5ad"))
    parser.add_argument("--meta-tsv", type=str,
                        default=str(OUT_GOVAERE / "cosmx_cell_metadata.tsv"))
    args = parser.parse_args()

    OUT_PREPROC.mkdir(parents=True, exist_ok=True)
    OUT_GOVAERE.mkdir(parents=True, exist_ok=True)

    selected = SAMPLES
    if args.samples:
        wanted = set(s.strip() for s in args.samples.split(","))
        selected = [s for s in SAMPLES if s["sample_id"] in wanted]
        if not selected:
            sys.exit(f"No samples match --samples={args.samples}")

    _log("=" * 70)
    _log(f"41_govaere2026_cosmx_load.py — loading {len(selected)} sample(s)")
    _log(f"Output AnnData: {args.output}")
    _log(f"Output meta TSV: {args.meta_tsv}")
    _log("=" * 70)

    adatas = []
    for samp in selected:
        adata = load_one_sample(samp)
        # Sanity check expected cell counts
        if adata.n_obs != samp["n_expected_cells"]:
            _log(f"  WARNING [{samp['sample_id']}]: n_cells={adata.n_obs} differs "
                 f"from expected {samp['n_expected_cells']}")
        adatas.append(adata)

    # Concatenate
    adata = concat_samples(adatas)
    del adatas
    gc.collect()

    # Validate every required marker survived the merge
    panel_summary = {
        "target_n": int(adata.n_vars),
        "samples": [s["sample_id"] for s in selected],
    }
    _log(f"  Merged target panel size: {adata.n_vars}")

    # QC
    adata = minimal_qc(adata,
                       min_counts=args.min_counts,
                       min_genes=args.min_genes,
                       min_cells_per_gene=args.min_cells_per_gene)

    # Clustering + cell-type labels
    adata = cluster_and_label(adata,
                              n_pcs=args.n_pcs,
                              n_neighbors=args.n_neighbors,
                              resolution=args.leiden_resolution,
                              hvg_n=args.hvg_n)

    # Restore raw counts at .X for downstream cell2location / squidpy
    _log("  Restoring raw counts to .X (preserved in .layers['counts'])")
    adata.X = adata.layers["counts"]
    # log-normalised view kept in .layers["lognorm"]
    # rebuild for downstream
    ln = adata.layers["counts"].copy().astype(np.float32)
    if sparse.issparse(ln):
        # sc.pp.normalize_total in-place mutates; do it on a temp copy of .X
        adata_norm = ad.AnnData(X=ln, obs=adata.obs[[]].copy(), var=adata.var[[]].copy())
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        adata.layers["lognorm"] = adata_norm.X
        del adata_norm
    gc.collect()

    # Final validation
    summary = validate(adata, selected)
    _log("=" * 70)
    _log("VALIDATION SUMMARY")
    for k, v in summary.items():
        _log(f"  {k}: {v}")
    _log("=" * 70)

    # Stash the panel summary in .uns
    adata.uns["panel_summary"] = panel_summary
    adata.uns["validation"] = {k: str(v) for k, v in summary.items()}

    # Write outputs
    _log(f"  Writing AnnData → {args.output}")
    adata.write_h5ad(args.output, compression="gzip")
    out_size = pathlib.Path(args.output).stat().st_size
    _log(f"  AnnData written ({out_size / 1e9:.2f} GB)")

    _log(f"  Writing cell metadata TSV → {args.meta_tsv}")
    meta_out = adata.obs.copy()
    meta_out.index.name = "obs_name"
    meta_out.to_csv(args.meta_tsv, sep="\t")
    _log(f"  Cell metadata TSV written ({pathlib.Path(args.meta_tsv).stat().st_size/1e6:.2f} MB)")

    _log("DONE")


if __name__ == "__main__":
    main()
