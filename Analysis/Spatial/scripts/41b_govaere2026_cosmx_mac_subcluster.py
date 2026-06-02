#!/usr/bin/env python3
"""
41b_govaere2026_cosmx_mac_subcluster.py — Sub-cluster the CosMx "KC" macrophage
bin from the Govaere 2026 deposit (GSE312698) to recover MetMac / TransMac /
KC / preMac / Monocyte. Paper's PRIMARY finding is the GPNMB+/MetMac axis;
this script directly replicates it by separating the merged macrophage cluster
the loader (41_*) produced at leiden res=0.5.

The loader already computed per-cell marker SCORES for each macrophage
lineage (sc.tl.score_genes on log-normalised counts):
  score_KC, score_MetMac, score_TransMac, score_preMac,
  score_Monocyte, score_cDC1, score_cDC2, score_migDC
and assigned `.obs["cell_type"]` via the cluster-mean argmax. The argmax merged
MetMac + TransMac into the "KC" bin because Leiden res=0.5 did not separate
them (resolution too coarse for the panel).

Strategy (per Govaere 2026 Methods page 14: "Cells were subclustered for
further analysis. They underwent new Louvain clustering (resolution = 1-2)
and UMAP reduction with 15-50 dimensions"):
  1. Subset to .obs.cell_type == "KC" (~100,525 cells).
  2. Re-find HVGs within the macrophage subset (~500 from 968 target genes).
  3. Scale + PCA (50 PCs).
  4. Compute neighbors (k=15, 50 PCs).
  5. Leiden res=1.5 (centre of paper's "1-2" range).
  6. Assign sub-cluster labels via score-based argmax + marker-based
     tie-breaker (GPNMB threshold for MetMac override).
  7. Run Wilcoxon MASH-vs-noMASH DE separately for each subtype.
  8. Write a NEW h5ad with the subset + mac_subtype col; DE CSVs;
     composition + marker summary TSVs.

Inputs (loader output, do NOT overwrite):
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad

F137/F138 CAVEAT (MASH contrast is NOT a clean case/control):
  The MASH-vs-no_MASH DE is built on a slide-level grouping with only 4 slides,
  and the single 'no_MASH' comparator slide (Leuven_2) is ~75% MASLD tissue
  (1 Normal + 1 MASL + 2 F3) — so the reference is mostly disease itself. The
  cell-level Wilcoxon p-values are also pseudoreplicated (cells within a slide are
  not independent). Therefore this is a NOISY MASH-ENRICHMENT descriptor, not a
  MASH-vs-control contrast, and the per-slide MetMac PROPORTIONS (cosmx_mac_subtype
  _composition.tsv) — which do NOT monotonically increase with MASH — are the
  primary caveated result, not the Wilcoxon p-values. Slide-level direction is
  surfaced as slide_logfc / slide_direction_concordant in each DE CSV.

Outputs:
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026_macsubclustered.h5ad
  Analysis/Spatial/results/govaere2026/cosmx_de_MetMac_MASH_vs_noMASH.csv  (KEY)
  Analysis/Spatial/results/govaere2026/cosmx_de_TransMac_MASH_vs_noMASH.csv
  Analysis/Spatial/results/govaere2026/cosmx_de_KC_post_subcluster_MASH_vs_noMASH.csv
  Analysis/Spatial/results/govaere2026/cosmx_de_preMac_MASH_vs_noMASH.csv      (if n>=50)
  Analysis/Spatial/results/govaere2026/cosmx_de_Monocyte_MASH_vs_noMASH.csv    (if n>=50)
  Analysis/Spatial/results/govaere2026/cosmx_mac_subtype_composition.tsv
  Analysis/Spatial/results/govaere2026/cosmx_mac_subtype_markers.tsv

Environment: micromamba activate spatial
"""

from __future__ import annotations

import argparse
import gc
import pathlib
import sys
import time
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Slide-aware helper (F138): slide-level pseudobulk for the MASH contrast.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from spatial_stats import pseudobulk_by_donor  # noqa: E402

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ADATA_IN = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad"
ADATA_OUT = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026_macsubclustered.h5ad"
OUT_DIR = PROJECT_ROOT / "Analysis/Spatial/results/govaere2026"

# ── Slide-level disease assignment (matches 42b_govaere2026_cosmx_de.py) ─────
SAMPLE_DISEASE = {
    "Leuven_1": {"sample_disease_str": "Explant_x2",        "stage": "MASH"},
    "Leuven_2": {"sample_disease_str": "F3_MASL_Normal_F3", "stage": "no_MASH"},
    "Leuven_3": {"sample_disease_str": "F1_F4",             "stage": "MASH"},
    "Leuven_4": {"sample_disease_str": "F3_F3_MASL",        "stage": "MASH"},
}

# Marker-based override panels for the score-based winner-take-all step.
# Mirrors CELLTYPE_MARKERS from 41_*.py but restricted to macrophage lineage
# subtypes we want to separate within the KC bin.
SUBTYPE_MARKER_PANELS = {
    "KC":       ["MARCO", "CD5L", "VSIG4", "CLEC4F", "TIMD4", "CD163"],
    "MetMac":   ["GPNMB", "HS3ST2", "LPL", "FABP5", "TREM2", "CD9"],
    "TransMac": ["CXCL10", "CXCL9", "CXCL11", "STAB1", "MGAT4A"],
    "preMac":   ["PCNX2", "ADAM28", "RUNX2", "PLTP", "CCL18"],
    "Monocyte": ["VCAN", "FCN1", "S100A8", "S100A9", "CD14"],
}

# Score columns produced by the loader (one per macrophage lineage label).
# Per task spec the loader populates these in .obs of the master h5ad.
SCORE_COLS = ["score_KC", "score_MetMac", "score_TransMac",
              "score_preMac", "score_Monocyte"]

# Subtypes we will run DE for (only if ≥ MIN_CELLS_PER_GROUP per disease group).
DE_TARGET_SUBTYPES = ["MetMac", "TransMac", "KC", "preMac", "Monocyte"]

# Required minimum cells per disease group for DE.
MIN_CELLS_PER_GROUP = 50

# GPNMB threshold for MetMac tie-breaker — if score_MetMac > score_KC AND
# mean log1p(GPNMB) in the cluster is above this fraction of the overall
# 80th percentile, force-MetMac.
GPNMB_TIEBREAK_QUANTILE = 0.80


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def assign_disease_stage(adata: ad.AnnData) -> ad.AnnData:
    """Map .obs.sample_id → MASH / no_MASH using the slide-level table."""
    stage_map = {sid: meta["stage"] for sid, meta in SAMPLE_DISEASE.items()}
    adata.obs["disease_stage"] = adata.obs["sample_id"].map(stage_map).astype("category")
    if adata.obs["disease_stage"].isna().any():
        bad = adata.obs.loc[adata.obs["disease_stage"].isna(), "sample_id"].unique()
        sys.exit(f"ERROR: unmapped sample_ids: {bad}")
    _log("Disease-stage assignment (slide-level):")
    for sid, meta in SAMPLE_DISEASE.items():
        n = int((adata.obs["sample_id"] == sid).sum())
        _log(f"  {sid:>9}  ({meta['sample_disease_str']:>20}) → {meta['stage']:>8}  "
             f"({n:,} cells)")
    # F137: the no_MASH comparator (Leuven_2) is ~75% MASLD tissue. Make the
    # confound explicit at runtime so nobody reads the DE as MASH-vs-control.
    _log("  WARNING (F137): 'no_MASH' = single slide Leuven_2, which is "
         "~75% MASLD tissue (1 Normal + 1 MASL + 2 F3). This is a NOISY "
         "MASH-ENRICHMENT contrast, NOT MASH-vs-control. Treat per-subtype DE "
         "p-values as descriptive; per-slide MetMac proportions are the primary "
         "caveated evidence.")
    return adata


def subcluster_macrophages(adata_mac: ad.AnnData,
                           n_pcs: int = 50,
                           n_neighbors: int = 15,
                           resolution: float = 1.5,
                           hvg_n: int = 500,
                           run_umap: bool = False) -> ad.AnnData:
    """Re-cluster the macrophage subset at finer resolution.

    Steps mirror Govaere 2026 paper Methods (page 14): per-subset normalize,
    HVG within subset, scale, PCA, k-NN, Louvain/Leiden at res=1-2, optional
    UMAP. We use Leiden (de-facto standard, deterministic with seed).
    """
    _log(f"  Subset shape: {adata_mac.shape} (cells × genes)")
    _log("  Normalising target_sum=1e4 + log1p (within macrophage subset)")
    # Restore raw counts to .X for re-normalisation (loader stashed raw at .X).
    # Defensive: if .X is already log-normalised (rare), re-derive from layers.
    if "counts" in adata_mac.layers:
        adata_mac.X = adata_mac.layers["counts"].copy()
    adata_mac.layers["counts_raw"] = adata_mac.X.copy()
    sc.pp.normalize_total(adata_mac, target_sum=1e4)
    sc.pp.log1p(adata_mac)
    adata_mac.layers["lognorm"] = adata_mac.X.copy()

    n_target = adata_mac.n_vars
    n_use = min(hvg_n, n_target)
    _log(f"  Selecting up to {n_use} HVGs (from {n_target} panel genes)")
    sc.pp.highly_variable_genes(adata_mac, n_top_genes=n_use, flavor="seurat")
    n_hvg = int(adata_mac.var["highly_variable"].sum())
    _log(f"  HVGs flagged: {n_hvg}")

    _log("  Scaling (max_value=10, zero_center=True)")
    sc.pp.scale(adata_mac, max_value=10, zero_center=True)
    n_pcs_use = min(n_pcs, adata_mac.n_vars - 1, adata_mac.n_obs - 1)
    _log(f"  PCA n_comps={n_pcs_use}")
    sc.tl.pca(adata_mac, n_comps=n_pcs_use)
    _log(f"  Neighbors n_pcs={n_pcs_use}, k={n_neighbors}")
    sc.pp.neighbors(adata_mac, n_pcs=n_pcs_use, n_neighbors=n_neighbors)
    _log(f"  Leiden resolution={resolution} (paper range 1-2)")
    sc.tl.leiden(adata_mac, resolution=resolution, key_added="mac_leiden",
                 random_state=0)
    n_clusters = adata_mac.obs["mac_leiden"].nunique()
    _log(f"  Found {n_clusters} sub-clusters")
    if run_umap:
        _log("  UMAP")
        sc.tl.umap(adata_mac, random_state=0)
    return adata_mac


def assign_subtype_labels(adata_mac: ad.AnnData) -> ad.AnnData:
    """Composite scoring with absolute marker thresholds.

    BACKGROUND (why pure score_genes argmax fails on this panel):
      The loader-computed score_genes values use a random-control comparison
      that is biased by the CosMx panel composition (968 liver-enriched
      genes). score_KC (MARCO/CD5L/VSIG4/CLEC4F/TIMD4) was the global argmax
      in every sub-cluster because these markers are panel-rich and
      shared across the macrophage lineage. score_MetMac (GPNMB/HS3ST2/LPL/
      FABP5/TREM2/CD9) never exceeded score_KC even in clusters with
      qualitatively high lipid-associated marker expression.

    REVISED STRATEGY (marker-rank based with priority cascade):
      Compute cluster-mean log1p of canonical markers per subtype. Then:
        1. MONOCYTE: cluster argmax of score_Monocyte (VCAN/FCN1/CD14 path
           is panel-distinctive and well separated; here cluster 13 had
           score_Monocyte = +2.157, > 4σ above any other cluster).
        2. METMAC: clusters whose mean GPNMB log1p ranks in the TOP-K
           (K = max(1, n_clusters // 6) ≈ 4 of 23) AND mean GPNMB >
           an absolute panel-rank threshold (≥ 75th-percentile of GPNMB
           across all sub-clusters). This isolates the lipid-associated
           macrophages even when score_KC dominates score_MetMac on
           absolute terms.
        3. TRANSMAC: cluster argmax of (CXCL10 + STAB1 + MGAT4A) /3 mean,
           gated by mean expression > 75th percentile of (CXCL10 + STAB1)/2
           across clusters (TransMac is hepatic-stellate/sinusoidal-niched).
        4. PREMAC: clusters where mean (PCNX2 + ADAM28 + RUNX2)/3 ranks top,
           OR fall-through assignment from score_preMac argmax when raw
           markers are too sparse.
        5. KC (residual): everything else — these are CD5L/MARCO/VSIG4-high
           resident Kupffer cells.

      This priority cascade matches the Govaere 2026 Fig 1d/Sup Table 11
      operational definition: MetMac is the GPNMB-high, lipid-loaded
      sub-cluster, regardless of where it lands in score-space.

    Decision audit stored in .uns["mac_subtype_decision"] as a flat
    DataFrame round-tripping through h5ad.
    """
    # Cluster-mean log1p of marker genes (from .layers["lognorm"])
    avail_markers: Dict[str, List[str]] = {}
    all_markers: List[str] = []
    for subtype, markers in SUBTYPE_MARKER_PANELS.items():
        present = [m for m in markers if m in adata_mac.var_names]
        avail_markers[subtype] = present
        all_markers.extend(present)
    all_markers = sorted(set(all_markers))
    _log(f"  Marker genes available on panel: {len(all_markers)}/"
         f"{sum(len(v) for v in SUBTYPE_MARKER_PANELS.values())}")
    _log(f"  Per-subtype marker availability: "
         f"{ {st: avail_markers[st] for st in avail_markers} }")

    gene_idx = [adata_mac.var_names.get_loc(g) for g in all_markers]
    ln = adata_mac.layers["lognorm"]
    if sparse.issparse(ln):
        marker_mat = ln[:, gene_idx].toarray()
    else:
        marker_mat = ln[:, gene_idx]
    marker_df = pd.DataFrame(marker_mat, columns=all_markers,
                             index=adata_mac.obs_names)
    marker_df["mac_leiden"] = adata_mac.obs["mac_leiden"].values
    cluster_marker = marker_df.groupby("mac_leiden", observed=True).mean()

    # Score columns (loader-computed; biased on this panel — used as a
    # secondary signal only).
    score_cols_use = [c for c in SCORE_COLS if c in adata_mac.obs.columns]
    cluster_scores = (adata_mac.obs
                      .groupby("mac_leiden", observed=True)[score_cols_use]
                      .mean()) if score_cols_use else pd.DataFrame()

    _log("  Per-sub-cluster mean scores (loader-computed):")
    for c in cluster_scores.index:
        row = "  ".join(f"{col.replace('score_', '')}={cluster_scores.loc[c, col]:+.3f}"
                        for col in score_cols_use)
        nc = int((adata_mac.obs["mac_leiden"] == c).sum())
        _log(f"    cluster {c:>3} ({nc:>6,} cells): {row}")

    _log("  Per-sub-cluster mean canonical marker log1p:")
    show_cols = [g for g in ["GPNMB", "LPL", "FABP5", "HS3ST2", "TREM2", "CD9",
                              "MARCO", "CD5L", "CD163",
                              "CXCL10", "CXCL9", "STAB1", "MGAT4A",
                              "VCAN", "FCN1", "CD14",
                              "PCNX2", "ADAM28", "RUNX2"]
                  if g in cluster_marker.columns]
    for c in cluster_marker.index:
        row = "  ".join(f"{g}={cluster_marker.loc[c, g]:.2f}" for g in show_cols)
        _log(f"    cluster {c:>3}: {row}")

    # ── Composite ranking ────────────────────────────────────────────────
    # MetMac composite = mean log1p of (GPNMB, LPL, FABP5, HS3ST2) — restricted
    # to present markers; absolute panel signal.
    mm_genes = [g for g in ["GPNMB", "LPL", "FABP5", "HS3ST2"]
                if g in cluster_marker.columns]
    tm_genes = [g for g in ["CXCL10", "CXCL9", "STAB1", "MGAT4A"]
                if g in cluster_marker.columns]
    pm_genes = [g for g in ["PCNX2", "ADAM28", "RUNX2"]
                if g in cluster_marker.columns]
    mo_genes = [g for g in ["VCAN", "FCN1", "CD14"]
                if g in cluster_marker.columns]
    kc_genes = [g for g in ["MARCO", "CD5L", "VSIG4", "TIMD4"]
                if g in cluster_marker.columns]
    # CD163 is shared between KC + MetMac in the literature; exclude from
    # composite_KC to avoid bias.

    cluster_marker["composite_MetMac"]   = cluster_marker[mm_genes].mean(axis=1) if mm_genes else np.nan
    cluster_marker["composite_TransMac"] = cluster_marker[tm_genes].mean(axis=1) if tm_genes else np.nan
    cluster_marker["composite_preMac"]   = cluster_marker[pm_genes].mean(axis=1) if pm_genes else np.nan
    cluster_marker["composite_Monocyte"] = cluster_marker[mo_genes].mean(axis=1) if mo_genes else np.nan
    cluster_marker["composite_KC"]       = cluster_marker[kc_genes].mean(axis=1) if kc_genes else np.nan
    # MetMac-discrimination index: contrast lipid-axis vs resident-KC axis.
    # High MetMac-discrim = clusters with high GPNMB/LPL/FABP5 AND low MARCO/CD5L.
    cluster_marker["metmac_discrim"] = (cluster_marker["composite_MetMac"].fillna(0.0)
                                        - cluster_marker["composite_KC"].fillna(0.0))

    # Quantile thresholds (75th percentile across all sub-clusters)
    metmac_q     = (cluster_marker["composite_MetMac"].quantile(0.75)
                    if mm_genes else np.inf)
    metmac_disc_q = cluster_marker["metmac_discrim"].quantile(0.65)
    transmac_q   = (cluster_marker["composite_TransMac"].quantile(0.75)
                    if tm_genes else np.inf)
    monocyte_q   = (cluster_marker["composite_Monocyte"].quantile(0.85)
                    if mo_genes else np.inf)
    premac_q     = (cluster_marker["composite_preMac"].quantile(0.75)
                    if pm_genes else np.inf)
    kc_q         = (cluster_marker["composite_KC"].quantile(0.75)
                    if kc_genes else np.inf)

    _log(f"  composite_MetMac   q75 = {metmac_q:.3f}  (genes used: {mm_genes})")
    _log(f"  metmac_discrim     q65 = {metmac_disc_q:.3f}  (MetMac − KC)")
    _log(f"  composite_TransMac q75 = {transmac_q:.3f}  (genes used: {tm_genes})")
    _log(f"  composite_preMac   q75 = {premac_q:.3f}  (genes used: {pm_genes})")
    _log(f"  composite_Monocyte q85 = {monocyte_q:.3f}  (genes used: {mo_genes})")
    _log(f"  composite_KC       q75 = {kc_q:.3f}  (genes used: {kc_genes})")

    # ── Assignment cascade ───────────────────────────────────────────────
    # Priority (Govaere 2026-aligned, score-anchored where the loader's
    # score is unambiguous — Monocyte and MetMac):
    #   1. Monocyte FIRST: rare type, loader score_Monocyte is unambiguous
    #      (cluster 13 here has +2.157, an order of magnitude above any other).
    #      Must beat KC priority because cluster 13 also has elevated MARCO/CD5L
    #      but is genuinely monocyte (VCAN/CD14 are dominant).
    #   2. MetMac (paper headline): high composite_MetMac AND metmac_discrim
    #      (rules out clusters with both KC + MetMac markers — those are
    #      resident KCs that happen to scavenge lipids).
    #   3. KC (resident, canonical): high MARCO/CD5L/VSIG4 AND NOT
    #      MetMac-discriminated.
    #   4. TransMac: CXCL10/CXCL9 axis.
    #   5. preMac (only if PCNX2/ADAM28/RUNX2 markers present on panel).
    #   6. KC residual (default for clusters that don't trigger any rule —
    #      these are intermediate macrophages, all lineage features moderate).
    decisions: Dict[str, Dict[str, object]] = {}
    final_label = pd.Series(index=cluster_marker.index, dtype="object")
    for c in cluster_marker.index:
        nc = int((adata_mac.obs["mac_leiden"] == c).sum())
        cmm = float(cluster_marker.loc[c, "composite_MetMac"]) if mm_genes else np.nan
        ctm = float(cluster_marker.loc[c, "composite_TransMac"]) if tm_genes else np.nan
        cpm = float(cluster_marker.loc[c, "composite_preMac"]) if pm_genes else np.nan
        cmo = float(cluster_marker.loc[c, "composite_Monocyte"]) if mo_genes else np.nan
        ckc = float(cluster_marker.loc[c, "composite_KC"]) if kc_genes else np.nan
        cd = float(cluster_marker.loc[c, "metmac_discrim"])
        decisions[c] = {
            "n_cells": nc,
            "composite_MetMac":   cmm,
            "composite_TransMac": ctm,
            "composite_preMac":   cpm,
            "composite_Monocyte": cmo,
            "composite_KC":       ckc,
            "metmac_discrim":     cd,
            "score_MetMac_loader":   float(cluster_scores.loc[c, "score_MetMac"])
                                       if "score_MetMac" in cluster_scores.columns else np.nan,
            "score_KC_loader":       float(cluster_scores.loc[c, "score_KC"])
                                       if "score_KC" in cluster_scores.columns else np.nan,
            "score_Monocyte_loader": float(cluster_scores.loc[c, "score_Monocyte"])
                                       if "score_Monocyte" in cluster_scores.columns else np.nan,
        }

        # 1. KC (resident): high MARCO/CD5L/VSIG4 (canonical resident KC)
        if (kc_genes and ckc >= kc_q
            and cd < metmac_disc_q):  # not MetMac-discriminated
            final_label[c] = "KC"
            decisions[c]["rule"] = f"KC (composite_KC>={kc_q:.3f} & not MetMac-discrim)"
            continue
        # 2. Monocyte
        if (mo_genes and cmo >= monocyte_q
            and "score_Monocyte" in cluster_scores.columns
            and cluster_scores.loc[c, "score_Monocyte"] > 1.0):
            final_label[c] = "Monocyte"
            decisions[c]["rule"] = "Monocyte (composite>=q85 & loader score>1.0)"
            continue
        # 3. MetMac — composite high AND MetMac-discrim
        if (mm_genes and cmm >= metmac_q
            and cd >= metmac_disc_q):
            final_label[c] = "MetMac"
            decisions[c]["rule"] = "MetMac (composite_MetMac>=q75 & discrim>=q65)"
            continue
        # 4. TransMac — CXCL10/CXCL9 axis
        if tm_genes and ctm >= transmac_q:
            final_label[c] = "TransMac"
            decisions[c]["rule"] = "TransMac (composite_TransMac>=q75)"
            continue
        # 5. preMac — only if markers actually exist on panel
        if pm_genes and cpm >= premac_q:
            final_label[c] = "preMac"
            decisions[c]["rule"] = "preMac (composite_preMac>=q75)"
            continue
        # 6. KC residual default
        final_label[c] = "KC"
        decisions[c]["rule"] = "KC (residual)"

    _log("  Final cluster → subtype assignment:")
    for c, lbl in final_label.items():
        nc = decisions[c]["n_cells"]
        rule = decisions[c]["rule"]
        _log(f"    cluster {c:>3} → {lbl:>9}   ({nc:>6,} cells)  [{rule}]")

    adata_mac.obs["mac_subtype"] = (adata_mac.obs["mac_leiden"]
                                    .map(final_label)
                                    .astype("category"))
    # Serialize decisions for h5ad — flat DataFrame
    dec_df = pd.DataFrame.from_dict(decisions, orient="index").reset_index().rename(
        columns={"index": "mac_leiden"})
    adata_mac.uns["mac_subtype_decision"] = dec_df.to_dict(orient="list")
    adata_mac.uns["mac_subtype_thresholds"] = {
        "metmac_q75": float(metmac_q),
        "transmac_q75": float(transmac_q),
        "premac_q75": float(premac_q),
        "monocyte_q85": float(monocyte_q),
    }
    return adata_mac


def per_subtype_de(adata_mac: ad.AnnData, subtype: str,
                   out_dir: pathlib.Path) -> Dict[str, object]:
    """Wilcoxon DE within one sub-type: MASH vs no_MASH."""
    mask = adata_mac.obs["mac_subtype"] == subtype
    n_total = int(mask.sum())
    if n_total == 0:
        _log(f"  [{subtype}] SKIP — 0 cells")
        return {"cell_type": subtype, "skipped": True, "reason": "absent",
                "n_cells_mash": 0, "n_cells_nomash": 0, "n_sig_padj05": 0}
    sub = adata_mac[mask].copy()
    n_mash = int((sub.obs["disease_stage"] == "MASH").sum())
    n_nomash = int((sub.obs["disease_stage"] == "no_MASH").sum())
    _log(f"  [{subtype}] subset {n_total:,} cells "
         f"(MASH={n_mash:,}; no_MASH={n_nomash:,})")
    if n_mash < MIN_CELLS_PER_GROUP or n_nomash < MIN_CELLS_PER_GROUP:
        _log(f"  [{subtype}] SKIP — group below {MIN_CELLS_PER_GROUP}")
        return {"cell_type": subtype, "skipped": True, "reason": "low_n",
                "n_cells_mash": n_mash, "n_cells_nomash": n_nomash,
                "n_sig_padj05": 0}

    # DE on log-normalised counts (loader convention + 42b convention)
    if "lognorm" in sub.layers:
        sub.X = sub.layers["lognorm"]
    sub.obs["disease_stage"] = (sub.obs["disease_stage"].astype("category")
                                  .cat.remove_unused_categories())
    sc.tl.rank_genes_groups(
        sub,
        groupby="disease_stage",
        groups=["MASH"],
        reference="no_MASH",
        method="wilcoxon",
        pts=True,
        use_raw=False,
        rankby_abs=False,
    )

    rg = sub.uns["rank_genes_groups"]
    de_df = pd.DataFrame({
        "gene":          rg["names"]["MASH"],
        "logfoldchange": rg["logfoldchanges"]["MASH"],
        "pval":          rg["pvals"]["MASH"],
        # NOTE (F138/F137): pval/pval_adj are CELL-LEVEL Wilcoxon p-values. The
        # disease label is slide-level (only 4 slides; the 'no_MASH' arm is a SINGLE
        # slide, Leuven_2). Cells within a slide are not independent, so these
        # p-values are pseudoreplicated/anti-conservative — NOT donor-level evidence.
        # Worse (F137), the no_MASH comparator slide (Leuven_2) is ~75% MASLD tissue,
        # so this is a noisy MASH-ENRICHMENT contrast, not MASH-vs-control. Column
        # names are kept for the atlas; the honest evidence is the slide-level
        # direction (slide_logfc) and the per-slide MetMac proportions, NOT these p's.
        "pval_adj":      rg["pvals_adj"]["MASH"],
        "score":         rg["scores"]["MASH"],
    })
    if "pts" in rg:
        pts = rg["pts"]
        de_df["pct_nz_MASH"]    = pts["MASH"].reindex(de_df["gene"]).values
        de_df["pct_nz_no_MASH"] = pts["no_MASH"].reindex(de_df["gene"]).values
    de_df["n_cells_MASH"]    = n_mash
    de_df["n_cells_no_MASH"] = n_nomash
    de_df["cell_type"]       = subtype

    # ── Slide-level pseudobulk direction (F138) ──────────────────────────────
    # Aggregate to one lognorm-mean profile per slide, then contrast MASH vs
    # no_MASH at the slide level. slide_logfc = mean(MASH slides) − mean(no_MASH
    # slide) on the log1p-CPM scale (a log-scale difference; agrees in sign with
    # scanpy's logFC but on a different scale). With 3 MASH vs 1 no_MASH slide this
    # supports only a direction, not a p-value.
    try:
        pb = pseudobulk_by_donor(
            sub, donor_col="sample_id", layer="lognorm",
            agg="mean", obs_cols=["disease_stage"])
        mash_slides = pb.index[pb["disease_stage"].astype(str) == "MASH"]
        nomash_slides = pb.index[pb["disease_stage"].astype(str) == "no_MASH"]
        gene_cols = [c for c in pb.columns if c != "disease_stage"]
        mash_mean = pb.loc[mash_slides, gene_cols].mean(axis=0)
        nomash_mean = pb.loc[nomash_slides, gene_cols].mean(axis=0)
        slide_logfc = (mash_mean - nomash_mean)
        de_df["slide_logfc"] = slide_logfc.reindex(de_df["gene"]).values
        de_df["n_slides_MASH"]    = int(len(mash_slides))
        de_df["n_slides_no_MASH"] = int(len(nomash_slides))
        de_df["slide_direction_concordant"] = (
            np.sign(de_df["slide_logfc"]) == np.sign(de_df["logfoldchange"]))
        _log(f"  [{subtype}] slide-level pseudobulk: "
             f"{len(mash_slides)} MASH vs {len(nomash_slides)} no_MASH slide(s); "
             f"direction-concordant genes = "
             f"{int(de_df['slide_direction_concordant'].sum())}/{de_df.shape[0]}")
    except Exception as e:
        _log(f"  [{subtype}] WARN — slide pseudobulk failed ({str(e)[:100]}); "
             f"slide_logfc set NaN")
        de_df["slide_logfc"] = np.nan
        de_df["n_slides_MASH"]    = np.nan
        de_df["n_slides_no_MASH"] = np.nan
        de_df["slide_direction_concordant"] = np.nan

    n_sig = int((de_df["pval_adj"] < 0.05).sum())
    n_sig_up = int(((de_df["pval_adj"] < 0.05) & (de_df["logfoldchange"] > 0)).sum())
    n_sig_down = int(((de_df["pval_adj"] < 0.05) & (de_df["logfoldchange"] < 0)).sum())
    _log(f"  [{subtype}] n_sig padj<0.05 = {n_sig} "
         f"(up={n_sig_up}, down={n_sig_down})")

    # File naming: MetMac/TransMac/preMac/Monocyte → standard;
    # post-subcluster KC residual gets a distinct suffix so 42b's output is preserved.
    if subtype == "KC":
        fname = "cosmx_de_KC_post_subcluster_MASH_vs_noMASH.csv"
    else:
        fname = f"cosmx_de_{subtype}_MASH_vs_noMASH.csv"
    out_path = out_dir / fname
    de_df.to_csv(out_path, index=False)
    _log(f"  [{subtype}] → {out_path.name} ({out_path.stat().st_size/1e3:.1f} KB)")

    return {
        "cell_type": subtype,
        "skipped": False,
        "n_cells_mash": n_mash,
        "n_cells_nomash": n_nomash,
        "n_genes_tested": int(de_df.shape[0]),
        "n_sig_padj05": n_sig,
        "n_sig_padj05_up": n_sig_up,
        "n_sig_padj05_down": n_sig_down,
        "out_path": str(out_path),
        "de_df": de_df,
    }


def write_composition_table(adata_mac: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    """Per-sample × per-subtype counts + proportion within macrophage bin."""
    comp = (adata_mac.obs
            .groupby(["sample_id", "mac_subtype"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index())
    sample_totals = comp.groupby("sample_id", observed=True)["n_cells"].transform("sum")
    comp["proportion_within_mac"] = comp["n_cells"] / sample_totals
    stage_map = {sid: meta["stage"] for sid, meta in SAMPLE_DISEASE.items()}
    comp["disease_stage"] = comp["sample_id"].map(stage_map)
    comp["sample_disease_str"] = comp["sample_id"].map(
        {sid: meta["sample_disease_str"] for sid, meta in SAMPLE_DISEASE.items()})
    out_path = out_dir / "cosmx_mac_subtype_composition.tsv"
    comp.to_csv(out_path, sep="\t", index=False)
    _log(f"Wrote composition → {out_path.name} ({comp.shape[0]} rows)")
    return comp


def write_marker_summary(adata_mac: ad.AnnData, out_dir: pathlib.Path,
                         top_n: int = 10) -> pd.DataFrame:
    """For each subtype, rank panel genes by mean log-norm expression and emit top-N."""
    if "lognorm" not in adata_mac.layers:
        _log("  WARNING: .layers['lognorm'] missing; skipping marker summary")
        return pd.DataFrame()
    rows = []
    subtypes = adata_mac.obs["mac_subtype"].cat.categories
    ln = adata_mac.layers["lognorm"]
    # Dense per-subtype mean (panel is 968 genes → cheap)
    for st in subtypes:
        m = (adata_mac.obs["mac_subtype"] == st).values
        if m.sum() == 0:
            continue
        ln_sub = ln[m]
        if sparse.issparse(ln_sub):
            mean_ln = np.asarray(ln_sub.mean(axis=0)).ravel()
        else:
            mean_ln = ln_sub.mean(axis=0)
        # Rank in descending order
        order = np.argsort(mean_ln)[::-1][:top_n]
        for rank, j in enumerate(order, start=1):
            rows.append({
                "subtype": st,
                "rank":    rank,
                "gene":    adata_mac.var_names[j],
                "mean_log1p_cpm": float(mean_ln[j]),
                "n_cells": int(m.sum()),
            })
    out_df = pd.DataFrame(rows)
    out_path = out_dir / "cosmx_mac_subtype_markers.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    _log(f"Wrote marker summary → {out_path.name} "
         f"({out_df.shape[0]} rows, top {top_n} per subtype)")
    return out_df


def metmac_vs_other_de(adata_mac: ad.AnnData, out_dir: pathlib.Path) -> Dict[str, object]:
    """Auxiliary DE: MetMac (all stages) vs all other macrophage subtypes.

    This is the contrast that surfaces the paper's MetMac-defining markers
    (GPNMB / LPL / FABP5 / HS3ST2). The primary MetMac-MASH-vs-noMASH
    contrast inside the MetMac subset cannot recover defining markers because
    every cell in MetMac shares the GPNMB-high phenotype.
    """
    if "mac_subtype" not in adata_mac.obs.columns:
        return {"skipped": True, "reason": "no mac_subtype"}
    if "MetMac" not in adata_mac.obs["mac_subtype"].cat.categories:
        return {"skipped": True, "reason": "MetMac not present"}
    sub = adata_mac.copy()
    sub.obs["mm_vs_other"] = pd.Categorical(
        np.where(sub.obs["mac_subtype"] == "MetMac", "MetMac", "Other"))
    n_mm = int((sub.obs["mm_vs_other"] == "MetMac").sum())
    n_oth = int((sub.obs["mm_vs_other"] == "Other").sum())
    _log(f"  [MetMac vs Other-mac] MetMac n={n_mm:,}; Other n={n_oth:,}")
    if n_mm < MIN_CELLS_PER_GROUP or n_oth < MIN_CELLS_PER_GROUP:
        return {"skipped": True, "reason": "low_n",
                "n_metmac": n_mm, "n_other": n_oth}
    if "lognorm" in sub.layers:
        sub.X = sub.layers["lognorm"]
    sub.obs["mm_vs_other"] = sub.obs["mm_vs_other"].cat.remove_unused_categories()
    sc.tl.rank_genes_groups(
        sub, groupby="mm_vs_other",
        groups=["MetMac"], reference="Other",
        method="wilcoxon", pts=True, use_raw=False,
    )
    rg = sub.uns["rank_genes_groups"]
    de_df = pd.DataFrame({
        "gene":          rg["names"]["MetMac"],
        "logfoldchange": rg["logfoldchanges"]["MetMac"],
        "pval":          rg["pvals"]["MetMac"],
        "pval_adj":      rg["pvals_adj"]["MetMac"],
        "score":         rg["scores"]["MetMac"],
    })
    if "pts" in rg:
        pts = rg["pts"]
        de_df["pct_nz_MetMac"] = pts["MetMac"].reindex(de_df["gene"]).values
        de_df["pct_nz_Other"]  = pts["Other"].reindex(de_df["gene"]).values
    de_df["n_cells_MetMac"] = n_mm
    de_df["n_cells_Other"]  = n_oth
    de_df["contrast"] = "MetMac_vs_other_mac"
    out_path = out_dir / "cosmx_de_MetMac_vs_other_mac.csv"
    de_df.to_csv(out_path, index=False)
    _log(f"  [MetMac vs Other-mac] → {out_path.name}")
    # Top markers
    up = de_df[de_df["logfoldchange"] > 0].sort_values("score", ascending=False)
    canon = ["GPNMB", "LPL", "FABP5", "HS3ST2", "TREM2", "CD9"]
    top50 = up.head(50)["gene"].tolist()
    found = [g for g in canon if g in top50]
    missing = [g for g in canon if g not in top50]
    _log(f"  [MetMac vs Other-mac] top 20 markers:")
    for _, row in up.head(20).iterrows():
        _log(f"    {row['gene']:>10}  lfc={row['logfoldchange']:+.3f}  "
             f"padj={row['pval_adj']:.2e}  pct_MetMac="
             f"{row.get('pct_nz_MetMac', float('nan')):.3f}")
    _log(f"  Canonical MetMac markers in top 50: {found}")
    _log(f"  Canonical MetMac markers absent:    {missing}")
    return {"skipped": False, "out_path": str(out_path),
            "n_metmac": n_mm, "n_other": n_oth,
            "canon_found_in_top50": found,
            "canon_missing_in_top50": missing,
            "de_df": de_df}


def report_metmac_validation(de_results: List[Dict], comp: pd.DataFrame) -> None:
    """Print MetMac validation: top 10/20 MASH-up genes and proportion shift."""
    _log("=" * 70)
    _log("VALIDATION: MetMac MASH-up axis")
    _log("=" * 70)
    metmac = next((r for r in de_results
                   if r["cell_type"] == "MetMac" and not r.get("skipped")), None)
    if metmac is None:
        _log("MetMac DE not available (subtype absent or low_n)")
    else:
        de_df = metmac["de_df"].copy()
        de_df_up = de_df[de_df["logfoldchange"] > 0].sort_values("score",
                                                                  ascending=False)
        _log("Top 20 MetMac MASH-up genes by score (paper Sup Table 11 anchors):")
        for _, row in de_df_up.head(20).iterrows():
            _log(f"  {row['gene']:>10}  lfc={row['logfoldchange']:+.3f}  "
                 f"padj={row['pval_adj']:.2e}  "
                 f"pct_MASH={row.get('pct_nz_MASH', float('nan')):.3f}")
        # MetMac canonical markers
        canon = ["GPNMB", "LPL", "FABP5", "HS3ST2", "TREM2", "CD9"]
        top50 = de_df_up.head(50)["gene"].tolist()
        found = [g for g in canon if g in top50]
        miss = [g for g in canon if g not in top50]
        _log(f"  Canonical MetMac markers in top 50: {found}")
        _log(f"  Canonical MetMac markers absent:    {miss}")

    # Per-sample MetMac proportion
    if comp is not None and "mac_subtype" in comp.columns:
        mm_rows = comp[comp["mac_subtype"] == "MetMac"]
        if mm_rows.empty:
            _log("MetMac proportion: no MetMac assignments in composition table")
        else:
            _log("Per-sample MetMac proportion (within macrophage bin):")
            for _, row in mm_rows.iterrows():
                _log(f"  {row['sample_id']:>9}  ({row['disease_stage']:>8})  "
                     f"MetMac n={int(row['n_cells']):,}  "
                     f"prop={row['proportion_within_mac']:.3f}")
            mash_props  = mm_rows.loc[mm_rows["disease_stage"] == "MASH",
                                       "proportion_within_mac"]
            nomash_props = mm_rows.loc[mm_rows["disease_stage"] == "no_MASH",
                                        "proportion_within_mac"]
            mash_mean   = mash_props.mean()  if len(mash_props)  else float("nan")
            nomash_mean = nomash_props.mean() if len(nomash_props) else float("nan")
            _log(f"  Mean MetMac proportion — MASH slides:    {mash_mean:.3f}")
            _log(f"  Mean MetMac proportion — no_MASH slides: {nomash_mean:.3f}")
            _log("  (Paper claim: MetMac 5.7% → 12% in MASH; gene-cluster proportion not panel-wide)")


def main():
    p = argparse.ArgumentParser(description="Sub-cluster CosMx KC macrophage bin")
    p.add_argument("--adata-in", type=str, default=str(ADATA_IN))
    p.add_argument("--adata-out", type=str, default=str(ADATA_OUT))
    p.add_argument("--out-dir", type=str, default=str(OUT_DIR))
    p.add_argument("--leiden-resolution", type=float, default=1.5)
    p.add_argument("--hvg-n", type=int, default=500)
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--run-umap", action="store_true",
                   help="Compute UMAP after Leiden (off by default for memory).")
    p.add_argument("--from-clustered", action="store_true",
                   help="Skip subcluster step and re-read the existing "
                        "cosmx_govaere2026_macsubclustered.h5ad (re-uses "
                        "Leiden + lognorm; useful when re-running the "
                        "subtype-assignment logic only).")
    args = p.parse_args()

    adata_in_path  = pathlib.Path(args.adata_in)
    adata_out_path = pathlib.Path(args.adata_out)
    out_dir        = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    adata_out_path.parent.mkdir(parents=True, exist_ok=True)

    _log("=" * 70)
    _log("41b_govaere2026_cosmx_mac_subcluster.py")
    _log(f"  AnnData in  : {adata_in_path}")
    _log(f"  AnnData out : {adata_out_path}")
    _log(f"  Out dir     : {out_dir}")
    _log(f"  Leiden res  : {args.leiden_resolution} (paper 1-2)")
    _log(f"  HVG n       : {args.hvg_n}")
    _log(f"  n_pcs       : {args.n_pcs}")
    _log(f"  n_neighbors : {args.n_neighbors}")
    _log("=" * 70)

    if args.from_clustered and adata_out_path.exists():
        _log(f"--from-clustered: re-reading "
             f"{adata_out_path.name} ({adata_out_path.stat().st_size/1e6:.1f} MB) — "
             f"skipping HVG/PCA/Leiden")
        adata_mac = sc.read_h5ad(adata_out_path)
        _log(f"  shape: {adata_mac.shape}")
        _log(f"  .obs columns: {list(adata_mac.obs.columns)}")
        _log(f"  .layers: {list(adata_mac.layers.keys())}")
        if "mac_leiden" not in adata_mac.obs.columns:
            sys.exit("ERROR: --from-clustered requires .obs.mac_leiden in the existing h5ad")
        if "lognorm" not in adata_mac.layers:
            sys.exit("ERROR: --from-clustered requires .layers['lognorm'] in the existing h5ad")
        # Make sure disease_stage is present (may already be in the h5ad)
        if "disease_stage" not in adata_mac.obs.columns:
            adata_mac = assign_disease_stage(adata_mac)
        # Drop any previously assigned mac_subtype so re-assignment is clean
        if "mac_subtype" in adata_mac.obs.columns:
            adata_mac.obs = adata_mac.obs.drop(columns=["mac_subtype"])
    else:
        _log(f"Reading master AnnData ({adata_in_path.stat().st_size/1e6:.1f} MB)")
        adata = sc.read_h5ad(adata_in_path)
        _log(f"  shape: {adata.shape}")
        _log(f"  .obs columns: {list(adata.obs.columns)}")
        _log(f"  .layers: {list(adata.layers.keys())}")
        _log(f"  cell_type tally: "
             f"{adata.obs['cell_type'].value_counts().to_dict()}")

        # Confirm score columns exist
        found_scores = [c for c in adata.obs.columns if c.startswith("score_")]
        _log(f"  score columns present: {found_scores}")

        # Subset to KC bin
        n_kc = int((adata.obs["cell_type"] == "KC").sum())
        _log(f"Subsetting to cell_type == 'KC'  ({n_kc:,} cells)")
        adata_mac = adata[adata.obs["cell_type"] == "KC"].copy()
        # Free master h5ad memory immediately — we will not write to it.
        del adata
        gc.collect()

        # Slide-level disease stage
        adata_mac = assign_disease_stage(adata_mac)

        # Sub-cluster (re-normalise within subset; HVG; PCA; neighbors; Leiden)
        adata_mac = subcluster_macrophages(
            adata_mac,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            resolution=args.leiden_resolution,
            hvg_n=args.hvg_n,
            run_umap=args.run_umap,
        )

    # Assign subtype labels via score-based winner-take-all + GPNMB tie-breaker
    adata_mac = assign_subtype_labels(adata_mac)

    # Show subtype tally
    tally = adata_mac.obs["mac_subtype"].value_counts().to_dict()
    _log(f"Final mac_subtype tally: {tally}")
    _log(f"  Number of distinct subtypes: {len(tally)}")

    # Composition table
    comp = write_composition_table(adata_mac, out_dir)

    # Marker summary
    marker_df = write_marker_summary(adata_mac, out_dir, top_n=10)

    # Per-subtype DE (MetMac is the KEY OUTPUT)
    de_results: List[Dict] = []
    for subtype in DE_TARGET_SUBTYPES:
        res = per_subtype_de(adata_mac, subtype, out_dir)
        de_results.append(res)

    # Auxiliary DE: MetMac vs other macrophage subtypes (surfaces defining
    # markers like GPNMB/LPL/FABP5 that primary MASH-vs-noMASH cannot).
    metmac_vs_other = metmac_vs_other_de(adata_mac, out_dir)

    # Validation reporting (MetMac MASH-up axis)
    report_metmac_validation(de_results, comp)

    # ── Write the macrophage-only h5ad (do NOT overwrite the master) ─────────
    # Convert tie-break dict to a JSON-able shape and stash on .uns so it
    # round-trips through h5ad.
    adata_mac.uns["mac_subcluster_args"] = {
        "leiden_resolution": float(args.leiden_resolution),
        "hvg_n": int(args.hvg_n),
        "n_pcs": int(args.n_pcs),
        "n_neighbors": int(args.n_neighbors),
        "gpnmb_tiebreak_quantile": float(GPNMB_TIEBREAK_QUANTILE),
        "min_cells_per_group": int(MIN_CELLS_PER_GROUP),
    }
    # h5ad can't store nested dicts of mixed dtypes inside .uns reliably; the
    # decision table is already a flat dict-of-lists (see assign_subtype_labels:
    # `dec_df.to_dict(orient="list")`). Re-flatten here was a bug — it produced
    # integer-keyed dicts that anndata cannot serialize. We leave the dict-of-
    # lists in place if already populated. If somehow legacy state contains
    # non-string keys, coerce them to strings.
    dec = adata_mac.uns.get("mac_subtype_decision", {})
    if isinstance(dec, dict) and dec and not all(isinstance(k, str) for k in dec.keys()):
        dec = {str(k): v for k, v in dec.items()}
        adata_mac.uns["mac_subtype_decision"] = dec
    # Cluster-mean tables → flattened lists for h5ad safety
    adata_mac.uns.pop("mac_subtype_cluster_scores", None)
    adata_mac.uns.pop("mac_subtype_cluster_markers", None)

    _log(f"Writing macrophage AnnData → {adata_out_path}")
    adata_mac.write_h5ad(adata_out_path, compression="gzip")
    out_size = adata_out_path.stat().st_size
    _log(f"  AnnData written ({out_size/1e6:.1f} MB)")

    # ── DE summary TSV ───────────────────────────────────────────────────────
    summary_rows = []
    for r in de_results:
        summary_rows.append({
            "subtype":             r["cell_type"],
            "contrast":            "MASH_vs_no_MASH",
            "n_cells_MASH":        r.get("n_cells_mash", 0),
            "n_cells_no_MASH":     r.get("n_cells_nomash", 0),
            "n_genes_tested":      r.get("n_genes_tested", 0),
            "n_sig_padj05":        r.get("n_sig_padj05", 0),
            "n_sig_padj05_up":     r.get("n_sig_padj05_up", 0),
            "n_sig_padj05_down":   r.get("n_sig_padj05_down", 0),
            "skipped":             r.get("skipped", False),
            "skip_reason":         r.get("reason", ""),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_dir / "cosmx_mac_subcluster_de_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    _log(f"Wrote DE summary → {summary_path.name}")

    _log("DONE")


if __name__ == "__main__":
    main()
