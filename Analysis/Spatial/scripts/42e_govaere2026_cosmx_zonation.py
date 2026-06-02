#!/usr/bin/env python3
"""
42e_govaere2026_cosmx_zonation.py — Score CosMx hepatocytes by liver zonation
(pericentral / midzonal / periportal) using canonical markers, then test for
MetMac proximity, IL32 expression, steatosis co-localisation, and disease-stage
composition.

Inputs:
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad

Disease-stage assignment (slide-level; documented in W2.B loader):
  Leuven_1 (GSM9351704)  "Explant_x2"          → 2 end-stage MASH explants  → MASH
  Leuven_2 (GSM9351705)  "F3_MASL_Normal_F3"   → 1 Normal + 1 MASL + 2 F3   → no_MASH (mixed-with-normal)
  Leuven_3 (GSM9351706)  "F1_F4"               → F1 + F4                    → MASH
  Leuven_4 (GSM9351707)  "F3_F3_MASL"          → 2 F3 + 1 MASL              → MASH

Zonation markers (canonical):
  Pericentral (PC): GLUL, CYP1A2, CYP2E1, CYP3A4, OAT, AXIN2
  Periportal  (PP): PCK1, ARG1, ALB, SDS, ASS1, HAL, HMGCS2
  Midzonal    (MZ): HAMP, CYP8B1  (transitional; less reliable)

Zone assignment (per hepatocyte):
  - PP   if score_periportal > score_pericentral + delta
  - PC   if score_pericentral > score_periportal + delta
  - MZ   otherwise (transitional)
  Default delta = 0.3; auto-tuned if it yields >80% in any single zone.

F150 CAVEAT: the CosMx 1000-plex panel lacks the canonical pericentral/periportal
zonation anchors (GLUL/CYP1A2/CYP2E1/CYP3A4/OAT/AXIN2 and ALB/PCK1/ARG1 are mostly
out of panel). The PC/PP axis is therefore driven by SURROGATE genes (LGR5/SERPINA1/
APOA1/HPD) and is a LOW-CONFIDENCE surrogate, NOT canonical zonation. The script
audits canonical-anchor survival, emits a zonation_confidence_flag.tsv, and writes a
caveat into zonation_summary.md; no canonical PC/PP concordance is claimed.

Outputs (Analysis/Spatial/results/govaere2026/zonation/):
  hep_zonation_scores.tsv        — per-cell zonation scores
  metmac_proximity_to_zone.tsv   — per-sample macrophage-to-PC vs PP distances
  il32_by_zone.tsv               — IL32 expression × zone × sample + KW test
  steatosis_by_zone.tsv          — lipid markers × zone × sample
  zone_disease_composition.tsv   — sample × zone × disease × cell counts/proportion
  zonation_confidence_flag.tsv   — F150: low-confidence flag + canonical-anchor count
  zonation_summary.md            — written interpretation (with F150 caveat header)

Environment:
  micromamba activate spatial

SLURM:
  --partition=io --qos=interactive --cpus-per-task=4 --mem=64G --time=4:00:00
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Dict, List, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.neighbors import KDTree

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ADATA_PATH   = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad"
OUT_DIR      = PROJECT_ROOT / "Analysis/Spatial/results/govaere2026/zonation"

# ── Disease stage assignment ─────────────────────────────────────────────────
SAMPLE_DISEASE = {
    "Leuven_1": {"sample_disease_str": "Explant_x2",        "stage": "MASH"},
    "Leuven_2": {"sample_disease_str": "F3_MASL_Normal_F3", "stage": "no_MASH"},
    "Leuven_3": {"sample_disease_str": "F1_F4",             "stage": "MASH"},
    "Leuven_4": {"sample_disease_str": "F3_F3_MASL",        "stage": "MASH"},
}

# ── Canonical zonation marker panels ────────────────────────────────────────
# Per-zone marker sets. Task spec markers come first; we extend with
# additional canonical Halpern 2017 / Aizarani 2019 / MacParland 2018 zonation
# markers because the CosMx 1000-plex panel has heavy CYP-family attrition
# (the panel ships only CYP1B1 + CYP2U1; ALB/CYP1A2/CYP2E1/CYP3A4/OAT/AXIN2 etc.
# are all out of panel). Extended markers below are documented zone-specific
# in primary single-cell atlases and survive the CosMx panel restriction.
ZONE_MARKERS: Dict[str, List[str]] = {
    # Task-spec PC: GLUL CYP1A2 CYP2E1 CYP3A4 OAT AXIN2  (5/6 absent from panel)
    # Extended PC: LGR5 (Wnt-target; PC stem; Halpern 2017)
    "pericentral": ["GLUL", "CYP1A2", "CYP2E1", "CYP3A4", "OAT", "AXIN2", "LGR5"],
    # Task-spec PP: PCK1 ARG1 ALB SDS ASS1 HAL HMGCS2  (6/7 absent from panel)
    # Extended PP: HPD (tyrosine catabolism; PP-enriched), SERPINA1 (alpha-1-antitrypsin; PP),
    #              APOA1 (apolipoprotein synthesis; PP-dominant)
    "periportal":  ["PCK1", "ARG1", "ALB", "SDS", "ASS1", "HAL", "HMGCS2",
                    "HPD", "SERPINA1", "APOA1"],
    # Task-spec MZ: HAMP CYP8B1  (both absent)
    # Extended MZ: IGFBP1 (zone 2 hepatocyte marker; Aizarani 2019)
    "midzonal":    ["HAMP", "CYP8B1", "IGFBP1"],
}

# F150: canonical zonation ANCHORS (the gold-standard task-spec markers, before
# the surrogate-gene extension). The PC/PP axis is only trustworthy if enough of
# THESE survive the CosMx 1000-plex panel. The extended markers (LGR5 / SERPINA1 /
# APOA1 / HPD / IGFBP1) are secretory/stem surrogates, not bulk PC/PP zonation
# genes, so an axis driven mostly by them does not reflect true zonation.
CANONICAL_ANCHORS: Dict[str, List[str]] = {
    "pericentral": ["GLUL", "CYP1A2", "CYP2E1", "CYP3A4", "OAT", "AXIN2"],
    "periportal":  ["PCK1", "ARG1", "ALB", "SDS", "ASS1", "HAL", "HMGCS2"],
    "midzonal":    ["HAMP", "CYP8B1"],
}
# Minimum canonical anchors per zone for the axis to be called confidence-bearing.
MIN_CANONICAL_ANCHORS = 2

# Lipid/steatosis markers (use those present in panel).
# Task spec: FABP1 PLIN2 PNPLA3. Extended: FASN CIDEA CIDEC + family members
# (FABP4/FABP5 — note FABP1 is not in panel, FABP4/5 are surrogate lipid-binding markers;
# FABP4 = adipocyte/macrophage so should be interpreted cautiously in hep DE);
# SREBF1 (master lipogenic TF) — included as expression-context marker.
LIPID_MARKERS = ["FABP1", "PLIN2", "PNPLA3", "FASN", "CIDEA", "CIDEC",
                 "FABP4", "FABP5", "SREBF1"]

# Cytokine of interest for spatial co-localisation analysis
IL32_GENE = "IL32"

# Macrophage cell-type label as set by 41_govaere2026_cosmx_load.py
MACROPHAGE_LABELS = ["KC", "MetMac", "TransMac"]


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def assign_disease_stage(adata: ad.AnnData) -> ad.AnnData:
    """Map .obs.sample_id → MASH / no_MASH at slide level."""
    stage_map = {sid: meta["stage"] for sid, meta in SAMPLE_DISEASE.items()}
    adata.obs["disease_stage"] = adata.obs["sample_id"].map(stage_map).astype("category")
    return adata


def verify_panel_markers(adata: ad.AnnData) -> Dict[str, List[str]]:
    """Verify which canonical zonation markers + lipid markers + IL32 are present."""
    available = set(adata.var_names)
    present_zone: Dict[str, List[str]] = {}
    _log("Panel marker availability:")
    for zone, markers in ZONE_MARKERS.items():
        present = [m for m in markers if m in available]
        missing = [m for m in markers if m not in available]
        present_zone[zone] = present
        _log(f"  {zone:>11}: present={present}  missing={missing}")
        if len(present) < 3:
            _log(f"    WARNING: only {len(present)} markers for {zone} (< 3) — proceeding with what's available")

    present_lipid = [m for m in LIPID_MARKERS if m in available]
    missing_lipid = [m for m in LIPID_MARKERS if m not in available]
    _log(f"  {'steatosis':>11}: present={present_lipid}  missing={missing_lipid}")
    present_zone["_lipid"] = present_lipid

    il32_present = IL32_GENE in available
    _log(f"  IL32 in panel: {il32_present}")
    present_zone["_il32"] = [IL32_GENE] if il32_present else []

    # ── F150: canonical-anchor audit + zonation-confidence flag ──────────────
    # Count how many GOLD-STANDARD task-spec anchors (not surrogates) survive the
    # panel per zone. If PC/PP fall below MIN_CANONICAL_ANCHORS, the PC/PP axis is
    # a low-confidence surrogate and must NOT be presented as canonical zonation.
    anchor_present: Dict[str, List[str]] = {}
    _log("  --- Canonical zonation anchor audit (F150) ---")
    for zone, anchors in CANONICAL_ANCHORS.items():
        present = [m for m in anchors if m in available]
        anchor_present[zone] = present
        _log(f"    {zone:>11}: canonical anchors present={present}  "
             f"({len(present)}/{len(anchors)})")
    present_zone["_canonical_anchors"] = anchor_present
    n_pc_anchor = len(anchor_present.get("pericentral", []))
    n_pp_anchor = len(anchor_present.get("periportal", []))
    low_conf = (n_pc_anchor < MIN_CANONICAL_ANCHORS
                or n_pp_anchor < MIN_CANONICAL_ANCHORS)
    present_zone["_zonation_low_confidence"] = [str(low_conf)]
    if low_conf:
        _log("    WARNING (F150): CosMx 1000-plex panel cannot resolve canonical "
             f"PC/PP zonation (PC canonical anchors={n_pc_anchor}, "
             f"PP canonical anchors={n_pp_anchor}; need >= {MIN_CANONICAL_ANCHORS} "
             "each). The PC/PP axis is driven by SURROGATE genes (LGR5/SERPINA1/"
             "APOA1/HPD), so it is a LOW-CONFIDENCE surrogate axis, NOT canonical "
             "zonation. Downstream PC/PP labels and zone-stratified results are "
             "flagged low-confidence; do not claim canonical pericentral/periportal "
             "concordance.")
    else:
        _log("    Canonical PC/PP anchors sufficient — zonation axis is "
             "confidence-bearing.")
    return present_zone


def score_zonation_per_cell(adata_hep: ad.AnnData,
                            present_zone: Dict[str, List[str]]) -> ad.AnnData:
    """Compute pericentral / periportal / midzonal scores via sc.tl.score_genes.

    Uses the log-normalised layer if available (matches loader settings).
    """
    if "lognorm" in adata_hep.layers:
        adata_hep.X = adata_hep.layers["lognorm"]
        _log("  Using .layers['lognorm'] for zone scoring")
    else:
        _log("  WARN: .layers['lognorm'] missing — re-normalising on hep subset")
        adata_hep.layers["counts"] = adata_hep.X.copy()
        sc.pp.normalize_total(adata_hep, target_sum=1e4)
        sc.pp.log1p(adata_hep)
        adata_hep.layers["lognorm"] = adata_hep.X.copy()

    for zone in ("pericentral", "periportal", "midzonal"):
        markers = present_zone[zone]
        if len(markers) == 0:
            _log(f"  {zone}: no markers available — setting score to 0")
            adata_hep.obs[f"score_{zone}"] = 0.0
            continue
        sc.tl.score_genes(
            adata_hep,
            gene_list=markers,
            score_name=f"score_{zone}",
            random_state=0,
            ctrl_size=min(50, adata_hep.n_vars - len(markers)),
            use_raw=False,
        )
        _log(f"  Scored {zone} (n={len(markers)} markers); "
             f"mean={float(adata_hep.obs[f'score_{zone}'].mean()):.4f}, "
             f"std={float(adata_hep.obs[f'score_{zone}'].std()):.4f}")

    return adata_hep


def assign_zone_labels(adata_hep: ad.AnnData,
                       delta: float = 0.3) -> Tuple[ad.AnnData, float, pd.Series]:
    """Categorical zone label per hepatocyte. Auto-tune delta if extreme bias."""

    def _classify(d: float) -> pd.Series:
        diff = adata_hep.obs["score_periportal"] - adata_hep.obs["score_pericentral"]
        zone = pd.Series("MZ", index=adata_hep.obs.index, dtype=object)
        zone.loc[diff > d] = "PP"
        zone.loc[diff < -d] = "PC"
        return zone

    zone = _classify(delta)
    frac = zone.value_counts(normalize=True)
    _log(f"  delta={delta}: PC={frac.get('PC', 0):.3f}  MZ={frac.get('MZ', 0):.3f}  "
         f"PP={frac.get('PP', 0):.3f}")

    # Auto-tune if any single zone > 80%
    if (frac.max() > 0.80) and (delta > 0.05):
        new_deltas = [0.2, 0.15, 0.1, 0.05]
        for d in new_deltas:
            zone_new = _classify(d)
            frac_new = zone_new.value_counts(normalize=True)
            _log(f"  retrying delta={d}: PC={frac_new.get('PC', 0):.3f}  "
                 f"MZ={frac_new.get('MZ', 0):.3f}  PP={frac_new.get('PP', 0):.3f}")
            if frac_new.max() <= 0.80:
                zone = zone_new
                delta = d
                frac = frac_new
                _log(f"  → using delta={delta}")
                break

    adata_hep.obs["hep_zone"] = pd.Categorical(zone, categories=["PC", "MZ", "PP"])
    return adata_hep, delta, frac


def write_per_cell_scores(adata_hep: ad.AnnData,
                          present_zone: Dict[str, List[str]],
                          out_dir: pathlib.Path) -> pathlib.Path:
    """Write per-hepatocyte zonation scores TSV."""
    cols = {
        "cell_id":   adata_hep.obs_names.values,
        "sample_id": adata_hep.obs["sample_id"].astype(str).values,
        "disease_stage": adata_hep.obs["disease_stage"].astype(str).values,
        "centerX_global_px": adata_hep.obs["centerX_global_px"].values,
        "centerY_global_px": adata_hep.obs["centerY_global_px"].values,
        "hep_zone":  adata_hep.obs["hep_zone"].astype(str).values,
        "score_PC":  adata_hep.obs["score_pericentral"].values,
        "score_PP":  adata_hep.obs["score_periportal"].values,
        "score_MZ":  adata_hep.obs["score_midzonal"].values,
    }

    # IL32 expression per cell (lognorm)
    if IL32_GENE in adata_hep.var_names:
        il32_idx = adata_hep.var_names.get_loc(IL32_GENE)
        X = adata_hep.layers["lognorm"]
        col = X[:, il32_idx]
        if sparse.issparse(col):
            col = np.asarray(col.todense()).ravel()
        else:
            col = np.asarray(col).ravel()
        cols["IL32_expr"] = col

    # Lipid marker expression
    for gene in present_zone["_lipid"]:
        gidx = adata_hep.var_names.get_loc(gene)
        col = adata_hep.layers["lognorm"][:, gidx]
        if sparse.issparse(col):
            col = np.asarray(col.todense()).ravel()
        else:
            col = np.asarray(col).ravel()
        cols[f"{gene}_expr"] = col

    df = pd.DataFrame(cols)
    out = out_dir / "hep_zonation_scores.tsv"
    df.to_csv(out, sep="\t", index=False)
    _log(f"  Wrote per-cell scores → {out.name} ({df.shape[0]:,} rows, "
         f"{out.stat().st_size/1e6:.1f} MB)")
    return out


def compute_metmac_proximity(adata: ad.AnnData,
                             hep_zone_idx: pd.Series,
                             out_dir: pathlib.Path) -> pathlib.Path:
    """For each macrophage cell, compute Euclidean distance to nearest PC vs PP hepatocyte.

    KDTree per (sample_id, zone). Returns a per-sample summary TSV.
    """
    summary_rows: List[Dict] = []
    cell_rows: List[Dict] = []

    mac_mask = adata.obs["cell_type"].isin(MACROPHAGE_LABELS)
    mac_cells = adata.obs.loc[mac_mask, ["sample_id", "cell_type",
                                          "centerX_global_px", "centerY_global_px"]].copy()
    # Filter out cells with NaN spatial coords (CosMx occasionally drops these)
    n_mac_before = mac_cells.shape[0]
    mac_cells = mac_cells.dropna(subset=["centerX_global_px", "centerY_global_px"])
    n_mac_after = mac_cells.shape[0]
    if n_mac_before != n_mac_after:
        _log(f"  Dropped {n_mac_before - n_mac_after} macrophage cells with NaN coords")
    _log(f"  Total macrophage cells (KC/MetMac/TransMac, valid coords): {n_mac_after:,}")

    # Build a hep-only DataFrame with zone labels; drop NaN coords too
    hep_df = pd.DataFrame({
        "sample_id":       adata.obs.loc[hep_zone_idx.index, "sample_id"].astype(str).values,
        "x":               adata.obs.loc[hep_zone_idx.index, "centerX_global_px"].values,
        "y":               adata.obs.loc[hep_zone_idx.index, "centerY_global_px"].values,
        "hep_zone":        hep_zone_idx.values,
    }, index=hep_zone_idx.index)
    n_hep_before = hep_df.shape[0]
    hep_df = hep_df.dropna(subset=["x", "y"])
    if hep_df.shape[0] != n_hep_before:
        _log(f"  Dropped {n_hep_before - hep_df.shape[0]} hepatocytes with NaN coords")

    for sid in sorted(adata.obs["sample_id"].unique()):
        mac_sub = mac_cells[mac_cells["sample_id"] == sid]
        hep_sub = hep_df[hep_df["sample_id"] == sid]
        if mac_sub.shape[0] < 10 or hep_sub.shape[0] < 100:
            _log(f"  [{sid}] SKIP — mac_n={mac_sub.shape[0]}, hep_n={hep_sub.shape[0]}")
            continue

        pc_pts = hep_sub.loc[hep_sub["hep_zone"] == "PC", ["x", "y"]].values
        pp_pts = hep_sub.loc[hep_sub["hep_zone"] == "PP", ["x", "y"]].values
        mac_pts = mac_sub[["centerX_global_px", "centerY_global_px"]].values

        # Final NaN safety net
        pc_pts = pc_pts[np.isfinite(pc_pts).all(axis=1)]
        pp_pts = pp_pts[np.isfinite(pp_pts).all(axis=1)]
        mac_finite_mask = np.isfinite(mac_pts).all(axis=1)
        mac_pts = mac_pts[mac_finite_mask]
        mac_sub = mac_sub.iloc[np.where(mac_finite_mask)[0]]

        if pc_pts.shape[0] < 10 or pp_pts.shape[0] < 10 or mac_pts.shape[0] < 10:
            _log(f"  [{sid}] SKIP — PC_n={pc_pts.shape[0]}, PP_n={pp_pts.shape[0]}, "
                 f"mac_n={mac_pts.shape[0]}")
            continue

        tree_pc = KDTree(pc_pts)
        tree_pp = KDTree(pp_pts)
        d_pc, _ = tree_pc.query(mac_pts, k=1)
        d_pp, _ = tree_pp.query(mac_pts, k=1)
        d_pc = d_pc.ravel()
        d_pp = d_pp.ravel()

        # Per-cell rows
        for i in range(mac_sub.shape[0]):
            cell_rows.append({
                "cell_id":        mac_sub.index[i],
                "sample_id":      sid,
                "cell_type":      mac_sub["cell_type"].iloc[i],
                "disease_stage":  SAMPLE_DISEASE[sid]["stage"],
                "dist_to_PC_px":  float(d_pc[i]),
                "dist_to_PP_px":  float(d_pp[i]),
                "closer_to":      "PC" if d_pc[i] < d_pp[i] else "PP",
            })

        # Wilcoxon paired test: distance to PC vs PP (per macrophage)
        try:
            w_stat, w_p = stats.wilcoxon(d_pc, d_pp, alternative="two-sided")
        except ValueError:
            w_stat, w_p = float("nan"), float("nan")

        median_pc = float(np.median(d_pc))
        median_pp = float(np.median(d_pp))
        mean_pc = float(np.mean(d_pc))
        mean_pp = float(np.mean(d_pp))

        summary_rows.append({
            "sample_id":        sid,
            "disease_stage":    SAMPLE_DISEASE[sid]["stage"],
            "n_macrophages":    int(mac_sub.shape[0]),
            "n_hep_PC":         int(pc_pts.shape[0]),
            "n_hep_PP":         int(pp_pts.shape[0]),
            "median_dist_PC_px": median_pc,
            "median_dist_PP_px": median_pp,
            "mean_dist_PC_px":   mean_pc,
            "mean_dist_PP_px":   mean_pp,
            "ratio_medPC_over_medPP": median_pc / median_pp if median_pp > 0 else float("nan"),
            "wilcoxon_stat":    float(w_stat),
            "wilcoxon_p":       float(w_p),
            "n_closer_to_PC":   int((d_pc < d_pp).sum()),
            "n_closer_to_PP":   int((d_pp < d_pc).sum()),
            "frac_closer_to_PC": float((d_pc < d_pp).mean()),
        })
        _log(f"  [{sid}] mac={mac_sub.shape[0]:,}  PC_hep={pc_pts.shape[0]:,}  "
             f"PP_hep={pp_pts.shape[0]:,}  median(d→PC)={median_pc:.0f}  "
             f"median(d→PP)={median_pp:.0f}  ratio={median_pc/median_pp if median_pp > 0 else float('nan'):.3f}  "
             f"frac_closer_PC={(d_pc < d_pp).mean():.3f}  Wp={w_p:.2e}")

    summary_df = pd.DataFrame(summary_rows)
    out = out_dir / "metmac_proximity_to_zone.tsv"
    summary_df.to_csv(out, sep="\t", index=False)
    _log(f"  Wrote macrophage proximity summary → {out.name} ({summary_df.shape[0]} rows)")

    # Write per-cell distances as well (useful for downstream plots)
    if cell_rows:
        cell_df = pd.DataFrame(cell_rows)
        cell_out = out_dir / "metmac_proximity_per_cell.tsv"
        cell_df.to_csv(cell_out, sep="\t", index=False)
        _log(f"  Wrote per-cell distances → {cell_out.name} ({cell_df.shape[0]:,} rows)")

    return out


def il32_by_zone(adata_hep: ad.AnnData,
                 out_dir: pathlib.Path) -> pathlib.Path:
    """Test whether IL32 expression differs across zones (KW + per-sample mean)."""
    if IL32_GENE not in adata_hep.var_names:
        _log(f"  IL32 absent from panel — skipping test 2")
        return None

    il32_idx = adata_hep.var_names.get_loc(IL32_GENE)
    X = adata_hep.layers["lognorm"]
    col = X[:, il32_idx]
    if sparse.issparse(col):
        il32 = np.asarray(col.todense()).ravel()
    else:
        il32 = np.asarray(col).ravel()

    obs = adata_hep.obs.copy()
    obs["IL32_expr"] = il32

    rows: List[Dict] = []
    for sid in sorted(obs["sample_id"].unique()):
        sub = obs[obs["sample_id"] == sid]
        pc_il32 = sub.loc[sub["hep_zone"] == "PC", "IL32_expr"].values
        mz_il32 = sub.loc[sub["hep_zone"] == "MZ", "IL32_expr"].values
        pp_il32 = sub.loc[sub["hep_zone"] == "PP", "IL32_expr"].values

        groups = [g for g in (pc_il32, mz_il32, pp_il32) if g.shape[0] >= 10]
        if len(groups) < 2:
            kw_p = float("nan")
            kw_stat = float("nan")
        else:
            try:
                kw_stat, kw_p = stats.kruskal(*groups)
            except ValueError:
                kw_stat, kw_p = float("nan"), float("nan")

        rows.append({
            "sample_id":     sid,
            "disease_stage": SAMPLE_DISEASE[sid]["stage"],
            "n_PC":          int(pc_il32.shape[0]),
            "n_MZ":          int(mz_il32.shape[0]),
            "n_PP":          int(pp_il32.shape[0]),
            "mean_IL32_PC":  float(np.mean(pc_il32)) if pc_il32.shape[0] else float("nan"),
            "mean_IL32_MZ":  float(np.mean(mz_il32)) if mz_il32.shape[0] else float("nan"),
            "mean_IL32_PP":  float(np.mean(pp_il32)) if pp_il32.shape[0] else float("nan"),
            "median_IL32_PC": float(np.median(pc_il32)) if pc_il32.shape[0] else float("nan"),
            "median_IL32_MZ": float(np.median(mz_il32)) if mz_il32.shape[0] else float("nan"),
            "median_IL32_PP": float(np.median(pp_il32)) if pp_il32.shape[0] else float("nan"),
            "pct_pos_PC":    float((pc_il32 > 0).mean()) if pc_il32.shape[0] else float("nan"),
            "pct_pos_MZ":    float((mz_il32 > 0).mean()) if mz_il32.shape[0] else float("nan"),
            "pct_pos_PP":    float((pp_il32 > 0).mean()) if pp_il32.shape[0] else float("nan"),
            "kruskal_stat":  float(kw_stat),
            "kruskal_p":     float(kw_p),
            "highest_zone":  ["PC", "MZ", "PP"][int(np.argmax([
                np.mean(pc_il32) if pc_il32.shape[0] else -np.inf,
                np.mean(mz_il32) if mz_il32.shape[0] else -np.inf,
                np.mean(pp_il32) if pp_il32.shape[0] else -np.inf,
            ]))],
        })
        _log(f"  [{sid}] IL32 mean: PC={rows[-1]['mean_IL32_PC']:.4f}  "
             f"MZ={rows[-1]['mean_IL32_MZ']:.4f}  PP={rows[-1]['mean_IL32_PP']:.4f}  "
             f"KW p={kw_p:.2e}  highest={rows[-1]['highest_zone']}")

    df = pd.DataFrame(rows)
    out = out_dir / "il32_by_zone.tsv"
    df.to_csv(out, sep="\t", index=False)
    _log(f"  Wrote IL32-by-zone → {out.name} ({df.shape[0]} rows)")
    return out


def steatosis_by_zone(adata_hep: ad.AnnData,
                      present_zone: Dict[str, List[str]],
                      out_dir: pathlib.Path) -> pathlib.Path:
    """Compute lipid-marker expression by zone (per sample). Optional steatosis score."""
    lipids = present_zone["_lipid"]
    if not lipids:
        _log("  No lipid markers in panel — skipping test 3")
        return None

    # Build a per-cell lipid_score = mean log-norm of available lipid markers
    lipid_idx = [adata_hep.var_names.get_loc(g) for g in lipids]
    X = adata_hep.layers["lognorm"]
    lip = X[:, lipid_idx]
    if sparse.issparse(lip):
        lip = np.asarray(lip.todense())
    else:
        lip = np.asarray(lip)
    lipid_score = lip.mean(axis=1)
    adata_hep.obs["steatosis_score"] = lipid_score

    rows: List[Dict] = []
    for sid in sorted(adata_hep.obs["sample_id"].unique()):
        sub_obs = adata_hep.obs[adata_hep.obs["sample_id"] == sid]
        sub_idx = sub_obs.index
        if sub_idx.shape[0] < 30:
            continue

        # Build per-gene + composite score arrays
        zone_vals = sub_obs["hep_zone"].astype(str).values
        sample_lipid = adata_hep[sub_idx].obs["steatosis_score"].values

        # Per-gene means by zone
        row: Dict[str, object] = {
            "sample_id":     sid,
            "disease_stage": SAMPLE_DISEASE[sid]["stage"],
            "n_PC":          int((zone_vals == "PC").sum()),
            "n_MZ":          int((zone_vals == "MZ").sum()),
            "n_PP":          int((zone_vals == "PP").sum()),
        }
        for z in ("PC", "MZ", "PP"):
            mask = zone_vals == z
            row[f"mean_steatosis_{z}"] = float(sample_lipid[mask].mean()) if mask.sum() else float("nan")
            row[f"median_steatosis_{z}"] = float(np.median(sample_lipid[mask])) if mask.sum() else float("nan")
        # Per-marker mean by zone
        for j, g in enumerate(lipids):
            gene_vals = lip[adata_hep.obs.index.get_indexer(sub_idx), j]
            for z in ("PC", "MZ", "PP"):
                mask = zone_vals == z
                row[f"mean_{g}_{z}"] = float(gene_vals[mask].mean()) if mask.sum() else float("nan")

        # Kruskal-Wallis on composite score across zones
        groups = []
        for z in ("PC", "MZ", "PP"):
            v = sample_lipid[zone_vals == z]
            if v.shape[0] >= 10:
                groups.append(v)
        if len(groups) >= 2:
            try:
                kw_stat, kw_p = stats.kruskal(*groups)
            except ValueError:
                kw_stat, kw_p = float("nan"), float("nan")
        else:
            kw_stat, kw_p = float("nan"), float("nan")
        row["kruskal_stat"] = float(kw_stat)
        row["kruskal_p"]    = float(kw_p)

        means = [row[f"mean_steatosis_{z}"] for z in ("PC", "MZ", "PP")]
        row["highest_zone"] = ["PC", "MZ", "PP"][int(np.nanargmax(means))]
        rows.append(row)
        _log(f"  [{sid}] steatosis mean: PC={row['mean_steatosis_PC']:.4f}  "
             f"MZ={row['mean_steatosis_MZ']:.4f}  PP={row['mean_steatosis_PP']:.4f}  "
             f"KW p={kw_p:.2e}  highest={row['highest_zone']}")

    df = pd.DataFrame(rows)
    out = out_dir / "steatosis_by_zone.tsv"
    df.to_csv(out, sep="\t", index=False)
    _log(f"  Wrote steatosis-by-zone → {out.name} ({df.shape[0]} rows; lipid markers used: {lipids})")
    return out


def zone_disease_composition(adata_hep: ad.AnnData,
                              out_dir: pathlib.Path) -> pathlib.Path:
    """Per-sample (× zone × disease) cell counts + proportions."""
    grp = (adata_hep.obs
           .groupby(["sample_id", "disease_stage", "hep_zone"], observed=True)
           .size()
           .rename("n_cells")
           .reset_index())
    totals = grp.groupby("sample_id", observed=True)["n_cells"].transform("sum")
    grp["proportion"] = grp["n_cells"] / totals

    # Also compute disease-stage aggregate
    disease_grp = (adata_hep.obs
                   .groupby(["disease_stage", "hep_zone"], observed=True)
                   .size()
                   .rename("n_cells")
                   .reset_index())
    d_totals = disease_grp.groupby("disease_stage", observed=True)["n_cells"].transform("sum")
    disease_grp["proportion"] = disease_grp["n_cells"] / d_totals
    disease_grp.insert(0, "sample_id", "ALL")  # tag aggregate rows

    out_df = pd.concat([grp, disease_grp], axis=0, ignore_index=True)
    out = out_dir / "zone_disease_composition.tsv"
    out_df.to_csv(out, sep="\t", index=False)
    _log(f"  Wrote zone-disease composition → {out.name} ({out_df.shape[0]} rows)")

    # Quick log of disease-aggregate proportions
    pivot = (disease_grp
             .pivot_table(index="disease_stage", columns="hep_zone",
                          values="proportion", fill_value=0.0))
    _log("Zone proportions by disease (aggregate):")
    for stage in pivot.index:
        row = " ".join(f"{z}={pivot.loc[stage, z]:.4f}" for z in pivot.columns)
        _log(f"    {stage:>8}: {row}")

    # Fisher's / chi-square test across disease groups for zone composition
    contingency = (grp.pivot_table(index="hep_zone", columns="disease_stage",
                                    values="n_cells", aggfunc="sum", fill_value=0))
    if contingency.shape[1] == 2 and contingency.shape[0] >= 2:
        chi2, p_chi, dof, _ = stats.chi2_contingency(contingency.values)
        _log(f"  Chi-square (zone × disease): chi2={chi2:.3f}, dof={dof}, p={p_chi:.2e}")
    return out


def write_summary(out_dir: pathlib.Path,
                  delta: float,
                  zone_frac: pd.Series,
                  present_zone: Dict[str, List[str]],
                  per_sample_zone: pd.DataFrame,
                  metmac_summary_path: pathlib.Path,
                  il32_path: pathlib.Path,
                  steatosis_path: pathlib.Path,
                  composition_path: pathlib.Path,
                  n_hep: int) -> None:
    """Write a markdown interpretation."""
    md_lines: List[str] = []
    md_lines.append("# CosMx Hepatocyte Zonation — Govaere 2026 (GSE312698)\n")
    md_lines.append(f"Script: `Analysis/Spatial/scripts/42e_govaere2026_cosmx_zonation.py`")
    md_lines.append(f"Input:  `Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad`")
    md_lines.append(f"Output: `Analysis/Spatial/results/govaere2026/zonation/`\n")

    # F150: zonation-confidence caveat at the very top of the report.
    anchor_present = present_zone.get("_canonical_anchors", {})
    low_conf = present_zone.get("_zonation_low_confidence", ["False"])[0] == "True"
    if low_conf:
        md_lines.append("> **CAVEAT (F150 — low-confidence zonation axis).** The CosMx "
                        "1000-plex panel has heavy CYP-family / ALB attrition: of the "
                        "canonical task-spec anchors, "
                        f"PC retains {anchor_present.get('pericentral', [])} and "
                        f"PP retains {anchor_present.get('periportal', [])}. The PC/PP "
                        "axis is therefore driven by SURROGATE genes (LGR5/SERPINA1/"
                        "APOA1/HPD), not bulk zonation markers, and does NOT reflect "
                        "true pericentral/periportal zonation. All PC/MZ/PP labels and "
                        "zone-stratified results below are a low-confidence surrogate; "
                        "no canonical zonation concordance is claimed.\n")

    md_lines.append("## Method\n")
    md_lines.append("Per-hepatocyte zone scores via `sc.tl.score_genes` on the log-normalised")
    md_lines.append("CosMx layer (target_sum=1e4 + log1p; loader defaults). Score panels:")
    md_lines.append("")
    for zone, markers in ZONE_MARKERS.items():
        present = present_zone[zone]
        canon = anchor_present.get(zone, [])
        md_lines.append(f"- {zone:>11}: panel={markers}  present={present}  "
                        f"canonical_anchors_present={canon}")
    md_lines.append(f"- lipid markers (steatosis): present={present_zone['_lipid']}")
    md_lines.append(f"- IL32 in panel: {IL32_GENE in present_zone['_il32']}")
    md_lines.append("")
    md_lines.append(f"Zone assignment (per hepatocyte):")
    md_lines.append(f"- PP if `score_PP - score_PC > {delta}`")
    md_lines.append(f"- PC if `score_PC - score_PP > {delta}`")
    md_lines.append(f"- MZ otherwise (transitional)\n")

    md_lines.append("## Results\n")
    md_lines.append(f"Total hepatocytes scored: {n_hep:,}.")
    md_lines.append(f"Aggregate zone proportions: "
                    f"PC={zone_frac.get('PC', 0):.3f}, "
                    f"MZ={zone_frac.get('MZ', 0):.3f}, "
                    f"PP={zone_frac.get('PP', 0):.3f}.\n")
    md_lines.append("Per-sample zone proportions:\n")
    md_lines.append("```")
    md_lines.append(per_sample_zone.to_string(index=False))
    md_lines.append("```\n")

    md_lines.append("## Outputs\n")
    md_lines.append(f"- `hep_zonation_scores.tsv` — per-cell scores + zone + IL32/lipid expression")
    md_lines.append(f"- `metmac_proximity_to_zone.tsv` — per-sample macrophage-to-PC vs PP distance")
    md_lines.append(f"- `metmac_proximity_per_cell.tsv` — per-macrophage distances")
    md_lines.append(f"- `il32_by_zone.tsv` — IL32 mean expression × zone × sample + KW")
    md_lines.append(f"- `steatosis_by_zone.tsv` — composite lipid score × zone × sample + KW")
    md_lines.append(f"- `zone_disease_composition.tsv` — zone × disease cell counts + proportions\n")

    md_lines.append("## Interpretation\n")
    md_lines.append("Across MASH-dominant slides we expect: (1) MASH-associated macrophages to ")
    md_lines.append("co-localise preferentially near pericentral hepatocytes early in disease, ")
    md_lines.append("then spread panlobularly with fibrosis progression; (2) IL32 — a MASH ")
    md_lines.append("hepatocyte-secreted cytokine — to track the pericentral metabolic zone; ")
    md_lines.append("(3) lipid-loading markers (FABP1/PLIN2/PNPLA3 when present) to be highest ")
    md_lines.append("in the same pericentral / mid-lobular zones in early disease. The ")
    md_lines.append("`metmac_proximity_to_zone.tsv` ratio (`median_dist_PC / median_dist_PP`) ")
    md_lines.append("and `frac_closer_to_PC` columns quantify this preference; per-sample KW ")
    md_lines.append("on IL32 and steatosis scores test the zone-specificity hypothesis ")
    md_lines.append("independently within each MASH/no-MASH slide. The Leuven_2 slide is the ")
    md_lines.append("noisy no-MASH comparator (mixed-with-normal; serves as a contrast).\n")

    out = out_dir / "zonation_summary.md"
    out.write_text("\n".join(md_lines))
    _log(f"  Wrote markdown summary → {out.name}")


def main():
    parser = argparse.ArgumentParser(description="Hepatocyte zonation scoring on Govaere 2026 CosMx")
    parser.add_argument("--adata", type=str, default=str(ADATA_PATH))
    parser.add_argument("--out-dir", type=str, default=str(OUT_DIR))
    parser.add_argument("--delta", type=float, default=0.3,
                        help="Score difference for PC/PP classification (auto-tunes if extreme bias).")
    args = parser.parse_args()

    adata_path = pathlib.Path(args.adata)
    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _log("=" * 70)
    _log("42e_govaere2026_cosmx_zonation.py")
    _log(f"  AnnData : {adata_path}")
    _log(f"  Out dir : {out_dir}")
    _log("=" * 70)

    _log(f"Reading AnnData ({adata_path.stat().st_size/1e6:.1f} MB)")
    adata = sc.read_h5ad(adata_path)
    _log(f"  shape: {adata.shape}")
    _log(f"  .obs columns: {list(adata.obs.columns)}")
    _log(f"  .layers: {list(adata.layers.keys())}")
    _log(f"  cell_type values: {sorted(adata.obs['cell_type'].unique().tolist())}")

    # Disease assignment (slide-level)
    adata = assign_disease_stage(adata)
    _log(f"  disease_stage counts: {adata.obs['disease_stage'].value_counts().to_dict()}")

    # Verify markers
    present_zone = verify_panel_markers(adata)

    # F150: persist a machine-readable zonation-confidence flag so downstream
    # consumers (atlas / figures) can gate on it instead of trusting PC/PP labels.
    anchor_present = present_zone.get("_canonical_anchors", {})
    low_conf = present_zone.get("_zonation_low_confidence", ["False"])[0] == "True"
    flag_df = pd.DataFrame([{
        "zonation_low_confidence": low_conf,
        "n_canonical_PC_anchors": len(anchor_present.get("pericentral", [])),
        "n_canonical_PP_anchors": len(anchor_present.get("periportal", [])),
        "canonical_PC_anchors_present": ";".join(anchor_present.get("pericentral", [])),
        "canonical_PP_anchors_present": ";".join(anchor_present.get("periportal", [])),
        "min_canonical_anchors_required": MIN_CANONICAL_ANCHORS,
        "note": ("PC/PP axis is a low-confidence surrogate (panel lacks canonical "
                 "zonation anchors); not canonical zonation" if low_conf
                 else "canonical anchors sufficient"),
    }])
    flag_path = out_dir / "zonation_confidence_flag.tsv"
    flag_df.to_csv(flag_path, sep="\t", index=False)
    _log(f"  Wrote zonation-confidence flag → {flag_path.name} "
         f"(low_confidence={low_conf})")

    # Filter to hepatocytes
    n_total = adata.n_obs
    hep_mask = adata.obs["cell_type"] == "Hepatocyte"
    n_hep = int(hep_mask.sum())
    _log(f"Hepatocytes: {n_hep:,} / {n_total:,} ({100*n_hep/n_total:.1f}%)")
    if n_hep < 1000:
        sys.exit(f"ERROR: only {n_hep} hepatocytes — aborting")

    adata_hep = adata[hep_mask].copy()
    _log(f"  hep AnnData shape: {adata_hep.shape}")

    # Score zonation per cell
    _log("Scoring zonation markers per hepatocyte")
    adata_hep = score_zonation_per_cell(adata_hep, present_zone)

    # Assign zone labels
    _log("Assigning zone labels (PC / MZ / PP)")
    adata_hep, delta, frac = assign_zone_labels(adata_hep, delta=args.delta)
    _log(f"  Final delta = {delta}")
    _log(f"  Final zone proportions: PC={frac.get('PC', 0):.3f}  "
         f"MZ={frac.get('MZ', 0):.3f}  PP={frac.get('PP', 0):.3f}")

    # Per-sample zone proportions for the summary table
    per_sample_zone = (adata_hep.obs
                       .groupby(["sample_id", "disease_stage", "hep_zone"], observed=True)
                       .size()
                       .rename("n_cells")
                       .reset_index())
    _log("Per-sample zone counts:")
    for sid in sorted(per_sample_zone["sample_id"].unique()):
        sub = per_sample_zone[per_sample_zone["sample_id"] == sid]
        tot = int(sub["n_cells"].sum())
        bits = []
        for _, row in sub.iterrows():
            bits.append(f"{row['hep_zone']}={int(row['n_cells']):>7,}({row['n_cells']/tot:.3f})")
        _log(f"  {sid:>9} ({SAMPLE_DISEASE[sid]['stage']:>8}): n={tot:>7,}  {'  '.join(bits)}")

    # Step 4 (write): per-cell scores TSV
    _log("Writing per-cell zonation scores TSV")
    write_per_cell_scores(adata_hep, present_zone, out_dir)

    # Step 6: MetMac/KC proximity to zonation
    _log("Test 1 — macrophage proximity to PC vs PP hepatocytes")
    hep_zone_series = adata_hep.obs["hep_zone"].astype(str)
    metmac_summary_path = compute_metmac_proximity(adata, hep_zone_series, out_dir)

    # Step 7: IL32 by zone
    _log("Test 2 — IL32 expression by zone")
    il32_path = il32_by_zone(adata_hep, out_dir)

    # Step 8: Steatosis by zone
    _log("Test 3 — Steatosis (lipid markers) by zone")
    steatosis_path = steatosis_by_zone(adata_hep, present_zone, out_dir)

    # Step 9: Zone × MASH disease composition
    _log("Test 4 — Zone × disease composition")
    composition_path = zone_disease_composition(adata_hep, out_dir)

    # Write markdown summary
    write_summary(out_dir, delta, frac, present_zone,
                  per_sample_zone, metmac_summary_path,
                  il32_path, steatosis_path, composition_path, n_hep)

    _log("DONE")


if __name__ == "__main__":
    main()
