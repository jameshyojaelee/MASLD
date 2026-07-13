#!/usr/bin/env python
# 42f_govaere2026_cosmx_il32_coloc.py
#
# SINGLE-CELL SPATIAL NOVELTY (a finding Visium's 55-um spots cannot resolve):
# in the Govaere et al. CosMx data (522,145 segmented cells, 968-gene panel,
# 4 slides), do IL32-high hepatocytes sit CLOSER to macrophages than IL32-low
# hepatocytes? IL32 is a bulk-DEG-replicated, slide-concordant disease gene; the
# question is whether its hepatocyte expression is spatially organised toward the
# macrophage compartment (a hepatocyte->macrophage cross-talk geometry).
#
# Design (per MASH slide; Leuven_1/3/4 = MASH, Leuven_2 = mixed-Normal, reported
# separately):
#   1. Within-slide, partition hepatocytes IL32-high vs IL32-low (lognorm layer).
#   2. Build a KD-tree on macrophage (KC) spatial coords; for each hepatocyte,
#      distance to the nearest macrophage (global px -> um, 0.12028 um/px).
#   3. Compare nearest-macrophage distance, IL32-high vs IL32-low hepatocytes.
#   4. IL32(hep) <-> CD74 of the k=5 nearest macrophages (ligand-receptor geometry).
#   5. Repeat nearest-distance to the MetMac subset (joined from the
#      mac-subclustered object) as a secondary readout.
#
# SIGNIFICANCE = slide-level DIRECTION CONCORDANCE across the 3 MASH slides, NOT
# cell-level p-values (cells within a slide are pseudoreplicated).
#
# CAVEATS (for legend): the recovered "MetMac" is MHC-II-high, NOT the canonical
# GPNMB+ lipid-associated macrophage (968-gene panel limit) -> do NOT claim LAM
# identity. Leuven_2 is a mixed-Normal slide. Distances are within-slide.
#
# Input : Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad  (full, lognorm layer)
#         Analysis/Spatial/results/preprocessed/cosmx_govaere2026_macsubclustered.h5ad (mac_subtype)
# Output: Analysis/Spatial/results/govaere2026/il32_colocalization/
# Env   : spatial   (compute node; 522k cells, KD-tree per slide; ~8-16 GB)

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5_FULL = os.path.join(BASE, "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad")
H5_MAC  = os.path.join(BASE, "Analysis/Spatial/results/preprocessed/cosmx_govaere2026_macsubclustered.h5ad")
OUT = os.path.join(BASE, "Analysis/Spatial/results/govaere2026/il32_colocalization")
os.makedirs(OUT, exist_ok=True)

PX_UM = 0.12028                       # CosMx SMI pixel size (um/px)
MASH_SLIDES = ["Leuven_1", "Leuven_3", "Leuven_4"]   # Leuven_2 = mixed-Normal
KNN = 5

print("[1] Loading full CosMx object ...", flush=True)
a = sc.read_h5ad(H5_FULL)
if "lognorm" in a.layers:
    a.X = a.layers["lognorm"].copy()
    print("    using lognorm layer", flush=True)
else:
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    print("    normalized counts -> log1p-CPM", flush=True)

# join macrophage subtype (KC / TransMac / MetMac) from the macsubclustered file
print("[2] Joining mac_subtype ...", flush=True)
mac = sc.read_h5ad(H5_MAC, backed="r")
sub_map = mac.obs["mac_subtype"].astype(str)
a.obs["mac_subtype"] = a.obs_names.map(sub_map).fillna("none").astype(str)
print("    mac_subtype:", a.obs["mac_subtype"].value_counts().to_dict(), flush=True)

xy = np.asarray(a.obsm["spatial"], dtype=float)
ct = a.obs["cell_type"].astype(str).values
il32 = np.asarray(a[:, "IL32"].X.todense()).ravel() if hasattr(a[:, "IL32"].X, "todense") \
       else np.asarray(a[:, "IL32"].X).ravel()
cd74 = np.asarray(a[:, "CD74"].X.todense()).ravel() if hasattr(a[:, "CD74"].X, "todense") \
       else np.asarray(a[:, "CD74"].X).ravel()
slide = a.obs["sample_id"].astype(str).values
mac_sub = a.obs["mac_subtype"].values

per_slide = []
hep_rows = []   # per-hepatocyte distances (for the figure violin), MASH only
for sl in sorted(set(slide)):
    m = slide == sl
    hep = m & (ct == "Hepatocyte")
    macm = m & (ct == "KC")
    metm = m & (ct == "KC") & (mac_sub == "MetMac")
    n_hep, n_mac, n_met = int(hep.sum()), int(macm.sum()), int(metm.sum())
    if n_hep < 50 or n_mac < 20:
        print(f"    {sl}: too few cells (hep={n_hep}, mac={n_mac}) -> skip", flush=True)
        continue

    hep_xy = xy[hep]; hep_il32 = il32[hep]
    # within-slide IL32-high / -low partition (sparse-aware)
    frac_pos = float((hep_il32 > 0).mean())
    if frac_pos < 0.5:
        hi = hep_il32 > 0; lo = ~hi; mode = "detected_vs_not"
    else:
        med = np.median(hep_il32); hi = hep_il32 > med; lo = ~hi; mode = "median_split"

    # nearest macrophage (any KC) distance, in um
    tree = cKDTree(xy[macm])
    d_mac, idx_mac = tree.query(hep_xy, k=min(KNN, n_mac))
    d_mac = np.atleast_2d(d_mac)
    nn_mac_um = d_mac[:, 0] * PX_UM
    # CD74 of the k nearest macrophages
    mac_cd74 = cd74[macm]
    idx_mac = np.atleast_2d(idx_mac)
    knn_cd74 = mac_cd74[idx_mac].mean(axis=1)

    # nearest MetMac distance (secondary)
    if n_met >= 20:
        tmet = cKDTree(xy[metm])
        nn_met_um = tmet.query(hep_xy, k=1)[0] * PX_UM
    else:
        nn_met_um = np.full(n_hep, np.nan)

    d_hi = nn_mac_um[hi]; d_lo = nn_mac_um[lo]
    med_hi, med_lo = float(np.median(d_hi)), float(np.median(d_lo))
    try:
        p_cell = mannwhitneyu(d_hi, d_lo, alternative="less").pvalue   # hi CLOSER
    except Exception:
        p_cell = np.nan
    rho_il32_cd74 = float(spearmanr(hep_il32, knn_cd74).correlation)

    is_mash = sl in MASH_SLIDES
    per_slide.append(dict(
        slide=sl, is_mash=is_mash, sample_disease=str(a.obs["sample_disease"][m][0]),
        n_hep=n_hep, n_mac=n_mac, n_metmac=n_met, il32_mode=mode, frac_il32_pos=frac_pos,
        median_dist_il32high_um=med_hi, median_dist_il32low_um=med_lo,
        delta_high_minus_low_um=med_hi - med_lo,     # <0 => IL32-high CLOSER
        il32high_closer=bool(med_hi < med_lo),
        median_dist_metmac_high_um=float(np.nanmedian(nn_met_um[hi])),
        median_dist_metmac_low_um=float(np.nanmedian(nn_met_um[lo])),
        rho_il32_vs_knn_cd74=rho_il32_cd74, cell_level_p_mwu=p_cell))
    print(f"    {sl} (MASH={is_mash}): IL32-high nearest-mac {med_hi:.1f}um vs low {med_lo:.1f}um "
          f"(delta={med_hi-med_lo:+.1f}, closer={med_hi<med_lo}); rho(IL32,knnCD74)={rho_il32_cd74:+.3f}",
          flush=True)

    if is_mash:
        df = pd.DataFrame(dict(slide=sl, nearest_mac_um=nn_mac_um,
                               il32_group=np.where(hi, "IL32-high", "IL32-low")))
        hep_rows.append(df)

ps = pd.DataFrame(per_slide)
ps.to_csv(os.path.join(OUT, "il32_coloc_per_slide.csv"), index=False)
if hep_rows:
    pd.concat(hep_rows, ignore_index=True).to_csv(
        os.path.join(OUT, "il32_hep_nearest_mac_distances_mash.csv"), index=False)

# slide-level direction concordance over the 3 MASH slides
mash = ps[ps["is_mash"]]
n_concord = int(mash["il32high_closer"].sum())
n_mash = int(mash.shape[0])
with open(os.path.join(OUT, "il32_coloc_summary.txt"), "w") as fh:
    fh.write("CosMx IL32 hepatocyte->macrophage spatial co-localization (Govaere CosMx)\n")
    fh.write(f"MASH slides where IL32-high hepatocytes are CLOSER to macrophages than "
             f"IL32-low: {n_concord}/{n_mash}\n")
    fh.write(f"Mean delta (high-low nearest-mac distance, um) over MASH slides: "
             f"{mash['delta_high_minus_low_um'].mean():.2f} (negative = high closer)\n")
    fh.write(f"Mean rho(IL32, kNN-CD74) over MASH slides: "
             f"{mash['rho_il32_vs_knn_cd74'].mean():+.3f}\n\n")
    fh.write(ps.to_string(index=False) + "\n\n")
    fh.write("SIGNIFICANCE = slide-level direction concordance (n=3 MASH), NOT cell-level p "
             "(cells pseudoreplicated). MetMac is MHC-II-high, not canonical GPNMB+ LAM "
             "(968-gene panel) -> no LAM identity claim. Leuven_2 = mixed-Normal (reported, not in concordance).\n")
print(f"[done] MASH slide concordance: {n_concord}/{n_mash} IL32-high-closer. Outputs in {OUT}", flush=True)
