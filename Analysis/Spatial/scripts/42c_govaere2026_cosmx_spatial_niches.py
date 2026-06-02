#!/usr/bin/env python3
"""
42c_govaere2026_cosmx_spatial_niches.py — Spatial neighborhood enrichment on
Govaere 2026 CosMx SMI 1000-plex deposit (GSE312698, Leuven_1..4 = 522,145
cells x 968 target genes). Tests paper Fig 5 spatial-paracrine claims:

  Test 1: GPNMB+ macrophages co-localize with IL32+ hepatocytes (kNN + perm)
  Test 2: Cell-type pair neighborhood enrichment (squidpy nhood_enrichment)
  Test 3: Portal-tract vs parenchymal segregation (cholangiocyte distance)
  Test 4: Lipogranuloma cluster detection (DBSCAN on KC cells)

Inputs:
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad
  data/external/govaere2026_natgenetics/cosmx/GSM935170{4,5,6,7}_Leuven_*_metadata_file.csv.gz

NOTE on coordinates: The h5ad .obs["centerX_global_px"], .obs["centerY_global_px"],
and .obsm["spatial"] are ALL NaN in the current deposit (loader 41_* aligned the
column but populated NaN — likely from a missing CenterX upstream join). The
RAW per-sample metadata CSVs (GSM*_metadata_file.csv.gz, columns 26 / 27
CenterX_global_px / CenterY_global_px) ARE populated. We rejoin per-cell
coordinates by (fov, cell_ID) from the raw CSVs and overwrite obsm["spatial"]
in-memory (the on-disk h5ad is not modified).

CosMx SMI pixel scale: 1 px = 0.18 microns (NanoString CosMx SMI specs).
  30 microns ~= 167 px  -> Test 2 inter-cell threshold
  200 microns ~= 1111 px -> Test 3 portal-tract radius
  500 microns ~= 2778 px -> Test 3 parenchymal threshold
  50 px ~= 9 microns      -> Test 4 DBSCAN eps (per task spec)
  100 px ~= 18 microns    -> Test 4 cluster composition radius

Outputs (Analysis/Spatial/results/govaere2026/niche/):
  gpnmb_il32_colocalization.tsv  — per sample: GPNMB+ count, IL32+ Hep count,
                                    mean_neighbor_il32_frac (F147 fix: now the
                                    fraction of the HEPATOCYTE neighbors of GPNMB+
                                    KC that are IL32+, under a composition-
                                    conditioned null that permutes IL32+ status
                                    within hepatocytes only), null mean, Z, p_perm.
                                    global_il32_hep_frac is now IL32+ fraction
                                    AMONG hepatocytes (denominator = Hep, not all
                                    cells); column names unchanged.
  nhood_enrichment_zscores.tsv   — per sample: cell_type1 x cell_type2 x Z.
  nhood_enrichment_summary.txt   — top 5 enriched / depleted pairs aggregated.
  portal_vs_parenchymal.tsv      — KC enrichment in portal regions per sample.
  lipogranuloma_clusters.tsv     — sample x cluster_id x n_cells x composition.
  niche_summary.md               — interpretation of all 4 tests.

Constraints:
  - DO NOT modify the master h5ad
  - DO NOT modify the multi-evidence atlas
  - Output dir is Analysis/Spatial/results/govaere2026/niche/

Environment: micromamba activate spatial
SLURM: --partition=io --qos=interactive --cpus-per-task=4 --mem=64G --time=4:00:00
"""
from __future__ import annotations

import argparse
import gc
import os
import pathlib
import sys
import time
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import fisher_exact
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors, BallTree

# ─── Paths ───────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ADATA_PATH   = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad"
COSMX_DIR    = PROJECT_ROOT / "data/external/govaere2026_natgenetics/cosmx"
OUT_DIR      = PROJECT_ROOT / "Analysis/Spatial/results/govaere2026/niche"

# Per-sample raw metadata files (for CenterX/Y_global_px which are NaN in h5ad)
RAW_META = {
    "Leuven_1": COSMX_DIR / "GSM9351704_Leuven_1_metadata_file.csv.gz",
    "Leuven_2": COSMX_DIR / "GSM9351705_Leuven_2_metadata_file.csv.gz",
    "Leuven_3": COSMX_DIR / "GSM9351706_Leuven_3_metadata_file.csv.gz",
    "Leuven_4": COSMX_DIR / "GSM9351707_Leuven_4_metadata_file.csv.gz",
}

# CosMx SMI px scale (NanoString): 1 px = 0.18 microns
PX_PER_UM = 1.0 / 0.18  # ≈ 5.556 px/μm
NEIGH_PX  = int(round(30.0 * PX_PER_UM))   # 30 μm ~= 167 px (Test 2)
PORTAL_PX = int(round(200.0 * PX_PER_UM))  # 200 μm ~= 1111 px (Test 3 inner)
PAREN_PX  = int(round(500.0 * PX_PER_UM))  # 500 μm ~= 2778 px (Test 3 outer)
DBSCAN_EPS_PX = 50                          # Test 4: spec'd 50 px
DBSCAN_MIN_SAMPLES = 5                      # spec
LIPOGR_CENTER_PX = 100                      # Test 4 composition radius (100 px)

# Test 1
KNN_K        = 10
N_PERM       = 200          # permutations per sample
SEED         = 42

# Sample-level disease assignment (from 42b loader)
SAMPLE_DISEASE = {
    "Leuven_1": "MASH",          # 2 end-stage MASH explants
    "Leuven_2": "no_MASH",       # 1 Normal + 1 MASL + 2 F3 (mixed)
    "Leuven_3": "MASH",          # F1 + F4
    "Leuven_4": "MASH",          # 2 F3 + 1 MASL
}

# ─── Helpers ─────────────────────────────────────────────────────────────────


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_raw_coords() -> pd.DataFrame:
    """Concatenate raw per-cell CenterX/Y_global_px from CosMx metadata CSVs.

    Returns a DataFrame indexed by obs_name (= f"{sample_id}_FOV{fov}_C{cell_ID}")
    with columns ['centerX_px', 'centerY_px'].
    """
    rows = []
    for sid, fp in RAW_META.items():
        _log(f"  reading raw coords {sid}: {fp.name}")
        meta = pd.read_csv(fp, usecols=[
            "fov", "cell_ID", "CenterX_global_px", "CenterY_global_px"
        ])
        # Build obs_name to match the h5ad index
        obs_name = (
            sid + "_FOV"
            + meta["fov"].astype(int).astype(str)
            + "_C"
            + meta["cell_ID"].astype(int).astype(str)
        )
        coords = pd.DataFrame({
            "obs_name":   obs_name.values,
            "centerX_px": meta["CenterX_global_px"].astype("float64").values,
            "centerY_px": meta["CenterY_global_px"].astype("float64").values,
        })
        rows.append(coords)
    out = pd.concat(rows, ignore_index=True)
    out = out.drop_duplicates(subset=["obs_name"]).set_index("obs_name")
    _log(f"  raw coord table: {out.shape[0]:,} cells; "
         f"X=[{out['centerX_px'].min():.1f},{out['centerX_px'].max():.1f}], "
         f"Y=[{out['centerY_px'].min():.1f},{out['centerY_px'].max():.1f}]")
    return out


def attach_coords(adata: ad.AnnData, coords: pd.DataFrame) -> ad.AnnData:
    """Rejoin per-cell CenterX/Y_global_px from raw metadata."""
    n_before = adata.n_obs
    matched = coords.reindex(adata.obs_names)
    n_match = int(matched.notna().all(axis=1).sum())
    _log(f"  matched coords: {n_match:,}/{n_before:,} cells "
         f"({100*n_match/n_before:.2f}%)")
    if n_match < 0.9 * n_before:
        sys.exit("ERROR: coord matching < 90%, check obs_name format")
    adata.obs["centerX_px"] = matched["centerX_px"].astype("float64").values
    adata.obs["centerY_px"] = matched["centerY_px"].astype("float64").values
    adata.obsm["spatial"]   = adata.obs[["centerX_px", "centerY_px"]].to_numpy(dtype="float64")
    return adata


def get_gene_lognorm(adata: ad.AnnData, gene: str) -> np.ndarray:
    """Return 1-D log-norm expression for `gene` (uses .layers['lognorm'])."""
    if gene not in adata.var_names:
        sys.exit(f"Gene {gene} missing from h5ad var_names")
    idx = adata.var_names.get_loc(gene)
    L = adata.layers["lognorm"]
    if sparse.issparse(L):
        col = L[:, idx].toarray().ravel()
    else:
        col = np.asarray(L[:, idx]).ravel()
    return col.astype("float64")


# ─── TEST 1: GPNMB+ macrophages near IL32+ hepatocytes ──────────────────────


def test1_gpnmb_il32(adata: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    _log("=" * 70)
    _log("TEST 1: GPNMB+ macrophages co-localize with IL32+ hepatocytes")
    _log("=" * 70)
    rng = np.random.default_rng(SEED)

    log_gpnmb = get_gene_lognorm(adata, "GPNMB")
    log_il32  = get_gene_lognorm(adata, "IL32")
    adata.obs["log_GPNMB"] = log_gpnmb
    adata.obs["log_IL32"]  = log_il32

    # GPNMB+ threshold: per-sample median + SD, AND cell_type==KC (macrophage lineage)
    # IL32+ threshold: per-sample median + SD AND cell_type==Hepatocyte
    rows = []
    for sid in sorted(SAMPLE_DISEASE.keys()):
        smask = (adata.obs["sample_id"].astype(str) == sid).values
        sub = adata[smask].copy()
        n_total = sub.n_obs

        # Drop cells without coordinates
        good = np.isfinite(sub.obs["centerX_px"].values) & np.isfinite(sub.obs["centerY_px"].values)
        sub = sub[good].copy()
        n_good = sub.n_obs
        _log(f"  [{sid}] kept {n_good:,}/{n_total:,} cells with coords")

        # Per-sample expression thresholds
        gpnmb_thr = float(np.median(sub.obs["log_GPNMB"].values) + np.std(sub.obs["log_GPNMB"].values))
        il32_thr  = float(np.median(sub.obs["log_IL32"].values)  + np.std(sub.obs["log_IL32"].values))

        is_gpnmb_pos = (
            (sub.obs["cell_type"].astype(str) == "KC").values
            & (sub.obs["log_GPNMB"].values > gpnmb_thr)
        )
        is_il32_pos_hep = (
            (sub.obs["cell_type"].astype(str) == "Hepatocyte").values
            & (sub.obs["log_IL32"].values > il32_thr)
        )

        n_gpnmb = int(is_gpnmb_pos.sum())
        n_il32  = int(is_il32_pos_hep.sum())

        _log(f"  [{sid}] GPNMB+_thr={gpnmb_thr:.3f}  IL32+_thr={il32_thr:.3f}")
        _log(f"  [{sid}] GPNMB+ KC = {n_gpnmb:,}; IL32+ Hep = {n_il32:,}")

        if n_gpnmb < 5 or n_il32 < 5:
            _log(f"  [{sid}] SKIP — group size too low")
            rows.append({
                "sample_id": sid,
                "disease_stage": SAMPLE_DISEASE[sid],
                "n_gpnmb_kc": n_gpnmb,
                "n_il32_hep": n_il32,
                "mean_neighbor_il32_frac": np.nan,
                "null_mean": np.nan,
                "null_sd": np.nan,
                "z_score": np.nan,
                "p_perm_one_tailed": np.nan,
                "global_il32_hep_frac": np.nan,
            })
            continue

        coords = sub.obs[["centerX_px", "centerY_px"]].to_numpy()
        # kNN on all sample cells
        nn = NearestNeighbors(n_neighbors=KNN_K + 1, metric="euclidean", algorithm="ball_tree", n_jobs=-1)
        nn.fit(coords)
        # Neighbors for GPNMB+ KC cells
        gp_idx = np.where(is_gpnmb_pos)[0]
        dist, ind = nn.kneighbors(coords[gp_idx], return_distance=True)
        # Drop self (first column is self)
        ind = ind[:, 1: KNN_K + 1]

        # F147 fix — composition-conditioned null. The previous null permuted the
        # IL32+Hep indicator across ALL cells, so null_mean collapsed to the global
        # IL32+Hep fraction. Because KC cells are strongly self-clustered, the
        # neighbors of a GPNMB+ KC are mostly other KCs (not hepatocytes), so the
        # observed neighbor-Hep fraction is structurally BELOW the global fraction
        # — the test reported DEPLETION as a pure cell-type-topology artifact.
        #
        # Conditioned statistic (option (a)/(c)): among the HEPATOCYTE neighbors of
        # each GPNMB+ KC, what fraction are IL32+? The null permutes the IL32+ label
        # ONLY within hepatocytes (cell-type field fixed), so the comparator is
        # "are nearby Hep more IL32+ than random Hep". This is the biologically
        # meaningful co-localization test.
        is_hep = (sub.obs["cell_type"].astype(str) == "Hepatocyte").values
        is_hep_int = is_hep.astype(np.int32)
        # IL32+ indicator restricted to hepatocytes (1 only for IL32+ Hep cells).
        il32_within_hep = is_il32_pos_hep.astype(np.int32)
        hep_idx = np.where(is_hep)[0]

        neigh_is_hep = is_hep_int[ind]                 # (n_gp, K) Hep neighbor mask
        neigh_is_il32hep = il32_within_hep[ind]        # (n_gp, K) IL32+Hep neighbor mask
        n_hep_neigh = neigh_is_hep.sum(axis=1)         # Hep neighbors per GPNMB+ KC
        n_il32hep_neigh = neigh_is_il32hep.sum(axis=1)
        # Pool across GPNMB+ KC: fraction of their Hep neighbors that are IL32+.
        total_hep_neigh = int(n_hep_neigh.sum())
        if total_hep_neigh > 0:
            obs_mean = float(n_il32hep_neigh.sum() / total_hep_neigh)
        else:
            obs_mean = np.nan

        # Null: fix the cell-type field; permute the IL32+ status among hepatocytes
        # only. Each permutation reassigns which hepatocytes are IL32+ (same count),
        # leaving Hep positions — and therefore which neighbors are Hep — unchanged.
        n_il32_hep_total = int(il32_within_hep.sum())
        null_means = np.empty(N_PERM, dtype="float64")
        n_cells_sub = is_hep_int.shape[0]
        for p in range(N_PERM):
            permuted = np.zeros(n_cells_sub, dtype=np.int32)
            if len(hep_idx) > 0 and n_il32_hep_total > 0:
                chosen = rng.choice(hep_idx, size=n_il32_hep_total, replace=False)
                permuted[chosen] = 1
            null_il32hep_neigh = permuted[ind].sum(axis=1)
            null_means[p] = (null_il32hep_neigh.sum() / total_hep_neigh
                             if total_hep_neigh > 0 else np.nan)

        if total_hep_neigh > 0:
            null_mean = float(np.nanmean(null_means))
            null_sd   = float(np.nanstd(null_means, ddof=1))
            z         = (obs_mean - null_mean) / null_sd if null_sd > 0 else np.nan
            # one-tailed test for ENRICHMENT (obs > null), matching co-loc hypothesis
            p_perm    = float((null_means >= obs_mean).sum() + 1) / (N_PERM + 1)
        else:
            # No GPNMB+ KC had any hepatocyte neighbor — statistic undefined.
            null_mean = null_sd = z = p_perm = np.nan
        # global IL32+ fraction *among hepatocytes* (denominator now Hep, not all cells)
        global_il32_hep_frac = (float(il32_within_hep.sum() / is_hep_int.sum())
                                if is_hep_int.sum() > 0 else np.nan)

        _log(f"  [{sid}] mean_neighbor_il32+frac = {obs_mean:.4f} "
             f"null = {null_mean:.4f} +/- {null_sd:.4f}  Z = {z:.2f}  p_perm = {p_perm:.4f}")
        _log(f"  [{sid}] global IL32+ Hep frac = {global_il32_hep_frac:.4f}")

        rows.append({
            "sample_id": sid,
            "disease_stage": SAMPLE_DISEASE[sid],
            "gpnmb_threshold_lognorm": gpnmb_thr,
            "il32_threshold_lognorm":  il32_thr,
            "n_gpnmb_kc": n_gpnmb,
            "n_il32_hep": n_il32,
            "global_il32_hep_frac": global_il32_hep_frac,
            "mean_neighbor_il32_frac": obs_mean,
            "null_mean": null_mean,
            "null_sd":   null_sd,
            "z_score":   z,
            "p_perm_one_tailed": p_perm,
        })

    df = pd.DataFrame(rows)
    out = out_dir / "gpnmb_il32_colocalization.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    _log(f"  wrote {out.name}")
    return df


# ─── TEST 2: neighborhood enrichment (squidpy) ──────────────────────────────


def test2_nhood(adata: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    """Pairwise neighborhood-enrichment Z-scores via squidpy.

    Build a radius-based neighbor graph per sample at 30 μm ~= 167 px and
    permutation-test cell-type pair adjacency (squidpy default: 1000 perms,
    we cap at 200 to stay within wall time).
    """
    import squidpy as sq
    _log("=" * 70)
    _log("TEST 2: squidpy nhood_enrichment per sample")
    _log("=" * 70)
    long_rows: List[Dict] = []

    for sid in sorted(SAMPLE_DISEASE.keys()):
        smask = (adata.obs["sample_id"].astype(str) == sid).values
        sub = adata[smask].copy()
        good = np.isfinite(sub.obs["centerX_px"].values) & np.isfinite(sub.obs["centerY_px"].values)
        sub = sub[good].copy()
        sub.obs["cell_type"] = sub.obs["cell_type"].astype("category")
        # Drop unused cats
        sub.obs["cell_type"] = sub.obs["cell_type"].cat.remove_unused_categories()
        sub.obsm["spatial"] = sub.obs[["centerX_px", "centerY_px"]].to_numpy(dtype="float64")
        _log(f"  [{sid}] n_cells={sub.n_obs:,}, cell types: {list(sub.obs['cell_type'].cat.categories)}")

        # Radius-based neighbor graph (cells within 30 μm = 167 px)
        sq.gr.spatial_neighbors(sub, coord_type="generic", radius=NEIGH_PX, delaunay=False)
        # Neighborhood enrichment z-scores (label permutation)
        sq.gr.nhood_enrichment(
            sub, cluster_key="cell_type",
            n_perms=200, seed=SEED, show_progress_bar=False,
        )

        zscores = sub.uns["cell_type_nhood_enrichment"]["zscore"]
        counts  = sub.uns["cell_type_nhood_enrichment"]["count"]
        cats    = list(sub.obs["cell_type"].cat.categories)
        for i, ct1 in enumerate(cats):
            for j, ct2 in enumerate(cats):
                long_rows.append({
                    "sample_id": sid,
                    "disease_stage": SAMPLE_DISEASE[sid],
                    "cell_type1": ct1,
                    "cell_type2": ct2,
                    "z_score": float(zscores[i, j]),
                    "count":   int(counts[i, j]),
                })
        _log(f"  [{sid}] computed {len(cats)*len(cats)} pair Z-scores")

    df = pd.DataFrame(long_rows)
    out = out_dir / "nhood_enrichment_zscores.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.4g")
    _log(f"  wrote {out.name}")

    # Aggregated summary text
    agg = (df.groupby(["cell_type1", "cell_type2"])
             .agg(mean_z=("z_score", "mean"),
                  median_z=("z_score", "median"),
                  n_samples=("z_score", "count"))
             .reset_index())
    # symmetric pairs — use sorted-pair key
    agg["pair_key"] = agg.apply(
        lambda r: " - ".join(sorted([r["cell_type1"], r["cell_type2"]])), axis=1
    )
    # For top/bottom, take MAX (more enriched of the two orderings) — or just keep all rows
    top_enr = agg.sort_values("mean_z", ascending=False).head(5)
    top_dep = agg.sort_values("mean_z", ascending=True).head(5)
    txt = ["Test 2 summary — squidpy nhood_enrichment mean Z across 4 samples",
           "(radius = 167 px ~= 30 microns; 200 label perms; seed=42)",
           "",
           "Top 5 ENRICHED ordered pairs (mean Z desc):"]
    for _, r in top_enr.iterrows():
        txt.append(f"  {r['cell_type1']:>14} ~ {r['cell_type2']:<14} mean Z = {r['mean_z']:+.2f} "
                   f"median = {r['median_z']:+.2f} (n={r['n_samples']})")
    txt += ["", "Top 5 DEPLETED ordered pairs (mean Z asc):"]
    for _, r in top_dep.iterrows():
        txt.append(f"  {r['cell_type1']:>14} ~ {r['cell_type2']:<14} mean Z = {r['mean_z']:+.2f} "
                   f"median = {r['median_z']:+.2f} (n={r['n_samples']})")
    txt += ["", "Per-sample KC×KC and KC×Hepatocyte Z-scores:"]
    for sid in sorted(SAMPLE_DISEASE.keys()):
        sub_df = df[df["sample_id"] == sid]
        kckc = sub_df[(sub_df["cell_type1"] == "KC") & (sub_df["cell_type2"] == "KC")]
        kchep = sub_df[(sub_df["cell_type1"] == "KC") & (sub_df["cell_type2"] == "Hepatocyte")]
        kckc_z = float(kckc["z_score"].iloc[0]) if not kckc.empty else np.nan
        kchep_z = float(kchep["z_score"].iloc[0]) if not kchep.empty else np.nan
        txt.append(f"  {sid}: KCxKC Z={kckc_z:+.2f}  KCxHep Z={kchep_z:+.2f}")

    out2 = out_dir / "nhood_enrichment_summary.txt"
    out2.write_text("\n".join(txt) + "\n")
    _log(f"  wrote {out2.name}")
    return df


# ─── TEST 3: portal-tract vs parenchymal segregation ────────────────────────


def test3_portal(adata: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    _log("=" * 70)
    _log("TEST 3: portal vs parenchymal (cholangiocyte distance, 200 / 500 μm)")
    _log("=" * 70)
    rows = []
    for sid in sorted(SAMPLE_DISEASE.keys()):
        smask = (adata.obs["sample_id"].astype(str) == sid).values
        sub = adata[smask].copy()
        good = np.isfinite(sub.obs["centerX_px"].values) & np.isfinite(sub.obs["centerY_px"].values)
        sub = sub[good].copy()
        coords = sub.obs[["centerX_px", "centerY_px"]].to_numpy()
        ctypes = sub.obs["cell_type"].astype(str).values
        chol_mask = (ctypes == "Cholangiocyte")
        n_chol = int(chol_mask.sum())
        if n_chol < 5:
            _log(f"  [{sid}] SKIP — only {n_chol} cholangiocytes")
            rows.append({
                "sample_id": sid, "disease_stage": SAMPLE_DISEASE[sid],
                "n_chol": n_chol, "skipped": True
            })
            continue
        tree = BallTree(coords[chol_mask], metric="euclidean")
        d_nearest_chol, _ = tree.query(coords, k=1, return_distance=True)
        d_nearest_chol = d_nearest_chol.ravel()  # px

        is_portal = (d_nearest_chol <= PORTAL_PX)
        is_paren  = (d_nearest_chol >  PAREN_PX)

        kc_mask = (ctypes == "KC")
        # 2x2: KC vs not-KC  ×  portal vs parenchymal (exclude in-between band)
        a = int((is_portal & kc_mask).sum())
        b = int((is_portal & ~kc_mask).sum())
        c = int((is_paren  & kc_mask).sum())
        d = int((is_paren  & ~kc_mask).sum())
        n_band = int((~is_portal & ~is_paren).sum())

        odds_ratio = np.nan
        p_val = np.nan
        if a + b > 0 and c + d > 0:
            table = np.array([[a, b], [c, d]])
            odds_ratio, p_val = fisher_exact(table, alternative="two-sided")

        frac_kc_portal = a / (a + b) if (a + b) else np.nan
        frac_kc_paren  = c / (c + d) if (c + d) else np.nan
        log2_or = float(np.log2(odds_ratio)) if (odds_ratio is not None and odds_ratio > 0 and not np.isnan(odds_ratio)) else np.nan

        _log(f"  [{sid}] n_chol={n_chol:,} portal_KC={a:,}/{a+b:,} ({100*frac_kc_portal:.2f}%) "
             f"parench_KC={c:,}/{c+d:,} ({100*frac_kc_paren:.2f}%)  band={n_band:,}  "
             f"OR={odds_ratio:.3f}  p={p_val:.2e}")
        rows.append({
            "sample_id": sid,
            "disease_stage": SAMPLE_DISEASE[sid],
            "n_chol": n_chol,
            "n_portal":  a + b,
            "n_parenchymal": c + d,
            "n_band":   n_band,
            "kc_portal":  a,
            "non_kc_portal": b,
            "kc_paren":  c,
            "non_kc_paren": d,
            "frac_kc_portal": frac_kc_portal,
            "frac_kc_parenchymal": frac_kc_paren,
            "odds_ratio_kc_portal_vs_paren": odds_ratio,
            "log2_odds_ratio": log2_or,
            "fisher_p":       p_val,
            "skipped":        False,
        })

    df = pd.DataFrame(rows)
    out = out_dir / "portal_vs_parenchymal.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    _log(f"  wrote {out.name}")
    return df


# ─── TEST 4: lipogranuloma cluster detection (DBSCAN on KC) ─────────────────


def test4_lipogranuloma(adata: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    _log("=" * 70)
    _log("TEST 4: lipogranuloma DBSCAN (eps=50 px, min_samples=5; KC cells)")
    _log("=" * 70)
    rows = []
    for sid in sorted(SAMPLE_DISEASE.keys()):
        smask = (adata.obs["sample_id"].astype(str) == sid).values
        sub = adata[smask].copy()
        good = np.isfinite(sub.obs["centerX_px"].values) & np.isfinite(sub.obs["centerY_px"].values)
        sub = sub[good].copy()
        ctypes = sub.obs["cell_type"].astype(str).values
        coords = sub.obs[["centerX_px", "centerY_px"]].to_numpy()
        kc_mask = (ctypes == "KC")
        kc_coords = coords[kc_mask]
        n_kc = int(kc_mask.sum())
        if n_kc < DBSCAN_MIN_SAMPLES:
            _log(f"  [{sid}] SKIP — only {n_kc} KC cells")
            continue

        db = DBSCAN(eps=DBSCAN_EPS_PX, min_samples=DBSCAN_MIN_SAMPLES, n_jobs=-1)
        labels = db.fit_predict(kc_coords)
        unique = np.unique(labels[labels >= 0])
        n_clust = unique.shape[0]
        n_noise = int((labels == -1).sum())
        sizes = []
        for cl in unique:
            sizes.append(int((labels == cl).sum()))
        sizes = np.array(sizes) if sizes else np.array([0])
        mean_size = float(sizes.mean()) if sizes.size > 0 else np.nan
        max_size  = int(sizes.max()) if sizes.size > 0 else 0
        _log(f"  [{sid}] n_KC={n_kc:,}  n_clusters={n_clust}  mean_size={mean_size:.2f}  "
             f"max_size={max_size}  noise={n_noise:,}")

        # Per-cluster composition (top-3 by size): cells within LIPOGR_CENTER_PX of centroid
        if n_clust > 0:
            top_clusters = unique[np.argsort(-sizes)][:3]
            tree_all = BallTree(coords)
            for rank, cl in enumerate(top_clusters):
                pts = kc_coords[labels == cl]
                centroid = pts.mean(axis=0)
                n_cells_clust = int((labels == cl).sum())
                idx = tree_all.query_radius(centroid.reshape(1, -1), r=LIPOGR_CENTER_PX, return_distance=False)[0]
                comp_ctypes = ctypes[idx]
                comp_counts = pd.Series(comp_ctypes).value_counts().to_dict()
                row = {
                    "sample_id": sid,
                    "disease_stage": SAMPLE_DISEASE[sid],
                    "cluster_id":   int(cl),
                    "cluster_rank": rank + 1,
                    "n_KC_cluster": n_cells_clust,
                    "centroid_X_px": float(centroid[0]),
                    "centroid_Y_px": float(centroid[1]),
                    "n_total_in_100px":   int(len(idx)),
                    "n_Hepatocyte":       int(comp_counts.get("Hepatocyte", 0)),
                    "n_KC":               int(comp_counts.get("KC", 0)),
                    "n_Mesenchymal":      int(comp_counts.get("Mesenchymal", 0)),
                    "n_Endothelial":      int(comp_counts.get("Endothelial", 0)),
                    "n_Cholangiocyte":    int(comp_counts.get("Cholangiocyte", 0)),
                    "n_Lymphocyte":       int(comp_counts.get("Lymphocyte", 0)),
                    "n_migDC":            int(comp_counts.get("migDC", 0)),
                }
                rows.append(row)
        # Always include a per-sample summary row (cluster_id = -1)
        rows.append({
            "sample_id":     sid,
            "disease_stage": SAMPLE_DISEASE[sid],
            "cluster_id":    -1,
            "cluster_rank":  0,
            "n_KC_cluster":  n_kc,
            "centroid_X_px": np.nan,
            "centroid_Y_px": np.nan,
            "n_total_in_100px":   n_clust,                # repurposed: cluster count
            "n_Hepatocyte":       int(np.round(mean_size)),
            "n_KC":               max_size,
            "n_Mesenchymal":      n_noise,
            "n_Endothelial":      0,
            "n_Cholangiocyte":    0,
            "n_Lymphocyte":       0,
            "n_migDC":            0,
        })

    df = pd.DataFrame(rows)
    out = out_dir / "lipogranuloma_clusters.tsv"
    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    _log(f"  wrote {out.name}")
    return df


# ─── Summary markdown ────────────────────────────────────────────────────────


def write_summary(test1, test2, test3, test4, out_dir: pathlib.Path) -> None:
    md = []
    md.append("# Govaere 2026 CosMx — Spatial Niche Validation (script 42c)")
    md.append("")
    md.append("**Dataset.** 522,145 CosMx SMI cells (968-gene panel) across 4 Leuven slides; "
              "cell-type labels from loader 41_* (7 categories: Cholangiocyte, Endothelial, "
              "Hepatocyte, KC, Lymphocyte, Mesenchymal, migDC).")
    md.append("")
    md.append("**Disease assignment (slide-level).** Leuven_1/3/4 = MASH; Leuven_2 = no_MASH "
              "(mixed-with-normal slide). See 42b loader docstring.")
    md.append("")
    md.append("**Coordinate fix.** h5ad `.obsm['spatial']` / `.obs[centerX/Y_global_px]` are NaN "
              "in the current loader output; per-cell CenterX_global_px / CenterY_global_px "
              "rejoined in-memory from raw `GSM*_metadata_file.csv.gz`. h5ad on disk untouched.")
    md.append("")
    md.append("**Pixel scale.** NanoString CosMx SMI: 1 px = 0.18 μm. 30 μm = 167 px (Test 2), "
              "200 / 500 μm = 1111 / 2778 px (Test 3), 50 px = 9 μm (Test 4 DBSCAN eps).")
    md.append("")
    # Test 1
    md.append("## Test 1 — GPNMB+ macrophages ~ IL32+ hepatocytes (kNN + label permutation)")
    md.append("")
    md.append("Threshold: per-sample (median + 1 SD) on `log_GPNMB` (KC subset) and `log_IL32` "
              "(Hepatocyte subset). k=10 nearest neighbors (Euclidean px); null = label "
              f"permutation over IL32+ Hep indicator, N={N_PERM} perms.")
    md.append("")
    md.append("| sample | stage | n GPNMB+ KC | n IL32+ Hep | obs frac | null mean | Z | p_perm |")
    md.append("|---|---|---|---|---|---|---|---|")
    for _, r in test1.iterrows():
        md.append(f"| {r['sample_id']} | {r['disease_stage']} | "
                  f"{int(r['n_gpnmb_kc']):,} | {int(r['n_il32_hep']):,} | "
                  f"{r['mean_neighbor_il32_frac']:.4f} | {r['null_mean']:.4f} | "
                  f"{r['z_score']:.2f} | {r['p_perm_one_tailed']:.4g} |")
    z_agg = float(test1["z_score"].mean(skipna=True))
    md.append("")
    md.append(f"**Aggregate.** Mean Z across {test1['z_score'].notna().sum()} samples = **{z_agg:.2f}**.")
    md.append("")
    # Test 2
    md.append("## Test 2 — Cell-type neighborhood enrichment (squidpy)")
    md.append("")
    md.append("Spatial neighbors built at radius 167 px (~30 μm). Z-scores from 200 label "
              "permutations on `cell_type`. Per-sample matrix in `nhood_enrichment_zscores.tsv`.")
    md.append("")
    # Top enriched / depleted (ordered pairs by mean Z)
    agg = (test2.groupby(["cell_type1", "cell_type2"])
                .agg(mean_z=("z_score", "mean"),
                     median_z=("z_score", "median"),
                     n_samples=("z_score", "count"))
                .reset_index())
    top_e = agg.sort_values("mean_z", ascending=False).head(5)
    top_d = agg.sort_values("mean_z", ascending=True).head(5)
    md.append("Top 5 enriched ordered pairs (mean Z across samples):")
    md.append("")
    md.append("| cell_type1 | cell_type2 | mean Z | median Z | n |")
    md.append("|---|---|---|---|---|")
    for _, r in top_e.iterrows():
        md.append(f"| {r['cell_type1']} | {r['cell_type2']} | {r['mean_z']:+.2f} | "
                  f"{r['median_z']:+.2f} | {int(r['n_samples'])} |")
    md.append("")
    md.append("Top 5 depleted ordered pairs (mean Z across samples):")
    md.append("")
    md.append("| cell_type1 | cell_type2 | mean Z | median Z | n |")
    md.append("|---|---|---|---|---|")
    for _, r in top_d.iterrows():
        md.append(f"| {r['cell_type1']} | {r['cell_type2']} | {r['mean_z']:+.2f} | "
                  f"{r['median_z']:+.2f} | {int(r['n_samples'])} |")
    md.append("")
    # Test 3
    md.append("## Test 3 — Portal vs parenchymal KC enrichment")
    md.append("")
    md.append("Portal = within 200 μm of nearest cholangiocyte; parenchymal = beyond 500 μm. "
              "Fisher exact: KC vs non-KC × portal vs parenchymal.")
    md.append("")
    md.append("| sample | stage | n_chol | KC portal | KC paren | OR | log2 OR | p |")
    md.append("|---|---|---|---|---|---|---|---|")
    for _, r in test3.iterrows():
        if r.get("skipped"):
            md.append(f"| {r['sample_id']} | {r['disease_stage']} | {int(r['n_chol'])} | skipped | | | | |")
            continue
        md.append(f"| {r['sample_id']} | {r['disease_stage']} | {int(r['n_chol']):,} | "
                  f"{int(r['kc_portal']):,}/{int(r['n_portal']):,} "
                  f"({100*r['frac_kc_portal']:.2f}%) | "
                  f"{int(r['kc_paren']):,}/{int(r['n_parenchymal']):,} "
                  f"({100*r['frac_kc_parenchymal']:.2f}%) | "
                  f"{r['odds_ratio_kc_portal_vs_paren']:.2f} | "
                  f"{r['log2_odds_ratio']:+.2f} | {r['fisher_p']:.2e} |")
    md.append("")
    # Test 4
    md.append("## Test 4 — KC DBSCAN clusters (lipogranuloma proxy)")
    md.append("")
    md.append("DBSCAN eps = 50 px (~9 μm), min_samples = 5. Top-3 clusters per sample: "
              "cell-type composition within 100 px of cluster centroid.")
    md.append("")
    md.append("| sample | stage | n_KC | n_clusters | mean size | max size | n_noise |")
    md.append("|---|---|---|---|---|---|---|")
    # We stored per-sample summary in cluster_id=-1 rows
    summ = test4[test4["cluster_id"] == -1].copy()
    for _, r in summ.iterrows():
        md.append(f"| {r['sample_id']} | {r['disease_stage']} | {int(r['n_KC_cluster']):,} | "
                  f"{int(r['n_total_in_100px']):,} | {r['n_Hepatocyte']:.2f} | "
                  f"{int(r['n_KC']):,} | {int(r['n_Mesenchymal']):,} |")
    md.append("")
    md.append("**Top-3 cluster composition** (Hep / KC / Mesenchymal counts within 100 px of centroid):")
    md.append("")
    md.append("| sample | rank | cluster | n_KC_cluster | n_total | Hep | KC | Mes |")
    md.append("|---|---|---|---|---|---|---|---|")
    top3 = test4[test4["cluster_id"] != -1].copy()
    for _, r in top3.iterrows():
        md.append(f"| {r['sample_id']} | {int(r['cluster_rank'])} | {int(r['cluster_id'])} | "
                  f"{int(r['n_KC_cluster']):,} | {int(r['n_total_in_100px']):,} | "
                  f"{int(r['n_Hepatocyte']):,} | {int(r['n_KC']):,} | {int(r['n_Mesenchymal']):,} |")
    md.append("")
    md.append("## Interpretation")
    md.append("")
    md.append("These results are spatially resolved cross-checks of the paper's Fig 5 claims "
              "(GPNMB+ MetMac-IL32+ hepatocyte juxtaposition, KC clustering as a "
              "lipogranuloma proxy, portal-tract macrophage bias). Z-scores and Fisher ORs "
              "encode the empirical excess of co-localization vs label-permutation null; "
              "magnitudes are conservative because cell-type labels are coarsened to 7 "
              "categories (MetMac and TransMac collapsed into KC by the loader's leiden + "
              "argmax annotation step).")
    (out_dir / "niche_summary.md").write_text("\n".join(md) + "\n")
    _log(f"  wrote niche_summary.md")


# ─── Main ────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="42c spatial niche tests on Govaere CosMx")
    parser.add_argument("--adata",   type=str, default=str(ADATA_PATH))
    parser.add_argument("--out-dir", type=str, default=str(OUT_DIR))
    args = parser.parse_args()

    adata_path = pathlib.Path(args.adata)
    out_dir    = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _log("=" * 70)
    _log("42c_govaere2026_cosmx_spatial_niches.py")
    _log(f"  adata    : {adata_path}")
    _log(f"  out_dir  : {out_dir}")
    _log(f"  pix scale: 1 px = 0.18 μm  (NEIGH_PX={NEIGH_PX}  PORTAL_PX={PORTAL_PX}  PAREN_PX={PAREN_PX})")
    _log("=" * 70)

    _log(f"Loading h5ad ({adata_path.stat().st_size/1e6:.1f} MB) into memory")
    t0 = time.time()
    adata = sc.read_h5ad(adata_path)
    _log(f"  shape: {adata.shape}  (load {time.time()-t0:.1f}s)")

    # Replace NaN obsm['spatial'] with raw per-cell CenterX/Y_global_px
    coords = load_raw_coords()
    adata = attach_coords(adata, coords)

    # Sanity per sample
    for sid in sorted(SAMPLE_DISEASE.keys()):
        m = (adata.obs["sample_id"].astype(str) == sid).values
        xs = adata.obs.loc[m, "centerX_px"].values
        ys = adata.obs.loc[m, "centerY_px"].values
        good = np.isfinite(xs) & np.isfinite(ys)
        _log(f"  [{sid}] coords good = {good.sum():,}/{m.sum():,}  "
             f"X=[{np.nanmin(xs):.0f},{np.nanmax(xs):.0f}]  "
             f"Y=[{np.nanmin(ys):.0f},{np.nanmax(ys):.0f}]")

    # Run tests
    test1_df = test1_gpnmb_il32(adata, out_dir)
    test2_df = test2_nhood(adata, out_dir)
    test3_df = test3_portal(adata, out_dir)
    test4_df = test4_lipogranuloma(adata, out_dir)

    # Summary markdown
    write_summary(test1_df, test2_df, test3_df, test4_df, out_dir)
    _log("DONE")


if __name__ == "__main__":
    main()
