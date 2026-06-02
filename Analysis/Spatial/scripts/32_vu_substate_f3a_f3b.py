#!/usr/bin/env python3
"""
32_vu_substate_f3a_f3b.py — F3 sub-state (spatial): per-spot F3a / F3b in Vu et al.

Tests whether the Li et al. F3a / F3b sub-state distinction is detectable at
spot resolution in the Vu et al. 2025 MASLD-spectrum Visium dataset (10 arrays,
5 individuals, 17,512 spots). Phase 1 bulk analysis ruled out discrete bulk-level
F3 sub-states. This script tests interpretation (3) — that sub-state biology
is real at spot resolution, averaged out at biopsy resolution.

The Vu et al. AnnData currently lacks per-individual fibrosis stage in the
pre-processed merge; this analysis is therefore stage-agnostic. We test whether
spots within MASLD-spectrum tissue separate on the F3a vs F3b axis in an
unbiased way (without conditioning on patient F3 status), and whether spatial
organization of F3a-rich vs F3b-rich regions exists.

Steps:
  1. Load Vu spatial AnnData (17,512 spots × 13,335 genes).
  2. Score every spot for Li F3a (UPR / SREBP / cholesterol / lipid) and
     Li F3b (ECM / senescence / Ig) signatures using scanpy.tl.score_genes.
  3. Test bimodality of F3a, F3b, and F3a-F3b mixing axis at spot level
     (Hartigan dip-MC + 2-Gaussian BIC).
  4. Per individual: compute spot-level F3a:F3b distribution.
  5. Spatial Moran's I on the F3a-F3b mixing axis (autocorrelation).
  6. Output spot-level CSV + per-individual summary + bimodality CSV.

Outputs (RNA-seq/results/granular_staging/):
  vu_spot_substate_scores.csv         — per-spot F3a/F3b scores + spatial coords
  vu_spot_bimodality.csv              — bimodality tests at spot resolution
  vu_per_individual_summary.csv       — per-patient F3a/F3b distribution stats
  vu_spatial_morans_i.csv             — spatial autocorrelation per signature

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00 (env: spatial)
"""
import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import stats

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INPUT_PATH = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/merged_spatial_vu.h5ad"
OUT_DIR = PROJECT_ROOT / "RNA-seq/results/granular_staging"
OUT_DIR.mkdir(parents=True, exist_ok=True)

print("== F3 sub-state (spatial) — Vu spot-level F3a/F3b ==")
print(f"Input: {INPUT_PATH}")
print(f"Output: {OUT_DIR}")

adata = sc.read_h5ad(INPUT_PATH)
print(f"Spots: {adata.n_obs}, Genes: {adata.n_vars}")
print(f"Individuals: {adata.obs['individual'].nunique()} ({list(adata.obs['individual'].unique())})")

# Verify spatial coordinates available
spatial_key = None
if "spatial" in adata.obsm:
    spatial_key = "spatial"
elif "X_spatial" in adata.obsm:
    spatial_key = "X_spatial"
print(f"Spatial coords key: {spatial_key}")

# Make sure we use log-normalized expression (already done in preprocessing)
print(f"Expression layer min/max: {adata.X.min():.2f} / {adata.X.max():.2f}")

# ----------------------------------------------------------------------------
# (1) Define Li et al. F3a / F3b signatures — same as 243_multi_anchor
# ----------------------------------------------------------------------------
li_f3a = [
    "HSPA5", "DDIT3", "ATF4", "ATF6", "EIF2AK3", "ERN1", "XBP1", "DNAJB9",
    "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "SQLE", "DHCR7", "DHCR24",
    "MVD", "MVK", "FASN", "ACACA", "SCD", "INSIG1",
    "FABP1", "ACOX1", "CPT1A", "PLIN2", "DGAT1",
]
li_f3b = [
    "IGFBP7", "BGN", "COL1A2", "COL3A1", "TIMP1",
    "COL1A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3",
    "FN1", "VCAN", "DCN", "LUM", "LOX", "LOXL1", "LOXL2",
    "MMP2", "MMP9", "TIMP2", "ACTA2", "TAGLN",
    "CDKN1A", "CDKN2A", "GLB1", "SERPINE1",
]

# Filter to genes present in adata
li_f3a_present = [g for g in li_f3a if g in adata.var_names]
li_f3b_present = [g for g in li_f3b if g in adata.var_names]
print(f"F3a genes present: {len(li_f3a_present)} / {len(li_f3a)}")
print(f"F3b genes present: {len(li_f3b_present)} / {len(li_f3b)}")

# ----------------------------------------------------------------------------
# (2) Score per spot
# ----------------------------------------------------------------------------
sc.tl.score_genes(adata, gene_list=li_f3a_present, score_name="li_f3a",
                  random_state=42)
sc.tl.score_genes(adata, gene_list=li_f3b_present, score_name="li_f3b",
                  random_state=42)

# Cohort-internal correction not needed (single dataset) but per-individual z
adata.obs["li_f3a_z"] = adata.obs.groupby("individual", observed=False)["li_f3a"].transform(
    lambda x: (x - x.mean()) / (x.std() + 1e-6)
)
adata.obs["li_f3b_z"] = adata.obs.groupby("individual", observed=False)["li_f3b"].transform(
    lambda x: (x - x.mean()) / (x.std() + 1e-6)
)
adata.obs["mixing_axis"] = adata.obs["li_f3b_z"] - adata.obs["li_f3a_z"]

# ----------------------------------------------------------------------------
# (3) Bimodality at spot level
# ----------------------------------------------------------------------------
def fit_gmm_2(x, n_iter=200, tol=1e-6):
    n = len(x)
    ord_ = np.sort(x)
    half = n // 2
    mu1, mu2 = ord_[:half].mean(), ord_[half:].mean()
    sd1, sd2 = max(ord_[:half].std(), 1e-3), max(ord_[half:].std(), 1e-3)
    p1 = 0.5
    ll_old = -np.inf
    for _ in range(n_iter):
        d1 = p1 * stats.norm.pdf(x, mu1, sd1)
        d2 = (1 - p1) * stats.norm.pdf(x, mu2, sd2)
        g1 = d1 / (d1 + d2 + 1e-300)
        p1 = g1.mean()
        mu1 = (g1 * x).sum() / g1.sum()
        mu2 = ((1 - g1) * x).sum() / (1 - g1).sum()
        sd1 = max(np.sqrt((g1 * (x - mu1) ** 2).sum() / g1.sum()), 1e-3)
        sd2 = max(np.sqrt(((1 - g1) * (x - mu2) ** 2).sum() / (1 - g1).sum()), 1e-3)
        ll = np.log(p1 * stats.norm.pdf(x, mu1, sd1) +
                    (1 - p1) * stats.norm.pdf(x, mu2, sd2) + 1e-300).sum()
        if abs(ll - ll_old) < tol:
            break
        ll_old = ll
    return dict(mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2, p1=p1, ll=ll)

def bic_compare(x):
    n = len(x)
    ll1 = stats.norm.logpdf(x, x.mean(), x.std()).sum()
    bic1 = -2 * ll1 + 2 * np.log(n)
    fit = fit_gmm_2(x)
    bic2 = -2 * fit["ll"] + 5 * np.log(n)
    sep = abs(fit["mu1"] - fit["mu2"]) / np.sqrt((fit["sd1"] ** 2 + fit["sd2"] ** 2) / 2)
    return dict(bic1=bic1, bic2=bic2, bic_diff=bic1 - bic2, sep_sigma=sep, p1=fit["p1"])

def hartigan_mc(x, n_mc=999, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    fx = stats.cumfreq(x, numbins=200)
    grid = np.linspace(x.min(), x.max(), 200)
    e = np.array([np.mean(x <= g) for g in grid])
    g_pdf = stats.norm.cdf(grid, x.mean(), x.std())
    obs_dip = float(np.max(np.abs(e - g_pdf)))
    null_dips = []
    for _ in range(n_mc):
        y = rng.normal(x.mean(), x.std(), len(x))
        ey = np.array([np.mean(y <= g) for g in grid])
        gy = stats.norm.cdf(grid, y.mean(), y.std())
        null_dips.append(np.max(np.abs(ey - gy)))
    p = (np.sum(np.array(null_dips) >= obs_dip) + 1) / (n_mc + 1)
    return obs_dip, p

bim_rows = []
for sig_name in ["li_f3a", "li_f3b", "li_f3a_z", "li_f3b_z", "mixing_axis"]:
    x = adata.obs[sig_name].values.astype(float)
    bic = bic_compare(x)
    obs_dip, dip_p = hartigan_mc(x, n_mc=499)  # smaller MC for speed at n=17k
    bim_rows.append({
        "signature": sig_name,
        "n_spots": len(x),
        "BIC_diff_2vs1": bic["bic_diff"],
        "GMM_sep_sigma": bic["sep_sigma"],
        "GMM_p1": bic["p1"],
        "dipMC_obs": obs_dip,
        "dipMC_p": dip_p,
    })
bim_df = pd.DataFrame(bim_rows)
bim_df.to_csv(OUT_DIR / "vu_spot_bimodality.csv", index=False)
print("\nBimodality at spot level:")
print(bim_df.to_string(index=False))

# ----------------------------------------------------------------------------
# (4) Per-individual summary
# ----------------------------------------------------------------------------
per_ind_summary = (
    adata.obs.groupby("individual", observed=False)
    .agg(
        n_spots=("li_f3a", "size"),
        f3a_mean=("li_f3a", "mean"), f3a_sd=("li_f3a", "std"),
        f3b_mean=("li_f3b", "mean"), f3b_sd=("li_f3b", "std"),
        f3a_f3b_corr=("li_f3a", lambda x: np.corrcoef(x, adata.obs.loc[x.index, "li_f3b"])[0, 1]),
        mixing_mean=("mixing_axis", "mean"), mixing_sd=("mixing_axis", "std"),
    )
    .reset_index()
)
per_ind_summary.to_csv(OUT_DIR / "vu_per_individual_summary.csv", index=False)
print("\nPer-individual F3a/F3b summary:")
print(per_ind_summary.to_string(index=False))

# ----------------------------------------------------------------------------
# (5) Spatial Moran's I per individual
# ----------------------------------------------------------------------------
moran_rows = []
for ind in adata.obs["individual"].unique():
    sub = adata[adata.obs["individual"] == ind].copy()
    if sub.n_obs < 50 or "spatial" not in sub.obsm:
        continue
    try:
        sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=6)
        # Score genes already done; need to compute autocorrelation on per-cell continuous scores.
        # squidpy spatial_autocorr operates on .X by default; for per-cell scalar scores, use
        # spatial_neighbors weight matrix and compute Moran's I manually.
        from libpysal.weights import KNN
        coords = sub.obsm["spatial"]
        from esda.moran import Moran
        kw = KNN.from_array(coords, k=6)
        mor_a = Moran(sub.obs["li_f3a"].values.astype(float), kw, permutations=99)
        mor_b = Moran(sub.obs["li_f3b"].values.astype(float), kw, permutations=99)
        mor_mix = Moran(sub.obs["mixing_axis"].values.astype(float), kw, permutations=99)
        moran_rows.extend([
            dict(individual=ind, signature="li_f3a", I=mor_a.I, p=mor_a.p_sim, n=sub.n_obs),
            dict(individual=ind, signature="li_f3b", I=mor_b.I, p=mor_b.p_sim, n=sub.n_obs),
            dict(individual=ind, signature="mixing", I=mor_mix.I, p=mor_mix.p_sim, n=sub.n_obs),
        ])
    except Exception as e:
        print(f"Moran's I skipped for {ind}: {e}")
moran_df = pd.DataFrame(moran_rows)
if len(moran_df) > 0:
    moran_df.to_csv(OUT_DIR / "vu_spatial_morans_i.csv", index=False)
    print("\nSpatial Moran's I:")
    print(moran_df.to_string(index=False))

# ----------------------------------------------------------------------------
# (6) Per-spot scatter table for downstream plotting
# ----------------------------------------------------------------------------
spot_df = adata.obs[["sample_id", "individual", "li_f3a", "li_f3b",
                     "li_f3a_z", "li_f3b_z", "mixing_axis"]].copy()
if "spatial" in adata.obsm:
    spot_df["x"] = adata.obsm["spatial"][:, 0]
    spot_df["y"] = adata.obsm["spatial"][:, 1]
spot_df.reset_index().rename(columns={"index": "spot_id"}).to_csv(
    OUT_DIR / "vu_spot_substate_scores.csv", index=False
)

print("\n== F3 sub-state (spatial) complete. ==")
