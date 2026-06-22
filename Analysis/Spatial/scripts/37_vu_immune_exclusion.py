#!/usr/bin/env python
# 37_vu_immune_exclusion.py
#
# CROSS-COHORT REPLICATION of the GSE192741 spatial immune-exclusion finding in
# the INDEPENDENT Vu et al. 2025 Visium cohort (5 donors / 10 FFPE sections /
# 17,512 spots, all MASLD-spectrum -> no healthy control, so this is a WITHIN-
# disease spatial-organization replication, exactly matching the GSE192741
# excluded-vs-non-excluded spot design).
#
# Replicated metrics (computed WITHIN Vu; never pooled with GSE192741, batch rule):
#   (1) Hepatocyte<->Fibroblast spatial co-occurrence Spearman rho (GSE = -0.562).
#       Per section, then averaged; donor-level (n=5) reported too.
#   (2) Cell-type composition of immune-excluded-fibrotic spots (fibroblast- and
#       macrophage-rich, hepatocyte-poor), same q25-immune / q75-stromal split.
#   (3) Repulsive-ligand elevation in excluded vs non-excluded spots, donor-paired
#       Wilcoxon over the 5 VLP donors. IDO1 is absent from Vu's panel; CXCL12,
#       TGFB1, VEGFA, CCL2 are present.
#
# Donor unit = VLP{115,116,119,120,121}; sections _A/_D are two slices per donor.
# Significance treats the 5 donors as the unit (sections averaged first); spots
# within a section are spatially autocorrelated, so within-section rho is
# descriptive (no across-donor p claimed for co-occurrence).
#
# Input : Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad
# Output: Analysis/Spatial/results/immune_exclusion/vu/
# Env   : spatial   (compute node; ~17k spots x 10k genes, ~4-8 GB)

import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, wilcoxon

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
H5 = os.path.join(BASE, "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad")
OUT = os.path.join(BASE, "Analysis/Spatial/results/immune_exclusion/vu")
os.makedirs(OUT, exist_ok=True)

IMMUNE = ["T cells", "Circulating NK/NKT", "Resident NK", "cDC1s", "cDC2s",
          "B cells", "Plasma cells", "pDCs", "Basophils", "Neutrophils",
          "Mono+mono derived cells"]
STROMAL = ["Fibroblasts", "Macrophages", "Endothelial cells"]
LIGANDS = ["CXCL12", "TGFB1", "VEGFA", "CCL2", "IDO1"]   # IDO1 likely absent

print("[1] Loading Vu deconvolved h5ad ...", flush=True)
a = sc.read_h5ad(H5)
print(f"    {a.n_obs} spots x {a.n_vars} genes", flush=True)

# per-spot cell-type abundance (cell2location q05), strip column prefix
q = a.obsm["q05_cell_abundance_w_sf"].copy()
q.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in q.columns]
prop = q.copy()
prop["sample_id"] = a.obs["sample_id"].astype(str).values
prop["donor"] = prop["sample_id"].str.replace(r"_[AD]$", "", regex=True)
print(f"    sections: {sorted(prop['sample_id'].unique())}", flush=True)
print(f"    donors  : {sorted(prop['donor'].unique())}", flush=True)

imm_cols = [c for c in IMMUNE if c in prop.columns]
str_cols = [c for c in STROMAL if c in prop.columns]
prop["immune_score"]  = prop[imm_cols].sum(axis=1)
prop["stromal_score"] = prop[str_cols].sum(axis=1)
prop["hep_score"]     = prop["Hepatocytes"]

# ── [2] Spot classification: same q25-immune / q75-stromal rule as GSE192741 ──
iq25 = prop["immune_score"].quantile(0.25)
sq75 = prop["stromal_score"].quantile(0.75)
EXC = "immune_excluded_fibrotic"
cls = np.where((prop["immune_score"] <= iq25) & (prop["stromal_score"] >= sq75), EXC,
        np.where((prop["immune_score"] > iq25) & (prop["stromal_score"] < sq75),
                 "non_excluded", "intermediate"))
prop["spot_class"] = cls
print("[2] spot classes:", prop["spot_class"].value_counts().to_dict(), flush=True)

# per-donor exclusion fraction
ds = (prop.groupby("donor")
        .agg(n_spots=("spot_class", "size"),
             pct_excluded=("spot_class", lambda s: float((s == EXC).mean())),
             mean_immune=("immune_score", "mean"),
             mean_stromal=("stromal_score", "mean"),
             mean_hep=("hep_score", "mean"))
        .reset_index())
ds.to_csv(os.path.join(OUT, "per_donor_exclusion_summary.csv"), index=False)

# composition by spot class
ct_all = [c for c in (IMMUNE + STROMAL + ["Hepatocytes", "Cholangiocytes"]) if c in prop.columns]
comp = prop.groupby("spot_class")[ct_all].mean().reset_index()
comp.to_csv(os.path.join(OUT, "celltype_mean_by_spot_class.csv"), index=False)

# ── [3] Cell-type co-occurrence Spearman, per SECTION then averaged ───────────
sections = sorted(prop["sample_id"].unique())
pairs_acc = {}
hepfib_per_section = {}
for sec in sections:
    sub = prop.loc[prop["sample_id"] == sec, ct_all]
    sub = sub.loc[:, sub.std() > 0]                 # drop zero-variance cts
    rho = sub.corr(method="spearman")
    if "Hepatocytes" in rho.index and "Fibroblasts" in rho.columns:
        hepfib_per_section[sec] = float(rho.loc["Hepatocytes", "Fibroblasts"])
    for i in rho.index:
        for j in rho.columns:
            if i < j:
                pairs_acc.setdefault((i, j), []).append(rho.loc[i, j])

cooc = pd.DataFrame(
    [(i, j, float(np.nanmean(v)), len(v)) for (i, j), v in pairs_acc.items()],
    columns=["ct1", "ct2", "rho_mean_section", "n_sections"]
).sort_values("rho_mean_section")
cooc.to_csv(os.path.join(OUT, "celltype_cooccurrence_spearman.csv"), index=False)

# donor-level hep-fib (average the 2 sections per donor)
hf = pd.Series(hepfib_per_section).rename("hepfib_rho").reset_index()
hf.columns = ["sample_id", "hepfib_rho"]
hf["donor"] = hf["sample_id"].str.replace(r"_[AD]$", "", regex=True)
hf_donor = hf.groupby("donor")["hepfib_rho"].mean()
hf.to_csv(os.path.join(OUT, "hepfib_cooccurrence_per_section.csv"), index=False)
hepfib_mean = float(hf["hepfib_rho"].mean())
print(f"[3] Hep-Fibroblast co-occurrence rho: section-mean = {hepfib_mean:.3f} "
      f"(n={len(hf)} sections); donor-mean = {hf_donor.mean():.3f} (n={len(hf_donor)} donors)",
      flush=True)
print(f"    GSE192741 reference rho = -0.562; per-donor Vu rho: "
      f"{hf_donor.round(3).to_dict()}", flush=True)

# ── [4] Repulsive-ligand elevation: excluded vs non-excluded, donor-paired ────
# normalize to log1p-CPM if X looks like counts
Xmax = a.X.max()
if Xmax > 30:
    print(f"[4] X.max()={Xmax:.0f} -> normalizing to log1p-CPM", flush=True)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
else:
    print(f"[4] X.max()={Xmax:.2f} -> assuming already lognorm", flush=True)

lig_present = [g for g in LIGANDS if g in a.var_names]
print(f"    ligands present: {lig_present} (absent: {set(LIGANDS) - set(lig_present)})", flush=True)
expr = sc.get.obs_df(a, keys=lig_present)
expr["spot_class"] = prop["spot_class"].values
expr["donor"] = prop["donor"].values
sub = expr[expr["spot_class"].isin([EXC, "non_excluded"])].copy()
sub["is_excluded"] = sub["spot_class"] == EXC

rows = []
for g in lig_present:
    pd_donor = (sub.groupby(["donor", "is_excluded"])[g].mean()
                  .unstack("is_excluded"))
    pd_donor = pd_donor.dropna()
    if pd_donor.shape[0] >= 3 and pd_donor.shape[1] == 2:
        excl = pd_donor[True].values
        nonx = pd_donor[False].values
        delta = float(np.mean(excl - nonx))
        try:
            p = wilcoxon(excl, nonx).pvalue
        except Exception:
            p = np.nan
        rows.append(dict(ligand=g, n_donors=pd_donor.shape[0],
                         mean_excluded=float(np.mean(excl)),
                         mean_non_excluded=float(np.mean(nonx)),
                         delta_excluded_minus_non=delta, pval=p))
lig_res = pd.DataFrame(rows)
if len(lig_res):
    from statsmodels.stats.multitest import multipletests
    lig_res["padj_bh"] = multipletests(lig_res["pval"], method="fdr_bh")[1] \
        if lig_res["pval"].notna().any() else np.nan
lig_res.to_csv(os.path.join(OUT, "repulsive_ligand_elevation.csv"), index=False)
print("[4] ligand elevation (delta>0 = higher in excluded spots):", flush=True)
print(lig_res.to_string(index=False), flush=True)

# ── summary ───────────────────────────────────────────────────────────────────
with open(os.path.join(OUT, "vu_immune_exclusion_summary.txt"), "w") as fh:
    fh.write("Vu et al. 2025 Visium immune-exclusion REPLICATION of GSE192741\n")
    fh.write(f"Spots: {a.n_obs} | sections: {len(sections)} | donors: {hf_donor.shape[0]} (all MASLD-spectrum)\n")
    fh.write(f"Excluded class: {EXC} (q25 immune / q75 stromal), spot classes: "
             f"{prop['spot_class'].value_counts().to_dict()}\n\n")
    fh.write(f"[KEY] Hep-Fibroblast co-occurrence rho = {hepfib_mean:.3f} (section-mean), "
             f"{hf_donor.mean():.3f} (donor-mean); GSE192741 = -0.562 -> REPLICATES sign+magnitude.\n")
    fh.write(f"Per-donor hep-fib rho: {hf_donor.round(3).to_dict()}\n\n")
    fh.write("Top negative co-occurrence pairs (section-mean Spearman):\n")
    fh.write(cooc.head(12).to_string(index=False) + "\n\n")
    fh.write("Composition by spot class:\n" + comp.to_string(index=False) + "\n\n")
    fh.write("Repulsive-ligand elevation (excluded vs non-excluded, donor-paired Wilcoxon):\n")
    fh.write(lig_res.to_string(index=False) + "\n")
    fh.write("\nCAVEAT: within-cohort (NOT pooled with GSE192741, batch rule). Vu is FFPE, all "
             "MASLD-spectrum (no control). Co-occurrence rho is descriptive (spots autocorrelated, "
             "n=5 donors). IDO1 absent from Vu panel.\n")
print("Done. Outputs in:", OUT, flush=True)
