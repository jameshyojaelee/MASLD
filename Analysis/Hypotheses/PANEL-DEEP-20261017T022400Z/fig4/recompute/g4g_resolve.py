"""Fig 4G: resolve nearest-coordinate ties from g4g_geometry.py. Drawn nodes that matched the same geometry node are
re-assigned (optimal assignment over the duplicated and unmatched geometry nodes, cost = coordinate distance + |colour value
difference|) and re-checked against my recomputed beta. Usage: python g4g_resolve.py <tag>"""
import sys
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
tag = sys.argv[1]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
n = pd.read_csv(f"{out}/g4g_{tag}.tsv.gz", sep="\t")
N = pd.read_csv(f"{out}/r4g_nodes_367_grid.tsv.gz", sep="\t")
N = pd.concat([N[(N.ymode == "zdis") & (N.axis == "axis_fp") & (N.collection != "hotspot")],
               N[(N.ymode == "zall") & (N.axis == "axis_fp") & (N.collection == "hotspot")]])
G = N[N.component_role == "main_system_component"][["feature_id", "collection", "umap_1", "umap_2", "cont_beta"]]
G = G.merge(pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz", sep="\t")
            .query("component_role == 'main_system_component'")[["feature_id", "collection"]], on=["feature_id", "collection"], how="right")
lo, hi = -0.7983, 0.7983
dup = n[n.feature_id.duplicated(keep=False)]
un = G[~G.feature_id.isin(n.feature_id)]
cand = pd.concat([G[G.feature_id.isin(dup.feature_id)], un]).reset_index(drop=True)
geo = pd.read_csv("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz", sep="\t")
cand = cand.drop(columns=["umap_1", "umap_2"], errors="ignore").merge(geo[["feature_id", "collection", "umap_1", "umap_2"]], on=["feature_id", "collection"])
C = np.hypot(dup.u1.values[:, None] - cand.umap_1.values[None], dup.u2.values[:, None] - cand.umap_2.values[None]) \
    + np.abs(dup.beta_drawn.values[:, None] - cand.cont_beta.clip(lo, hi).fillna(99).values[None])
ri, ci = linear_sum_assignment(C)
fixed = n.copy()
for i, j in zip(ri, ci):
    k = dup.index[i]
    fixed.loc[k, ["feature_id", "beta_exp"]] = [cand.feature_id.iloc[j], cand.cont_beta.iloc[j]]
fixed["beta_exp_clip"] = fixed.beta_exp.clip(lo, hi)
fixed["off"] = (fixed.beta_drawn - fixed.beta_exp_clip).abs() > 0.008
print(tag, "| re-assigned", len(dup), "| unique features", fixed.feature_id.nunique(), "| checked", len(fixed), "| off after resolution", int(fixed.off.sum()),
      "| max|d| (testable)", round((fixed.beta_drawn - fixed.beta_exp_clip).abs().max(), 4), "| untestable drawn as", fixed[fixed.beta_exp.isna()].fill.unique().tolist())
fixed.to_csv(f"{out}/g4g_{tag}_resolved.tsv.gz", sep="\t", index=False)
