"""Independent recomputation of Figure 5D from the ATAC v3 COLOC-r2 replay output.

Upstream: Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/
          genetics/context/genetic_lineage_context_primary_pairs.tsv (primary signal pair per gene-study pair x lineage)
State per pair-lineage from per-source normalized posterior mass in peaks:
  both  : both sources >= 0.5
  one   : max source >= 0.5, other < 0.5
  partial: 0.05 <= max source < 0.5
  under5: max source < 0.05 (includes no overlap)
Gene-level: best state over the gene's study pairs within the trait class, per lineage.
"""
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = BASE / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute"
f = BASE / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-coloc-r2-20260924T053244Z/genetics/context/genetic_lineage_context_primary_pairs.tsv"
d = pd.read_csv(f, sep="\t")
print("rows", len(d), "primary flag", d.primary_signal_pair.value_counts().to_dict())
print("trait_class", d.trait_class.value_counts().to_dict())
print("lineages", d.lineage.value_counts().to_dict())
print("evidence_state", d.evidence_state.value_counts().to_dict())
d = d[d.primary_signal_pair.astype(str).str.upper() == "TRUE"].copy()
pairs = d[["gwas_name", "gene"]].drop_duplicates()
print("pairs", len(pairs), "genes", d.gene.nunique(), "ensembl", d.ensembl.nunique())
print("mapped mass min", d.mapped_posterior_mass.min())
for col in ["gse244832_accessible_mass", "gse281367_accessible_mass"]:
    d[col + "_n"] = d[col] / d.mapped_posterior_mass
a, b = d.gse244832_accessible_mass_n, d.gse281367_accessible_mass_n
mx, mn = np.maximum(a, b), np.minimum(a, b)
d["state"] = np.select([mn >= 0.5, mx >= 0.5, mx >= 0.05], ["both", "one", "partial"], "under5")
d.loc[d.mapped_posterior_mass < 0.95, "state"] = "untestable"
rank = {"both": 3, "one": 2, "partial": 1, "under5": 0, "untestable": -1}
d["rank"] = d.state.map(rank)
cls_map = {c: ("direct" if "direct" in c.lower() else "enzyme") for c in d.trait_class.unique()}
print("class map", cls_map)
d["cls"] = d.trait_class.map(cls_map)
pc = d.groupby(["cls", "lineage"]).apply(lambda x: x.drop_duplicates(["gwas_name", "gene"]).shape[0]).unstack()
print("pairs per class (per lineage)", pc.to_dict())
gl = d.groupby(["cls", "lineage", "ensembl"]).agg(gene=("gene", "first"), best=("rank", "max")).reset_index()
inv = {v: k for k, v in rank.items()}
gl["state"] = gl.best.map(inv)
cnt = gl.groupby(["cls", "lineage", "state"]).size().unstack(fill_value=0)
for s in ["both", "one", "partial", "under5"]:
    if s not in cnt:
        cnt[s] = 0
cnt = cnt[["under5", "partial", "one", "both"]]
cnt["n_genes"] = cnt.sum(axis=1)
print(cnt.to_string())
cnt.to_csv(OUT / "5d_recomputed_gene_level_counts.tsv", sep="\t")
gl.to_csv(OUT / "5d_recomputed_gene_level_states.tsv", sep="\t", index=False)
genes_by_cls = d.groupby("cls").ensembl.nunique().to_dict()
pairs_by_cls = d.drop_duplicates(["gwas_name", "gene"]).groupby("cls").size().to_dict()
both_cls = set(d[d.cls == "direct"].ensembl) & set(d[d.cls == "enzyme"].ensembl)
print("genes by class", genes_by_cls, "pairs by class", pairs_by_cls, "genes in both", len(both_cls), "union", d.ensembl.nunique())
# genes under 5% in every lineage
u = gl.groupby(["cls", "ensembl"]).best.max()
print("genes under 5% in every lineage:", (u == 0).groupby(level=0).sum().to_dict())

# ---- sensitivity: pair-level counts (no best-state collapse) and no 5% floor
d["state_nofloor"] = np.select([mn >= 0.5, mx >= 0.5, mx > 0], ["both", "one", "partial"], "none")
alt = d.groupby(["cls", "lineage"]).state_nofloor.value_counts().unstack(fill_value=0)
print("pair-level, no floor\n", alt.to_string())
gl2 = d.assign(r2=d.state_nofloor.map({"both": 3, "one": 2, "partial": 1, "none": 0})).groupby(["cls", "lineage", "ensembl"]).r2.max().reset_index()
alt2 = gl2.groupby(["cls", "lineage"]).r2.value_counts().unstack(fill_value=0)
print("gene-level, no floor (3 both, 2 one, 1 >0 partial, 0 none)\n", alt2.to_string())
# per-source thresholds swapped to 'joint' columns for comparison
ja, jb = d.joint_gse244832_accessible_mass / d.mapped_posterior_mass, d.joint_gse281367_accessible_mass / d.mapped_posterior_mass
print("joint-column max mass median", np.median(np.maximum(ja, jb)), "per-source max mass median", np.median(mx))
