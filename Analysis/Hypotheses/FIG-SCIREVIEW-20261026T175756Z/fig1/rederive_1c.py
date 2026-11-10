"""Independent re-derivation of Figure 1C trait-class cross-tab from COLOC r2 per-study table."""
import json
import pandas as pd, numpy as np
A = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/"
R = A + "results/susie_coloc_r2_20260924T053244Z/"
tier = pd.read_csv(A + "config/gwas_trait_tier.tsv", sep="\t")
main = tier[(tier.placement == "main") & tier.tier.isin([1, 2])]
cls = dict(zip(main.study_name, np.where(main.tier == 1, "direct", "enzyme")))
print("studies per class", pd.Series(cls).value_counts().to_dict())
a = pd.read_csv(R + "susie_coloc_all_gwas.csv", usecols=["gwas_name", "ensembl", "PP.H4.susie", "method", "n_cs_pairs_eligible", "eligibility_source"], low_memory=False)
print("method values", a.method.value_counts(dropna=False).head(8).to_dict())
a = a[a.gwas_name.isin(cls)].copy()
print("main studies present", a.gwas_name.nunique())
a["cls"] = a.gwas_name.map(cls)
a["tested"] = a["PP.H4.susie"].notna()
a["pos"] = a["PP.H4.susie"] > 0.5
gl = pd.read_csv(R + "gene_level_coloc_tier12.csv")
pos436 = set(gl.loc[gl.susie_state_t12 == "positive", "ensembl"])
g = a.groupby(["ensembl", "cls"]).agg(tested=("tested", "any"), pos=("pos", "any")).unstack("cls", fill_value=False)
def state(c):
    return np.where(g[("pos", c)], "positive", np.where(g[("tested", c)], "tested_le", "none"))
st = pd.DataFrame({"enzyme": state("enzyme"), "direct": state("direct")}, index=g.index)
anypos = set(st.index[(st.enzyme == "positive") | (st.direct == "positive")])
out = {"n_anypos_recomputed": len(anypos), "n_pos436": len(pos436), "anypos_minus_436": sorted(anypos - pos436)[:20], "436_minus_anypos": sorted(pos436 - anypos)[:20]}
s4 = st.loc[sorted(pos436 & set(st.index))]
out["crosstab_on_436"] = pd.crosstab(s4.enzyme, s4.direct).to_dict()
out["enzyme_pos"] = int((s4.enzyme == "positive").sum()); out["direct_pos"] = int((s4.direct == "positive").sum())
out["both_pos"] = int(((s4.enzyme == "positive") & (s4.direct == "positive")).sum())
out["tested_enzyme_any"] = int(g[("tested", "enzyme")].sum()); out["tested_direct_any"] = int(g[("tested", "direct")].sum())
both = st[(st.enzyme != "none") & (st.direct != "none")]
out["tested_both"] = len(both); out["tested_both_enzyme_pos"] = int((both.enzyme == "positive").sum()); out["tested_both_direct_pos"] = int((both.direct == "positive").sum())
out["tested_both_both_pos"] = int(((both.enzyme == "positive") & (both.direct == "positive")).sum())
# genes with a positive PP.H4 in a study where the gene-level state is untestable
out["gene_level_untestable_ids"] = sorted(gl.loc[gl.susie_state_t12 == "untestable", "ensembl"])
print(json.dumps(out, indent=1, default=str))
json.dump(out, open("rederive_1c.json", "w"), indent=1, default=str)
