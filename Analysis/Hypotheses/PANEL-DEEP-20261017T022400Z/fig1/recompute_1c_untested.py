"""For the 1C classes, how many 'only' genes were SuSiE-tested in the other trait class."""
import pandas as pd, json
from pathlib import Path
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
REL = ROOT / "GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
OUT = ROOT / "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1"
tier = pd.read_csv(ROOT / "GWAS/finemapping/config/gwas_trait_tier.tsv", sep="\t")
main = tier[(tier.placement == "main") & tier.tier.isin([1, 2])]
a = pd.read_csv(REL / "susie_coloc_all_gwas.csv", usecols=["gwas_name", "gene", "ensembl", "PP.H4.susie", "method"], low_memory=False)
a = a[a.gwas_name.isin(main.study_name)].merge(main[["study_name", "tier_label", "trait"]], left_on="gwas_name", right_on="study_name")
cls = pd.read_csv(OUT / "recompute_1c_gene_classes.tsv", sep="\t", index_col=0)["class"]
res = {}
for c, other in (("enzyme_only", "direct_MASLD"), ("direct_only", "liver_enzyme")):
    genes = cls[cls == c].index
    sub = a[a.ensembl.isin(genes) & (a.tier_label == other)]
    tested = sub[sub["PP.H4.susie"].notna()].groupby("ensembl")["PP.H4.susie"].max()
    unt = sub[sub.method.str.startswith("susie_untestable")].ensembl.unique()
    res[c] = dict(n=len(genes), tested_in_other_class=int(tested.size),
                  untested_in_other_class=int(len(genes) - tested.size),
                  untestable_pairs_in_other=int(len(unt)),
                  max_pp4_other_among_tested=float(tested.max()) if tested.size else None,
                  n_other_pp4_025_to_05=int(((tested > 0.25) & (tested <= 0.5)).sum()) if tested.size else 0)
# by trait for enzyme-only 
genes = cls[cls == "enzyme_only"].index
sub = a[a.ensembl.isin(genes) & (a.tier_label == "direct_MASLD") & a["PP.H4.susie"].notna()]
res["enzyme_only_direct_tested_by_trait"] = sub.groupby("trait").ensembl.nunique().to_dict()
pos = a[a["PP.H4.susie"] > 0.5]
res["positive_genes_by_trait"] = pos.groupby("trait").ensembl.nunique().to_dict()
res["tested_genes_by_trait"] = a[a["PP.H4.susie"].notna()].groupby("trait").ensembl.nunique().to_dict()
# genes with direct positive: which traits (PDFF vs NAFLD vs NASH)
d = pos[pos.tier_label == "direct_MASLD"]
res["direct_positive_genes_trait_combo"] = d.groupby("ensembl").trait.agg(lambda s: "+".join(sorted(set(s)))).value_counts().to_dict()
res["direct_positive_gene_symbols"] = sorted(d.gene.unique().tolist())
json.dump(res, open(OUT / "recompute_1c_untested.json", "w"), indent=1)
print(json.dumps(res, indent=1))
