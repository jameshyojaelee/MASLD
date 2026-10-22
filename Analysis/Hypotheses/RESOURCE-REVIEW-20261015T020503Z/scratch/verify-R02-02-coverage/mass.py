import pandas as pd
D="GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
m=pd.read_csv(D+"credible_set_members_annotated.tsv",sep="\t",dtype=str,keep_default_na=False)
g=m[m.pip_sum_gate_passed=="true"].copy()
g["w"]=g.normalized_susie_pip.astype(float)
print(g.consequence_category.value_counts().to_dict())
nc=g[~g.consequence_category.isin(["protein_altering","canonical_splice"])].copy()
nc["unres"]=nc.promoter_context_status=="unresolved_liftover"
for sc,d in nc.groupby("trait_scope"):
    print(sc,"noncoding mass",round(d.w.sum(),4),"unresolved mass",round(d[d.unres].w.sum(),4),"instances",d.credible_set_uid.nunique())
# instance counts by ancestry, per scope
inst=g.drop_duplicates("credible_set_uid")
print(pd.crosstab(inst.trait_scope,inst.ancestry,margins=True))
# unresolved share of noncoding mass by ancestry x scope
t=nc.groupby(["trait_scope","ancestry"]).apply(lambda d: pd.Series({"nc_mass":d.w.sum(),"unres_frac":d[d.unres].w.sum()/d.w.sum(),"n_cs":d.credible_set_uid.nunique()}))
print(t.round(4))
# site-level: non-EUR rows at sites resolved in a EUR study
g["site"]=g.chromosome+":"+g.position
eur_res=set(g[(g.ancestry=="EUR")&(g.promoter_context_status=="resolved")].site)
ne=g[(g.ancestry!="EUR")&(g.site.isin(eur_res))]
print("nonEUR rows at EUR-resolved sites",len(ne),"unresolved among them",(ne.promoter_context_status=="unresolved_liftover").sum())
ne2=g[(g.ancestry!="EUR")&(g.site.isin(eur_res))].drop_duplicates(["study","site"])
print("unique study-site", len(ne2), (ne2.promoter_context_status=="unresolved_liftover").sum())
# enzyme non-EUR share
enz=g[g.trait_scope=="tier2_liver_enzyme"]
print("enzyme rows nonEUR",(enz.ancestry!="EUR").mean(),"enzyme mass nonEUR",enz[enz.ancestry!="EUR"].w.sum()/enz.w.sum(), "enzyme instances nonEUR", (enz.drop_duplicates('credible_set_uid').ancestry!="EUR").mean())
