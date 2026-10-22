import pandas as pd, collections
D="GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19783321/"
m=pd.read_csv(D+"credible_set_members_annotated.tsv",sep="\t",dtype=str,keep_default_na=False)
L=pd.read_csv("GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv",dtype=str,keep_default_na=False)
print("liftover rows",len(L), "empty hg38:", (L.pos_hg38=="").sum())
studies=collections.Counter(s for x in L.study_list for s in x.split(";"))
print("studies in liftover table:",len(studies)); print(sorted(studies.items()))
L["vid"]=L.chromosome+":"+L.position+":"+L.allele1+":"+L.allele2
L["site"]=L.chromosome+":"+L.position
L["vid_sw"]=L.chromosome+":"+L.position+":"+L.allele2+":"+L.allele1
vids=set(L.vid); sites=set(L.site); sw=set(L.vid_sw)
g=m[m.pip_sum_gate_passed=="true"].copy()
g["site"]=g.chromosome+":"+g.position
u=g[g.promoter_context_status=="unresolved_liftover"]
print("unresolved rows",len(u))
u=u.assign(in_vid=u.variant_id.isin(vids), swapped=u.variant_id.isin(sw), site_in=u.site.isin(sites))
print(u.groupby("ancestry")[["in_vid","swapped","site_in"]].sum().assign(n=u.groupby("ancestry").size()))
print("unresolved studies absent from liftover study_list:")
print(u.groupby("study").size().to_frame("n").assign(in_list=lambda d: [s in studies for s in d.index]).to_string())
# Is MVP/BBJ allele order REF-first? check vs liftover table for swapped
