import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
x=pd.read_pickle("t12.pkl")
b=pd.read_pickle("b_m.pkl")
s=x[x["PP.H4.susie"].notna()]
top=s.sort_values(["PP.H4.susie","gwas_name"],ascending=[False,True]).groupby("ensembl").head(1)
pos=top[top["PP.H4.susie"]>0.5]
print("436?",len(pos), "symbols unique",pos.gene.nunique())
inB=pos.gene.isin(b.gene)
print("multi-signal genes drawn in 2B:",inB.sum()," absent:",(~inB).sum())
ab=pos[~inB]
mvp=ab.gwas_name.isin(["MVP_ALT_EUR","MVP_AST_EUR","MVP_NAFLD_EUR"])
print("absent with max in MVP EUR:",mvp.sum()," other:",(~mvp).sum(), ab[~mvp].gwas_name.value_counts().to_dict())
# check 2B y vs our max for drawn
mm=b.merge(top[["gene","PP.H4.susie","gwas_name","top_snp"]],on="gene",how="left")
print("2B y == max SuSiE PP.H4 over 35:", (abs(mm.susie_pp4-mm["PP.H4.susie"])<1e-9).mean(), " gwas same:",(mm.gwas==mm.gwas_name).mean(), " lead same:",(mm.lead_variant==mm.top_snp).mean())
print(mm[(mm.gwas!=mm.gwas_name)][["gene","gwas","gwas_name","susie_pp4","PP.H4.susie"]].head(10).to_string())
print("2B points with PP.H4>0.5:",(b.susie_pp4>0.5).sum())
# For the 59 others, is the lead in ct?
ct=pd.read_pickle("ct.pkl")
ab2=ab[~mvp].copy()
ab2["chr"]=ab2.top_snp.str.split(":").str[0]; ab2["pos"]=ab2.top_snp.str.split(":").str[1]
ab2=ab2[ab2.chr.str.isnumeric()]
ab2["chr"]=ab2.chr.astype(int); ab2["pos"]=ab2.pos.astype(int)
k=set(zip(ct.study,ct.chromosome,ct.position))
print("others: lead present in coloc-targeted run for same study:",sum((r.gwas_name,r.chr,r.pos) in k for r in ab2.itertuples()), "of",len(ab2))
hit=[r for r in ab2.itertuples() if (r.gwas_name,r.chr,r.pos) in k]
for r in hit:
    print(r.gene,r.gwas_name,r.top_snp, round(r._3 if False else 0,1), ct[(ct.study==r.gwas_name)&(ct.chromosome==r.chr)&(ct.position==r.pos)][["locus","susie_pip","susie_cs","recommended_pip"]].to_dict("records"))
xr=pd.read_csv(R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F2abce/out/fig2b_same_row_x_by_rule.tsv",sep="\t")
print(xr.columns.tolist())
print(xr[xr.gene.isin([r.gene for r in hit])].to_string())
