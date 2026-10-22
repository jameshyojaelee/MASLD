import pandas as pd, numpy as np, sys
S=sys.argv[1]; R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
m=pd.read_csv(f"{R}/figures/main/fig2_genetics/panels/Fig2E_ancestry_unique_coloc_GWS_membership.tsv",sep="\t")
print("Fig2E membership rows",len(m)); print(m.display_category.value_counts().to_string())
eo=m[m.display_category=="EUR_only"].gene
c=pd.read_csv(f"{R}/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",usecols=["gwas_name","gene","PP.H4.abf","PP.H4.susie","top_snp","ancestry"])
print("coloc rows",len(c),"GWAS",c.gwas_name.nunique()); print(c.groupby("ancestry").gwas_name.nunique().to_string())
e=c[(c.ancestry=="EUR")&c.gene.isin(eo)&c["PP.H4.susie"].notna()]
best=e.sort_values("PP.H4.susie",ascending=False).drop_duplicates("gene")
print("EUR-only genes",len(eo),"with an EUR SuSiE row",best.gene.nunique())
vc=best.gwas_name.value_counts(); print(vc.to_string())
alt_ast={"UKBB_ALT","UKBB_AST","MVP_ALT_EUR","MVP_AST_EUR"}
print("best EUR GWAS is ALT/AST (fits the UKBB-discovery, MVP-estimation design):",best.gwas_name.isin(alt_ast).sum())
print("best EUR GWAS is UKBB_GGT (no MVP GGT stratum):",(best.gwas_name=="UKBB_GGT").sum())
best.to_csv(f"{S}/eur_only_best_eur_row.tsv",sep="\t",index=False)
# map to nearest clumped UKBB lead (same trait) within 500 kb
for T in ["ALT","AST"]:
    cl=pd.read_csv(f"{S}/ukbb_{T}_clumped500kb.tsv",sep="\t",dtype={"chromosome":str})
    sub=best[best.gwas_name.isin({f"UKBB_{T}",f"MVP_{T}_EUR"})].copy()
    sp=sub.top_snp.str.split(":",expand=True); sub["chr"]=sp[0]; sub["pos"]=pd.to_numeric(sp[1],errors="coerce")
    hit=[]
    for r in sub.itertuples():
        cc=cl[cl.chromosome==r.chr]
        d=(cc.position-r.pos).abs()
        hit.append(d.min() if len(d) else np.nan)
    sub["dist_nearest_lead"]=hit
    print(T,"genes",len(sub),"top SNP within 500 kb of a clumped non/palindromic UKBB lead:",(sub.dist_nearest_lead<=500000).sum(),"exact lead:",(sub.dist_nearest_lead==0).sum())
