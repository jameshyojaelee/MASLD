# H71 prespec: margins only. Counts list sizes in each universe and unmatched symbols.
# Does not compute the overlap of any competitor list with the genetic set inside
# U_full or U_eval, any expectation, or any test statistic.
import pandas as pd, numpy as np
R="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
G4=R+"Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG10/"
C=R+"GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
B=R+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
P=R+"figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/"
X=R+"data/public_r4/biorxiv_2026.08.26.747338/media-1.xlsx"
E=R+"data/external/elison2025_sc_chromatin/elison2025_sc_chromatin_calls.csv"
gi=pd.read_parquet(P+"passport_gene_index.parquet")
print("catalog genes",len(gi),"symbol_collision True",int(gi.symbol_collision.fillna(False).sum()),"dup symbols",int(gi.symbol.duplicated(keep=False).sum()))
gi["g"]=gi.ensembl_id.str.split(".").str[0]
lists={}
for s in ["NAS","Fibrosis"]:
    k=pd.read_csv(G4+f"kumak_{s}.tsv",sep="\t")
    print(s,"tsv rows",len(k),"unique genes",k.Gene.nunique())
    k["dir"]=np.where(k.total_upregulated>k.total_downregulated,"up",np.where(k.total_upregulated<k.total_downregulated,"down","tie"))
    print(s,"dir",k.dir.value_counts().to_dict(),"GTEX",k.GTEX.value_counts().to_dict(),"pubmed_mentioned",k.pubmed_mentioned.value_counts().to_dict())
    m=k.merge(gi[["symbol","g"]],left_on="Gene",right_on="symbol",how="left")
    print(s,"unmatched",sorted(m.loc[m.g.isna(),"Gene"].tolist()))
    lists[s]=set(m.g.dropna())
print("NAS & Fibrosis matched overlap",len(lists["NAS"]&lists["Fibrosis"]),"union",len(lists["NAS"]|lists["Fibrosis"]))
c=pd.read_csv(C); c["g"]=c.ensembl.astype(str).str.split(".").str[0]
c=c.sort_values(["coloc_best_susie_pp4","coloc_best_abf_pp4"],ascending=False).drop_duplicates("g")
d=pd.read_csv(B); d["g"]=d.gene.astype(str).str.split(".").str[0]; d=d.drop_duplicates("g")
jt=d.merge(c[["g","coloc_best_susie_pp4","susie_state_t12"]],on="g",how="inner")
jt["susie"]=pd.to_numeric(jt.coloc_best_susie_pp4,errors="coerce")
ue=jt[jt.susie_state_t12.isin(["positive","tested_le_0.5"])]
print("U_full",len(jt),"G",int((jt.susie>=0.5).sum()),"| U_eval",len(ue),"G",int((ue.susie>=0.5).sum()))
print("catalog genes in bulk table",int(gi.g.isin(d.g).sum()),"in U_full",int(gi.g.isin(jt.g).sum()))
for s,L in lists.items():
    inb=d.g.isin(L).sum(); inf=jt.g.isin(L).sum(); inv=ue.g.isin(L).sum()
    st=jt[jt.g.isin(L)].susie_state_t12.value_counts().to_dict()
    print(s,"matched",len(L),"in bulk table",int(inb),"in U_full",int(inf),"in U_eval",int(inv),"SuSiE state in U_full",st)
    print(s,"AveExpr quantiles in U_eval",ue[ue.g.isin(L)].AveExpr.quantile([0,.1,.5,.9,1]).round(3).to_dict())
print("U_eval AveExpr quantiles",ue.AveExpr.quantile([0,.1,.5,.9,1]).round(3).to_dict())
el=pd.read_csv(E); print("Elison rows",len(el),"types",el.elison_finding_type.value_counts().to_dict())
