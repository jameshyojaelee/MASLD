import pandas as pd, numpy as np
from scipy import stats
B="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
prot=pd.read_csv(B+"Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv",sep="\t")[["gene","logFC","padj"]].rename(columns={"logFC":"p_lfc","padj":"p_padj"})
atlas=pd.read_csv(B+"RNA-seq/results/multi_evidence/multi_evidence_atlas.csv",usecols=["human_symbol","ensembl_id"])
deg=pd.read_csv(B+"RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv",usecols=["gene","symbol","logFC","padj"])
deg["ensembl_id"]=deg.gene.str.replace(r"\..*$","",regex=True)
m=atlas.rename(columns={"human_symbol":"g"}).merge(deg[["ensembl_id","symbol"]],on="ensembl_id")
m=m[m.g==m.symbol][["g","ensembl_id"]]
up=deg[deg.symbol.notna()&(deg.symbol!="")].symbol.value_counts(); up=set(up[up==1].index)
rep=deg[deg.symbol.isin(up)&~deg.symbol.isin(m.g)][["symbol","ensembl_id"]].rename(columns={"symbol":"g"})
m=pd.concat([m,rep]); assert not m.g.duplicated().any()
bulk=m.merge(deg[["ensembl_id","logFC","padj"]],on="ensembl_id")
c=prot.merge(bulk,left_on="gene",right_on="g")
c=c[np.isfinite(c.logFC)&np.isfinite(c.p_lfc)]
uni=set(x.strip() for x in open(B+"Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt") if x.strip())
tx=set(x.strip() for x in open(B+"Analysis/Spatial/results/universe_validation/universe_transcriptomic.txt") if x.strip())
ge=set(x.strip() for x in open(B+"Analysis/Spatial/results/universe_validation/universe_genetic.txt") if x.strip())
print("full measured",len(c),"rho_full",stats.spearmanr(c.logFC,c.p_lfc)[0])
u=c[c.gene.isin(uni)]
print("universe measured",len(u),"rho",stats.spearmanr(u.logFC,u.p_lfc)[0],"dir",np.mean(np.sign(u.logFC)==np.sign(u.p_lfc)))
t=u[u.gene.isin(tx)]; go=u[~u.gene.isin(tx)]; gg=u[u.gene.isin(ge)&~u.gene.isin(tx)]
print("tx arm",len(t),stats.spearmanr(t.logFC,t.p_lfc)[0],"not-tx",len(go),stats.spearmanr(go.logFC,go.p_lfc)[0],"genetic&not-tx",len(gg), "in neither",len(u[~u.gene.isin(tx)&~u.gene.isin(ge)]))
s=u[u.p_padj<0.05]; s=s.assign(agree=np.sign(s.logFC)==np.sign(s.p_lfc))
print("prot-sig",len(s),"agree",s.agree.sum(),"notsig",len(u)-len(s))
for lab,x in [("rna<.05",s[s.padj<0.05]),("rna>=.05",s[~(s.padj<0.05)])]:
    print(lab,len(x),x.agree.sum(),round(x.agree.mean()*100,2),stats.binomtest(int(x.agree.sum()),len(x),0.5).pvalue)
print("nan rna padj among sig:",s.padj.isna().sum())
