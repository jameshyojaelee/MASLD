"""Recompute Fig 3F: 436 primary SuSiE genes (tier-1/2, PP.H4>0.5, COLOC r2) x four adjacent fibrosis contrasts
(774-unit stage_extension_all_gene_results.tsv); selection BH FDR<0.05 in >=1 contrast."""
import pandas as pd, numpy as np
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
S="figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv"
st=pd.read_csv(S,sep="\t")
print(st.groupby(['axis','contrast']).agg(n=('gene_id_base','size'),nref=('n_reference','first'),ncmp=('n_comparison','first'),fam=('bh_family_size','first'),sig=('FDR',lambda s:(s<0.05).sum())))
adj=st[st.axis=='fibrosis_adjacent']
t12=pd.read_csv("GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv")
prim=t12[t12.coloc_best_susie_pp4>0.5]
a=adj[adj.gene_id_base.isin(prim.ensembl)]
print("primary genes in stage table:",a.gene_id_base.nunique())
w=a.pivot_table(index=['gene_id_base','gene_name'],columns='contrast',values='logFC')
q=a.pivot_table(index=['gene_id_base','gene_name'],columns='contrast',values='FDR')
sel=q.index[(q<0.05).any(axis=1)]
print("selected genes:",len(sel), "marks per contrast:",(q.loc[sel]<0.05).sum().to_dict(),"total",(q.loc[sel]<0.05).values.sum())
qs=q.loc[sel]; ws=w.loc[sel]
only01=((qs['F0_to_F1']<0.05)&((qs.drop(columns='F0_to_F1')>=0.05).all(axis=1))).sum(); print("genes only via F0->F1:",only01)
split=(qs<0.0125); print("q<0.0125 genes:",split.any(axis=1).sum(),"marks",split.sum().to_dict(), sorted(qs.index[split.any(axis=1)].get_level_values(1)))
absmax=ws.abs().stack().sort_values(ascending=False); print("top |lfc|:",absmax.head(8).round(3).to_dict())
print("cells |lfc|>0.688:",(ws.abs()>0.688).values.sum())
pp=prim.set_index('ensembl').loc[sel.get_level_values(0)]
out=ws.join(qs,lsuffix='_lfc',rsuffix='_q'); out['pp4_best_susie']=pp.coloc_best_susie_pp4.values; out['susie_driving_gwas']=pp.susie_driving_gwas.values
out.to_csv(O+"fig3f_recomputed_matrix.tsv",sep="\t")
print("PP4 range",out.pp4_best_susie.min(),out.pp4_best_susie.max(), "below 0.8:",(out.pp4_best_susie<0.8).sum())
print("driving gwas:",out.susie_driving_gwas.value_counts().to_dict())
