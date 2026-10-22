"""Recompute Fig 3E bar heights (F1-F4 vs F0, BH q<0.05 split by sign, no fold-change floor) and the
alternative under the canonical |log2FC|>0.5 floor; also positive-control genes."""
import pandas as pd, numpy as np
O="Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig3/out/"
st=pd.read_csv("figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv",sep="\t")
f=st[st.axis=='fibrosis']
rows=[]
for c,g in f.groupby('contrast'):
    s=g.FDR<0.05
    rows.append(dict(contrast=c,n_ref=g.n_reference.iloc[0],n_cmp=g.n_comparison.iloc[0],up=(s&(g.logFC>0)).sum(),down=(s&(g.logFC<0)).sum(),
       total=s.sum(),up_lfc05=(s&(g.logFC>0.5)).sum(),down_lfc05=(s&(g.logFC<-0.5)).sum()))
r=pd.DataFrame(rows); r['total_lfc05']=r.up_lfc05+r.down_lfc05; print(r.to_string()); r.to_csv(O+"fig3e_recomputed_counts.tsv",sep="\t",index=False)
pc=['COL1A1','COL1A2','COL3A1','LUM','EPCAM','AKR1B10','GPNMB','TREM2','SPP1','CYP2C19','GNMT','MAT1A','CYP7A1','FABP4','THBS2','PLIN2']
x=f[f.gene_name.isin(pc)].pivot_table(index='gene_name',columns='contrast',values='logFC').round(2)
qq=f[f.gene_name.isin(pc)].pivot_table(index='gene_name',columns='contrast',values='FDR')
print(x.join(qq.map(lambda v: f"{v:.1e}"),rsuffix='_q').to_string())
