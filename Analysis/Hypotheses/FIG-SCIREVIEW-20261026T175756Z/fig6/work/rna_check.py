import pandas as pd, numpy as np
W='/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/FIG-SCIREVIEW-20261026T175756Z/fig6/work/'
d=pd.read_csv('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv')
d['ens']=d.gene.str.split('.').str[0]
print('genes in fit',len(d),'unique ens',d.ens.nunique())
d['sup']=(d.padj<0.05)&(d.logFC.abs()>0.5)
print('supported (canonical)',d.sup.sum(),' q<0.05 & |lfc|<=0.5',((d.padj<0.05)&(d.logFC.abs()<=0.5)).sum())
m=pd.read_csv(W+'catalog_gene_calls_with_coloc.tsv',sep='\t')
j=m.merge(d[['ens','logFC','SE','padj','sup']],left_on='ensembl_id',right_on='ens',how='left')
print(pd.crosstab(j.rna_call,j.sup.fillna('absent')))
ti=j[j.rna_call=='indeterminate']; print('catalog RNA indeterminate with q<0.05',(ti.padj<0.05).sum())
print('max |lfc diff| cat vs fit', (j.lfc-j.logFC).abs().max(), ' max q diff',(j.q-j.padj).abs().max())
for s in ['GLP1R','HKDC1','MTARC1','THRB']:
    r=d[d.symbol==s]
    for _,x in r.iterrows():
        print(s,x.ens,round(x.logFC,4),'SE',round(x.SE,4),'CI',round(x.logFC-1.96*x.SE,3),round(x.logFC+1.96*x.SE,3),'padj',f'{x.padj:.3g}')
    if r.empty: print(s,'absent from 774 fit')
