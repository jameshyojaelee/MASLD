import pandas as pd, numpy as np
from scipy.stats import false_discovery_control, spearmanr
m=pd.read_csv('data/PXD051911/meta_data.txt',sep='\t',keep_default_na=False,na_values=['','NA'])
l=m[m.liver_proteomics_filename.notna()]
print('liver rows',len(l),'unique patients',l.patient_name.nunique())
u=l.drop_duplicates('patient_name')
print('saf_diagnosis',u.saf_diagnosis.value_counts(dropna=False).to_dict())
print('kleiner',u.kleiner_fibrosis_grade.value_counts(dropna=False).to_dict())
print('sample_group',l.sample_group.value_counts().to_dict())
print('bmi missing',u.bmi.isna().sum(),'age missing',u.alder.isna().sum())
p=pd.read_csv('figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/protein_de_adjusted.tsv',sep='\t',keep_default_na=False,na_values=[''])
print('proteins',len(p),'unique gene',p.gene.nunique(),'padj<0.05',int((p.padj<0.05).sum()))
q=false_discovery_control(p.pvalue.values,method='bh'); print('BH recompute maxdiff',np.nanmax(np.abs(q-p.padj.values)))
print('genes with ; or empty',int(p.gene.astype(str).str.contains(';').sum()),int(p.gene.isna().sum()))
deg=pd.read_csv('RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv',keep_default_na=False,na_values=[''])
print('deg cols',list(deg.columns)[:20],'rows',len(deg))
sym=deg[deg.symbol.notna()]; cnt=sym.groupby('symbol').size(); uniq=cnt[cnt==1].index
j=p.merge(sym[sym.symbol.isin(uniq)][['symbol','logFC']].rename(columns={'symbol':'gene','logFC':'rna'}),on='gene')
print('joined',len(j),'rho',round(spearmanr(j.logFC,j.rna)[0],4),'dir',round((np.sign(j.logFC)==np.sign(j.rna)).mean(),4))
