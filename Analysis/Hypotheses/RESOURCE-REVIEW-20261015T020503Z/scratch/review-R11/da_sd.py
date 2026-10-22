import pandas as pd, numpy as np
d=pd.read_csv('da/da_peak_results.tsv.gz',sep='\t',keep_default_na=False,na_values=[''])
sd=d[d.evidence_state=='source_dependent'].copy()
sd['which']=np.where(sd.qvalue_gse244832<0.05,'sig_in_244832','sig_in_281367')
sd['same_sign']=np.sign(sd.logFC_gse244832)==np.sign(sd.logFC_gse281367)
print(sd.groupby(['lineage','which']).agg(n=('same_sign','size'),same_sign=('same_sign','sum')))
print('total source_dependent same sign', sd.same_sign.sum(),'/',len(sd))
ind=d[d.evidence_state=='indeterminate']
print('indeterminate: min q per cohort <0.1', ((ind.qvalue_gse244832<0.1)|(ind.qvalue_gse281367<0.1)).sum())
print('indeterminate: |logFC|>1 in both same sign', ((ind.logFC_gse244832.abs()>1)&(ind.logFC_gse281367.abs()>1)&(np.sign(ind.logFC_gse244832)==np.sign(ind.logFC_gse281367))).sum())
sup=d[d.evidence_state=='supported']; print(sup[['lineage','peak_coordinate','logFC_gse244832','qvalue_gse244832','logFC_gse281367','qvalue_gse281367','promoter_genes']].to_string())
# expected overlap of sig sets under independence
for L in ['hepatocyte','stellate']:
    x=d[d.lineage==L]; a=(x.qvalue_gse244832<0.05).sum(); b=(x.qvalue_gse281367<0.05).sum(); n=len(x)
    print(L,'n',n,'sig244832',a,'sig281367',b,'expected both under independence',round(a*b/n,2),'observed both',((x.qvalue_gse244832<0.05)&(x.qvalue_gse281367<0.05)).sum())
