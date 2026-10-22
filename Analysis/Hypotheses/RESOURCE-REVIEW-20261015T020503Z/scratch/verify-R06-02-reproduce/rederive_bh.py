import pandas as pd, numpy as np
reg=pd.read_csv('program_registry_v2.tsv',sep='\t')
print([c for c in reg.columns if 'select' in c.lower() or 'state' in c.lower() or 'pass' in c.lower()])
sub=reg[(reg.cell_type=='hepatocytes')&(reg.module.isin([8,20]))]
print(sub[['cell_type','module','module_name','n_source_genes','primary_n_donors','primary_beta','primary_pvalue','primary_qvalue','primary_hc3_qvalue']].to_string())
pe=pd.read_csv('primary_effects.tsv',sep='\t')
print('rows',len(pe),'estimable',pe.estimable.sum())
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p)
    q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; q=np.minimum(q,1); out=np.empty(n); out[o]=q; return out
m=pe.estimable==True
pe.loc[m,'q_re']=bh(pe.loc[m,'pvalue'])
pe.loc[m,'hc3q_re']=bh(pe.loc[m,'hc3_pvalue'])
print('family size',m.sum())
for mod in [8,20]:
    r=pe[(pe.cell_type=='hepatocytes')&(pe.module==mod)].iloc[0]
    print('module',mod,'beta',r.beta,'p',r.pvalue,'q_file',r.qvalue,'q_re',r.q_re,'hc3p',r.hc3_pvalue,'hc3q_file',r.hc3_qvalue,'hc3q_re',r.hc3q_re,'n',r.n_donors)
print('max abs diff q', (pe.q_re-pe.qvalue).abs().max(), 'hc3', (pe.hc3q_re-pe.hc3_qvalue).abs().max())
print('q<0.05:'); print(pe.loc[pe.qvalue<0.05,['cell_type','module','beta','qvalue']].to_string())
print('hc3 q<0.05 count', (pe.hc3_qvalue<0.05).sum())
