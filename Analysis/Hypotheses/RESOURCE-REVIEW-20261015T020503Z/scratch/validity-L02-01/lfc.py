import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
deg=pd.read_csv(B+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
print(list(deg.columns)[:20])
can=(deg.padj<0.05)&(deg.logFC.abs()>0.5)
print('canonical n',can.sum(),'logFC quantiles',np.round(deg.logFC[can].quantile([0,.05,.25,.5,.75,.95,1]).values,3))
print('|logFC| canonical quantiles',np.round(deg.logFC[can].abs().quantile([.05,.25,.5,.75,.95]).values,3))
nc=~can
print('non-call n',nc.sum(),'logFC sd',round(deg.logFC[nc].std(),3),'IQR',np.round(deg.logFC[nc].quantile([.25,.75]).values,3), 'canonical sd', round(deg.logFC[can].std(),3))
for g in ['SAA1','SAA2','CRP','LBP','ORM1','HP','FGB','SERPINA3','LCN2','MT1X','MT2A','HAMP','CYP7A1','HMOX1','IL6','S100A8','S100A9','FOS','JUN','CXCL8','SOCS3']:
    r=deg[deg.symbol==g]
    if len(r): print(g, round(float(r.logFC.iloc[0]),3), '%.2g'%float(r.padj.iloc[0]), 'canonical' if bool(can[r.index[0]]) else '')
