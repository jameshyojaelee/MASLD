import pandas as pd, numpy as np
from statsmodels.stats.power import TTestPower
E='Analysis/MASLD_Model_Benchmark/executions/'
d=pd.read_csv(E+'gse83452-paired-intervention-design-20260830/results_20260830T134500/participant_deltas.tsv',sep='\t')
print(d.columns.tolist())
arm=[c for c in d.columns if 'arm' in c.lower() or 'interv' in c.lower()][0]
dc=[c for c in d.columns if 'delta' in c.lower()][0]
for a,g in d.groupby(arm):
    x=g[dc]; print('prior130',a,len(x),round(x.mean(),3),round(x.std(),3),'d=',round(x.mean()/x.std(),3),'n_down',(x<0).sum())
r=pd.read_csv(E+'model-data-072-21080651-participants/gse83452_records.tsv',sep='\t')
p=r[r.paired_expression_stress_evaluable==True]
# DPT ENSG00000143196 (from Catalog r8 index); COL1A1 ENSG00000108821
want={'ENSG00000143196':'DPT','ENSG00000108821':'COL1A1'}
rows=[]
with open(E+'model-data-088-21112137-summarization/summary/gse83452_gpl16686_gene_matrix.tsv') as f:
    hdr=f.readline().rstrip('\n').split('\t')
    for line in f:
        k=line.split('\t',1)[0].split('.')[0]
        if k in want: rows.append(line.rstrip('\n').split('\t'))
m=pd.DataFrame([x[1:] for x in rows],index=[want[x[0].split('.')[0]] for x in rows],columns=hdr[1:]).astype(float)
for a,g in p.groupby('intervention'):
    dd=m[g.sample_accession].values-m[g.paired_baseline_accession].values
    for i,gn in enumerate(m.index):
        x=dd[i]; print('PC',gn,a,len(x),'mean',round(x.mean(),3),'sd',round(x.std(ddof=1),3),'d',round(x.mean()/x.std(ddof=1),3),'n_down',(x<0).sum())
tp=TTestPower()
for al in [0.05,0.025]:
  print('MDE n=28 alpha',al,round(tp.solve_power(nobs=28,alpha=al,power=0.8,alternative='two-sided'),3))
print('MDE n=28 one-sided 0.05',round(tp.solve_power(nobs=28,alpha=0.05,power=0.8,alternative='larger'),3))
print('pairs',p.groupby('intervention').size().to_dict(),'scan same year', (p.groupby('intervention').scan_date.apply(lambda s: s.str[-4:].value_counts().to_dict())).to_dict())
