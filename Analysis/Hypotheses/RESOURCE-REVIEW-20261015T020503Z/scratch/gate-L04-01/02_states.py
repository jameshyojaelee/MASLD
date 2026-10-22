import pandas as pd, numpy as np
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
d=pd.read_pickle(S+'/four.pkl'); d=d[d.ensembl.notna()].copy()
def state(r):
    h4=r['PP.H4.susie']; h3=r['PP.H3.susie']
    if pd.isna(h4): return 'no_susie'
    if h4>0.5: return 'H4'
    if h3>0.8: return 'H3'
    return 'other'
d['state']=d.apply(state,axis=1)
print(pd.crosstab(d.gwas_name,d.state))

W=d.pivot(index='ensembl',columns='gwas_name',values='state')
P4=d.pivot(index='ensembl',columns='gwas_name',values='PP.H4.susie')
rows=[]
for a,b in [('UKBB_ALT','MVP_ALT_EUR'),('MVP_ALT_EUR','UKBB_ALT'),('UKBB_AST','MVP_AST_EUR'),('MVP_AST_EUR','UKBB_AST')]:
    for s in ['H4','H3']:
        sub=W[W[a]==s]
        vc=sub[b].value_counts().to_dict()
        print(f'{a} {s} (n={len(sub)}) -> {b}: {vc}')
    pos=W[W[a]=='H4']; tested=pos[pos[b]!='no_susie']
    rep=(tested[b]=='H4').sum()
    print(f'  directed {a}->{b}: discovery H4 {len(pos)}, tested in B {len(tested)}, replicated {rep}, rate {rep/len(tested):.3f}')
    x=P4.loc[tested.index,a]
    print('  discovery PP.H4 quantiles among tested', x.quantile([0,.1,.25,.5,.75,.9,1]).round(3).to_dict())
    print('  bins: 0.5-0.8', ((x>0.5)&(x<0.8)).sum(), ' 0.8-0.95', ((x>=0.8)&(x<0.95)).sum(), ' >=0.95', (x>=0.95).sum())
    for lo,hi in [(0.5,0.8),(0.8,0.95),(0.95,1.01)]:
        m=(x>lo)&(x<hi) if lo==0.5 else (x>=lo)&(x<hi)
        idx=x[m].index
        print(f'   bin {lo}-{hi}: n={m.sum()} replicated={(W.loc[idx,b]=="H4").sum()}')
