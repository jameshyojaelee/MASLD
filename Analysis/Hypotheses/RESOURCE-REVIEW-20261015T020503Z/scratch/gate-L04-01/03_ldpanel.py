import pandas as pd, numpy as np
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
f='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/results/ld_panel_comparison/ld_panel_concordance.csv'
L=pd.read_csv(f,low_memory=False)
print('MVP studies in LD file:', [g for g in L.gwas_name.unique() if 'MVP' in g])
L=L[L.gwas_name.isin(['UKBB_ALT','UKBB_AST'])&L.ensembl.notna()&(L.ensembl!='')]
print('dup gene x study', L.duplicated(['gwas_name','ensembl']).sum())
d=pd.read_pickle(S+'/four.pkl'); d=d[d.ensembl.notna()]
for st in ['UKBB_ALT','UKBB_AST']:
    l=L[L.gwas_name==st].set_index('ensembl'); r=d[d.gwas_name==st].set_index('ensembl')
    cols=['PP.H4.susie_1KG','PP.H4.susie_TOP-LD','PP.H4.susie_UKBB_v1']
    print(st, 'non-NA susie per panel', {c:int(l[c].notna().sum()) for c in cols}, ' H4>0.5 per panel', {c:int((l[c]>0.5).sum()) for c in cols})
    j=r[['PP.H4.susie']].join(l[cols],how='inner')
    both=j.dropna(subset=['PP.H4.susie','PP.H4.susie_UKBB_v1'])
    print('  r2 vs LD-file UKBB_v1: both tested', len(both), ' |diff|<1e-6', int((abs(both['PP.H4.susie']-both['PP.H4.susie_UKBB_v1'])<1e-6).sum()), ' spearman', round(both['PP.H4.susie'].corr(both['PP.H4.susie_UKBB_v1'],method='spearman'),3))
    print('  r2 H4 n', int((j['PP.H4.susie']>0.5).sum()), ' r2 H4 & UKBB_v1 H4', int(((j['PP.H4.susie']>0.5)&(j['PP.H4.susie_UKBB_v1']>0.5)).sum()))
    # LD-panel replicate: discovery = UKBB_v1 positives, replicate = 1KG / TOP-LD
    for rep in ['PP.H4.susie_1KG','PP.H4.susie_TOP-LD']:
        pos=l[l['PP.H4.susie_UKBB_v1']>0.5]; t=pos[pos[rep].notna()]
        print(f'  UKBB_v1 H4 {len(pos)} -> {rep}: tested {len(t)}, H4 {(t[rep]>0.5).sum()}, rate {(t[rep]>0.5).mean():.3f}')
        x=t['PP.H4.susie_UKBB_v1']
        for lo,hi in [(0.5,0.8),(0.8,0.95),(0.95,1.01)]:
            m=(x>lo)&(x<hi) if lo==0.5 else (x>=lo)&(x<hi)
            print(f'     bin {lo}-{hi}: n={m.sum()} rep={(t.loc[m,rep]>0.5).sum()}')
    # also r2-discovery -> 1KG
    pos=j[j['PP.H4.susie']>0.5]
    for rep in ['PP.H4.susie_1KG','PP.H4.susie_TOP-LD']:
        t=pos[pos[rep].notna()]
        print(f'  r2 H4 {len(pos)} -> {rep}: tested {len(t)}, H4 {(t[rep]>0.5).sum()}')
