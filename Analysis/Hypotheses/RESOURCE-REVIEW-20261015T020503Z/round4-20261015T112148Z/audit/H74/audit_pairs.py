# H74 audit, part 2: coverage of the r2 SuSiE signal-pair file (hit1 = GWAS-side CS lead) for
# tier-1/2 positive rows and for locus anchors. Structure only: no gene distance, no nearest gene.
import sys, numpy as np, pandas as pd
np.random.seed(20260923)
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D=R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); t12=set(tier.study_name[tier.tier.isin([1,2])])
d=pd.read_csv(D+'susie_coloc_all_gwas.csv',usecols=['gwas_name','ensembl','PP.H4.susie','PP.H4.abf','top_snp','top_snp_PP'])
s=d[d.gwas_name.isin(t12)&d['PP.H4.susie'].notna()].copy()
s['chr']=s.top_snp.str.split(':').str[0]; s['pos']=s.top_snp.str.split(':').str[1].astype(int)
pos=s[s['PP.H4.susie']>0.5].copy()
sp=pd.read_csv(D+'susie_coloc_signal_pairs.tsv.gz',sep='\t')
sp=sp[sp.eligible.astype(str).str.upper()=='TRUE']
best=sp.sort_values('PP.H4',ascending=False).drop_duplicates(['ensembl','gwas'])
m=pos.merge(best[['ensembl','gwas','hit1','hit2','PP.H4']],left_on=['ensembl','gwas_name'],right_on=['ensembl','gwas'],how='left')
print('tier12 positive rows',len(pos),'with an eligible signal pair',m.hit1.notna().sum())
mm=m[m.hit1.notna()]
print('best pair PP.H4 equals PP.H4.susie (abs diff<1e-9):',(abs(mm['PP.H4']-mm['PP.H4.susie'])<1e-9).sum(),'of',len(mm))
print('hit1 == top_snp',(mm.hit1==mm.top_snp).sum(),'hit2 == top_snp',(mm.hit2==mm.top_snp).sum())
# loci and anchors (rule in PRESPEC_H74)
locs=[]
for c,g in pos.groupby('chr'):
    ps=sorted(g.pos.unique()); cur=[ps[0]]
    for p in ps[1:]:
        if p-cur[-1]<=1_000_000: cur.append(p)
        else: locs.append((c,cur[0],cur[-1])); cur=[p]
    locs.append((c,cur[0],cur[-1]))
cov=0; abf_low=0
for c,a,b in locs:
    w=m[(m.chr==c)&(m.pos>=a)&(m.pos<=b)].sort_values(['PP.H4.susie','top_snp_PP','gwas_name','ensembl'],ascending=[False,False,True,True])
    r=w.iloc[0]; cov+=pd.notna(r.hit1); abf_low+=r['PP.H4.abf']<=0.5
print('loci',len(locs),'anchor row has eligible signal pair',cov,'anchor row PP.H4.abf<=0.5',abf_low)
