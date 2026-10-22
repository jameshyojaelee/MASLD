# Independent re-derivation of E1, E2, untested count, T_A and T_B from raw inputs with separate code
# (check only; written after the final run; does not feed the result).
import re, numpy as np, pandas as pd
from scipy.stats import rankdata
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
F=R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t'); t12=set(tier.study_name[tier.tier<=2])
d=pd.read_csv(F+'susie_coloc_all_gwas.csv',usecols=['gwas_name','ensembl','PP.H3.susie','PP.H4.susie','top_snp','top_snp_PP'])
s=d[d.gwas_name.isin(t12)&d['PP.H4.susie'].notna()].copy()
s['c']=s.top_snp.str.split(':').str[0]; s['p']=s.top_snp.str.split(':').str[1].astype(int)
pc={}; tssg={}
for l in open(R+'data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'):
    if l.startswith('#'): continue
    f=l.split('\t')
    if f[2]!='gene' or '_PAR_Y' in f[8]: continue
    g=re.search(r'gene_id "(ENSG\d+)',f[8]).group(1); c=f[0][3:]; t=int(f[3]) if f[6]=='+' else int(f[4])
    tssg[(g,c)]=t
    if 'gene_type "protein_coding"' in f[8]: pc.setdefault(c,[]).append((g,t))
sp=pd.read_csv(F+'susie_coloc_signal_pairs.tsv.gz',sep='\t'); sp=sp[sp.eligible==True]
sp=sp.sort_values(['ensembl','gwas','PP.H4','idx1','idx2'],ascending=[True,True,False,True,True]).groupby(['ensembl','gwas']).head(1)
h1={(e,g):h for e,g,h in zip(sp.ensembl,sp.gwas,sp.hit1)}
pos=s[s['PP.H4.susie']>0.5]
n=k=kb=ku=0; TA=0; mu=0; TB=[]; ninf=0
for c in sorted(pos.c.unique()):
    ps=np.sort(pos[pos.c==c].p.unique()); brk=np.flatnonzero(np.diff(ps)>1_000_000)
    for a,b in zip(np.r_[ps[0],ps[brk+1]], np.r_[ps[brk],ps[-1]]):
        n+=1
        w=s[(s.c==c)&s.p.between(a-1_000_000,b+1_000_000)]
        gs=w.groupby('ensembl')[['PP.H4.susie','PP.H3.susie']].max()
        P=pos[(pos.c==c)&pos.p.between(a,b)].sort_values(['PP.H4.susie','top_snp_PP','gwas_name','ensembl'],ascending=[False,False,True,True]).iloc[0]
        x=int(h1[(P.ensembl,P.gwas_name)].split(':')[1])
        arr=pc[c]; dd=np.array([abs(t-x) for _,t in arr]); near=[arr[i][0] for i in np.flatnonzero(dd==dd.min())]
        inC=[g for g in near if g in gs.index and gs.at[g,'PP.H4.susie']>0.5]
        if not inC:
            k+=1
            if any(g in gs.index and gs.at[g,'PP.H3.susie']>0.8 for g in near): kb+=1
            elif all(g not in gs.index for g in near): ku+=1
        m=[(g,tssg[(g,c)],gs.at[g,'PP.H4.susie']>0.5) for g in gs.index if (g,c) in tssg]
        lab=np.array([z for _,_,z in m]); 
        if lab.sum()==0 or lab.all(): continue
        ninf+=1; dist=np.array([abs(t-x) for _,t,_ in m]); nn=dist==dist.min()
        TA+=lab[nn].mean(); mu+=lab.mean(); r=rankdata(dist); TB.append(((r[lab]-0.5)/len(m)).mean())
print('loci',n,'E1',k,'bystander',kb,'untested(all tied untested)',ku,'informative',ninf,'T_A',TA,'null mean',round(mu,3),'T_B',round(np.mean(TB),6))
