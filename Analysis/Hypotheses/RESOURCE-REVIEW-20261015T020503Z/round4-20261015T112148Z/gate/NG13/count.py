# NG13 gate: counts structure only (n loci, tested genes per locus, null nearest rate, MDE, SORT1 control). Does not compute the observed nearest-gene contradiction rate.
import pandas as pd, numpy as np, re, json, sys
from scipy.stats import norm, binom
np.random.seed(20260923)
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
D=R+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/NG13/'
tier=pd.read_csv(R+'GWAS/finemapping/config/gwas_trait_tier.tsv',sep='\t')
t12=set(tier.study_name[tier.tier.isin([1,2])])
cols=['gwas_name','gene','ensembl','chr','PP.H3.susie','PP.H4.susie','method','top_snp','top_snp_PP']
d=pd.read_csv(R+'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',usecols=cols)
print('rows',len(d),'methods',d.method.value_counts().to_dict())
s=d[d.gwas_name.isin(t12)&d['PP.H4.susie'].notna()].copy()
print('tier12 susie rows',len(s),'studies',s.gwas_name.nunique(),'genes',s.ensembl.nunique())
s['pos']=s.top_snp.str.split(':').str[1].astype(int); s['chr']=s.top_snp.str.split(':').str[0]
# gene TSS hg19 protein-coding
tss={}
for l in open(R+'data/gsmap_resource/genome_annotation/gtf/gencode.v46lift37.basic.annotation.gtf'):
    if l[0]=='#': continue
    f=l.split('\t')
    if f[2]!='gene': continue
    gid=re.search('gene_id "([^".]+)',f[8]).group(1); gt=re.search('gene_type "([^"]+)',f[8]).group(1); gn=re.search('gene_name "([^"]+)',f[8]).group(1)
    t=int(f[3]) if f[6]=='+' else int(f[4])
    tss[gid]=(f[0].replace('chr',''),t,gt,gn)
s['tss']=s.ensembl.map(lambda g: tss.get(g,(None,None,None,None))[1])
s['gtype']=s.ensembl.map(lambda g: tss.get(g,(None,None,None,None))[2])
print('susie rows with TSS',s.tss.notna().mean(), 'protein_coding frac',(s.gtype=='protein_coding').mean())
pos=s[s['PP.H4.susie']>0.5]
# 1-Mb clustering of positive top SNPs (pooled across tier12 studies)
locs=[]
for c,g in pos.groupby('chr'):
    ps=sorted(g.pos.unique()); cur=[ps[0]]
    for p in ps[1:]:
        if p-cur[-1]<=1_000_000: cur.append(p)
        else: locs.append((c,cur[0],cur[-1])); cur=[p]
    locs.append((c,cur[0],cur[-1]))
print('tier12 positive loci (1Mb clusters):',len(locs))
rows=[]
for i,(c,a,b) in enumerate(locs):
    w=s[(s.chr==c)&(s.pos>=a-1_000_000)&(s.pos<=b+1_000_000)]
    gs=w.groupby('ensembl').agg(h4=('PP.H4.susie','max'),h3=('PP.H3.susie','max'),gtype=('gtype','first'),gene=('gene','first'))
    rows.append(dict(locus=i,chr=c,start=a,end=b,n_tested=len(gs),n_tested_pc=(gs.gtype=='protein_coding').sum(),n_h4=(gs.h4>0.5).sum(),n_h3=((gs.h3>0.8)&(gs.h4<=0.5)).sum()))
L=pd.DataFrame(rows); L.to_csv(D+'locus_structure.tsv',sep='\t',index=False)
print(L[['n_tested','n_tested_pc','n_h4','n_h3']].describe().round(2).to_string())
# null nearest rate: P(a colocalized gene is the nearest among tested genes) under random labels = n_h4/n_tested per locus
p0=(L.n_h4/L.n_tested).clip(upper=1)
print('loci with n_tested==n_h4 (uninformative for permutation):',(L.n_tested==L.n_h4).sum())
inf=L[L.n_tested>L.n_h4]; p0i=(inf.n_h4/inf.n_tested)
print('informative loci',len(inf),'mean null P(nearest tested gene is H4)',round(p0i.mean(),3))
# MDE: one-sided alpha 0.05, sum of Bernoullis; normal approx
n=len(inf); mu=p0i.sum(); var=(p0i*(1-p0i)).sum(); se=np.sqrt(var)/n
mde=(norm.ppf(0.95)+norm.ppf(0.8))*se
print('perm test MDE (excess proportion, one-sided 0.05, 80%):',round(mde,3))
# CI half-width for estimand (1) at p=0.5 over all loci
print('95% CI half-width at p=0.5 over',len(L),'loci:',round(1.96*np.sqrt(0.25/len(L)),3))
# SORT1 control
sub=s[(s.chr=='1')&(s.pos.between(109_700_000,109_950_000))]
print(sub.groupby(['gene']).agg(h4=('PP.H4.susie','max'),h3=('PP.H3.susie','max'),n=('gwas_name','size')).to_string())
snp=109817590
for gid,(c,t,gt,gn) in tss.items():
    if c=='1' and gt=='protein_coding' and abs(t-snp)<150000: print('pc gene near rs12740374',gn,t,abs(t-snp))
