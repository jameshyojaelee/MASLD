# H72 audit: power on the 281 applicable testable gene SNPs (PP.H4 >= 0.5), from per-sex SE and both_sexes beta only.
import math, numpy as np, pandas as pd
from scipy.stats import norm
A='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H72'
SEED=20260923
d=pd.read_csv(A+'/audit_power_rows.tsv',sep='\t')
g=d[(d.kind=='gene')&(d.pp>=0.5)]
out=open(A+'/audit_power_281.out','w')
def p(*a): print(*a,file=out,flush=True)
p('applicable testable genes',len(g),'PP>=0.8',int((g.pp>=0.8).sum()))
u=g.drop_duplicates(['trait','snp']); p('unique trait x SNP',len(u),u.trait.value_counts().to_dict(),'PP>=0.8 unique',len(g[g.pp>=0.8].drop_duplicates(['trait','snp'])))
u=u.assign(c=u.snp.str.split(':').str[0],pos=u.snp.str.split(':').str[1].astype(int))
n=0; per={}
for (t,c),x in u.groupby(['trait','c']):
    ps=sorted(x.pos); k=1+sum(1 for a,b in zip(ps,ps[1:]) if b-a>500000); n+=k; per[t]=per.get(t,0)+k
p('500 kb single-linkage clusters of the testable SNPs',n,per)
def power(x,label,nsim=2000):
    m=len(x); sD=np.sqrt(x.sF**2+x.sM**2); zc=norm.isf(0.05/m/2); mde=(zc+norm.isf(0.2))*sD.median()
    def pe(r):
        b=abs(r.bB); se=math.sqrt((r.sF/b)**2+(r.sM/b)**2)
        if b/r.sM<1.959964: return 0.0
        h=math.log(1.5)-1.959964*se
        return max(0.0,2*norm.cdf(h/se)-1) if h>0 else 0.0
    pev=x.apply(pe,axis=1)
    p(f'{label}: n {m}; median SE(diff) {sD.median():.5f}; BH single-hit z {zc:.3f}; MDE80 {mde:.4f}; median |b_both| {x.bB.abs().median():.4f}; n |b_both|>=MDE80 {int((x.bB.abs()>=(zc+norm.isf(0.2))*sD).sum())}; prespec expected bounded {pev.sum():.1f}; n P>=0.8 {(pev>=0.8).sum()}')
    rng=np.random.default_rng(SEED); b=x.bB.abs().values; sF=x.sF.values; sM=x.sM.values; z=1.959964
    nf=[]; nb=[]; nboth=[]
    for _ in range(nsim):
        bF=b+rng.standard_normal(m)*sF; bM=b+rng.standard_normal(m)*sM
        za=(bF-(2/3)*bM)/np.sqrt(sF**2+(4/9)*sM**2); zb=(1.5*bM-bF)/np.sqrt(sF**2+2.25*sM**2)
        fiel=(np.abs(bM/sM)>z)&(za>z)&(zb>z); peq=np.maximum(norm.sf(za),norm.sf(zb))
        o=np.argsort(peq); q=np.minimum.accumulate((peq[o]*m/np.arange(1,m+1))[::-1])[::-1]; qq=np.empty(m); qq[o]=np.minimum(q,1)
        bh=qq<0.05; nf.append(fiel.sum()); nb.append(bh.sum()); nboth.append((fiel&bh).sum())
    p(f'   simulated equal effects, one scale ({nsim} draws, seed {SEED}): Fieller-95% bounded mean {np.mean(nf):.1f}; BH(TOST) q<0.05 mean {np.mean(nb):.1f}; both mean {np.mean(nboth):.1f} (5th-95th pct {np.percentile(nboth,5):.0f}-{np.percentile(nboth,95):.0f})')
power(u,'gene arm 281')
power(u[u.pp>=0.8],'gene arm PP>=0.8 subset (one SNP per trait x SNP; pp of first gene)')
out.close()
