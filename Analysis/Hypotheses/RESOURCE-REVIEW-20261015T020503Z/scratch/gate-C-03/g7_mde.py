import numpy as np
from scipy.stats import binom, binomtest, norm
def region(n,a=0.05):
    return [k for k in range(n+1) if binomtest(k,n,0.5).pvalue<=a]
def power(n,p,a=0.05):
    R=region(n,a); return sum(binom.pmf(k,n,p) for k in R)
def mde(n,a=0.05):
    lo,hi=0.5,1.0
    for _ in range(50):
        mid=(lo+hi)/2
        if power(n,mid,a)>=0.8: hi=mid
        else: lo=mid
    return hi
for n in [83,86,89,90,112,122,123]:
    R=[k for k in region(n) if k>n/2]
    print(f'n={n} reject upper at k>={min(R)} ({min(R)/n:.3f}); MDE80={mde(n):.3f}; power@0.60={power(n,0.60):.2f}; power@0.69={power(n,0.69):.2f}')
# Bounds on the locus statistic implied by the disclosed gene-level count 58/96 (outcome-blind derivation):
# 96 eligible genes = 83 single-gene loci + 6 two-gene loci (12 genes) + 1 MHC gene.
res=[]
for mhc_agree in [0,1]:
    tot=58-mhc_agree
    for k2 in range(7):
        for k1 in range(7-k2):
            k0=6-k2-k1; b=2*k2+k1; a=tot-b
            if 0<=a<=83:
                n=83+k2+k0; agree=a+k2
                res.append((agree/n,agree,n,binomtest(agree,n,0.5).pvalue))
res.sort()
print('locus agreement range',f'{res[0][1]}/{res[0][2]}={res[0][0]:.3f} p={res[0][3]:.3f}','to',f'{res[-1][1]}/{res[-1][2]}={res[-1][0]:.3f} p={res[-1][3]:.4f}')
print('share of configurations with exact p<0.05:',round(np.mean([r[3]<0.05 for r in res]),2),'of',len(res))
print('gene-level 58/96 exact two-sided p',round(binomtest(58,96,0.5).pvalue,4))
d=(norm.ppf(0.975)+norm.ppf(0.8))*np.sqrt(1/12+1/46); d3=(norm.ppf(1-0.0005)+norm.ppf(0.8))*np.sqrt(1/12+1/46)
print('protein per-protein MDE SD units: alpha .05',round(d,2),'alpha .001',round(d3,2))
