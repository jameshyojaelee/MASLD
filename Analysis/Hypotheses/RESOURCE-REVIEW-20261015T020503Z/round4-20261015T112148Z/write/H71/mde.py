# H71 prespec (copied from plan-H57/mde.py, margins changed): minimum detectable odds ratio from margins only (unstratified approximation).
# Equal-tail exact critical values at two-sided alpha 0.05; power under Fisher noncentral hypergeometric.
import numpy as np
from scipy.stats import nchypergeom_fisher as NCH, hypergeom
from scipy.optimize import brentq
def mde(N,G,D,alpha=0.05,power=0.8):
    H=hypergeom(N,G,D); lo=max(0,G+D-N); hi=min(G,D); xs=np.arange(lo,hi+1)
    cdf=H.cdf(xs); sf=H.sf(xs-1)
    cL=xs[cdf<=alpha/2].max(); cU=xs[sf<=alpha/2].min()
    pw=lambda o: NCH(N,G,D,o).cdf(cL)+NCH(N,G,D,o).sf(cU-1)
    up=brentq(lambda o: pw(o)-power,1.0001,20); dn=brentq(lambda o: pw(o)-power,0.01,0.9999)
    e0=G*D/N
    return dict(N=N,G=G,D=D,E0=round(e0,3),cL=int(cL),cU=int(cU),size=round(pw(1.0),4),
                OR_up=round(up,3),E_up=round(NCH(N,G,D,up).mean(),2),OR_dn=round(dn,3),E_dn=round(NCH(N,G,D,dn).mean(),2))
for args in [(3769,426,241),(3769,426,228),(14704,426,902),(14704,426,943)]:
    print(mde(*args))
