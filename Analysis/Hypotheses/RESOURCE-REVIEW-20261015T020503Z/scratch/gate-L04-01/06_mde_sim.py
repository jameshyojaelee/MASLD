import pandas as pd, numpy as np, warnings
import statsmodels.api as sm
from scipy.stats import norm
warnings.filterwarnings('ignore')
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L04-01'
O=pd.read_csv(S+'/pooled_obs.tsv',sep='\t')
rng=np.random.default_rng(20261015)
zcrit=norm.ppf(1-0.025)  # one-sided, Holm smallest of 2 at alpha 0.05
res=[]
for t,o in O.groupby('trait'):
    x=o.x.values; xc=x-x.mean(); loc=pd.factorize(o.locus)[0]; nl=loc.max()+1
    p0=o.y.mean(); a0=np.log(p0/(1-p0))
    X=sm.add_constant(xc)
    for sig in [1.0,2.0]:
        for beta in [0,2,4,6,8]:
            rej=0; nsim=400
            for s in range(nsim):
                u=rng.normal(0,sig,nl)[loc]
                pr=1/(1+np.exp(-(a0+beta*xc+u)))
                y=rng.binomial(1,pr)
                try:
                    m=sm.GLM(y,X,family=sm.families.Binomial()).fit(cov_type='cluster',cov_kwds={'groups':loc})
                    if m.params[1]/m.bse[1]>zcrit: rej+=1
                except Exception: pass
            res.append(dict(trait=t,sigma_locus=sig,beta_per_unit=beta,power=rej/nsim))
            print(t,sig,beta,rej/nsim,flush=True)
R=pd.DataFrame(res); R.to_csv(S+'/mde_sim.tsv',sep='\t',index=False)
