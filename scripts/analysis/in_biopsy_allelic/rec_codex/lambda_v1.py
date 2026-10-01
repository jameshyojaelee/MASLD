"""P2 reliability callback using fixed public-pair draws and the same F resample.

No allele counts enter this estimator. Error draws and Gate A information
weights must be supplied; the callback never substitutes outcome-based weights.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr


X,W=np.polynomial.hermite.hermgauss(60)
W=W/np.sqrt(np.pi)


def category_probabilities(z,cuts,sd):
    cumulative=ndtr((np.asarray(cuts)[None,:]-np.asarray(z)[:,None])/sd)
    return np.diff(np.column_stack((np.zeros(len(z)),cumulative,np.ones(len(z)))),axis=1)


def target_distribution(labels,collapsed=False):
    cuts=(1.5,2.5) if collapsed else (.5,1.5,2.5,3.5)
    labels=np.asarray(labels,int)
    if not len(labels) or not np.isin(labels,np.arange(len(cuts)+1)).all():
        raise ValueError('invalid target labels')
    counts=np.bincount(labels,minlength=len(cuts)+1)
    def loss(v):
        prob=category_probabilities(np.array([v[0]]),cuts,np.exp(v[1]))[0].clip(1e-12)
        return -counts@np.log(prob)
    fit=minimize(loss,[2.,0.],method='L-BFGS-B',bounds=[(-2,6),(np.log(.05),np.log(6))])
    if not fit.success:
        raise ValueError('target label fit failed: '+str(fit.message))
    return float(fit.x[0]),float(np.exp(fit.x[1]))


def marginal_reliability(mu,sd,var_e):
    if not np.isfinite(var_e) or var_e<=0:
        return float('nan')
    latent=mu+np.sqrt(2*max(0.,sd*sd-var_e))*X
    prob=category_probabilities(latent,(1.5,2.5),np.sqrt(var_e))
    conditional_mean=prob@np.arange(3)
    marginal=W@prob
    mean=W@conditional_mean
    between=W@(conditional_mean-mean)**2
    total=marginal@np.arange(3)**2-mean*mean
    return float(between/total) if total>0 else float('nan')


def conditional_reliability(marginal,r_squared):
    if not np.isfinite(r_squared) or r_squared>=1:
        return float('nan')
    return float((marginal-r_squared)/(1-r_squared))


class PairedReliability:
    """source_specs is ordered source names with a fixed carried Y5 or S3.

    error_draws[source][b] must correspond to public-pair bootstrap draw b.
    target_y5 maps original participant IDs to original Y5 labels; GSE213621
    is fitted from its sampled Participant.stage S3 labels instead.
    """
    def __init__(self,source_specs,error_draws,target_y5,information,
                 observed_r_squared=None,recompute_r_squared=False):
        self.source_specs=dict(source_specs)
        self.error_draws={k:np.asarray(v,float) for k,v in error_draws.items()}
        self.target_y5=dict(target_y5)
        self.information=dict(information)
        self.observed_r_squared=None if observed_r_squared is None else dict(observed_r_squared)
        self.recompute_r_squared=bool(recompute_r_squared)
        if not self.recompute_r_squared and self.observed_r_squared is None:
            raise ValueError('v1.3 requires observed adjusted R² fixed before resampling')
        if self.observed_r_squared is not None and set(self.observed_r_squared)!=set(self.information):
            raise ValueError('observed R² must cover the information-weighted cohorts')
        if any(v not in ('Y5','S3') for v in self.source_specs.values()):
            raise ValueError('carried source specification must be fixed')
        if set(self.error_draws)!=set(self.source_specs):
            raise ValueError('source draw/spec mismatch')
        weights=np.asarray(list(self.information.values()),float)
        if not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0:
            raise ValueError('fixed outcome-free information weights required')

    def evaluate(self,design,source_variances):
        r_squared=design.stage_r_squared() if self.recompute_r_squared else self.observed_r_squared
        if set(r_squared)!=set(self.information):
            raise ValueError('information weights must cover the sampled cohorts')
        target={}
        for cohort in sorted(set(design.cohort)):
            people=[p for p in design.participants if p.cohort==cohort]
            collapsed=cohort=='GSE213621'
            labels=[p.stage if collapsed else self.target_y5[p.individual.split('__copy')[0]] for p in people]
            target[cohort]=target_distribution(labels,collapsed)
        result={}
        denominator=sum(self.information.values())
        for source in self.source_specs:
            per={c:conditional_reliability(marginal_reliability(*params,source_variances[source]),r_squared[c])
                 for c,params in target.items()}
            # Zero-weight cohorts must not poison a pool with an undefined lambda.
            pool=sum(w*per[c] for c,w in self.information.items() if w>0)/denominator
            result[source]={'per_cohort':per,'lambda':float(pool)}
        return {'sources':result,'r_squared':r_squared}

    def __call__(self,design,b):
        variances={source:self.error_draws[source][b] for source in self.source_specs}
        result=self.evaluate(design,variances)
        return np.asarray([result['sources'][source]['lambda'] for source in self.source_specs])
