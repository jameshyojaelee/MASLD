"""v1 test drivers; factory inputs must recenter on the complete sampled F.

Fixtures may supply S* and cluster-resample factories directly. Lambda draws
must be paired to the same resample; this module never estimates them from
allelic data. Failed refits count as exceedances.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm
from likelihood_v1 import Model,conditional_draw


def rng_for(test,b):
    return np.random.default_rng(np.random.SeedSequence([20260923,test,b]))


def mc_result(exceed):
    b=len(exceed)
    p=(1+sum(exceed))/(b+1)
    return {'p':p,'mc_se':float(np.sqrt(p*(1-p)/b)),
            'exceedance':list(map(bool,exceed))}


def ordinal_stage_probabilities(stage,cohort,x):
    """Cohort-specific ordered thresholds and common slopes; no intercept in x."""
    stage=np.asarray(stage,int);x=np.asarray(x,float)
    groups=sorted(set(cohort));ci=np.array([groups.index(c) for c in cohort])
    c=len(groups)
    def probabilities(v):
        low=v[:c];high=low+np.exp(v[c:2*c])
        predictor=x@v[2*c:]
        p0=expit(low[ci]-predictor);p01=expit(high[ci]-predictor)
        return np.column_stack((p0,p01-p0,1-p01))
    def loss(v):
        p=probabilities(v)
        chosen=np.maximum(p[np.arange(len(stage)),stage],1e-14)
        qlow=p[:,0];qhigh=p[:,:2].sum(axis=1)
        dl=qlow*(1-qlow)*np.where(stage==0,1,np.where(stage==1,-1,0))/chosen
        dh=qhigh*(1-qhigh)*np.where(stage==1,1,np.where(stage==2,-1,0))/chosen
        gradient=np.r_[-np.bincount(ci,weights=dl+dh,minlength=c),
                       -np.bincount(ci,weights=dh*np.exp(v[c:2*c])[ci],minlength=c),
                       x.T@(dl+dh)]
        return -np.log(chosen).sum(),gradient
    initial=np.r_[np.full(c,-.7),np.full(c,np.log(1.4)),np.zeros(x.shape[1])]
    result=minimize(loss,initial,jac=True,method='BFGS',options={'gtol':1e-7,'maxiter':2000})
    return {'probabilities':probabilities(result.x),'parameters':result.x,
            'converged':bool(np.max(np.abs(result.jac),initial=0)<1e-6),'message':str(result.message)}


def p1(observed_model,draw_model_factories,clusters=None):
    observed=observed_model.fit()
    if not observed['converged']:return {'computed':False,'reason':'observed fit failed'}
    se=observed_model.sandwich(observed['parameters'],clusters)
    if not np.isfinite(se) or se<=0:return {'computed':False,'reason':'observed sandwich SE invalid'}
    z=observed['parameters'][2*observed_model.G]/se
    exceed=[];failed=[];statistics=[]
    for b,factory in enumerate(draw_model_factories):
        try:
            model=factory();fit=model.fit(observed['parameters'])
            sd=model.sandwich(fit['parameters'],clusters)
            value=fit['parameters'][2*model.G]/sd
            okay=fit['converged'] and np.isfinite(value)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):okay=False
        if not okay:failed.append(b);value=float('nan')
        statistics.append(value);exceed.append(not okay or abs(value)>=abs(z)-1e-6)
    return {**mc_result(exceed),'computed':True,'z':float(z),'failed':failed,'draw_z':statistics}


def p2_kappa(kappa_draws,lambda_draws,delta=.125):
    """Rows of lambda_draws are bootstrap draws, columns are label sources."""
    kappa=np.asarray(kappa_draws,float);lam=np.asarray(lambda_draws,float)
    if lam.ndim==1:lam=lam[:,None]
    if len(kappa)!=len(lam):raise ValueError('unpaired bootstrap draws')
    result=[]
    for col in lam.T:
        invalid=(col<=.05)|~np.isfinite(col)|~np.isfinite(kappa)
        ratio=np.divide(kappa,col,out=np.full_like(kappa,np.nan),where=~invalid)
        up=mc_result(invalid|(ratio>=delta));lo=mc_result(invalid|(ratio<=-delta))
        result.append({'p_up':up['p'],'p_lo':lo['p'],
                       'mc_se_up':up['mc_se'],'mc_se_lo':lo['mc_se'],
                       'invalid':int(invalid.sum()),'invalid_draws':np.flatnonzero(invalid).tolist(),
                       'interval90':np.nanquantile(ratio,[.05,.95]).tolist() if (~invalid).any() else [None,None]})
    return {'sources':result,'p_kappa':max(max(x['p_up'],x['p_lo']) for x in result)}


def p2_tau(model,marginal_fit,tau_max,lambda_draws,progress=None):
    lam=np.asarray(lambda_draws,float)
    if lam.ndim==1:lam=lam[:,None]
    results=[]
    for source,col in enumerate(lam.T):
        # v1.2: retain all draws, entering invalid reliability as zero.
        cleaned=np.where(np.isfinite(col)&(col>.05),col,0.)
        fifth=float(np.quantile(cleaned,.05))
        if fifth<=.05:
            results.append({'p_tau':1.,'lambda_fifth':fifth,'reason':'lambda fifth percentile <= 0.05'})
            continue
        tau0=tau_max*fifth
        if progress is not None:progress('tau_profile_start',{'source':source,'tau0':tau0})
        fit=model.fit_marginal(tau_max,marginal_fit['parameters'],fixed_v=tau0**2,
                               null_fit=marginal_fit.get('null'),profile_start=marginal_fit['parameters'])
        if progress is not None:progress('tau_profile_fit',{'source':source,'fit':fit})
        if not fit['converged'] or not marginal_fit['converged']:
            results.append({'p_tau':1.,'reason':'profile fit failed'});continue
        maximum=max(marginal_fit['loglik'],marginal_fit.get('null',{}).get('loglik',-np.inf))
        lr=max(0.,2*(maximum-fit['loglik']))
        root=np.sign(tau0-np.sqrt(marginal_fit['v']))*np.sqrt(lr)
        results.append({'p_tau':float(norm.sf(root)),'tau0':tau0,'signed_root':float(root)})
    return {'sources':results,'p_tau':max(r['p_tau'] for r in results)}


def p3(model,tau_max,uniforms):
    observed=model.fit_marginal(tau_max)
    if not observed['converged'] or not observed['null']['converged']:
        return {'computed':False,'reason':'observed fit failed'}
    null=observed['null'];lr=max(0.,2*(observed['loglik']-null['loglik']))
    exceed=[];failed=[];statistics=[]
    for b,draw in enumerate(uniforms):
        if np.asarray(draw).shape!=(len(model.rows),2):raise ValueError('two uniforms per ordered row required')
        try:
            rows=[]
            for row,(u1,u2) in zip(model.rows,draw):
                eta,r,_,_=model.predictors(row,null['parameters'])
                rows.append(conditional_draw(row,eta,r,u1,u2))
            fit=Model(rows).fit_marginal(tau_max,null['parameters'],section4_from_null=True)
            okay=fit['converged'] and fit['null']['converged']
            value=max(0.,2*(fit['loglik']-fit['null']['loglik']))
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):okay=False
        if not okay:failed.append(b);value=float('nan')
        statistics.append(value);exceed.append(not okay or value>=lr-1e-6)
    return {**mc_result(exceed),'computed':True,'lr':lr,'draw_lr':statistics,'failed':failed}


def holm_three(p):
    p=np.asarray(p,float)
    if p.shape!=(3,):raise ValueError('family retains three hypotheses')
    order=np.argsort(p);adjusted=np.empty(3)
    adjusted[order]=np.minimum(1,np.maximum.accumulate(p[order]*np.array([3,2,1])))
    return adjusted


def p1_design(design, stage_draws=None, draws=1000):
    """Complete-F recentering for provided B-FIX draws or seeded stage draws."""
    if stage_draws is None:
        stage_draws=design.stage_draws(draws)
    factories=(lambda s=np.asarray(stage):design.with_stage(s).model() for stage in stage_draws)
    return p1(design.model(),factories,design.clusters)


def p2(design,lambda_factory,resamples=None,draws=2000,delta=.125,tau_max=None,progress=None):
    """Lambda callback takes (sampled complete-F Design, draw b).

    It must use the corresponding public-pair error draw b, refit each target
    label distribution on this sample, and use fixed Gate A I_c and the
    observed adjusted R² (v1.3). Recomputed R² is a sensitivity only.
    Returning different target indices would invalidate the pairing.
    """
    model=design.model();observed=model.fit()
    if progress is not None:progress('observed_fit',observed)
    if not observed['converged']:
        return {'computed':False,'reason':'observed fit failed'}
    if tau_max is None:
        tau_max=delta*float(np.median(np.abs(observed['parameters'][:model.G])))
    if resamples is None:
        resamples=(design.cluster_indices(b) for b in range(draws))
    kappas=[];lambdas=[];failed=[];lambda_failed=[]
    genes={r.gene for r in design.rows}
    for b,indices in enumerate(resamples):
        sampled=design.resample(indices)
        fit=None
        try:
            lam=np.atleast_1d(lambda_factory(sampled,b)).astype(float)
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
            lam=np.full(len(getattr(lambda_factory,'source_specs',{'GSE193066':0,'PXD051911_A1':0})),np.nan)
            lambda_failed.append(b)
        try:
            fit=sampled.model().fit(observed['parameters'])
            okay=fit['converged']
            k=fit['parameters'][2*model.G]
        except (ValueError,FloatingPointError,np.linalg.LinAlgError):
            okay=False
        if not okay:
            k=float('nan');failed.append(b)
        kappas.append(float(k));lambdas.append(lam)
        if progress is not None:progress('kappa_draw',{'b':b,'fit':fit,'lambda':lam})
    kappa_result=p2_kappa(kappas,lambdas,delta)
    marginal=model.fit_marginal(tau_max,observed['parameters'])
    if progress is not None:progress('marginal_fit',marginal)
    tau_result=p2_tau(model,marginal,tau_max,lambdas,progress)
    return {'computed':True,'p':max(kappa_result['p_kappa'],tau_result['p_tau']),
            'kappa':kappa_result,'tau':tau_result,'tau_max':tau_max,
            'failed':failed,'lambda_failed':lambda_failed,
            'kappa_draws':kappas,'lambda_draws':np.asarray(lambdas).tolist()}


def p3_design(design,tau_max=None,uniforms=None,draws=1000):
    model=design.model()
    if tau_max is None:
        fit=model.fit()
        if not fit['converged']:
            return {'computed':False,'reason':'observed fit failed'}
        tau_max=.125*float(np.median(np.abs(fit['parameters'][:model.G])))
    if uniforms is None:
        uniforms=(rng_for(3,b).random((len(model.rows),2)) for b in range(draws))
    return p3(model,tau_max,uniforms)
