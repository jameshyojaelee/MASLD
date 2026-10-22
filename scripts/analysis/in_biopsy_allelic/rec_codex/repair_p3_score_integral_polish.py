"""Isolated null-point polish using directional-score likelihood increments.

The exact identity is Delta ell=integral_0^1 grad ell(p+t*d).d dt.
16/32-rule agreement is empirical, not a certified sign/error bound.
This diagnostic changes neither an existing fit nor its failure classification.
"""
import hashlib
import math
import os
from pathlib import Path
import threading

ORIGINAL_SHA='41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
GAMMA_SHA='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def polish(model,start,bounds,original,gamma,maxiter=20):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and threading.active_count()==1,'Individual single-thread compute allocation required')
    import numpy as np
    for module,digest in ((original,ORIGINAL_SHA),(gamma,GAMMA_SHA)):
        require(hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()==digest,'Numerical source differs')
    require(model.G==8 and model.Z==7 and model.size==28 and 0<len(model.rows)<=256
        and all(type(r.n) is int and 0<=r.n<=5000 for r in model.rows),
        'Small synthetic null-point scope exceeded')
    require(model.evaluate.__func__ is original.Model.evaluate
        and model.predictors.__func__ is original.Model.predictors,'Original nuisance evaluator required')
    require(len(bounds)==28 and type(maxiter) is int and 0<maxiter<=20,'Diagnostic iteration/bound scope differs')
    input_p=np.asarray(start);before=input_p.tobytes();p=np.asarray(start,dtype=float).copy()
    require(p.shape==(28,) and np.isfinite(p).all(),'Finite fixed synthetic start required')
    require(all((lo is None or p[j]>=lo) and (hi is None or p[j]<=hi) for j,(lo,hi) in enumerate(bounds)),
        'Start outside original bounds')
    rows_before=repr(model.rows);original_bb=original.bb_terms
    require(original_bb.__name__=='bb_terms','Unpatched original BB required')
    rules={order:np.polynomial.legendre.leggauss(order) for order in (16,32)}
    history=[];evaluations=0;original_integral_evaluations=0
    def evaluate(point,with_hessian=True):
        nonlocal evaluations
        require(original.bb_terms is original_bb,'Original BB changed between evaluations')
        try:
            original.bb_terms=gamma.stable_bb_terms
            result=model.evaluate(point,need_hessian=with_hessian);evaluations+=1
        finally:original.bb_terms=original_bb
        require(np.isfinite(result[0]) and np.isfinite(result[1]).all()
            and (not with_hessian or np.isfinite(result[2]).all()),'Nonfinite nuisance evaluation')
        return result
    def projected(point,g):
        answer=g.copy()
        for j,(lo,hi) in enumerate(bounds):
            if ((lo is not None and lo==hi) or (lo is not None and point[j]<=lo+1e-8 and g[j]<0)
                or (hi is not None and point[j]>=hi-1e-8 and g[j]>0)):answer[j]=0.
        return answer
    def integral(point,delta,order,original_arithmetic=False):
        nonlocal original_integral_evaluations
        nodes,weights=rules[order];parts=[]
        for node,weight in zip(nodes,weights):
            trial=point+(.5+.5*node)*delta
            if original_arithmetic:
                _,g,_=model.evaluate(trial,need_hessian=False);original_integral_evaluations+=1
                require(np.isfinite(g).all(),'Nonfinite original integral score')
            else:_,g,_=evaluate(trial,False)
            parts.append(float(.5*weight)*float(np.dot(g,delta)))
        return math.fsum(parts),math.fsum(abs(x) for x in parts)
    for iteration in range(maxiter):
        ll,g,h=evaluate(p);old_ll,old_g,_=model.evaluate(p,need_hessian=False)
        maximum=float(np.max(np.abs(projected(p,g))))
        original_maximum=float(np.max(np.abs(projected(p,old_g))))
        record=dict(iteration=iteration,parameters=p.tolist(),Gamma_loglik=float(ll),original_loglik=float(old_ll),
            Gamma_max_projected_gradient=maximum,original_max_projected_gradient=original_maximum,
            trials=[]);history.append(record)
        if max(maximum,original_maximum)<1e-6:record['already_converged']=True;break
        # The active set depends on score direction, exactly as in the original.
        active=np.array([(lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0)
            or (hi is not None and p[j]>=hi-1e-8 and g[j]>0) for j,(lo,hi) in enumerate(bounds)])
        free=np.flatnonzero(~active)
        if not len(free):break
        curvature=-h[np.ix_(free,free)];values,vectors=np.linalg.eigh((curvature+curvature.T)/2)
        floor=max(1e-9,float(np.max(np.abs(values)))*1e-10)
        regularized=(vectors*np.maximum(values,floor))@vectors.T
        step=np.zeros(28);step[free]=np.linalg.solve(regularized,g[free]);accepted=False
        for damping in 2.**-np.arange(24):
            trial=p+damping*step
            for j,(lo,hi) in enumerate(bounds):
                trial[j]=np.clip(trial[j],-np.inf if lo is None else lo,np.inf if hi is None else hi)
            delta=trial-p
            if not np.any(delta):continue
            tl,tg,_=evaluate(trial,False);ol,og,_=model.evaluate(trial,need_hessian=False)
            new_max=float(np.max(np.abs(projected(trial,tg))))
            new_original_max=float(np.max(np.abs(projected(trial,og))))
            evidence=dict(damping=float(damping),parameters=trial.tolist(),Gamma_max_projected_gradient=new_max,
                original_max_projected_gradient=new_original_max,Gamma_LL_difference=float(tl-ll),
                original_LL_difference=float(ol-old_ll),step_max=float(np.max(np.abs(delta))))
            record['trials'].append(evidence)
            # A strict decrease alone can be rounding noise at tiny damping.
            if (new_max>=(1-.5*damping)*maximum
                    or new_original_max>=(1-.5*damping)*original_maximum):continue
            coarse,_=integral(p,delta,16);fine,absolute=integral(p,delta,32)
            difference=abs(fine-coarse);rounding=64*np.finfo(float).eps*absolute
            lower=fine-8*difference-rounding
            old_coarse,_=integral(p,delta,16,True);old_fine,old_absolute=integral(p,delta,32,True)
            old_lower=old_fine-8*abs(old_fine-old_coarse)-64*np.finfo(float).eps*old_absolute
            noise_proxy=abs(old_fine-fine)
            evidence.update(score_integral16=coarse,score_integral32=fine,
                quadrature_difference=difference,floating_sum_allowance=rounding,
                empirical_gain_lower_estimate=lower,original_score_integral16=old_coarse,
                original_score_integral32=old_fine,original_empirical_gain_lower_estimate=old_lower,
                arithmetic_noise_proxy=noise_proxy,noise_multiplier=32,
                floating_sum_allowance_excludes_score_error=True,empirical_not_certified=True)
            if min(lower,old_lower)>32*noise_proxy:
                p=trial;accepted=True;evidence['accepted']=True;break
        record['accepted_step']=accepted
        if not accepted:break
    ll,g,_=evaluate(p);ol,og,_=model.evaluate(p,need_hessian=False)
    require(original.bb_terms is original_bb and input_p.tobytes()==before
        and repr(model.rows)==rows_before,'Fixed inputs or original BB changed')
    maximum=float(np.max(np.abs(projected(p,g))));original_maximum=float(np.max(np.abs(projected(p,og))))
    return p,dict(history=history,Gamma_evaluations=evaluations,
        original_integral_evaluations=original_integral_evaluations,Gamma_loglik=float(ll),original_loglik=float(ol),
        Gamma_max_projected_gradient=maximum,original_max_projected_gradient=original_maximum,
        converged=max(maximum,original_maximum)<1e-6,scientific_convergence_cutoff=1e-6,
        no_original_result_reclassification=True,empirical_gain_checks_not_certified=True)
