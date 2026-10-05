#!/usr/bin/env python3
"""Compute-only fifteen-point synthetic integration comparison; no fits."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback

HERE = Path(__file__).resolve().parent
REC = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z')
PINS = {'check_quadrature_v1.py':'40a51fd0090d797b631c897f1c35c046b750e040a780432f0f614960b9a17735',
        'likelihood_v1.py':'41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'}
COUNTS = (5,25,75,200,500)
MULTIPLES = (.25,1.,4.)
CAP = 32*1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-kernel-sha',required=True)
    parser.add_argument('--expected-gl-sha',required=True)
    parser.add_argument('--gene',type=int,choices=range(5))
    parser.add_argument('--tau-multiple',type=float,choices=MULTIPLES)
    args = parser.parse_args()
    require((args.gene is None)==(args.tau_multiple is None),'Gene and tau-multiple must be supplied together')
    selected_points = [(multiple,gene) for multiple in MULTIPLES for gene in range(5)
                       if args.gene is None or (gene==args.gene and multiple==args.tau_multiple)]
    pins = {HERE/name:digest for name,digest in PINS.items()}
    pins.update({Path(__file__):args.expected_checker_sha,
                 HERE/'repair_p3_integral.py':args.expected_kernel_sha,
                 HERE/'repair_p3_composite_gl.py':args.expected_gl_sha})
    for path,digest in pins.items():
        require(sha(path)==digest,'Reviewed source bytes differ: '+path.name)
        ast.parse(path.read_text())
    metadata = dict(seed=20260923,rows_per_gene=COUNTS,tau_max=.05,
                    tau_multiples=MULTIPLES,required_points=15,
                    requested_points=selected_points,requested_point_count=len(selected_points),
                    combined_family_all15_remains_required=True,
                    source_sha256={str(p):h for p,h in pins.items()},
                    no_model_fits=True,no_calibration=True,output_cap_bytes=CAP)
    if args.inspect:
        print(json.dumps(metadata,sort_keys=True))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required before numerical imports')
    output = REC/('repair_p3_five_gene_accuracy_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for source,digest in pins.items():
        target=archive/source.name
        target.write_bytes(source.read_bytes())
        require(sha(target)==digest,'Archived source differs')
    import numpy as np
    import scipy
    from scipy.integrate import quad
    from scipy.special import betaln,gammaln,logsumexp,expit
    from scipy.stats import norm
    import warnings
    original = load('five_gene_original',archive/'likelihood_v1.py')
    private = load('five_gene_private',archive/'likelihood_v1.py')
    kernel = load('five_gene_kernel',archive/'repair_p3_integral.py')
    gl = load('five_gene_gl',archive/'repair_p3_composite_gl.py')
    private.bb_terms = kernel.stable_bb_terms

    def serial(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)

    def save(name,value):
        body=json.dumps(value,default=serial,allow_nan=False,indent=2)+'\n'
        used=sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
        require(used+len(body.encode())<=CAP,'Synthetic output cap exhausted')
        path=output/(name+'.json')
        with path.open('x') as handle:handle.write(body)

    def diagnostic(value):
        """Preserve nonfinite diagnostic allowances without inventing numbers."""
        if isinstance(value,np.ndarray):return diagnostic(value.tolist())
        if isinstance(value,np.generic):return diagnostic(value.item())
        if isinstance(value,dict):return {k:diagnostic(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [diagnostic(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):
            return {'nonfinite_float':'nan' if np.isnan(value) else 'positive_infinity' if value>0 else 'negative_infinity'}
        return value

    metadata['environment']=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__)
    save('protocol',metadata)
    # Exact construction and RNG call order from check_quadrature_v1.py.
    rows=[]
    for gene,count in enumerate(COUNTS):
        for j in range(count):
            rows.append(original.Row(gene,f'{j:04d}','fixture',20+j%17,10,.6,.3,.1,.01,.02,.05,
                 (-1)**(gene+j),float(j%3-1),(float(j%5-2)/10,),float(j%7-3)/10,float(j%5-2)/10,10,2,10))
    model=original.Model(rows)
    p=model.initial();p[:5]=[-1.2,-.7,.4,.9,1.3];p[10]=.08;p[11]=.04
    rng=np.random.default_rng(20260923)
    drawn=[]
    for row in model.rows:
        eta,r,_,_=model.predictors(row,p)
        drawn.append(original.conditional_draw(row,eta,r,*rng.random(2)))
    model=original.Model(drawn)
    save('fixture',dict(parameters=p,drawn_counts=[r.a for r in drawn],rows_per_gene=COUNTS,
                        native_units='synthetic oriented RNA allele read counts',inferential_analysis=False))

    def direct(rs,u):
        total=0.
        for row in rs:
            # Explicit predictors, independent of conditional_terms/analytic scores.
            eta=p[row.gene]*(1+p[10]*row.stage+np.dot(p[11:12],row.z))
            eta+=u*row.stage-row.orientation*p[12]*row.stage-row.orientation*row.omega
            rr=p[5+row.gene]+np.dot(p[-3:],(row.stage,row.background,row.duplication))
            m=expit(eta);rho=1e-6+(1-2e-6)*expit(rr)
            lo=max(row.minor,(row.percent*row.n+99)//100)
            k=np.arange(lo,row.n-lo+1)
            require(row.n>=row.depth and row.a in k,'Synthetic selected support differs')
            prior=np.array([2*row.p_ref*row.p_alt,row.p_ref**2+2*row.p_ref*row.p_x,
                            row.p_alt**2+2*row.p_alt*row.p_x])
            if row.orientation==-1:prior=prior[[0,2,1]]
            parts=[]
            for mean,disp in ((m,rho),(row.e,row.phi),(1-row.e,row.phi)):
                t=(1-disp)/disp
                parts.append(gammaln(row.n+1)-gammaln(k+1)-gammaln(row.n-k+1)
                    +betaln(k+mean*t,row.n-k+(1-mean)*t)-betaln(mean*t,(1-mean)*t))
            mass=logsumexp(np.asarray(parts)+np.log(prior)[:,None],axis=0)
            total+=mass[row.a-lo]-logsumexp(mass)
        require(np.isfinite(total) and total<=1e-8,'Independent selected probability invalid')
        return float(total)

    def scalar(rs,tau,bound):
        grid=np.linspace(-bound,bound,int(np.ceil(2*bound/.25))+1)
        values=[direct(rs,tau*z)-z*z/2-.5*np.log(2*np.pi) for z in grid]
        shift=max(values)
        with warnings.catch_warnings(record=True) as caught:
            result=quad(lambda z:np.exp(direct(rs,tau*z)-z*z/2-.5*np.log(2*np.pi)-shift),
                        -bound,bound,points=grid[1:-1],epsabs=1e-12,epsrel=1e-12,limit=6000,full_output=1)
        mass,error,info=result[:3]
        require(len(result)==3 and mass>error>=0 and np.isfinite(mass),'Independent scalar integral failed')
        lower=shift+np.log(mass-error)
        require(np.isfinite(shift+np.log(mass)),'Independent scalar logmass nonfinite')
        tail=float(np.exp(np.log(2.)+norm.logsf(bound)-lower))
        return dict(logmass=float(shift+np.log(mass)),z_bound=bound,evaluations=info['neval'],
                    estimated_logmass_error=float(-np.log1p(-error/mass)+np.log1p(tail)),
                    relative_prior_mass_tail_allowance=tail,warnings=[str(w.message) for w in caught],
                    error_estimate_not_certified=True,probability_arithmetic_error_unbounded=True)

    records=[]
    for multiple in MULTIPLES:
        tau=multiple*.05
        for gene in range(5):
            if (multiple,gene) not in selected_points:continue
            rs=[r for r in drawn if r.gene==gene]
            record=dict(gene=gene,rows=len(rs),tau=tau,variance=tau*tau)
            for order in (20,40):
                try:
                    gm=original.Model(rs,gene_count=5)
                    ll,gradient=gm.marginal(p,tau*tau,nodes=order)
                    require(np.isfinite(ll) and np.isfinite(gradient).all(),'Original GH result nonfinite')
                    record['original_GH'+str(order)]=dict(logmass=ll,gradient=gradient)
                except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                    record['original_GH'+str(order)]=dict(refused=True,traceback=traceback.format_exc())
            try:
                candidate=private.Model(rs,gene_count=5)
                ll,g,h,report=gl.integrate_gene(candidate,p,tau*tau,rs,private.conditional_terms,
                    with_hessian=True,expected_kernel_sha=args.expected_kernel_sha,kernel_path=archive/'repair_p3_integral.py')
                require(np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),'Candidate output nonfinite')
                record['candidate']=dict(logmass=ll,gradient=g,hessian=h,integration=diagnostic(report))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError) as failure:
                record['candidate']=dict(refused=True,traceback=traceback.format_exc(),
                    partial_integration_reports=diagnostic(getattr(failure,'integration_reports',None)))
            try:
                bound=12.
                for attempt in range(11):
                    first=scalar(rs,tau,bound)
                    if first['relative_prior_mass_tail_allowance']<=1e-10:break
                    bound+=8
                require(first['relative_prior_mass_tail_allowance']<=1e-10,'Scalar relative tail unresolved')
                expanded=scalar(rs,tau,bound+2)
                record['scalar']=dict(first=first,expanded=expanded,
                    refinement_difference=expanded['logmass']-first['logmass'])
                require(abs(expanded['logmass']-first['logmass'])+expanded['estimated_logmass_error']<=1e-9,
                        'Scalar refinement unresolved')
                if 'logmass' in record['candidate']:
                    delta=record['candidate']['logmass']-expanded['logmass']
                    record['candidate_scalar_difference']=delta
                    record['candidate_scalar_relative_LL_pass']=bool(abs(delta)<=1e-8*abs(expanded['logmass']))
                for order in (20,40):
                    gh=record['original_GH'+str(order)]
                    if 'logmass' in gh:
                        delta=gh['logmass']-expanded['logmass']
                        gh['scalar_difference']=delta
                        gh['scalar_relative_LL_pass']=bool(abs(delta)<=1e-8*abs(expanded['logmass']))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError):
                record['scalar_refusal']=traceback.format_exc()
            if all('logmass' in record['original_GH'+str(n)] for n in (20,40)):
                delta=record['original_GH20']['logmass']-record['original_GH40']['logmass']
                record['original_GH20_GH40_difference']=delta
                record['original_GH20_GH40_absolute_pass']=bool(abs(delta)<=1e-8)
            save(f'gene{gene}_tau{multiple}',record)
            records.append(record)
            print(json.dumps(dict(gene=gene,tau=tau,completed_points=len(records),required_points=15)),flush=True)
    for path,digest in pins.items():require(sha(path)==digest,'Source changed during synthetic comparison')
    save('summary',dict(all15points_recorded=len(records)==15,
        requested_point_count=len(selected_points),all_requested_points_recorded=len(records)==len(selected_points),
        combined_family_all15_remains_required=True,
        candidate_scalar_relative_LL_pass_requested_points=all(r.get('candidate_scalar_relative_LL_pass') is True for r in records),
        candidate_scalar_relative_LL_pass_all15=len(records)==15 and all(r.get('candidate_scalar_relative_LL_pass') is True for r in records),
        original_GH20_GH40_pass_all15=len(records)==15 and all(r.get('original_GH20_GH40_absolute_pass') is True for r in records),
        original_GH20_scalar_pass_all15=len(records)==15 and all(r['original_GH20'].get('scalar_relative_LL_pass') is True for r in records),
        independent_derivative_accuracy_verified=False,no_model_fits=True,no_calibration=True,
        quadrature_and_floating_error_bounds_not_certified=True))


if __name__=='__main__':main()
