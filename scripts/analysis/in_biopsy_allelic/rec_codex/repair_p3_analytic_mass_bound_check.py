#!/usr/bin/env python3
"""Test a proposed analytic mass initializer; no integrator or fit changes."""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
COLLECTION=REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z/five_gene_accuracy_collection.json'
COUNTS=(5,25,75,200,500)


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec)
    sys.modules[name]=m;spec.loader.exec_module(m);return m


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins changed')
    pins=json.loads(a.pins.read_text())
    for p,h in pins.items():
        p=Path(p);require(p.is_absolute() and p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Input changed')
        if p.suffix=='.py':ast.parse(p.read_text())
    sources=('likelihood_v1.py','repair_p3_bb_gamma_difference_v3.py')
    require(all(str(HERE/n) in pins for n in sources) and str(COLLECTION) in pins,'Required source/reference absent')
    collection=json.loads(COLLECTION.read_text());require(collection['execution_complete'] is True,'Scalar family incomplete')
    refs=collection['results'];require(len(refs)==15 and {(r['gene'],r['tau_multiple']) for r in refs}
        =={(g,t) for g in range(5) for t in (.25,1.,4.)},'Scalar family differs')
    fixtures=[]
    for r in refs:
        d=REC/('repair_p3_five_gene_accuracy_'+r['job_id']);p=d/'fixture.json';rp=d/f"gene{r['gene']}_tau{r['tau_multiple']}.json"
        require(str(p) in pins and str(rp) in pins and json.loads(rp.read_text())==r['record'],'Scalar evidence differs')
        fixtures.append(p)
    require(len({sha(p) for p in fixtures})==1 and sha(fixtures[0])==collection['shared_fixture_sha256'],'Fixtures differ')
    if a.inspect:print(json.dumps(dict(points=15,source_checks=True,no_integration=True,no_fit=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Threads not fixed')
    out=REC/('repair_p3_analytic_mass_bound_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in sources)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),
            'Archived source changed')
    import numpy as np
    import scipy
    from scipy.special import erfcx
    original=load('analytic_bound_original',archive/'likelihood_v1.py')
    gamma=load('analytic_bound_gamma',archive/'repair_p3_bb_gamma_difference_v3.py')
    rows=[]
    for gene,count in enumerate(COUNTS):
        for j in range(count):
            rows.append(original.Row(gene,f'{j:04d}','fixture',20+j%17,10,.6,.3,.1,.01,.02,.05,
                (-1)**(gene+j),float(j%3-1),(float(j%5-2)/10,),float(j%7-3)/10,float(j%5-2)/10,10,2,10))
    model=original.Model(rows);p=model.initial();p[:5]=[-1.2,-.7,.4,.9,1.3];p[10]=.08;p[11]=.04
    rng=np.random.default_rng(20260923);drawn=[]
    for row in model.rows:
        eta,r,_,_=model.predictors(row,p);drawn.append(original.conditional_draw(row,eta,r,*rng.random(2)))
    fixture=json.loads(fixtures[0].read_text())
    require(np.array_equal(p,fixture['parameters']) and [r.a for r in drawn]==fixture['drawn_counts'],'Fixture reconstruction differs')
    p.setflags(write=False);bb=original.bb_terms;start=time.perf_counter()
    def save(name,data):
        with (out/(name+'.json')).open('x') as f:json.dump(data,f,allow_nan=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    save('protocol',dict(source_and_input_sha256=pins,seed=20260923,
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        formula='ell(0)+log(erfcx(L*sqrt(v)/sqrt(2))); L=2*sum(n*abs(stage))',
        computed_ell0_not_certified=True,no_integrator_change=True,no_fits=True,no_calibration=True,no_human_data=True))
    results=[];score_probes=[]
    try:
        original.bb_terms=gamma.stable_bb_terms
        for ref in refs:
            gene=ref['gene'];tau=.05*ref['tau_multiple'];rs=[r for r in drawn if r.gene==gene]
            ell0=sum(original.conditional_terms(row,*model.predictors(row,p,0.,False)[:2],False)[0] for row in rs)
            L=2*math.fsum(row.n*abs(row.stage) for row in rs);c=L*tau
            lower=ell0+float(np.log(erfcx(c/np.sqrt(2.))))
            scalar=ref['record']['scalar']['expanded'];require(math.isfinite(lower),'Nonfinite analytic lower estimate')
            gap=scalar['logmass']-scalar['estimated_logmass_error']-lower
            checks=[]
            for u in (-8*tau,-tau,tau,8*tau):
                ell=sum(original.conditional_terms(row,*model.predictors(row,p,u,False)[:2],False)[0] for row in rs)
                excess=abs(ell-ell0)-L*abs(u);require(excess<=1e-8,'Pointwise Lipschitz estimate violated')
                checks.append(dict(u=u,loglik=ell,excess_over_bound=excess))
            require(gap>=-1e-8,'Analytic lower estimate exceeds scalar reference')
            require(float(erfcx(0.))==1.,'Zero-L constant-law identity differs')
            results.append(dict(gene=gene,tau=tau,L=L,ell0=ell0,logmass_lower_estimate=lower,
                scalar_logmass=scalar['logmass'],gap_after_scalar_error=gap,pointwise_checks=checks))
        save('mass_results',dict(points=results))
        # Assay-native adversarial score checks include selection-window edges,
        # zero mixture priors and the original deep-count allowance.
        for n in (20,36,3000):
            template=replace(drawn[0],n=n);lo=max(template.minor,(template.percent*n+99)//100)
            for priors in ((.6,.3,.1),(1.,0.,0.),(0.,1.,0.),(.5,.5,0.)):
                for count in (lo,n-lo):
                    row=replace(template,a=count,p_ref=priors[0],p_alt=priors[1],p_x=priors[2])
                    require(np.all(row.prior>=0) and row.prior.sum()>0,'Nonzero genotype support required')
                    for eta in (-40.,-10.,0.,10.,40.):
                        _,score,_=original.conditional_terms(row,eta,-3.,False);step=1e-4
                        fd=(original.conditional_terms(row,eta+step,-3.,False)[0]
                            -original.conditional_terms(row,eta-step,-3.,False)[0])/(2*step)
                        record=dict(n=n,source_probabilities=priors,genotype_weights=row.prior.tolist(),
                            count=count,eta=eta,score=float(score[0]),finite_difference=float(fd),difference=float(fd-score[0]))
                        if not (np.isfinite(score).all() and np.isfinite(fd)
                            and abs(score[0])<=n+1e-8*n and abs(fd-score[0])<=1e-5+1e-7*abs(score[0])):
                            save('score_probe_refusal',dict(current=record,completed=score_probes))
                            raise RuntimeError('Finite eta-score bound or finite-difference check failed')
                        score_probes.append(record)
    finally:original.bb_terms=bb
    require(original.bb_terms is bb and np.array_equal(p,fixture['parameters']),'Original state changed')
    for path,h in pins.items():require(sha(Path(path))==h,'Input changed during check')
    save('results',dict(points=results,score_probes=score_probes))
    save('summary',dict(completed=True,points=len(results),score_probes=len(score_probes),
        max_score_FD_difference=max(abs(x['difference']) for x in score_probes),
        min_scalar_lower_gap=min(x['gap_after_scalar_error'] for x in results),seconds=time.perf_counter()-start,
        only_initializer_hypothesis_checked=True,no_integrator_change=True,no_fits=True,no_calibration=True))
    print(json.dumps(dict(completed=True,points=len(results),score_probes=len(score_probes))))


if __name__=='__main__':main()
