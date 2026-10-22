#!/usr/bin/env python3
"""Independent scalar check of a Gaussian variance-score identity.

For I(v)=integral exp(ell(u))*Normal(u;0,v)du, integration by parts gives
I'(v)/I(v)=.5 E_posterior[ell_uu+ell_u**2]. This avoids subtracting nearly
equal Gaussian second moments. Only two saved synthetic gene points; no fit.
Quadrature error estimates and mass-based allowances are not certificates.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
POINT=REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_bb_gamma_difference_v3.py',
    'repair_p3_batched_conditional.py','repair_p3_prepared_conditional_v2.py')
GL_POINTS={2:REC/'GL_initial64_refusal_gene2_22084480/summary.json',
    5:REC/'GL_initial64_refusal_gene5_22084481/summary.json'}


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gene',type=int,choices=(2,5),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    for path,h in pins.items():
        p=Path(path);require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),*(POINT/n for n in ('protocol.json','recipe_summary.json','draw_counts.json')),
        *GL_POINTS.values(),*(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Missing source/input pin')
    read=lambda p:json.loads(p.read_text());failure=read(POINT/'recipe_summary.json')['fitting_integration_failure']
    gl=read(GL_POINTS[a.gene]);protocol=read(POINT/'protocol.json')
    require(failure['call']==14 and gl['parameters']==failure['parameters'] and gl['variance']==failure['variance']
        and gl['gene']==a.gene,'Saved point differs')
    require(len(gl['parameters'])==28 and all(math.isfinite(x) for x in gl['parameters'])
        and len(gl['gradient'])==29 and all(math.isfinite(x) for x in gl['gradient'])
        and gl['report']['full_parameter_dimension']==28
        and gl['report']['augmented_global_indices']==sorted([a.gene,8+a.gene,*range(16,28)])+[28],
        'Expected gene-local nuisance and final variance coordinates')
    require(gl['source_and_input_sha256'][str(POINT/'draw_counts.json')]==pins[str(POINT/'draw_counts.json')],
        'GL comparison counts differ')
    # Bind the GL final coordinate to the reviewed d/dv moment formula,
    # rather than trusting an unqualified final gradient entry.
    for name,expected in (
        ('repair_p3_local_moments.py','95ec9f270d18bfbbc8f31f0857308f5f2acad053661a30b208049b7a3f59f31f'),
        ('repair_p3_integral.py','461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'),
        ('repair_p3_local_gl_initial64.py','cea90f5189767b99de7171d8506cad2f583e4cd6c7eac533b8f31cf6b8bbde29')):
        matches=[h for path,h in gl['source_and_input_sha256'].items() if Path(path).name==name]
        require(matches==[expected],'GL derivative formula/source differs')
    require(all(protocol['source_and_input_sha256'][str(HERE/n)]==pins[str(HERE/n)] for n in SOURCES),
        'Original fixed-point numerical sources differ')
    if a.inspect:print(json.dumps(dict(gene=a.gene,backends=['Gamma_scalar','count_prefix'],bounds=[24,26],no_fit=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Single numerical thread required')
    out=REC/('heat_variance_score_gene'+str(a.gene)+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),'Archive differs')
    import numpy as np
    import scipy
    from scipy.integrate import quad_vec
    from scipy.special import log_ndtr
    sys.path.insert(0,str(archive));original=importlib.import_module('likelihood_v1');bb=original.bb_terms
    require(bb.__name__=='bb_terms' and Path(bb.__code__.co_filename).resolve()==(archive/'likelihood_v1.py').resolve(),
        'Original unpatched BB required')
    loader=importlib.import_module('bfix_v1')
    gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3');prefix=importlib.import_module('repair_p3_batched_conditional')
    require(loader.FIXTURE.resolve()==FIXTURE.resolve(),'Fixture loader path differs')
    require(original.bb_terms is bb,'Loader changed original BB')
    for name in ('likelihood_v1','bfix_v1','design_v1','tests_v1','repair_p3_bb_gamma_difference_v3','repair_p3_batched_conditional'):
        require(Path(sys.modules[name].__file__).resolve().parent==archive.resolve(),'Import escaped archive')
    design,settings=loader.load();template=design.model();counts=read(POINT/'draw_counts.json')
    require(counts['full_row_keys']==[[r.individual,settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts'])==len(template.rows)==173,'Fixed count row order differs')
    model=original.Model([replace(r,a=c) for r,c in zip(template.rows,counts['counts'])],gene_count=8)
    rows=tuple(r for r in model.rows if r.gene==a.gene);p=np.asarray(failure['parameters']);p.setflags(write=False)
    require(model.G==8 and model.Z==7 and model.size==28 and 0<len(rows)<=32 and all(r.n<=5000 for r in rows),'Fixed synthetic scope differs')
    require(p.shape==(model.size,) and gl['rows']==len(rows),'Parameter/row scope differs')
    v=failure['variance'];sqrtv=math.sqrt(v);before=p.tobytes();row_before=repr(rows)
    context=prefix.BatchedConditional(model,p,rows);records=[];began=time.perf_counter()
    require(context.original is original and original.bb_terms is bb,'Prefix context changed original arithmetic')
    # Same conservative native selected-BB eta derivative envelopes as the
    # existing integration kernel, applied to the scalar u direction.
    cp=2*sum(r.n*abs(r.stage) for r in rows)
    hp=sum((16*r.n+8*r.n*r.n)*r.stage*r.stage for r in rows)
    heat_bound=.5*(hp+cp*cp)
    def safe(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:json.dump(safe(value),f,indent=2,allow_nan=False);f.flush();os.fsync(f.fileno())
    save('protocol',dict(gene=a.gene,parameters=p,variance=v,rows=len(rows),source_and_input_sha256=pins,
        node=os.environ.get('SLURMD_NODENAME'),environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        conditional_score_u_envelope=cp,conditional_curvature_u_envelope=hp,heat_numerator_envelope=heat_bound,
        quadrature_epsabs=1e-12,quadrature_epsrel=1e-12,quadrature_limit=512,quad_cache_bytes=4*1024**2,
        no_fit=True,no_calibration=True,no_human_data=True,estimated_errors_not_certified=True))
    for backend in ('Gamma_scalar','count_prefix'):
        calls=0
        def point(u):
            nonlocal calls
            values=[];scores=[];curvatures=[];calls+=1
            for row in rows:
                eta,r,_,_=model.predictors(row,p,u,need_hessian=False)
                if backend=='Gamma_scalar':value,score,hess=original.conditional_terms(row,eta,r,need_hessian=True)
                else:
                    value,score,hess=context.terms_batch(row,np.asarray([eta]),r,p,need_hessian=True)
                    value,score,hess=float(value[0]),score[0],hess[0]
                values.append(float(value));scores.append(float(score[0])*row.stage)
                curvatures.append(float(hess[0,0])*row.stage*row.stage)
            ll=math.fsum(values);su=math.fsum(scores);hu=math.fsum(curvatures);heat=.5*(hu+su*su)
            require(math.isfinite(ll) and math.isfinite(heat),'Nonfinite heat integrand')
            return ll,heat
        try:
            original.bb_terms=gamma.stable_bb_terms if backend=='Gamma_scalar' else bb
            ell0,heat0=point(0.);scale=max(1.,abs(heat0));domain=[]
            def integrand(z):
                ll,heat=point(float(z)*sqrtv)
                weight=math.exp(ll-ell0-.5*z*z)/math.sqrt(2*math.pi)
                return np.asarray([weight,weight*heat/scale])
            for bound in (24.,26.):
                start=time.perf_counter()
                result,error,info=quad_vec(integrand,-bound,bound,epsabs=1e-12,epsrel=1e-12,norm='max',
                    cache_size=4*1024**2,limit=512,workers=1,full_output=True)
                require(info.success and np.isfinite(result).all() and math.isfinite(error) and error>=0,'Scalar quadrature refused')
                mass,numerator=float(result[0]),float(result[1])*scale
                tail=math.exp(math.log(2.)+float(log_ndtr(-bound))-ell0)
                lower_mass=mass-error-tail;require(lower_mass>0,'Estimated mass allowance unresolved')
                score=numerator/mass
                allowance=(abs(score)*(error+tail)+scale*error+heat_bound*tail)/lower_mass
                record=dict(z_bound=bound,logmass=ell0+math.log(mass),variance_score=score,
                    estimated_logmass_allowance=(error+tail)/lower_mass,estimated_score_allowance=allowance,
                    shifted_mass=mass,shifted_heat_numerator=numerator,heat_scale=scale,
                    quadrature_error_estimate=float(error),absolute_shifted_mass_tail_bound=tail,
                    heat_numerator_tail_bound=heat_bound*tail,quad_evaluations=info.neval,
                    quad_status=info.status,quad_message=info.message,seconds=time.perf_counter()-start)
                domain.append(record);save(backend+'_bound'+str(int(bound)),record)
        finally:original.bb_terms=bb
        first,last=domain;dl=abs(first['logmass']-last['logmass']);dg=abs(first['variance_score']-last['variance_score'])
        stable=(dl+first['estimated_logmass_allowance']+last['estimated_logmass_allowance']<=1e-9
            and dg+first['estimated_score_allowance']+last['estimated_score_allowance']<=1e-8)
        records.append(dict(backend=backend,domain=domain,stable=stable,domain_LL_difference=dl,domain_score_difference=dg,
            conditional_point_calls=calls,GL_variance_score_difference=gl['gradient'][-1]-last['variance_score'],
            GL_LL_difference=gl['logmass']-last['logmass'],native40_80_differences=[dict(order=n['order'],
                LL_difference=n['logmass']-last['logmass'],variance_score_difference=n['gradient'][-1]-last['variance_score'])
                for n in gl['native40_80']]))
    left,right=(r['domain'][-1] for r in records)
    agreement=(abs(left['logmass']-right['logmass'])<=1e-9 and abs(left['variance_score']-right['variance_score'])<=1e-8)
    require(p.tobytes()==before and repr(rows)==row_before and original.bb_terms is bb,'Fixed inputs changed')
    for path,h in pins.items():require(sha(Path(path))==h,'Input changed during check')
    passed=all(r['stable'] for r in records) and agreement
    save('summary',dict(completed=True,gene=a.gene,records=records,both_backends_stable=all(r['stable'] for r in records),
        backend_agreement=agreement,reference_check_passed=passed,seconds=time.perf_counter()-began,
        empirical_not_certified=True,no_Hessian_validation=True,no_fit=True,no_calibration=True,no_human_validation=True))
    print(json.dumps(dict(completed=True,gene=a.gene,reference_check_passed=passed)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
