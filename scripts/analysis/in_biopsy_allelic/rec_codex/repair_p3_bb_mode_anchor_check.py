#!/usr/bin/env python3
"""Bounded synthetic anchored-BB regression; no integration, fit or adoption."""
import argparse
import ast
import csv
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = Path(__file__).resolve().parent
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
CAMPAIGN = REC/'p3_parallel_20261003T130323Z'
KERNEL = REC/'repair_p3_launch_20261005T121740Z/repair_p3_integral.py'
DECIMAL_SOURCE = HERE/'repair_p3_boundary_decimal_fast.py'
PROTOTYPE = HERE/'repair_p3_bb_mode_anchor.py'
GENES = (1,2,3)
STEPS = ('1e-15','1e-20')
DIGITS = (80,110)
COUNT_SHA = '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
CAP = 32*1024**2


def require(ok,message):
    if not ok: raise ValueError(message)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name,path):
    spec = importlib.util.spec_from_file_location(name,path)
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def decimal_fd(n,k,eta,r,bb_vector):
    """Probability-only closed rising-prefix expression; no analytic scores."""
    records = []
    for digits in DIGITS:
        with localcontext() as context:
            context.prec = digits
            one,two = Decimal(1),Decimal(2)
            e,rr = Decimal.from_float(float(eta)),Decimal.from_float(float(r))
            eps,scale = Decimal.from_float(1e-6),Decimal.from_float(1-2e-6)
            logchoose = [Decimal(math.comb(n,j)).ln() for j in range(n+1)]
            def probability_log(x,y):
                rho = eps+scale/(one+(-y).exp())
                t = (one-rho)/rho
                a,b = t/(one+(-x).exp()),t/(one+x.exp())
                return bb_vector((k,),n,a,b,logchoose)[0]
            f0 = probability_log(e,rr)
            for step in STEPS:
                h = Decimal(step)
                pe,me = probability_log(e+h,rr),probability_log(e-h,rr)
                pr,mr = probability_log(e,rr+h),probability_log(e,rr-h)
                mixed = (probability_log(e+h,rr+h)-probability_log(e+h,rr-h)
                         -probability_log(e-h,rr+h)+probability_log(e-h,rr-h))/(4*h*h)
                g = [(pe-me)/(two*h),(pr-mr)/(two*h)]
                H = [[(pe-two*f0+me)/(h*h),mixed],[mixed,(pr-two*f0+mr)/(h*h)]]
                records.append(dict(digits=digits,step=step,logpmf=str(f0),
                                    score=list(map(str,g)),hessian=[[str(x) for x in row] for row in H]))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-prototype-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args = parser.parse_args()
    pins = {Path(__file__).resolve():args.expected_checker_sha,PROTOTYPE:args.expected_prototype_sha,
        KERNEL:'461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3',
        DECIMAL_SOURCE:'d192f1aa311b54d3fe94af3744cf14a8c02856a95ffc5e3e590c7efa1432c017',
        CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
        FIXTURE/'SHA256SUMS':'0d3d06fc32a8e588d5e1c54c943a582e5054b6d2ab565290a4d57ced545d2925'}
    for path,digest in pins.items(): require(sha(path)==digest,'Source/input pin differs: '+str(path))
    campaign = json.loads((CAMPAIGN/'campaign.json').read_text())
    for name in ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py'):
        path = CAMPAIGN/'executed_sources'/name
        pins[path] = campaign['files'][str(path)]
    published = {name.removeprefix('./'):digest for digest,name in
        (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','genes.tsv','tags.tsv',
                 'params/theta4_perturbed.json'):
        pins[FIXTURE/name] = published[name]
    for path,digest in pins.items():
        require(sha(path)==digest,'Source/input pin differs: '+str(path))
        if path.suffix=='.py': ast.parse(path.read_text())
    # Explicit bounded families, fixed before any numerical evaluation.
    moderate = [(n,n//3,eta,r) for n in (5,20,100) for eta,r in ((-2.,-3.),(0.,0.),(2.,2.))]
    broad = [(0,e,r) for e in (-1200.,0.,1200.) for r in (-15.,0.,15.)]
    broad += [(n,e,r) for n in (5,20,100) for e,r in ((-2.,-3.),(0.,0.),(2.,2.))]
    broad += [(n,e,r) for n in (2000,2466) for e,r in
              ((-1200.,0.),(0.,0.),(1200.,0.),(0.,-15.),(0.,15.))]
    fixture_settings=json.loads((FIXTURE/'globals.json').read_text())
    with (FIXTURE/'rows.tsv').open() as handle:
        metadata_rows=list(csv.DictReader(handle,delimiter='\t'))
    representative_depths={gene:max(int(row['n']) for row in metadata_rows
        if row['gene_id']==fixture_settings['gene_order'][gene]) for gene in GENES}
    planned_work=108*(sum(n for n,_,_,_ in moderate)+3*sum(representative_depths.values()))
    protocol = dict(source_input_pins={str(p):h for p,h in pins.items()},genes=GENES,
        theta='theta4_perturbed',seed=20260930,expected_count_sha=COUNT_SHA,
        moderate_cases=moderate,broad_cases=broad,decimal_digits=DIGITS,decimal_steps=STEPS,
        independent_family_points=18,representative_rule='maximum n row per target gene, original row order breaks ties; u=0,+sqrt(v),-sqrt(v)',
        native_conditional_comparisons=9,expected_conditional_rows=210,
        maximum_decimal_n=5000,maximum_decimal_prefix_log_evaluations=20000000,
        representative_depths=representative_depths,
        planned_decimal_prefix_log_evaluation_upper_bound=planned_work,
        independent_criterion=dict(logpmf_absolute=1e-9,score_max_absolute=1e-8,hessian_max_absolute=1e-5,
            decimal_step_and_precision_stability=1e-25,scope='18 selected logPMF points only; not all tails, conditional mixtures or integration'),
        old_helper_is_not_ground_truth=True,output_cap_bytes=CAP,
        prototype_peak_memory_estimate_bytes=48*1024**2,memory_estimate_not_measurement=True,
        no_integration=True,no_fits=True,no_calibration=True,no_adoption=True)
    require(max(representative_depths.values())<=protocol['maximum_decimal_n']
            and planned_work<=protocol['maximum_decimal_prefix_log_evaluations'],
            'Pinned regression metadata exceeds prespecified compute bounds')
    if args.inspect:
        print(json.dumps(protocol)); return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output = REC/('repair_p3_bb_mode_anchor_check_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False); archive=output/'executed_sources'; archive.mkdir()
    for index,(path,digest) in enumerate(pins.items()):
        target = archive/path.name if path.suffix=='.py' else output/('input_'+str(index)+'_'+path.name)
        require(not target.exists(),'Duplicate archive basename')
        target.write_bytes(path.read_bytes()); require(sha(target)==digest,'Archive differs')
    sys.path.insert(0,str(archive))
    require(not any(name in sys.modules for name in ('bfix_v1','likelihood_v1','design_v1','tests_v1')),
            'Original numerical fixture modules already imported')
    import numpy as np
    import scipy
    from scipy.special import expit
    versions = dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    require(versions==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Numerical environment differs')
    def safe(value):
        if isinstance(value,np.ndarray): return safe(value.tolist())
        if isinstance(value,np.generic): return safe(value.item())
        if isinstance(value,float) and not math.isfinite(value): return {'nonfinite_float':str(value)}
        if isinstance(value,dict): return {str(k):safe(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)): return [safe(v) for v in value]
        return value
    def save(name,value):
        body=(json.dumps(safe(value),allow_nan=False,separators=(',',':'))+'\n').encode()
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body)<=CAP,
                'Saved-output32MiB cap exceeded')
        with (output/(name+'.json')).open('xb') as handle: handle.write(body)
    def attempt(function):
        began=time.perf_counter()
        try: return dict(status='returned',value=function(),seconds=time.perf_counter()-began)
        except Exception as exc:
            return dict(status='refused',error_type=type(exc).__name__,message=str(exc),
                        seconds=time.perf_counter()-began)
    kernel=load('_mode_check_kernel',archive/KERNEL.name)
    prototype=load('_mode_check_prototype',archive/PROTOTYPE.name)
    decimal_source=load('_mode_check_decimal',archive/DECIMAL_SOURCE.name)
    loader=importlib.import_module('bfix_v1'); importlib.import_module('tests_v1')
    for name in ('bfix_v1','likelihood_v1','design_v1','tests_v1'):
        require(Path(sys.modules[name].__file__).resolve().parent==archive.resolve(),
                'Fixture import escaped pinned archive')
    design,settings=loader.load(); model=design.model()
    require(len(design.participants)==60 and len(model.rows)==173 and model.G==8 and model.Z==7,
            'Exact fixture dimensions differ')
    p,v=loader.parameters(FIXTURE/'params/theta4_perturbed.json')
    counts=hashlib.sha256(np.asarray([row.a for row in model.rows],dtype='<i8').tobytes()).hexdigest()
    require(counts==COUNT_SHA and v==.0004 and settings['seed']==20260930,'Exact fixture count/variance/seed differs')
    protocol.update(environment=versions,p=p,v=v,counts_sha=counts)
    target_rows={gene:[row for row in model.rows if row.gene==gene] for gene in GENES}
    require([len(target_rows[g]) for g in GENES]==[22,24,24],'Target gene row counts differ')
    points=[dict(n=n,k=k,eta=e,r=r,origin='moderate') for n,k,e,r in moderate]
    for gene in GENES:
        row=max(target_rows[gene],key=lambda row:row.n)
        for u in (0.,math.sqrt(v),-math.sqrt(v)):
            eta,r=model.predictors(row,p,u)[:2]
            points.append(dict(n=row.n,k=row.a,eta=eta-row.orientation*row.omega,r=r,
                               origin='theta4_representative',gene=gene,individual=row.individual,u=u))
    require(max(point['n'] for point in points)<=protocol['maximum_decimal_n'],'Decimal depth bound exceeded')
    work=sum(3*point['n']*36 for point in points)
    require(work<=protocol['maximum_decimal_prefix_log_evaluations'],'Decimal work bound exceeded')
    protocol.update(actual_independent_points=points,decimal_prefix_log_evaluation_upper_bound=work)
    save('protocol',protocol)
    independent=[]
    for index,point in enumerate(points):
        n,k,eta,r=(point[key] for key in ('n','k','eta','r'))
        new=attempt(lambda:prototype.stable_bb_terms([k],n,eta,r))
        reference=attempt(lambda:decimal_fd(n,k,eta,r,decimal_source.bb_vector))
        record=dict(point=point,candidate=new,reference=reference,selected_family_criterion_met=False)
        if new['status']==reference['status']=='returned':
            lp,g,H=new['value']; ladder=reference['value']; final=ladder[-1]
            with localcontext() as context:
                context.prec=110
                flatten=lambda a:[Decimal(a['logpmf']),*map(Decimal,a['score']),
                                  *(Decimal(x) for row in a['hessian'] for x in row)]
                differences=[[str(x-y) for x,y in zip(flatten(a),flatten(final))] for a in ladder]
                stability=max(abs(Decimal(x)) for row in differences for x in row)
                raw=[Decimal.from_float(float(lp[0]))-Decimal(final['logpmf']),
                     *(Decimal.from_float(float(x))-Decimal(y) for x,y in zip(g[0],final['score'])),
                     *(Decimal.from_float(float(x))-Decimal(y) for x,y in zip(H[0].flat,[x for row in final['hessian'] for x in row]))]
                record.update(raw_candidate_minus_decimal=list(map(str,raw)),
                    decimal_ladder_minus_final=differences,max_decimal_stability_difference=str(stability))
                record['selected_family_criterion_met']=bool(np.isfinite(lp).all() and np.isfinite(g).all()
                    and np.isfinite(H).all() and abs(raw[0])<=Decimal('1e-9')
                    and max(map(abs,raw[1:3]))<=Decimal('1e-8')
                    and max(map(abs,raw[3:]))<=Decimal('1e-5') and stability<=Decimal('1e-25'))
        save('independent_'+str(index),record); independent.append(record)
    for index,(n,eta,r) in enumerate(broad):
        ks=np.arange(n+1)
        old=attempt(lambda:kernel.stable_bb_terms(ks,n,eta,r))
        new=attempt(lambda:prototype.stable_bb_terms(ks,n,eta,r))
        record=dict(n=n,eta=eta,r=r,old=old,new=new,descriptive_only=True)
        if new['status']=='returned':
            lp,g,H=new['value']; probability=np.exp(lp)
            m,q=expit(eta),expit(-eta)
            rho=1e-6+(1-2e-6)*expit(r)
            expected_mean=n*m
            expected_variance=n*m*q*(1+(n-1)*rho)
            empirical_mean=math.fsum(map(float,probability*ks))
            variance_about_native_mean=math.fsum(map(float,probability*(ks-expected_mean)**2))
            record['invariants']=dict(total_probability=math.fsum(map(float,probability)),
                weighted_score=np.sum(probability[:,None]*g,axis=0),
                weighted_H_plus_score_outer=np.sum(probability[:,None,None]*(H+g[:,:,None]*g[:,None,:]),axis=0),
                native_expected_mean=expected_mean,native_expected_variance=expected_variance,
                empirical_mean=empirical_mean,variance_about_native_mean=variance_about_native_mean,
                mean_difference=empirical_mean-expected_mean,
                variance_difference=variance_about_native_mean-expected_variance,
                variance_relative_difference=None if expected_variance==0 else
                    (variance_about_native_mean-expected_variance)/expected_variance,
                zero_variance_float_extreme_or_n0=bool(expected_variance==0),
                identities_not_independent_accuracy_certificate=True,
                finite=bool(np.isfinite(lp).all() and np.isfinite(g).all() and np.isfinite(H).all()))
        if old['status']==new['status']=='returned':
            record['raw_new_minus_old']=[a-b for a,b in zip(new['value'],old['value'])]
        save('broad_'+str(index),record)
    oldmodule=load('_mode_check_old_native',archive/'likelihood_v1.py')
    newmodule=load('_mode_check_new_native',archive/'likelihood_v1.py')
    oldmodule.bb_terms=kernel.stable_bb_terms; newmodule.bb_terms=prototype.stable_bb_terms
    for gene in GENES:
        for position,u in enumerate((0.,math.sqrt(v),-math.sqrt(v))):
            def evaluate(module):
                total,gradient,hessian=0.,np.zeros(model.size),np.zeros((model.size,model.size))
                rows=[]
                for row in target_rows[gene]:
                    eta,r,jac,etaH=model.predictors(row,p,u,need_hessian=True)
                    ll,g,H=module.conditional_terms(row,eta,r,need_hessian=True)
                    total+=ll; gradient+=jac.T@g; hessian+=jac.T@H@jac+g[0]*etaH
                    rows.append(dict(individual=row.individual,n=row.n,a=row.a,eta=eta,r=r,ll=ll,g=g,H=H))
                return dict(ll=total,g=gradient,H=hessian,rows=rows)
            old,new=attempt(lambda:evaluate(oldmodule)),attempt(lambda:evaluate(newmodule))
            record=dict(gene=gene,u=u,old=old,new=new,descriptive_only=True)
            if old['status']==new['status']=='returned':
                record['raw_new_minus_old']={key:new['value'][key]-old['value'][key] for key in ('ll','g','H')}
                record['row_raw_new_minus_old']=[{key:b[key]-a[key] for key in ('ll','g','H')}
                    for a,b in zip(old['value']['rows'],new['value']['rows'])]
            save('conditional_gene'+str(gene)+'_u'+str(position),record)
    for path,digest in pins.items():
        require(sha(path)==digest,'Source/input changed during regression: '+str(path))
    save('summary',dict(required_independent_points=18,completed_independent_points=len(independent),
        selected_decimal_family_criterion_met=all(r['selected_family_criterion_met'] for r in independent),
        criterion_scope=protocol['independent_criterion']['scope'],broad_points=len(broad),
        no_universal_allpass=True,no_integration_or_model_accuracy_certificate=True,
        no_adoption=True,no_calibration=True))


if __name__=='__main__': main()
