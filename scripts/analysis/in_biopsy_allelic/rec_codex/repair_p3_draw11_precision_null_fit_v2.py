#!/usr/bin/env python3
"""V2 synthetic null pilot: independent precision budget for each original start.

Only the final constrained-polish proposals/acceptance are revised. The original
two starts, tie, bounds, optimizer and strict gradient cutoff remain unchanged.
This runner does not fit positive variance, replace a family member or calibrate.
"""
import argparse
from dataclasses import replace
from decimal import Decimal,localcontext
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
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT=REC/'repair_p3_gamma_draw_11_22081580'
REFERENCE=REC/'repair_p3_draw11_decimal_gradient_22087909/executed_checker.py'
REFERENCE_SHA='8e7e7914db125007d7fe99d81144975f17d33d9d7f98c39bb71e61a7795819fd'
WARM=REC/'repair_p3_gamma_observed_native_22080548/observed_null.json'


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def metadata():
    require(sha(REFERENCE)==REFERENCE_SHA,'Independent gradient source changed')
    spec=importlib.util.spec_from_file_location('draw11_precision_reference',REFERENCE)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    base,pins=module.metadata();pins[str(REFERENCE)]=REFERENCE_SHA
    require(str(WARM) in pins and sha(WARM)==pins[str(WARM)],'Original observed warm start missing')
    required=[POINT/'draw_counts.json',WARM,*[POINT/'executed_sources'/(n+'.py')
        for n in ('likelihood_v1','design_v1','bfix_v1','tests_v1')]]
    require(all(str(p) in pins and sha(p)==pins[str(p)] for p in required),'Required source/count pin missing')
    return module,base,pins


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();reference,base,pins=metadata()
    if args.inspect:
        print(json.dumps({'pins_verified':len(pins),'numeric_imports':False,'two_starts_unchanged':True,
            'only_constrained_polish_revised':True,'positive_variance_not_fitted':True}));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS')),'One numerical thread required')
    out=REC/('repair_p3_draw11_precision_null_fit_v2_'+os.environ['SLURM_JOB_ID']);out.mkdir()
    (out/'executed_runner.py').write_bytes(Path(__file__).read_bytes());began=time.perf_counter()
    recording_errors=[]
    try:
        import numpy as np
        import scipy
        from scipy.optimize import minimize
        require((platform.python_version(),np.__version__,scipy.__version__)==('3.10.18','2.2.6','1.15.2'),
            'Recorded environment required')
        def safe(x):
            if isinstance(x,np.ndarray):return safe(x.tolist())
            if isinstance(x,np.generic):return safe(x.item())
            if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
            if isinstance(x,dict):return {str(k):safe(v) for k,v in x.items()}
            if isinstance(x,(tuple,list)):return [safe(v) for v in x]
            return x
        def save(name,x):
            body=json.dumps(safe(x),indent=2,allow_nan=False)+'\n'
            require(len(body)<8*1024**2,'Bounded diagnostic output exceeded')
            (out/(name+'.json')).write_text(body)
        archive=POINT/'executed_sources';sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve()==(archive/(name+'.py')).resolve(),'Import escaped archive')
        design,settings=loader.load();template=design.model();counts=json.loads((POINT/'draw_counts.json').read_text())
        require(len(template.rows)==len(counts['counts'])==173 and len(design.participants)==60,'Dimensions differ')
        require(counts['full_row_keys']==[[x.individual,settings['gene_order'][x.gene]] for x in template.rows],'Row order differs')
        require(hashlib.sha256(np.asarray(counts['counts'],dtype='<i8').tobytes()).hexdigest()==counts['count_vector_sha256'],
            'Synthetic counts differ')
        rows=[replace(row,a=int(a)) for row,a in zip(template.rows,counts['counts'])]
        private=sys.modules['likelihood_v1'];cache={};precision_receipts=[];history=[];starts=[]
        context={'start_index':None,'bounds':None};precision_counts=[]
        def high_precision(p):
            key=np.asarray(p,dtype='<f8').tobytes()
            if key not in cache:
                require(len(cache)<64,'Per-start precision budget exhausted; preserve incomplete fit')
                values={digits:base.likelihood(rows,list(map(float,p)),digits)['loglik'] for digits in (80,110)}
                with localcontext() as ctx:
                    ctx.prec=130
                    require(abs(Decimal(values[80])-Decimal(values[110]))<Decimal('1e-50'),
                        'Independent precision refinement disagrees')
                cache[key]=values;precision_receipts.append({'start_index':context['start_index'],
                    'parameters':list(map(float,p)),'logliks':values})
                save('precision_receipts',precision_receipts)
            return cache[key]
        class PrecisionNullModel(private.Model):
            def constrained_polish(self,start,bounds,maxiter=100):
                require(context['bounds']==bounds,'Passed polish bounds differ from original single-fit bounds')
                p=np.asarray(start,float).copy();local_history=[]
                for iteration in range(maxiter):
                    ll,g,h=self.evaluate(p)
                    active=np.array([(lo is not None and lo==hi) or
                        (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or
                        (hi is not None and p[j]>=hi-1e-8 and g[j]>0)
                        for j,(lo,hi) in enumerate(bounds)])
                    free=np.flatnonzero(~active);maximum=float(np.max(np.abs(g[free]),initial=0))
                    record={'start_index':context['start_index'],'iteration':iteration,
                        'parameters':p.copy(),'bounds':bounds,'loglik':ll,
                        'max_projected_gradient':maximum,'proposals':[]}
                    local_history.append(record)
                    history.append(record);save('polish_history',history)
                    if maximum<1e-6:break
                    require(np.isfinite(ll) and np.isfinite(h).all() and np.isfinite(g).all(),
                        'Nonfinite polish likelihood/derivatives')
                    curvature=-h[np.ix_(free,free)];values,vectors=np.linalg.eigh((curvature+curvature.T)/2)
                    floor=max(1e-9,float(np.max(np.abs(values)))*1e-10)
                    regularized=(vectors*np.maximum(values,floor))@vectors.T
                    limits=[(None if bounds[j][0] is None else bounds[j][0]-p[j],
                        None if bounds[j][1] is None else bounds[j][1]-p[j]) for j in free]
                    def quadratic(s):return .5*s@regularized@s-g[free]@s,regularized@s-g[free]
                    bounded=minimize(quadratic,np.zeros(len(free)),jac=True,method='SLSQP',bounds=limits,
                        options={'ftol':1e-14,'maxiter':1000})
                    directions=(('lstsq',np.linalg.lstsq(curvature,g[free],rcond=1e-10)[0]),
                        ('original_regularized_quadratic',bounded.x))
                    before=high_precision(p);accepted=False
                    for label,direction in directions:
                        if accepted:break
                        step=np.zeros(self.size);step[free]=direction
                        for damping in 2.**-np.arange(15):
                            trial=p+damping*step
                            for j,(lo,hi) in enumerate(bounds):trial[j]=np.clip(trial[j],
                                -np.inf if lo is None else lo,np.inf if hi is None else hi)
                            after=high_precision(trial)
                            with localcontext() as ctx:
                                ctx.prec=130
                                changes={d:Decimal(after[d])-Decimal(before[d]) for d in (80,110)}
                                require(abs(changes[80]-changes[110])<Decimal('1e-50'),'Precision gain disagrees')
                                accepted=all(x>0 for x in changes.values())
                            record['proposals'].append({'direction':label,'damping':float(damping),
                                'parameters':trial.copy(),'independent_LL_changes':{d:str(v) for d,v in changes.items()},
                                'accepted':accepted})
                            save('polish_history',history)
                            if accepted:p=trial;break
                    if not accepted:break
                return p,local_history
        for name in ('fit','fit_single','bounds','initial','cold_start','informative','predictors','evaluate'):
            require(getattr(PrecisionNullModel,name) is getattr(private.Model,name),'Unexpected numeric method revision')
        model=PrecisionNullModel(rows,gene_count=8);require((model.G,model.Z,model.size)==(8,7,28),'Coordinates differ')
        original_single=model.fit_single
        def recorded_single(start):
            nonlocal cache
            cache={};context['start_index']=len(starts)
            start=np.asarray(start,float)
            context['bounds']=[b if model.informative()[j] else (start[j],start[j])
                for j,b in enumerate(model.bounds())]
            record={'start':start.copy(),'bounds_used_by_original_fit_single':context['bounds']}
            try:result=original_single(start)
            except BaseException as exc:
                record.update(status='raised',exception_type=type(exc).__name__,message=str(exc),
                    traceback=traceback.format_exc());starts.append(record)
                try:save('start_receipts',starts)
                except Exception as error:recording_errors.append({'type':type(error).__name__,'message':str(error)})
                raise
            precision_counts.append(len(cache));record.update(status='returned',terminal=result)
            starts.append(record);save('start_receipts',starts)
            return result
        model.fit_single=recorded_single
        warm=json.loads(WARM.read_text())['fit']['parameters'];fit=model.fit(np.asarray(warm,dtype=float))
        save('null_fit',{'fit':fit,'only_constrained_polish_numeric_revision':True,
            'fit_single_instance_instrumented':True,'original_fit_selection_method_retained':True})
        independent=[]
        for index,receipt in enumerate(starts):
            terminal=receipt['terminal'];p=terminal['parameters'].tolist()
            checks=[]
            for digits in (80,110):
                gref=reference.gradient(base,rows,p,digits);gs=list(map(Decimal,gref['gradient']))
                bounds=receipt['bounds_used_by_original_fit_single']
                require(all((lo is None or p[j]>=lo-1e-10) and (hi is None or p[j]<=hi+1e-10)
                    for j,(lo,hi) in enumerate(bounds)),'Terminal point escaped original bounds')
                with localcontext() as ctx:
                    ctx.prec=130
                    projected=[Decimal(0) if ((lo is not None and lo==hi) or
                        (lo is not None and p[j]<=lo+1e-8 and gs[j]<0) or
                        (hi is not None and p[j]>=hi-1e-8 and gs[j]>0)) else gs[j]
                        for j,(lo,hi) in enumerate(bounds)]
                    checks.append({'digits':digits,'loglik':gref['loglik'],'gradient':gref['gradient'],
                        'bounds':bounds,'projected_gradient':list(map(str,projected)),
                        'max_projected_gradient':str(max(map(abs,projected))),
                        'strict_original_cutoff_met':max(map(abs,projected))<Decimal('0.000001')})
            independent.append({'start_index':index,'terminal_parameters':p,'checks':checks})
            save('independent_start_checks',independent)
        require(len(starts)==2,'Original two-start recipe differs')
        for absolute,digest in pins.items():require(sha(Path(absolute))==digest,'Input changed during pilot')
        save('summary',{'source_and_input_sha256':pins,'runner_sha256':sha(Path(__file__)),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'elapsed_seconds':time.perf_counter()-began,'native_converged':fit['converged'],
            'start_converged':fit['start_converged'],'kept_start':fit['kept_start'],
            'all_independent_start_scores_pass':all(c['strict_original_cutoff_met'] for x in independent for c in x['checks']),
            'qualified_two_start_null_pilot':bool(fit['converged'] and all(
                starts[i]['terminal']['converged'] and all(c['strict_original_cutoff_met'] for c in x['checks'])
                for i,x in enumerate(independent))),
            'precision_points_per_start':precision_counts,'independent_budget_per_start':True,
            'actual_single_fit_bounds_recorded':True,'original_family_member_not_replaced':True,
            'section4_not_fitted':True,'complete_fifty_draw_agreement':False,'calibration':False})
    except BaseException:
        (out/'failure.json').write_text(json.dumps({'traceback':traceback.format_exc(),
            'recording_errors':recording_errors,
            'original_family_member_not_replaced':True},indent=2)+'\n');raise


if __name__=='__main__':main()
