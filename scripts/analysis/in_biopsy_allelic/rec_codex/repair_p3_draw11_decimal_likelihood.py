#!/usr/bin/env python3
"""Independent Decimal likelihood at two fixed synthetic null points; no fits."""
import argparse
from dataclasses import replace
from decimal import Decimal,localcontext
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT=REC/'repair_p3_gamma_draw_11_22081580'
PROBE=REC/'repair_p3_draw11_stationarity_probe_22087829'
EXTRA={'summary.json':'93c51252ce5213ee1e7dd2478ff3ea865e561d07f2ae7ba24f23d15e4f0c6aae',
 'probe_01.json':'6e6651887c73a7df7a8874de940624f32b38ed64307406f70841c6dc1fe671a3',
 'executed_checker.py':'a6458089995c4c0c02a1ec437e2dab8fa1215734d28b7709031eab34bbaeedb0'}


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def metadata():
    for name,digest in EXTRA.items():require(sha(PROBE/name)==digest,'Fixed probe changed')
    spec=importlib.util.spec_from_file_location('draw11_fixed_probe_metadata',PROBE/'executed_checker.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    pins=module.metadata();pins.update({str(PROBE/n):h for n,h in EXTRA.items()})
    return pins


def dec(x):return Decimal.from_float(float(x))


def logistic(x):return Decimal(1)/(Decimal(1)+(-x).exp())


def bb_support(n,mean,rho):
    """Direct probability products/recurrence, independent of gamma functions."""
    one=Decimal(1);require(0<mean<1 and 0<=rho<1,'Invalid Decimal BB inputs')
    if rho==0:
        value=(one-mean)**n;mass=[value]
        for k in range(n):
            value=value*Decimal(n-k)/Decimal(k+1)*mean/(one-mean);mass.append(value)
    else:
        t=(one-rho)/rho;a=mean*t;b=(one-mean)*t
        value=one
        for j in range(n):value*= (b+j)/(t+j)
        mass=[value]
        for k in range(n):
            value=value*Decimal(n-k)/Decimal(k+1)*(a+k)/(b+n-k-1);mass.append(value)
    require(all(x>=0 for x in mass),'Negative Decimal mass')
    error=abs(sum(mass,Decimal(0))-one)
    return mass,error


def likelihood(rows,parameters,digits):
    with localcontext() as ctx:
        ctx.prec=digits
        p=list(map(dec,parameters));G,Z=8,7;eps=Decimal('0.000001')
        total=Decimal(0);records=[];maximum=Decimal(0);hom_cache={}
        for row in rows:
            stage=dec(row.stage);orientation=Decimal(row.orientation)
            factor=1+p[2*G]*stage+sum((p[2*G+1+j]*dec(z) for j,z in enumerate(row.z)),Decimal(0))
            eta=p[row.gene]*factor-orientation*p[2*G+1+Z]*stage-orientation*dec(row.omega)
            r=p[G+row.gene]+sum((p[-3+j]*dec(z) for j,z in enumerate(
                (row.stage,row.background,row.duplication))),Decimal(0))
            require(abs(eta)<200 and abs(r)<1000 and 0<=row.n<=5000,'Fixed synthetic scope exceeded')
            rho=eps+(1-2*eps)*logistic(r)
            het,error=bb_support(row.n,logistic(eta),rho);maximum=max(maximum,error)
            key=(row.n,float(row.e),float(row.phi))
            if key not in hom_cache:
                ref,re=bb_support(row.n,dec(row.e),dec(row.phi))
                alt,ae=bb_support(row.n,1-dec(row.e),dec(row.phi))
                hom_cache[key]=(ref,alt);maximum=max(maximum,re,ae)
            ref,alt=hom_cache[key]
            pr,pa,px=map(dec,(row.p_ref,row.p_alt,row.p_x))
            weights=[2*pr*pa,pr*pr+2*pr*px,pa*pa+2*pa*px]
            if row.orientation==-1:weights[1],weights[2]=weights[2],weights[1]
            lo=max(row.minor,(row.percent*row.n+99)//100)
            require(row.n>=row.depth and lo<=row.a<=row.n-lo,'Integer selection predicate differs')
            masses=(het,ref,alt)
            numerator=sum((w*m[row.a] for w,m in zip(weights,masses)),Decimal(0))
            denominator=sum((w*sum(m[lo:row.n-lo+1],Decimal(0)) for w,m in zip(weights,masses)),Decimal(0))
            ll=(numerator/denominator).ln();total+=ll
            records.append(str(ll))
        require(maximum<Decimal(10)**(-(digits-15)),'Probability normalization precision check failed')
        return {'digits':digits,'loglik':str(total),'row_logliks':records,
            'maximum_support_normalization_error':str(maximum),'parameters':parameters,
            'predictors_evaluated_in_Decimal':True,'raw_three_class_weights_not_renormalized':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();pins=metadata()
    if args.inspect:
        print(json.dumps({'pins_verified':len(pins),'numeric_imports':False,'no_fits':True}));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS')),'One numerical thread required')
    out=REC/('repair_p3_draw11_decimal_likelihood_'+os.environ['SLURM_JOB_ID']);out.mkdir()
    (out/'executed_checker.py').write_bytes(Path(__file__).read_bytes())
    try:
        import numpy as np
        archive=POINT/'executed_sources';sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve()==(archive/(name+'.py')).resolve(),
                'Fixture import escaped archive')
        design,settings=loader.load();template=design.model()
        counts=json.loads((POINT/'draw_counts.json').read_text())
        require(len(template.rows)==len(counts['counts'])==173 and len(design.participants)==60,
            'Exact synthetic dimensions required')
        require(counts['full_row_keys']==[[x.individual,settings['gene_order'][x.gene]] for x in template.rows],
            'Row order differs')
        require(hashlib.sha256(np.asarray(counts['counts'],dtype='<i8').tobytes()).hexdigest()==
            counts['count_vector_sha256'],'Synthetic count hash differs')
        rows=[replace(row,a=int(a)) for row,a in zip(template.rows,counts['counts'])]
        baseline=json.loads((POINT/'p3_011.json').read_text())['fit']['null']['parameters']
        trial=json.loads((PROBE/'probe_01.json').read_text())['parameters']
        records=[]
        for digits in (80,110):
            for label,p in (('baseline',baseline),('full_original_Newton_probe',trial)):
                record={'point':label,**likelihood(rows,p,digits)};records.append(record)
                (out/(label+'_'+str(digits)+'.json')).write_text(json.dumps(record,indent=2)+'\n')
        with localcontext() as ctx:
            ctx.prec=130
            differences=[{'digits':digits,'trial_minus_baseline':str(
                Decimal(records[i+1]['loglik'])-Decimal(records[i]['loglik']))}
                for i,digits in ((0,80),(2,110))]
        for absolute,digest in pins.items():require(sha(Path(absolute))==digest,'Input changed during diagnostic')
        (out/'summary.json').write_text(json.dumps({'source_and_input_sha256':pins,
            'checker_sha256':sha(Path(__file__)),'differences':differences,
            'no_fit_or_proposal_accepted':True,'gradient_not_independently_evaluated':True,
            'precision_refinement_not_a_formal_error_bound':True,'calibration':False},indent=2)+'\n')
    except BaseException:
        (out/'failure.json').write_text(json.dumps({'traceback':traceback.format_exc(),
            'no_fit_or_proposal_accepted':True},indent=2)+'\n');raise


if __name__=='__main__':main()
