#!/usr/bin/env python3
"""Independent finite-product Decimal scores at the fixed draw11 point pair."""
import argparse
from dataclasses import replace
from decimal import Decimal,localcontext
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT=REC/'repair_p3_gamma_draw_11_22081580'
PROBE=REC/'repair_p3_draw11_stationarity_probe_22087829'
DECIMAL=REC/'repair_p3_draw11_decimal_likelihood_22087880/executed_checker.py'
DECIMAL_SHA='38b3539bf30798b1ab34d738a45936ed00c62dc88e47170ac4644b1ea9a7b99b'


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def metadata():
    require(sha(DECIMAL)==DECIMAL_SHA,'Independent Decimal source changed')
    spec=importlib.util.spec_from_file_location('draw11_decimal_probability_reference',DECIMAL)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    pins=module.metadata();pins[str(DECIMAL)]=DECIMAL_SHA
    return module,pins


def gradient(base,rows,parameters,digits):
    with localcontext() as ctx:
        ctx.prec=digits;p=list(map(base.dec,parameters));G,Z=8,7;eps=Decimal('0.000001')
        total=[Decimal(0)]*28;row_records=[];ll=Decimal(0);hom_cache={}
        for row in rows:
            stage=base.dec(row.stage);orientation=Decimal(row.orientation)
            zs=list(map(base.dec,row.z))
            factor=1+p[16]*stage+sum((p[17+j]*z for j,z in enumerate(zs)),Decimal(0))
            eta=p[row.gene]*factor-orientation*p[24]*stage-orientation*base.dec(row.omega)
            r=p[G+row.gene]+sum((p[25+j]*base.dec(z) for j,z in enumerate(
                (row.stage,row.background,row.duplication))),Decimal(0))
            require(abs(eta)<200 and abs(r)<1000 and 0<=row.n<=5000,'Synthetic scope exceeded')
            mean=base.logistic(eta);q=1-mean;vr=base.logistic(r)
            rho=eps+(1-2*eps)*vr;t=(1-rho)/rho
            tp=-(1-2*eps)*vr*(1-vr)/(rho*rho)
            a,b=mean*t,q*t
            het,normalization=base.bb_support(row.n,mean,rho)
            sa=[Decimal(0)];sb=[Decimal(0)];st=Decimal(0)
            for j in range(row.n):
                sa.append(sa[-1]+1/(a+j));sb.append(sb[-1]+1/(b+j));st+=1/(t+j)
            # Derivatives of the finite rising products; no candidate scores,
            # digamma, gamma functions, finite differences or complex steps.
            seta=[t*mean*q*(sa[k]-sb[row.n-k]) for k in range(row.n+1)]
            sr=[tp*(mean*sa[k]+q*sb[row.n-k]-st) for k in range(row.n+1)]
            key=(row.n,float(row.e),float(row.phi))
            if key not in hom_cache:
                ref,re=base.bb_support(row.n,base.dec(row.e),base.dec(row.phi))
                alt,ae=base.bb_support(row.n,1-base.dec(row.e),base.dec(row.phi))
                hom_cache[key]=(ref,alt)
                require(max(re,ae)<Decimal(10)**(-(digits-15)),'Homozygote normalization differs')
            require(normalization<Decimal(10)**(-(digits-15)),'Het normalization differs')
            ref,alt=hom_cache[key];pr,pa,px=map(base.dec,(row.p_ref,row.p_alt,row.p_x))
            weights=[2*pr*pa,pr*pr+2*pr*px,pa*pa+2*pa*px]
            if row.orientation==-1:weights[1],weights[2]=weights[2],weights[1]
            lo=max(row.minor,(row.percent*row.n+99)//100);hi=row.n-lo
            require(row.n>=row.depth and lo<=row.a<=hi,'Integer selection differs')
            numerator=sum((w*m[row.a] for w,m in zip(weights,(het,ref,alt))),Decimal(0))
            denominator=sum((w*sum(m[lo:hi+1],Decimal(0)) for w,m in zip(weights,(het,ref,alt))),Decimal(0))
            row_ll=(numerator/denominator).ln();ll+=row_ll
            scores=[]
            for s in (seta,sr):
                scores.append(weights[0]*het[row.a]*s[row.a]/numerator-
                    weights[0]*sum((het[k]*s[k] for k in range(lo,hi+1)),Decimal(0))/denominator)
            ge,gr=scores;contribution=[Decimal(0)]*28
            contribution[row.gene]=factor*ge;contribution[G+row.gene]=gr
            contribution[16]=p[row.gene]*stage*ge
            for j,z in enumerate(zs):contribution[17+j]=p[row.gene]*z*ge
            contribution[24]=-orientation*stage*ge
            for j,z in enumerate((row.stage,row.background,row.duplication)):contribution[25+j]=base.dec(z)*gr
            total=[x+y for x,y in zip(total,contribution)]
            row_records.append({'gene':row.gene,'individual':row.individual,'loglik':str(row_ll),
                'gradient_contribution':list(map(str,contribution))})
        return {'digits':digits,'loglik':str(ll),'gradient':list(map(str,total)),'rows':row_records}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args();base,pins=metadata()
    if args.inspect:
        print(json.dumps({'pins_verified':len(pins),'numeric_imports':False,'no_fits':True}));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS')),'One numerical thread required')
    out=REC/('repair_p3_draw11_decimal_gradient_'+os.environ['SLURM_JOB_ID']);out.mkdir()
    (out/'executed_checker.py').write_bytes(Path(__file__).read_bytes())
    try:
        import numpy as np
        import scipy
        require((platform.python_version(),np.__version__,scipy.__version__)==('3.10.18','2.2.6','1.15.2'),
            'Recorded environment required')
        archive=POINT/'executed_sources';sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve()==(archive/(name+'.py')).resolve(),'Import escaped archive')
        design,settings=loader.load();template=design.model();counts=json.loads((POINT/'draw_counts.json').read_text())
        require(len(template.rows)==len(counts['counts'])==173 and len(design.participants)==60,'Dimensions differ')
        require(counts['full_row_keys']==[[x.individual,settings['gene_order'][x.gene]] for x in template.rows],'Row order differs')
        require(hashlib.sha256(np.asarray(counts['counts'],dtype='<i8').tobytes()).hexdigest()==counts['count_vector_sha256'],
            'Count hash differs')
        rows=[replace(row,a=int(a)) for row,a in zip(template.rows,counts['counts'])]
        model=sys.modules['likelihood_v1'].Model(rows,gene_count=8)
        require((model.G,model.Z,model.size)==(8,7,28),'Parameter coordinates differ')
        original=json.loads((POINT/'p3_011.json').read_text())['fit']['null']
        trial=json.loads((PROBE/'probe_01.json').read_text());records=[]
        for digits in (80,110):
            for label,saved in (('baseline',original),('full_original_Newton_probe',trial)):
                require(len(saved['parameters'])==len(saved['gradient'])==28,'Saved coordinate count differs')
                p=saved['parameters'];record={'point':label,**gradient(base,rows,p,digits)}
                bounds=[b if model.informative()[j] else (p[j],p[j]) for j,b in enumerate(model.bounds())]
                with localcontext() as ctx:
                    ctx.prec=digits;gs=list(map(Decimal,record['gradient']));projected=[]
                    for j,(lo,hi) in enumerate(bounds):
                        active=((lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-10 and gs[j]<0)
                            or (hi is not None and p[j]>=hi-1e-10 and gs[j]>0))
                        projected.append(Decimal(0) if active else gs[j])
                    record['bounds']=bounds;record['parameters']=p
                    record['projected_gradient']=list(map(str,projected))
                    record['max_projected_gradient']=str(max(map(abs,projected)))
                    record['strict_original_cutoff_met']=max(map(abs,projected))<Decimal('0.000001')
                    record['saved_double_gradient_differences']=[str(x-base.dec(y)) for x,y in zip(gs,saved['gradient'])]
                records.append(record);body=json.dumps(record,indent=2)+'\n'
                require(len(body)<2*1024**2,'Record cap exceeded')
                (out/(label+'_'+str(digits)+'.json')).write_text(body)
        for absolute,digest in pins.items():require(sha(Path(absolute))==digest,'Input changed during diagnostic')
        (out/'summary.json').write_text(json.dumps({'source_and_input_sha256':pins,'checker_sha256':sha(Path(__file__)),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'records':[{k:v for k,v in x.items() if k!='rows'} for x in records],
            'no_fit_or_proposal_accepted':True,'independent_Hessian_not_evaluated':True,'calibration':False},indent=2)+'\n')
    except BaseException:
        (out/'failure.json').write_text(json.dumps({'traceback':traceback.format_exc(),
            'no_fit_or_proposal_accepted':True},indent=2)+'\n');raise


if __name__=='__main__':main()
