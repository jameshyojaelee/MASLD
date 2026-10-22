#!/usr/bin/env python3
"""Regenerate supplied synthetic P3 counts by independent Decimal inverse CDFs."""
import argparse
from bisect import bisect_left
from decimal import Decimal,localcontext
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import struct
import sys
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT=REC/'repair_p3_gamma_draw_11_22081580'
REFERENCE=REC/'repair_p3_draw11_decimal_likelihood_22087880/executed_checker.py'
REFERENCE_SHA='38b3539bf30798b1ab34d738a45936ed00c62dc88e47170ac4644b1ea9a7b99b'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
NULL=REC/'p3_parallel_20261003T130323Z/baseline/observed_null.json'


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata(args):
    require(args.manifest.resolve().is_relative_to(REC.resolve()) and args.manifest.stat().st_size<128*1024,
        'Bounded synthetic manifest required')
    require(sha(args.manifest)==args.expected_manifest_sha and sha(REFERENCE)==REFERENCE_SHA,'Manifest/reference differs')
    spec=importlib.util.spec_from_file_location('Decimal_P3_probability_reference',REFERENCE)
    base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
    pins=base.metadata();pins[str(REFERENCE)]=REFERENCE_SHA;pins[str(args.manifest)]=args.expected_manifest_sha
    manifest=json.loads(args.manifest.read_text());records=manifest['records']
    require(sorted(x['draw'] for x in records)==list(range(50)),'All50 unique supplied draws required')
    for record in records:
        path=Path(record['counts_path'])
        require(path==REC/('repair_p3_gamma_draw_%d_%s'%(record['draw'],record['job_id']))/'draw_counts.json',
            'Unexpected synthetic count path')
        require(path.stat().st_size<64*1024 and sha(path)==record['sha256'],'Count receipt differs')
        pins[str(path)]=record['sha256']
    require(str(NULL) in pins and str(FIXTURE/'draws/p3_uniforms.npy') in pins,'Generating null/uniform pins missing')
    return base,pins,records


def cdfs(base,rows,parameters,digits):
    with localcontext() as ctx:
        ctx.prec=digits;p=list(map(base.dec,parameters));tables=[];maximum=Decimal(0)
        for row in rows:
            stage=base.dec(row.stage);orientation=Decimal(row.orientation)
            factor=1+p[16]*stage+sum((p[17+j]*base.dec(z) for j,z in enumerate(row.z)),Decimal(0))
            eta=p[row.gene]*factor-orientation*p[24]*stage-orientation*base.dec(row.omega)
            r=p[8+row.gene]+sum((p[25+j]*base.dec(z) for j,z in enumerate(
                (row.stage,row.background,row.duplication))),Decimal(0))
            require(abs(eta)<200 and abs(r)<1000 and 0<=row.n<=5000,'Synthetic probability scope exceeded')
            rho=Decimal('0.000001')+(1-2*Decimal('0.000001'))*base.logistic(r)
            masses=[]
            for mean,dispersion in ((base.logistic(eta),rho),(base.dec(row.e),base.dec(row.phi)),
                    (1-base.dec(row.e),base.dec(row.phi))):
                mass,error=base.bb_support(row.n,mean,dispersion);masses.append(mass);maximum=max(maximum,error)
            lo=max(row.minor,(row.percent*row.n+99)//100);hi=row.n-lo
            require(row.n>=row.depth and lo<=hi,'Integer selection predicate differs')
            pr,pa,px=map(base.dec,(row.p_ref,row.p_alt,row.p_x))
            weights=[2*pr*pa,pr*pr+2*pr*px,pa*pa+2*pa*px]
            if row.orientation==-1:weights[1],weights[2]=weights[2],weights[1]
            inclusions=[sum(m[lo:hi+1],Decimal(0)) for m in masses]
            den=sum((w*c for w,c in zip(weights,inclusions)),Decimal(0));require(den>0,'No selected probability')
            acc=Decimal(0);classes=[]
            for weight,inclusion in zip(weights,inclusions):
                acc+=weight*inclusion/den;classes.append(acc)
            conditional=[]
            for mass,inclusion in zip(masses,inclusions):
                require(inclusion>0,'Undefined selected class')
                acc=Decimal(0);cdf=[]
                for value in mass[lo:hi+1]:acc+=value/inclusion;cdf.append(acc)
                conditional.append(cdf)
            tables.append((lo,classes,conditional))
        require(maximum<Decimal(10)**(-(digits-15)),'Full-support normalization refinement failed')
        return tables,str(maximum)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--expected-manifest-sha',required=True)
    parser.add_argument('--inspect',action='store_true');args=parser.parse_args()
    base,pins,records=metadata(args)
    if args.inspect:
        print(json.dumps({'pins_verified':len(pins),'draws':len(records),'numeric_imports':False}));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
        'NUMEXPR_NUM_THREADS')),'One numerical thread required')
    out=REC/('repair_p3_decimal_count_regeneration_'+os.environ['SLURM_JOB_ID']);out.mkdir()
    (out/'executed_checker.py').write_bytes(Path(__file__).read_bytes())
    try:
        import numpy as np
        import scipy
        require((platform.python_version(),np.__version__,scipy.__version__)==('3.10.18','2.2.6','1.15.2'),
            'Recorded numerical environment required')
        archive=POINT/'executed_sources';sys.path.insert(0,str(archive))
        loader=importlib.import_module('bfix_v1');importlib.import_module('tests_v1')
        for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
            require(Path(sys.modules[name].__file__).resolve()==(archive/(name+'.py')).resolve(),'Fixture loader escaped archive')
        design,settings=loader.load();model=design.model();rows=model.rows
        require((model.G,model.Z,model.size,len(rows),len(design.participants))==(8,7,28,173,60),
            'Synthetic dimensions differ')
        null=json.loads(NULL.read_text());require(null['converged'] is True,'Canonical generating null incomplete')
        uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy',mmap_mode='r',allow_pickle=False)
        require(uniforms.shape==(50,173,2) and np.isfinite(uniforms).all() and
            np.all((uniforms>=0)&(uniforms<1)),'Supplied uniform stream differs')
        expected={}
        for record in records:
            raw=json.loads(Path(record['counts_path']).read_text());b=record['draw']
            require(raw['full_row_keys']==[[x.individual,settings['gene_order'][x.gene]] for x in rows],
                'Saved row order differs')
            require(len(raw['counts'])==173 and all(isinstance(x,int) for x in raw['counts']),'Integer count shape differs')
            digest=hashlib.sha256(struct.pack('<173q',*raw['counts'])).hexdigest()
            require(digest==raw['count_vector_sha256']==raw['original_baseline_count_vector_sha256'],
                'Saved native/original count hashes differ')
            expected[b]=raw['counts']
        results=[]
        for digits in (80,110):
            tables,normalization=cdfs(base,rows,null['parameters'],digits);draws=[]
            for b in range(50):
                values=[]
                for (lo,classes,conditional),(uclass,ucount) in zip(tables,uniforms[b]):
                    genotype=min(2,bisect_left(classes,base.dec(uclass)))
                    selected=conditional[genotype]
                    values.append(lo+min(len(selected)-1,bisect_left(selected,base.dec(ucount))))
                differences=[j for j,(a,c) in enumerate(zip(values,expected[b])) if a!=c]
                draws.append({'draw':b,'counts':values,'count_vector_sha256':hashlib.sha256(
                    struct.pack('<173q',*values)).hexdigest(),'differing_row_indices':differences})
            record={'digits':digits,'maximum_normalization_error':normalization,'draws':draws};results.append(record)
            (out/('regenerated_'+str(digits)+'.json')).write_text(json.dumps(record,indent=2)+'\n')
        for path,digest in pins.items():require(sha(Path(path))==digest,'Input changed during regeneration')
        require(all(a['counts']==b['counts'] for a,b in zip(results[0]['draws'],results[1]['draws'])),
            'Independent count streams differ between precisions')
        (out/'summary.json').write_text(json.dumps({'source_and_input_sha256':pins,'checker_sha256':sha(Path(__file__)),
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
            'draws_per_precision':50,'precision_levels':[80,110],'rows_per_draw':173,
            'mismatched_draws':{x['digits']:[d['draw'] for d in x['draws'] if d['differing_row_indices']] for x in results},
            'sampler_and_probabilities_independently_implemented':True,'shared_fixture_row_loader':True,
            'no_new_random_draws':True,'no_fit_result_changed':True,'complete_fit_agreement':False,'calibration':False},indent=2)+'\n')
    except BaseException:
        (out/'failure.json').write_text(json.dumps({'traceback':traceback.format_exc(),
            'no_fit_result_changed':True},indent=2)+'\n');raise


if __name__=='__main__':main()
