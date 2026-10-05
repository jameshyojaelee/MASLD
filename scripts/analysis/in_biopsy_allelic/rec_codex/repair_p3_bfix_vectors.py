#!/usr/bin/env python3
"""All-gene shared synthetic fixed-vector baselines; no fits/calibration."""
import argparse
import ast
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
HERE=Path(__file__).resolve().parent
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
CAMPAIGN=REC/'p3_parallel_20261003T130323Z'
OLD=REC/'bfix_fixed_vectors_v12.json'
THETAS=('theta1_start_like','theta2_truth_v0','theta3_truth_vpos','theta4_perturbed')
MODULES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py')
CAP=32*1024**2


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--theta',choices=THETAS)
    parser.add_argument('--gene',type=int,choices=range(8))
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-kernel-sha',required=True)
    parser.add_argument('--expected-gl-sha',required=True)
    args=parser.parse_args()
    require((args.theta is None)==(args.gene is None),'Theta and gene must be paired')
    requested=[(theta,gene) for theta in THETAS for gene in range(8)
               if args.theta is None or (theta==args.theta and gene==args.gene)]
    pins={FIXTURE/'SHA256SUMS':'0d3d06fc32a8e588d5e1c54c943a582e5054b6d2ab565290a4d57ced545d2925',
          CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
          OLD:'25270e3294800c7f49a5124b8b4982b58a1c6df8338f36f5a05a4e72c56e0f98',
          HERE/'repair_p3_complex_reference.py':'f722eb76a2fb288f4f4e697dcd149d714b89221ec791ba6a43db8cfebb2f7ec3',
          Path(__file__):args.expected_checker_sha,HERE/'repair_p3_integral.py':args.expected_kernel_sha,
          HERE/'repair_p3_composite_gl.py':args.expected_gl_sha}
    for path,digest in pins.items():require(sha(path)==digest,'Reviewed source/metadata differs: '+path.name)
    campaign=json.loads((CAMPAIGN/'campaign.json').read_text())
    for name in MODULES:
        pins[HERE/name]=campaign['files'][str(CAMPAIGN/'executed_sources'/name)]
        require(sha(HERE/name)==pins[HERE/name],'Original module differs')
    published={name.removeprefix('./'):digest for digest,name in
               (line.split(maxsplit=1) for line in (FIXTURE/'SHA256SUMS').read_text().splitlines())}
    inputs=('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','genes.tsv','tags.tsv')+tuple('params/'+t+'.json' for t in THETAS)
    for name in inputs:
        pins[FIXTURE/name]=published[name]
        require(sha(FIXTURE/name)==published[name],'Original synthetic fixture differs')
    for path in pins:
        if path.suffix=='.py':ast.parse(path.read_text())
    protocol=dict(requested_points=requested,requested_point_count=len(requested),required_combined_points=32,
                  all32_remain_required=True,participants=60,genes=8,rows=173,nuisance_coordinates=28,
                  positive_variance_coordinates=29,no_model_fits=True,no_calibration=True,
                  independent_gradient_accuracy_verified=False,independent_Hessian_accuracy_verified=False,
                  source_and_input_sha256={str(p):h for p,h in pins.items()},output_cap_bytes=CAP)
    if args.inspect:
        print(json.dumps(protocol,sort_keys=True));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output=REC/('repair_p3_bfix_vectors_'+os.environ['SLURM_JOB_ID']);output.mkdir(exist_ok=False)
    archive=output/'executed_sources';archive.mkdir()
    for source,digest in pins.items():
        if source.suffix=='.py':
            target=archive/source.name;target.write_bytes(source.read_bytes())
            require(sha(target)==digest,'Archived source differs')
    sys.path.insert(0,str(archive))
    require(not any(n[:-3] in sys.modules for n in MODULES),'Original numerical module already imported')
    import numpy as np
    import scipy
    from scipy.stats import norm
    original=importlib.import_module('likelihood_v1')
    loader=importlib.import_module('bfix_v1')
    importlib.import_module('tests_v1')
    for name in MODULES:
        module=sys.modules[name[:-3]]
        require(Path(module.__file__).resolve().parent==archive.resolve() and sha(module.__file__)==pins[HERE/name],
                'Import escaped pinned archive')
    versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    require(versions==dict(python='3.10.18',numpy='2.2.6',scipy='1.15.2'),'Reviewed numerical environment differs')
    protocol['environment']=versions
    private=load('bfix_private_original',archive/'likelihood_v1.py')
    kernel=load('bfix_private_kernel',archive/'repair_p3_integral.py')
    gl=load('bfix_private_gl',archive/'repair_p3_composite_gl.py')
    reference=load('bfix_independent_probability',archive/'repair_p3_complex_reference.py')
    private.bb_terms=kernel.stable_bb_terms
    Candidate=kernel.candidate_class(private)

    def serial(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)

    def diagnostic(value):
        if isinstance(value,np.ndarray):return diagnostic(value.tolist())
        if isinstance(value,np.generic):return diagnostic(value.item())
        if isinstance(value,dict):return {k:diagnostic(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)):return [diagnostic(v) for v in value]
        if isinstance(value,float) and not np.isfinite(value):
            return {'nonfinite_float':'nan' if np.isnan(value) else 'positive_infinity' if value>0 else 'negative_infinity'}
        return value

    def save(name,value):
        body=json.dumps(value,default=serial,allow_nan=False,indent=2)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())<=CAP,
                'Synthetic output cap exhausted')
        with (output/(name+'.json')).open('x') as handle:handle.write(body)

    save('protocol',protocol)
    design,settings=loader.load();model=design.model()
    require(len(design.participants)==60 and len(model.rows)==173 and model.G==8 and model.Z==7 and model.size==28,
            'Shared synthetic dimensions differ')
    old=json.loads(OLD.read_text())['vectors']
    records=[]
    for theta,gene in requested:
        p,v=loader.parameters(FIXTURE/'params'/(theta+'.json'))
        require(p.shape==(28,) and np.isfinite(p).all() and v>=0 and np.isfinite(v),'Fixed vector differs')
        rows=[r for r in model.rows if r.gene==gene]
        candidate=Candidate(rows,gene_count=8)
        record=dict(theta=theta,gene=gene,rows=len(rows),p=p,v=v,
                    independent_gradient_accuracy_verified=False,independent_Hessian_accuracy_verified=False)
        try:
            if v>0:
                ll,g,h,report=gl.integrate_gene(candidate,p,v,rows,private.conditional_terms,with_hessian=True,
                    expected_kernel_sha=args.expected_kernel_sha,kernel_path=archive/'repair_p3_integral.py')
            else:
                ll,g,h=candidate.marginal(p,0.,with_hessian=True)
                # Original boundary Hessian deliberately has undefined variance
                # row/column. Retain only the exact nuisance block, never fill 0.
                h=h[:-1,:-1]
                report=dict(exact_variance_zero_branch=True,positive_variance_Hessian_not_applicable=True)
            require(np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),'Candidate result nonfinite')
            record['candidate']=dict(logmass=ll,gradient=g,hessian=h,integration=diagnostic(report))
            record['historical_original_per_gene_logmass']=old[theta]['per_gene_loglik'][str(gene)]
            record['candidate_minus_historical_original']=ll-record['historical_original_per_gene_logmass']
            record['historical_difference_descriptive_not_revised_approximation_agreement']=True
        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as failure:
            record['candidate_refusal']=dict(traceback=traceback.format_exc(),
                    partial_reports=diagnostic(getattr(failure,'integration_reports',None)))
        try:
            # Each gene has its own independent random effect. Never integrate
            # all genes jointly using a single shared u.
            ref=reference.Reference(rows,8,7)
            if v==0:
                z=ref.loglik(p,0.,ref.prepare(p))
                require(np.isfinite(z) and abs(z.imag)<1e-12,'Independent boundary probability invalid')
                scalar=float(z.real)
                record['independent_scalar']=dict(logmass=scalar,exact_variance_zero_probability=True,
                    arithmetic_error_unbounded=True)
            else:
                first=None
                for bound in range(24,95,2):
                    geometry=ref.geometry(p,v,float(bound),.5)
                    result=ref.mass(p,v,geometry)
                    ratio=float(np.exp(np.log(2.)+norm.logsf(bound)-result['logmass_lower_estimate']))
                    if ratio<=1e-10:
                        first=result;first['relative_prior_mass_tail_allowance']=ratio;break
                require(first is not None,'Independent scalar relative-tail domain cap exhausted')
                geometry=ref.geometry(p,v,float(bound+2),.25)
                expanded=ref.mass(p,v,geometry)
                ratio=float(np.exp(np.log(2.)+norm.logsf(bound+2)-expanded['logmass_lower_estimate']))
                delta=expanded['logmass']-first['logmass']
                error=-np.log1p(-expanded['estimated_error']/expanded['value'])+np.log1p(ratio)
                require(np.isfinite(expanded['logmass']) and abs(delta)+error<=1e-9,
                        'Independent expanded scalar refinement unresolved')
                scalar=float(expanded['logmass'])
                record['independent_scalar']=dict(first=first,expanded=expanded,logmass=scalar,
                    refinement_difference=delta,estimated_quad_plus_tail_log_error=float(error),
                    expanded_relative_prior_mass_tail_allowance=ratio,certified_error_bound=False)
            if 'candidate' in record:
                difference=record['candidate']['logmass']-scalar
                record['candidate_scalar_logmass_difference']=difference
                record['candidate_scalar_per_gene_budget_pass']=abs(difference)<=1e-6/8
                record['per_gene_budget_not_full_vector_agreement']=True
        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as failure:
            record['independent_scalar_refusal']=dict(traceback=traceback.format_exc(),
                     partial_reference_record=diagnostic(getattr(failure,'reference_record',None)))
        save(theta+'_gene'+str(gene),record);records.append(record)
    for path,digest in pins.items():require(sha(path)==digest,'Pinned input/source changed')
    totals={}
    for theta in THETAS:
        points=[r for r in records if r['theta']==theta]
        if len(points)==8 and all('candidate' in r and 'independent_scalar' in r for r in points):
            ll=sum(r['candidate']['logmass'] for r in points)
            scalar=sum(r['independent_scalar']['logmass'] for r in points)
            totals[theta]=dict(candidate_logmass=ll,independent_scalar_logmass=scalar,
                gradient=sum((np.asarray(r['candidate']['gradient']) for r in points),np.zeros(29)),
                hessian=sum((np.asarray(r['candidate']['hessian']) for r in points),
                            np.zeros_like(np.asarray(points[0]['candidate']['hessian']))),
                hessian_coordinate_scope='nuisance28_only' if points[0]['v']==0 else 'all29_positive_variance',
                scalar_logmass_difference=ll-scalar,scalar_likelihood_only_pass=abs(ll-scalar)<=1e-6)
    save('summary',dict(all32points_recorded=len(records)==32,all_requested_points_recorded=len(records)==len(requested),
         required_combined_points=32,requested_point_count=len(requested),all32_remain_required=True,vector_totals=totals,
         independent_gradient_accuracy_verified=False,independent_Hessian_accuracy_verified=False,
         frozen20_node_approximation_agreement_claimed=False,no_model_fits=True,no_calibration=True))


if __name__=='__main__':main()
