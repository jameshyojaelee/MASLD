#!/usr/bin/env python3
"""Full eight-gene fixed-point independent probability scores; compute only.

One baseline job saves all eight independent scalar masses and candidate scores.
Separate single-axis jobs reuse pinned real geometry/masses, never candidate
scores in the reference calculation. No fits, acceptance or calibration.
"""
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
import time
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
CAMPAIGN=REC/'p3_parallel_20261003T130323Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
COMPLEX=REC/'repair_p3_complex_launch_20261005T141156Z'
BATCHED=REC/'repair_p3_batched_gl_positive18_22040259'
CALLBACK_CHECK=REC/'repair_p3_batched_conditional_positive18_22040052'
SCALAR=REC/'repair_p3_scalar_prefix_positive18_22039119'
STEPS=(1e-10,1e-15,1e-20)
DOMAINS=((24.,.5),(26.,.25))
GENES=8
SIZE=28
CAP=32*1024**2
EVAL_CAP=100000
REFERENCE_SHA='f722eb76a2fb288f4f4e697dcd149d714b89221ec791ba6a43db8cfebb2f7ec3'
KERNEL_SHA='461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
SOURCE_PINS={
 COMPLEX/'repair_p3_complex_reference.py':REFERENCE_SHA,
 BATCHED/'executed_sources/repair_p3_integral.py':KERNEL_SHA,
 BATCHED/'executed_sources/repair_p3_batched_gl.py':'8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52',
 BATCHED/'executed_sources/repair_p3_batched_conditional.py':'184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76',
 BATCHED/'executed_sources/repair_p3_prepared_conditional_v2.py':'289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b',
}
METADATA_PINS={
 CAMPAIGN/'campaign.json':'bb3e82833476f67b0b11d2388e2dc00358eb987367b7a6a666d11938c6baca94',
 CAMPAIGN/'baseline/observed_null.json':'a617d331fd3fb93576e5bd8588581aaab8f77946936536401905490438660181',
 CAMPAIGN/'baseline/p3_018.json':'7b0161543b8d9ceeabdc4acbe62ed202b7e0cf2e325e0a82941149525a2a0f97',
 SCALAR/'protocol.json':'26f061f8e481a1b9d1e11f2a5dfc0984c82010d3b1dd9c8f270eb1418d701001',
 CALLBACK_CHECK/'protocol.json':'f5ca433498f1f1c5d9cd556a86210e76d533a736704273ba33871e5d8de34c3b',
 CALLBACK_CHECK/'summary.json':'b40d8b577dd35fe1ee814cce8f32fb03f7ad761c98aa0e22be5d0e9f95543aeb',
 CALLBACK_CHECK/'batched_conditional_comparison.json':'a5ed572af944fd9688da7990236b1ddd5baca128e33e858ceb565f66d9858007',
 BATCHED/'protocol.json':'4fff9a95d228ce2f12ab0ff08350df24d225cabda53044d00ab2a96582e43b56',
 BATCHED/'summary.json':'3daebbcc620794a84aaa75ccdcd28d5817d027e4f5fa3f9795d47dca40cb4982',
}


def require(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1<<20),b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    require(path.suffix=='.json' and path.stat().st_size<=8*1024**2,'Only bounded JSON metadata allowed')
    def invalid(value):
        raise ValueError('Nonfinite JSON constant '+value)
    return json.loads(path.read_text(),parse_constant=invalid)


def metadata():
    guards={**SOURCE_PINS,**METADATA_PINS}
    campaign=read_json(CAMPAIGN/'campaign.json')
    for path,digest in guards.items():
        require(sha(path)==digest,'Pinned source/metadata differs: '+str(path))
    for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
        path=CAMPAIGN/'executed_sources'/(name+'.py')
        guards[path]=campaign['files'][str(path)]
        require(sha(path)==guards[path],'Immutable original source differs')
        ast.parse(path.read_text())
    for source in SOURCE_PINS:
        ast.parse(source.read_text())
    summary=read_json(CALLBACK_CHECK/'summary.json')
    check=read_json(CALLBACK_CHECK/'batched_conditional_comparison.json')
    require(summary.get('batched_conditional_comparison_completed') is True and summary['case']=='positive18' and
        check['all_comparisons_within_allowance'] and check['outside_allowance_count']==0 and
        check['comparison_records_completed'] and check['changed_parameters_rejected'] and
        check['same_row_field_mutation_rejected'] and check['too_many_nodes_rejected'] and
        check['original_callables_unchanged'] and check['persistent_cache_layout_unchanged'],
        'Pinned callback checks are not passing')
    prior=read_json(BATCHED/'summary.json')
    require(prior['fixed_point_check_completed'] and prior['candidate_integral_returned'] and prior['case']=='positive18',
        'Reviewed batched integration was not complete')
    return campaign,guards


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=('baseline','axis'),required=True)
    parser.add_argument('--axis',type=int,choices=range(SIZE+1))
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--expected-baseline-sha')
    parser.add_argument('--expected-helper-sha')
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    require((args.mode=='axis')==(args.axis is not None),'Exactly one axis required in axis mode')
    require(args.mode=='axis' or (args.baseline is None and args.expected_baseline_sha is None),
        'Baseline producer does not consume a previous baseline')
    campaign,guards=metadata()
    if args.inspect:
        print(json.dumps(dict(mode=args.mode,axis=args.axis,required_combined_axes=list(range(29)),draw=18,
            full_rows=173,genes=8,parameters=28,z_columns=7,fixed_v=.00711885901713953,
            steps=STEPS,domains=DOMAINS,per_gene_axis_budget=EVAL_CAP,
            candidate_baseline_separate_job=True,shared_independent_mass_reused=True,
            proposed_header=dict(partition='cpu',account='nslab',qos='nslab',cpus=1,memory='4G',time='90:00:00'),
            estimates=dict(axis_gene_specific_hours='0.2–2',axis_shared_hours='1–12',baseline_hours='1–12',
                peak_RSS_GiB='less than1',not_measured_for_full_model=True),
            output_cap=CAP,numerical_execution=False,accuracy_accepted=False,fit=False,calibration=False)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required before numerical imports')
    require(args.expected_helper_sha==sha(Path(__file__)),'Reviewed helper changed')
    baseline=None
    if args.mode=='axis':
        require(args.baseline is not None and args.baseline.name=='shared_baseline.json' and
            args.baseline.parent.parent.resolve()==REC.resolve() and
            args.baseline.parent.name.startswith('repair_p3_full_gradient_baseline_'),
            'Axis requires an owned full-model baseline receipt')
        require(args.expected_baseline_sha is not None and sha(args.baseline)==args.expected_baseline_sha,
            'Shared baseline pin differs')
        baseline=read_json(args.baseline)
    suffix='baseline' if args.mode=='baseline' else 'axis%02d'%args.axis
    output=REC/('repair_p3_full_gradient_'+suffix+'_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    archive=output/'executed_sources';archive.mkdir()
    began=time.monotonic();phase='archive_before_import'
    def serial(value):
        if isinstance(value,np.ndarray):return value.tolist()
        if isinstance(value,np.generic):return value.item()
        raise TypeError(type(value).__name__)
    def save(name,value):
        body=json.dumps(value,indent=2,default=serial,allow_nan=False)+'\n'
        used=sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
        require(used+len(body.encode())<=CAP-65536,'Output cap reached; failure reserve retained')
        final=output/(name+'.json');partial=output/(name+'.json.partial')
        require(not final.exists(),'Result overwrite refused')
        with partial.open('x') as handle:
            handle.write(body);handle.flush();os.fsync(handle.fileno())
        partial.rename(final)
    def refusal(error,**fields):
        return dict(**fields,completed=False,exception_type=type(error).__name__,traceback=traceback.format_exc(),
            reference_record=getattr(error,'reference_record',None),
            integration_reports=getattr(error,'integration_reports',None),
            likelihood_evaluations=getattr(error,'likelihood_evaluations',None))
    try:
        sources={p:d for p,d in guards.items() if p.suffix=='.py'}
        sources[Path(__file__).resolve()]=args.expected_helper_sha
        require(len({p.name for p in sources})==len(sources),'Archive basename collision')
        for source,digest in sources.items():
            target=archive/source.name
            with target.open('xb') as handle:handle.write(source.read_bytes())
            require(sha(target)==digest,'Executed source archive differs')
        for name in ('globals.json','individuals.tsv','rows.tsv','cohorts.tsv','draws/p3_uniforms.npy'):
            path=FIXTURE/name
            require(path.resolve().is_relative_to(FIXTURE.resolve()),'Synthetic input escaped root')
            require(sha(path)==campaign['files'][str(path)],'Frozen synthetic fixture differs')
            guards[path]=campaign['files'][str(path)]
        phase='immutable_draw18_reconstruction'
        require(not any(n in sys.modules for n in ('bfix_v1','design_v1','likelihood_v1','tests_v1')),
            'Original numerical modules already imported')
        sys.path.insert(0,str(archive))
        global np
        import numpy as np
        import scipy
        from scipy.stats import norm
        loader=importlib.import_module('bfix_v1')
        original=importlib.import_module('likelihood_v1')
        importlib.import_module('tests_v1')
        for name in ('bfix_v1','design_v1','likelihood_v1','tests_v1'):
            path=Path(sys.modules[name].__file__)
            require(path.resolve().parent==archive.resolve() and sha(path)==guards[CAMPAIGN/'executed_sources'/path.name],
                'Original numerical import escaped archive')
        reference=load_module('full_gradient_independent_reference',archive/'repair_p3_complex_reference.py')
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        scalar_protocol=read_json(SCALAR/'protocol.json')
        require(environment==scalar_protocol['environment'] and reference.STEPS==STEPS,'Environment/step family differs')
        design,settings=loader.load();model=design.model()
        require(model.G==GENES and model.Z==7 and model.size==SIZE and len(model.rows)==173 and settings['B']==50,
            'Full synthetic model shape differs')
        null=read_json(CAMPAIGN/'baseline/observed_null.json')
        require(null['converged'],'Original generating null failed')
        null_p=np.asarray(null['parameters'])
        uniforms=np.load(FIXTURE/'draws/p3_uniforms.npy',mmap_mode='r',allow_pickle=False)
        require(uniforms.shape==(50,173,2),'Supplied uniform shape differs')
        require([(r['individual_id'],r['gene_id']) for r in loader.table('rows.tsv')]==
            [(r.individual,settings['gene_order'][r.gene]) for r in model.rows],'Uniform row join differs')
        rows=[]
        for row,(u1,u2) in zip(model.rows,uniforms[18]):
            eta,dispersion,_,_=model.predictors(row,null_p)
            rows.append(original.conditional_draw(row,eta,dispersion,u1,u2))
        model=original.Model(rows,gene_count=GENES)
        saved=read_json(CAMPAIGN/'baseline/p3_018.json')
        require(saved['b']==18 and saved['fit']['converged'],'Saved positive18 fit differs')
        p=np.asarray(saved['fit']['parameters']);v=float(saved['fit']['v'])
        count_sha=hashlib.sha256(np.asarray([r.a for r in rows],dtype='<i8').tobytes()).hexdigest()
        require(p.shape==(SIZE,) and np.isfinite(p).all() and v==.00711885901713953 and
            scalar_protocol['fixed_p']==p.tolist() and scalar_protocol['fixed_v']==v and
            scalar_protocol['count_vector_sha256']==count_sha,'Saved full-model point/counts differ')
        for x,(lo,hi) in zip(p,model.bounds()):
            require((lo is None or x>=lo) and (hi is None or x<=hi),'Saved parameter outside original bounds')
        identity=dict(case='positive18',draw=18,fixed_p=p.tolist(),fixed_v=v,count_vector_sha256=count_sha,
            rows=173,genes=GENES,parameter_count=SIZE,z_columns=7,gene_order=settings['gene_order'],environment=environment)
        source_declarations={str(k):d for k,d in sources.items()}
        save('protocol',dict(**identity,mode=args.mode,axis=args.axis,steps=STEPS,domains=DOMAINS,
            helper_sha256=args.expected_helper_sha,source_pins=source_declarations,
            metadata_pins={str(k):d for k,d in guards.items() if k.suffix!='.py'},
            per_gene_axis_cap=EVAL_CAP,required_combined_axes=list(range(29)),output_cap=CAP,
            fixed_point_scope_not_fit_accuracy=True,accuracy_accepted=False,no_fits=True,no_calibration=True))
        gene_rows={g:[r for r in rows if r.gene==g] for g in range(GENES)}
        require(all(gene_rows.values()),'Every fixture gene must be present')
        if args.mode=='baseline':
            phase='independent_all8_scalar_mass'
            masses={};mass_refusals=[]
            for gene,rs in gene_rows.items():
                engine=reference.Reference(rs,GENES,7,maximum_evaluations=EVAL_CAP)
                gene_mass={}
                for bound,spacing in DOMAINS:
                    try:
                        geometry=engine.geometry(p,v,bound,spacing)
                        mass=engine.mass(p,v,geometry)
                        tail=float(np.exp(np.log(2.)+norm.logsf(bound)-mass['logmass_lower_estimate']))
                        allowance=float(mass['estimated_error']/(mass['value']-mass['estimated_error'])+tail)
                        gene_mass[str(bound)]=dict(geometry=geometry,mass=mass,
                            logmass_quad_plus_prior_tail_estimate=allowance)
                        save('mass_gene%02d_B%02d'%(gene,bound),gene_mass[str(bound)])
                    except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                        item=refusal(error,gene=gene,bound=bound,probability_evaluations=engine.evaluations)
                        mass_refusals.append(item);save('mass_refusal_gene%02d_B%02d'%(gene,bound),item)
                masses[str(gene)]=dict(domains=gene_mass,shared_geometry_mass_counter=engine.evaluations,
                    source_row_count=len(rs))
            phase='single_candidate_full_sum'
            gl=load_module('full_gradient_batched_GL',archive/'repair_p3_batched_gl.py')
            batch=load_module('full_gradient_batched_conditional',archive/'repair_p3_batched_conditional.py')
            kernel=load_module('full_gradient_base_kernel',archive/'repair_p3_integral.py')
            target=kernel.Settings(logmass_tolerance=1e-9/GENES,gradient_tolerance=1e-8/GENES)
            candidates=[];candidate_refusals=[]
            for gene,rs in gene_rows.items():
                try:
                    context=batch.BatchedConditional(model,p,rs)
                    estimate=(gl.GL_workspace_estimate_bytes(model,rs,32)+context.memory_estimate['estimated_static_bytes']+
                        context.memory_estimate['estimated_maximum_batch_working_bytes'])
                    require(estimate<=3*1024**3,'Candidate allocation estimate exceeds3GiB')
                    def callback(row,etas,r,need_hessian=True):
                        return context.terms_batch(row,etas,r,p,need_hessian)
                    ll,g,_,report=gl.integrate_gene(model,p,v,rs,callback,with_hessian=False,
                        expected_kernel_sha=KERNEL_SHA,kernel_path=archive/'repair_p3_integral.py',settings=target,batch_size=32)
                    item=dict(gene=gene,logmass=ll,gradient=g,integration=report,allocation_estimate=estimate)
                    candidates.append(item);save('candidate_gene%02d'%gene,item)
                except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                    item=refusal(error,gene=gene);candidate_refusals.append(item);save('candidate_refusal_gene%02d'%gene,item)
            totals=[]
            for bound,_ in DOMAINS:
                key=str(bound)
                if all(key in masses[str(g)]['domains'] for g in range(GENES)):
                    totals.append(dict(bound=bound,logmass=sum(masses[str(g)]['domains'][key]['mass']['logmass'] for g in range(GENES)),
                        uncertified_logmass_error=sum(masses[str(g)]['domains'][key]['logmass_quad_plus_prior_tail_estimate'] for g in range(GENES))))
            candidate_total=(dict(logmass=sum(c['logmass'] for c in candidates),
                gradient=np.sum([c['gradient'] for c in candidates],axis=0),
                logmass_uncertified_error=sum(c['integration']['empirical_refinement_plus_tail_allowances'][0] for c in candidates),
                gradient_uncertified_error=np.sum([c['integration']['empirical_refinement_plus_tail_allowances'][1] for c in candidates],axis=0))
                if len(candidates)==GENES else None)
            comparisons=[dict(bound=r['bound'],difference=candidate_total['logmass']-r['logmass'],
                formal_absolute_tolerance=1e-6,within_nominal_tolerance=abs(candidate_total['logmass']-r['logmass'])<=1e-6,
                summed_uncertified_error=candidate_total['logmass_uncertified_error']+r['uncertified_logmass_error'],accuracy_accepted=False)
                for r in totals] if candidate_total is not None else []
            save('shared_baseline',dict(identity=identity,helper_sha256=args.expected_helper_sha,source_pins=source_declarations,
                masses=masses,independent_total_scalar_masses=totals,candidate_total=candidate_total,
                candidate_gene_count=len(candidates),candidate_per_gene_targets=dict(logmass=1e-9/GENES,gradient=1e-8/GENES),
                likelihood_comparisons=comparisons,mass_refusals=mass_refusals,candidate_refusals=candidate_refusals,
                reference_mass_complete=len(totals)==2,candidate_complete=len(candidates)==GENES,
                shared_probability_evaluations=sum(r['shared_geometry_mass_counter'] for r in masses.values()),
                no_fits=True,no_calibration=True,accuracy_accepted=False))
            complete=not mass_refusals and not candidate_refusals
        else:
            phase='validate_shared_baseline'
            require(baseline['identity']==identity and baseline['source_pins']==source_declarations and
                baseline['helper_sha256']==args.expected_helper_sha,'Baseline fixedpoint/source identity differs')
            require(baseline['reference_mass_complete'],'Baseline scalar full sum incomplete')
            active=(args.axis if args.axis<GENES else args.axis-GENES) if args.axis<2*GENES else None
            selected=[active] if active is not None else list(range(GENES))
            engines={};records=[];refusals=[]
            for gene in selected:
                engine=reference.Reference(gene_rows[gene],GENES,7,maximum_evaluations=EVAL_CAP)
                shared=baseline['masses'][str(gene)]['shared_geometry_mass_counter']
                require(isinstance(shared,int) and 0<=shared<=EVAL_CAP,'Shared probability counter invalid')
                engine.evaluations=shared;engines[gene]=engine
                phase='independent_complex_axis'
                for step in STEPS:
                    for bound,_ in DOMAINS:
                        key=str(bound);fields=dict(gene=gene,axis=args.axis,step=step,bound=bound)
                        try:
                            saved_mass=baseline['masses'][str(gene)]['domains'][key]
                            result=engine.derivative(p,v,args.axis,step,saved_mass['geometry'],saved_mass['mass'])
                            item=dict(**fields,completed=True,result=result,axis_budget_counter=engine.evaluations)
                        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as error:
                            item=refusal(error,**fields,axis_budget_counter=engine.evaluations)
                            refusals.append(item)
                        records.append(item);save('derivative_gene%02d_%s_B%02d'%(gene,str(step),bound),item)
            sums=[]
            candidate=baseline['candidate_total']
            for step in STEPS:
                for bound,_ in DOMAINS:
                    selected_records=[r for r in records if r['step']==step and r['bound']==bound]
                    item=dict(axis=args.axis,step=step,bound=bound,required_genes=selected,
                        other_gene_contributions_exact_zero=([g for g in range(GENES) if g!=active] if active is not None else []),
                        other_genes_zero_by_disjoint_gene_parameters=True if active is not None else False,
                        completed=len(selected_records)==len(selected) and all(r['completed'] for r in selected_records))
                    if item['completed']:
                        gradient=sum(r['result']['gradient'] for r in selected_records)
                        uncertainty=sum(r['result']['estimated_combined_gradient_uncertainty'] for r in selected_records)
                        item.update(gradient=gradient,summed_uncertified_error=uncertainty,
                            formal_gradient_tolerance=1e-8+1e-6*abs(gradient))
                        if candidate is not None:
                            difference=gradient-candidate['gradient'][args.axis]
                            combined=uncertainty+candidate['gradient_uncertified_error'][args.axis]
                            item.update(independent_minus_candidate=difference,candidate_gradient=candidate['gradient'][args.axis],
                                nominal_tolerance_within=abs(difference)<=item['formal_gradient_tolerance'],
                                combined_uncertified_error=combined,uncertainty_overlap=abs(difference)<=combined)
                    sums.append(item)
            pairs=[]
            for i,left in enumerate(sums):
                for right in sums[i+1:]:
                    if left['completed'] and right['completed']:
                        difference=right['gradient']-left['gradient']
                        error=right['summed_uncertified_error']+left['summed_uncertified_error']
                        pairs.append(dict(first=dict(step=left['step'],bound=left['bound']),second=dict(step=right['step'],bound=right['bound']),
                            difference=difference,summed_uncertified_error=error,disagreement_outside_estimate=abs(difference)>error))
            save('axis_results',dict(axis=args.axis,records=records,full_sum_results=sums,all_step_domain_comparisons=pairs,
                refusals=refusals,baseline_sha256=args.expected_baseline_sha,
                independent_total_scalar_masses_reused=baseline['independent_total_scalar_masses'],
                likelihood_comparisons_reused=baseline['likelihood_comparisons'],
                probability_counters=[dict(gene=g,shared_counter=baseline['masses'][str(g)]['shared_geometry_mass_counter'],
                    total_counter_including_refused_attempt=e.evaluations,
                    successful_new_derivative_calls=max(0,min(e.evaluations,EVAL_CAP)-baseline['masses'][str(g)]['shared_geometry_mass_counter']))
                    for g,e in engines.items()],per_gene_axis_cap=EVAL_CAP,aggregate_axis_cap=len(selected)*EVAL_CAP,
                baseline_evaluations_physically_performed_in_separate_job=True,
                uncertainty_not_certified=True,complex_truncation_and_arithmetic_unbounded=True,
                accuracy_accepted=False,no_fits=True,no_calibration=True))
            complete=not refusals and all(r['completed'] for r in sums)
        phase='source_postcheck'
        for source,digest in sources.items():
            require(sha(source)==digest and sha(archive/source.name)==digest,'Source changed during execution')
        save('summary',dict(mode=args.mode,axis=args.axis,numerical_records_completed=complete,
            required_combined_axes=list(range(29)),all29_accuracy_accepted=False,
            source_hashes_checked_before_after=True,fixed_point_not_keptfit_accuracy=True,
            elapsed_seconds=time.monotonic()-began,accuracy_accepted=False,no_fits=True,no_calibration=True,
            errors_uncertified=True,output_cap=CAP))
    except BaseException as error:
        body=json.dumps(dict(phase=phase,completed=False,exception_type=type(error).__name__,traceback=traceback.format_exc(),
            partials_preserved=True,no_fits=True,no_calibration=True),indent=2)+'\n'
        require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body.encode())<=CAP,'Failure cap exhausted')
        with (output/'failure.json').open('x') as handle:handle.write(body)
        raise


if __name__=='__main__':
    main()
