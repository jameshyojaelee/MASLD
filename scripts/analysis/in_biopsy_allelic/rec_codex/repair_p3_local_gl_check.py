#!/usr/bin/env python3
"""Immutable synthetic reconstruction and local GL checks; no fits/calibration."""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCE = REC/'repair_p3_fixed_point_failed0gene2_22039132/executed_sources/repair_p3_fixed_point_check.py'
SOURCE_SHA = 'd3f529ad20666c58dc8a84d01785d3f62c42f18a430a3a572a055c180b9b4a9e'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
HERE = Path(__file__).resolve().parent
LOCAL_GL = HERE/'repair_p3_local_gl.py'
MOMENTS = HERE/'repair_p3_local_moments.py'
MOMENTS_SHA = '95ec9f270d18bfbbc8f31f0857308f5f2acad053661a30b208049b7a3f59f31f'
BATCHED = HERE/'repair_p3_batched_conditional.py'
BATCHED_SHA = '184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76'
V2 = HERE/'repair_p3_prepared_conditional_v2.py'
V2_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
FROZEN_GL = HERE/'repair_p3_batched_gl.py'
FROZEN_GL_SHA = '8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52'
# protocol, summary, numerical-output hashes from completed synthetic checks.
CASES = {
 'failed0gene2': (
  'repair_p3_batched_gl_failed0gene2_22040257',
  ('1173cd333901e57ef4690922650beec0c1d9567c5df25c62e6a7e8ca173cee7b',
   '47eb4d2838a3dc96165acdf1b1c8a73ef9f429504025e6f4099cc247081e2535',
   'b2207a47757ed50d63d59802bc6d81566c23529337d6e5465cee6c35554bf6c5'),
  'repair_p3_local_moments_failed0gene2_22040916',
  ('87ef1ff104b617db4c72b492fadf8a9d93d1beebe4245a2d10a56343215a02da',
   '492c8f193e2ec5214d8820a99845fa9d17814e2e0df6d68fa87acf4a8b8d3858',
   '73f76d9d2788fe10e21798fb457b8ffc08119a9901885d74f91c93db3dfdc7fb')),
 'positive18': (
  'repair_p3_batched_gl_positive18_22040259',
  ('4fff9a95d228ce2f12ab0ff08350df24d225cabda53044d00ab2a96582e43b56',
   '3daebbcc620794a84aaa75ccdcd28d5817d027e4f5fa3f9795d47dca40cb4982',
   '8e84dcabf557dd71824936f9a581ed8a1405068330602359f6431bb3038425bf'),
  'repair_p3_local_moments_positive18_22040917',
  ('e895bad7d6333d1f4140f7eeb3bcf80f05a8f20114be31f73b275d41f7a914a1',
   'e3ca6fa6c3976bc539a793abad0b0c589a26f1c797bda52b1c5a9471a4a17483',
   'cdfa812c12ae3e3b4137c82531aadb5fd16a98a7e2782a1c70c2723ea15022a4'))}
IDENTITY = ('case','draw','fixed_gene','fixed_gene_rows','fixed_p','fixed_v',
            'count_vector_sha256','environment')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(ok,message):
    if not ok: raise ValueError(message)


def comparisons(module, model, p, rows, kernel, batch_module, moments_module,
                baseline, variance):
    """Compute only: both parameter dimensions, same rows and all nine ladders."""
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Individual compute required')
    import numpy as np
    original = sys.modules[model.evaluate.__func__.__module__]
    require(model.G==8 and model.Z==7 and model.size==28 and len(rows)==24
            and {r.gene for r in rows}=={2},'Exact archived synthetic context required')
    reference_active = moments_module.LocalGeneMoments(model,p,rows,kernel).augmented
    prior = baseline['integration']
    require(len(prior['ladder'])==9 and prior['requested_node_batch_size']==32,
            'Completed batched32 nine-ladder baseline required')
    qmodel = original.Model(rows,gene_count=379)
    q = np.zeros(qmodel.size)
    q[:model.G] = p[:model.G]
    q[qmodel.G:qmodel.G+model.G] = p[model.G:2*model.G]
    q[2*qmodel.G:] = p[2*model.G:]
    require(qmodel.size==770,'Dimension-only coordinate mapping differs')
    original_callables=(original.conditional_terms,original.bb_terms,original.Model.predictors,original.Model.evaluate)
    records=[]
    for current,parameters,label in ((model,p,'G8_P28'),(qmodel,q,'dimension_only_G379_P770')):
        phase='preallocation'
        record=dict(shape=label,G=current.G,P=current.size,Z=current.Z,row_count=len(rows),
                    fullscale_validation=False,no_fit=True,no_calibration=True)
        started=time.perf_counter()
        try:
            gl_bytes=module.GL_workspace_estimate_bytes(current,rows,32)
            interface_bytes=8*((current.size+1)**2+current.size+1)
            # Compare output matrices sequentially. Eight interface-sized arrays
            # cover full comparison intermediates, without full ladder storage.
            comparison_bytes=8*interface_bytes
            require(gl_bytes<256*1024**2,'Local integration workspace exceeds unchanged cap')
            phase='callback_preparation'
            context=batch_module.BatchedConditional(current,parameters,rows)
            memory=context.memory_estimate
            estimate=gl_bytes+memory['estimated_static_bytes']+memory['estimated_maximum_batch_working_bytes']+comparison_bytes
            require(estimate<1024**3,'Declared synthetic combined allocation estimate exceeds1GiB')
            record.update(GL_workspace_estimate_bytes=gl_bytes,interface_output_bytes=interface_bytes,
                comparison_workspace_estimate_bytes=comparison_bytes,callback_memory_estimate=memory,
                combined_declared_allocation_estimate_bytes=estimate,not_measured_process_peak=True,
                callback_setup_seconds=time.perf_counter()-started)
            print(json.dumps(dict(event='local_GL_preallocation',shape=label,
                                 combined_declared_allocation_estimate_bytes=estimate)),flush=True)
            def callback(row,etas,r,need_hessian=True):
                return context.terms_batch(row,etas,r,parameters,need_hessian)
            phase='full_nine_ladder_integration'
            integration_started=time.perf_counter()
            result=module.integrate_gene(current,parameters,variance,rows,callback,True,
                expected_kernel_sha=KERNEL_SHA,kernel_path=Path(kernel.__file__),batch_size=32)
            record['integration_seconds']=time.perf_counter()-integration_started
            local=moments_module.LocalGeneMoments(current,parameters,rows,kernel)
            active=local.augmented
            oldgradient=np.asarray(baseline['gradient'])
            oldhessian=np.asarray(baseline['hessian'])
            # All baseline inactive coordinates must be exactly zero; no masked
            # disagreement may be hidden by the subsequent local projection.
            projected_g=oldgradient[reference_active]
            projected_h=oldhessian[np.ix_(reference_active,reference_active)]
            oldlocal=moments_module.LocalGeneMoments(model,p,rows,kernel)
            require(np.array_equal(oldlocal.embed_vector(projected_g),oldgradient)
                and np.array_equal(oldlocal.embed_matrix(projected_h),oldhessian),
                'Full batched baseline has nonzero inactive moments')
            dg=result[1]-local.embed_vector(projected_g)
            dh=result[2]-local.embed_matrix(projected_h)
            primary=dict(logmass_absolute_difference=abs(result[0]-baseline['logmass']),
                         gradient_max_absolute_difference=float(np.max(np.abs(dg))),
                         hessian_max_absolute_difference=float(np.max(np.abs(dh))))
            def sparse_matrix(a):
                i,j=np.nonzero(a)
                return [[int(x),int(y),float(a[x,y])] for x,y in zip(i,j)]
            phase='all_ladder_and_error_estimate_comparison'
            newreport=result[3]
            require(len(newreport['ladder'])==9,'Changed nine-ladder count')
            ladder_records=[]
            array_fields={'gradient':1,'hessian':2,'gradient_error_allowance':1,
                'hessian_error_allowance':2,'relative_gradient_tail_allowance':1,
                'relative_hessian_tail_allowance':2}
            scalar_fields=('logmass','logmass_error_allowance','relative_mass_tail_allowance',
                           'log_absolute_prior_mass_tail_bound')
            grid_fields=('z_bound','grid_spacing','panels','order','node_evaluations','cumulative_evaluations')
            def compare_mixed(actual,old):
                deltas=[]
                for index,(a,b) in enumerate(zip(actual,old)):
                    if index==0: deltas.append(float(abs(a-b)))
                    else:
                        b=np.asarray(b)
                        expected=b[reference_active] if index==1 else b[np.ix_(reference_active,reference_active)]
                        deltas.append(float(np.max(np.abs(np.asarray(a)-expected))))
                return deltas
            for ordinal,(a,b) in enumerate(zip(newreport['ladder'],prior['ladder'])):
                grid_equal=all(a[k]==b[k] for k in grid_fields)
                differences={k:float(abs(a[k]-b[k])) for k in scalar_fields}
                array_differences={}
                for key,dimension in array_fields.items():
                    oldarray=np.asarray(b[key]); expected=(oldarray[reference_active] if dimension==1
                        else oldarray[np.ix_(reference_active,reference_active)])
                    embedded_old=(oldlocal.embed_vector(expected) if dimension==1 else oldlocal.embed_matrix(expected))
                    require(np.array_equal(embedded_old,oldarray),'Saved ladder has nonzero inactive '+key)
                    delta=np.asarray(a[key])-expected
                    array_differences[key]=dict(max_absolute_difference=float(np.max(np.abs(delta))),
                        all_nonzero_differences=(sparse_matrix(delta) if dimension==2 else
                            [[int(j),float(delta[j])] for j in np.nonzero(delta)[0]]))
                order_deltas={key:compare_mixed(a[key],b[key]) for key in
                    ('order_refinement_16_32_allowances','order_refinement_32_64_allowances') if key in a or key in b}
                passed=(grid_equal and max(differences.values())<=1e-9 and
                    all(v['max_absolute_difference']<=(1e-8 if array_fields[k]==1 else 1e-5)
                        for k,v in array_differences.items()) and
                    all(d[0]<=1e-9 and d[1]<=1e-8 and d[2]<=1e-5 for d in order_deltas.values()))
                ladder_records.append(dict(ordinal=ordinal,grid_fields_equal=grid_equal,
                    scalar_differences=differences,array_differences=array_differences,
                    order_refinement_allowance_differences=order_deltas,within_unchanged_allowances=bool(passed)))
            outer_mixed=compare_mixed(newreport['empirical_refinement_plus_tail_allowances'],
                                      prior['empirical_refinement_plus_tail_allowances'])
            outer_differences={}
            for key,dimension in (('gradient_refinement_difference',1),('hessian_refinement_difference',2)):
                oldarray=np.asarray(prior[key]); expected=(oldarray[reference_active] if dimension==1
                    else oldarray[np.ix_(reference_active,reference_active)])
                outer_differences[key]=float(np.max(np.abs(np.asarray(newreport[key])-expected)))
            outer_differences['logmass_refinement_difference']=abs(newreport['logmass_refinement_difference']-prior['logmass_refinement_difference'])
            passed=(primary['logmass_absolute_difference']<=1e-9 and primary['gradient_max_absolute_difference']<=1e-8
                and primary['hessian_max_absolute_difference']<=1e-5 and all(r['within_unchanged_allowances'] for r in ladder_records)
                and newreport['likelihood_evaluations']==prior['likelihood_evaluations']
                and outer_mixed[0]<=1e-9 and outer_mixed[1]<=1e-8 and outer_mixed[2]<=1e-5
                and outer_differences['logmass_refinement_difference']<=1e-9
                and outer_differences['gradient_refinement_difference']<=1e-8
                and outer_differences['hessian_refinement_difference']<=1e-5)
            require(original_callables==(original.conditional_terms,original.bb_terms,original.Model.predictors,original.Model.evaluate),
                'Archived original conditional/predictor/evaluator callable changed')
            local.validate(parameters)
            record.update(original_callables_unchanged=True,integration_returned=True,all_comparisons_within_unchanged_allowances=bool(passed),
                primary_differences=primary,all_nonzero_interface_gradient_differences=[[int(j),float(dg[j])] for j in np.nonzero(dg)[0]],
                all_nonzero_interface_hessian_differences=sparse_matrix(dh),
                active_indices=local.active,augmented_indices=active,full_interface_shape=list(result[2].shape),
                ladder_comparisons=ladder_records,full_local_integration_report=newreport,
                outer_error_estimate_differences=outer_mixed,outer_refinement_differences=outer_differences,
                adaptive_chunk_ledger=context.chunk_ledger())
        except Exception as exc:
            record.update(integration_returned=False,all_comparisons_within_unchanged_allowances=False,
                failure_phase=phase,exception_type=type(exc).__name__,exception_message=str(exc),
                traceback=traceback.format_exc(),partial_integration_reports=getattr(exc,'integration_reports',None),
                partial_likelihood_evaluations=getattr(exc,'likelihood_evaluations',None))
        record['elapsed_seconds']=time.perf_counter()-started
        records.append(record)
    return dict(comparison_records_completed=True,records=records,record_count=len(records),
        all_comparisons_within_unchanged_allowances=all(r['all_comparisons_within_unchanged_allowances'] for r in records),
        all_failures_retained=True,unchanged24_synthetic_rows=True,
        dimension_only_G379_not_fullscale_validation=True,no_fits=True,no_calibration=True,
        no_independent_probability_reference=True,source_refinement_errors_not_certified=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=tuple(CASES),required=True)
    parser.add_argument('--expected-checker-sha',required=True)
    parser.add_argument('--expected-local-gl-sha',required=True)
    parser.add_argument('--inspect',action='store_true')
    args=parser.parse_args()
    bd,bhash,md,mhash=CASES[args.case]
    guards={SOURCE:SOURCE_SHA,Path(__file__).resolve():args.expected_checker_sha,
        LOCAL_GL:args.expected_local_gl_sha,MOMENTS:MOMENTS_SHA,BATCHED:BATCHED_SHA,V2:V2_SHA,FROZEN_GL:FROZEN_GL_SHA}
    for directory,names,digests in ((bd,('protocol.json','summary.json','candidate_baseline.json'),bhash),
                                  (md,('protocol.json','summary.json','local_moment_comparison.json'),mhash)):
        guards.update({REC/directory/n:d for n,d in zip(names,digests)})
    for path,digest in guards.items():
        require(sha(path)==digest,'Reviewed source/completed synthetic receipt changed: '+str(path))
    baseline_protocol=json.loads((REC/bd/'protocol.json').read_text())
    moment_protocol=json.loads((REC/md/'protocol.json').read_text())
    require(all(baseline_protocol[k]==moment_protocol[k] for k in IDENTITY),'Saved synthetic identities differ')
    bs=json.loads((REC/bd/'summary.json').read_text()); ms=json.loads((REC/md/'summary.json').read_text())
    mc=json.loads((REC/md/'local_moment_comparison.json').read_text())
    require(bs['candidate_integral_returned'] and bs['fixed_point_check_completed']
        and ms['local_moment_comparison_completed'] and mc['comparison_records_completed']
        and mc['all_comparisons_within_unchanged_allowances'] and mc['record_count']==12,
        'Required batched integral/local moment checks are incomplete or failing')
    require(MOMENTS_SHA in moment_protocol['source_guards'].values(), 'Passed moments used a different source')
    source=SOURCE.read_text()
    changes={'repair_p3_fixed_point_':'repair_p3_local_gl_',
        "save('predictor_cache_comparison', report)":"save('local_gl_comparison', report)",
        'predictor_comparison_completed=True':'local_gl_comparison_completed=True',
        'no_integration=True':'no_integration=False'}
    for old,new in changes.items():
        require(source.count(old)==1,'Immutable reconstruction output label differs')
        source=source.replace(old,new)
    namespace={'__file__':str(SOURCE),'__name__':'private_local_GL_reconstruction'}
    exec(compile(ast.parse(source),str(SOURCE),'exec'),namespace)
    namespace['PREDICTOR_CACHE']=LOCAL_GL
    namespace['GUARDS'].update(guards)
    original_loader=namespace['load_module']
    def local_loader(name,path):
        module=original_loader(name,path)
        if Path(path).name==LOCAL_GL.name:
            output=REC/('repair_p3_local_gl_'+args.case+'_'+os.environ['SLURM_JOB_ID'])
            protocol=json.loads((output/'protocol.json').read_text())
            require(all(protocol[k]==baseline_protocol[k] for k in IDENTITY),'Exact frozen reconstructed context differs')
            batch_module=original_loader('repair_p3_local_GL_batch',Path(path).with_name(BATCHED.name))
            moments_module=original_loader('repair_p3_local_GL_moments',Path(path).with_name(MOMENTS.name))
            baseline=json.loads((REC/bd/'candidate_baseline.json').read_text())
            def hook(model,p,rows):
                return comparisons(module,model,p,rows,sys.modules['repair_p3_private_kernel'],
                    batch_module,moments_module,baseline,protocol['fixed_v'])
            module.compare_cached_predictors=hook
        return module
    namespace['load_module']=local_loader
    sys.argv=[str(SOURCE),'--case',args.case,'--predictor-points-sha',args.expected_local_gl_sha,
              '--expected-helper-sha',SOURCE_SHA,'--expected-kernel-sha',KERNEL_SHA]
    if args.inspect:sys.argv.append('--inspect')
    namespace['main']()


if __name__=='__main__':main()
