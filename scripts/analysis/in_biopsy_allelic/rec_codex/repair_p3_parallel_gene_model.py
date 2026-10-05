"""Isolated fork-based per-gene evaluation of the frozen batched GL rule.

Small G8 synthetic scope only. No fitter, likelihood, boundary, numerical
budget or integration-node change. Workers never aggregate the model sum.
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import hashlib
import multiprocessing
import os
from pathlib import Path
import sys
import threading
import time
import traceback

KERNEL_SHA='461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
GL_SHA='8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52'
CONDITIONAL_SHA='184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76'
SEQUENTIAL_SHA='2f64fdaffb44e5d78fd327ca6caf8672242ff53e17e53ad89b7e249bd77c0c3c'
_WORK_CONTEXT=None


def _worker_gene(gene):
    """Top-level callable pickled by name; state inherited via explicit fork."""
    model,p,v,with_hessian,private,kernel,gl,conditional_module,kernel_path,batch_size,settings=_WORK_CONTEXT
    rows=tuple(row for row in model.rows if row.gene==gene)
    if not rows:return dict(gene=gene,present=False,success=True)
    signatures=tuple(tuple(vars(r).items()) for r in rows)
    began=time.perf_counter()
    phase='per_gene_preparation'
    try:
        kernel.require(not p.flags.writeable,'Worker parameter vector must be read-only')
        kernel.require(sum(r.n+1 for r in rows)<=100000,'Per-gene prepared support scope exceeded')
        setup=time.perf_counter()
        context=conditional_module.BatchedConditional(model,p,rows)
        setup_seconds=time.perf_counter()-setup
        kernel.require(context.original is private and context.original.conditional_terms is private.conditional_terms,
            'Prepared conditional escaped private original')
        combined=(gl.GL_workspace_estimate_bytes(model,rows,batch_size)+context.memory_estimate['estimated_static_bytes']
                  +context.memory_estimate['estimated_maximum_batch_working_bytes'])
        kernel.require(combined<=3*1024**3,'Existing declared3GiB workspace allowance exceeded')
        def conditional(row,etas,r,need_hessian=True):
            return context.terms_batch(row,etas,r,p,need_hessian)
        phase='unchanged_per_gene_GL'
        ll,g,h,report=gl.integrate_gene(model,p,v,rows,conditional,with_hessian=with_hessian,
            expected_kernel_sha=KERNEL_SHA,kernel_path=kernel_path,settings=settings,batch_size=batch_size)
        kernel.require(tuple(tuple(vars(r).items()) for r in rows)==signatures and not p.flags.writeable,
            'Worker fixed row/parameter state changed')
        context.prepared.validate(p)
        report['batched_callback']=dict(setup_seconds=setup_seconds,callback_memory_estimate=context.memory_estimate,
            combined_declared_allocation_estimate_bytes=combined,estimate_not_measured_process_peak=True,
            adaptive_chunk_ledger=context.chunk_ledger())
        return dict(gene=gene,present=True,success=True,logmass=ll,gradient=g,hessian=h,
                    report=dict(gene=gene,**report),elapsed_seconds=time.perf_counter()-began,
                    worker_pid=os.getpid(),fixed_inputs_unchanged=True)
    except Exception as exc:
        return dict(gene=gene,present=True,success=False,failure_phase=phase,
            error_type=type(exc).__name__,message=str(exc),value_error=isinstance(exc,ValueError),
            linalg_error=isinstance(exc,kernel.np.linalg.LinAlgError),
            refused_gene_ladders=getattr(exc,'integration_reports',None),
            refused_gene_evaluations=getattr(exc,'likelihood_evaluations',None),
            traceback=traceback.format_exc(),elapsed_seconds=time.perf_counter()-began,worker_pid=os.getpid())


def parallel_model_class(private,kernel,gl,conditional_module,sequential_module,*,kernel_path,batch_size=32,event_sink=None):
    if not os.environ.get('SLURM_JOB_ID','').isdigit():raise ValueError('Individual compute allocation required')
    import numpy as np
    for module,expected in ((kernel,KERNEL_SHA),(gl,GL_SHA),(conditional_module,CONDITIONAL_SHA),(sequential_module,SEQUENTIAL_SHA)):
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()!=expected:
            raise ValueError('Frozen parallel evaluation dependency differs')
    Sequential=sequential_module.batched_model_class(private,kernel,gl,conditional_module,
        kernel_path=kernel_path,batch_size=batch_size)
    settings=kernel.Settings()
    kernel_path=Path(kernel_path)
    class ParallelGeneModel(Sequential):
        def marginal(self,p,v,nodes=20,with_hessian=False,gene_terms=None):
            global _WORK_CONTEXT
            kernel.require(len(self.rows)<=256 and 0<self.G<=8,'Small synthetic model scope exceeded')
            if v<=0:return super().marginal(p,v,nodes,with_hessian,gene_terms)
            kernel.require(sys.platform.startswith('linux') and threading.active_count()==1,
                'Explicit Linux fork requires a single Python-thread parent')
            kernel.require(int(os.environ.get('SLURM_CPUS_PER_TASK','0'))>=8
                and len(os.sched_getaffinity(0))>=8,'Eight allocated and available CPUs required')
            kernel.require(all(os.environ.get(name)=='1' for name in
                ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
                'Inherited numerical worker threads must be fixed to1 before imports')
            kernel.require(_WORK_CONTEXT is None,'Concurrent/reentrant fork evaluator is unsupported')
            kernel.require(np.asarray(p).shape==(self.size,) and np.isfinite(p).all() and np.isfinite(v) and v>0,
                'Finite model parameters and variance required')
            began=time.perf_counter()
            self.marginal_call_count=getattr(self,'marginal_call_count',0)+1
            call=self.marginal_call_count
            self.last_integration_failure=None
            self.last_integration_reports=[]
            self.last_total_error_allowance=None
            self.last_evaluated_point=dict(call=call,parameters=np.asarray(p).copy(),variance=float(v),with_hessian=with_hessian)
            present=sum(any(r.gene==gene for r in self.rows) for gene in range(self.G))
            gene_settings=replace(settings,logmass_tolerance=min(settings.logmass_tolerance,1e-8/max(1,present)),
                gradient_tolerance=min(settings.gradient_tolerance,1e-7/max(1,present)))
            estimates=[]
            for gene in range(self.G):
                rows=[r for r in self.rows if r.gene==gene]
                if not rows:continue
                entries=sum(r.n+1 for r in rows)
                kernel.require(entries<=100000,'Per-gene prepared support scope exceeded')
                glbytes=gl.GL_workspace_estimate_bytes(self,rows,batch_size)
                kernel.require(glbytes<=256*1024**2,'Unchanged GL workspace cap exceeded')
                estimates.append(dict(gene=gene,GL_workspace=glbytes,
                    callback_static=512*entries+1048576,callback_maximum_batch=512*100000+1048576))
            aggregate=sum(e['GL_workspace']+e['callback_static']+e['callback_maximum_batch'] for e in estimates)
            # Explicit generous interpreter/library/result/parent capacity estimate,
            # not measured RSS or a certified maximum. G8/P28 only is admitted.
            overhead=9*256*1024**2
            kernel.require(self.size==28 and self.Z==7 and aggregate+overhead<=6*1024**3,
                'Declared small-synthetic sixGiB capacity estimate exceeded')
            readonly=np.asarray(p).copy();readonly.setflags(write=False)
            parent_rows=tuple(tuple(vars(r).items()) for r in self.rows)
            original_callables=(private.conditional_terms,private.bb_terms,private.Model.predictors,private.Model.evaluate)
            _WORK_CONTEXT=(self,readonly,float(v),with_hessian,private,kernel,gl,conditional_module,kernel_path,batch_size,gene_settings)
            results=[]
            try:
                with ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('fork')) as pool:
                    for record in pool.map(_worker_gene,range(self.G)):
                        results.append(record)
            except Exception as exc:
                # Preserve the collected ordered prefix; unknown worker results
                # are explicitly unavailable, never interpreted as zero/success.
                self.last_parallel_worker_records=results
                self.last_integration_reports=[r['report'] for r in results if r['present'] and r['success']]
                self.last_integration_failure=dict(call=call,variance=float(v),parameters=np.asarray(p).copy(),
                    failure_phase='fork_pool_infrastructure',error_type=type(exc).__name__,message=str(exc),
                    traceback=traceback.format_exc(),collected_worker_prefix=results,
                    unavailable_worker_genes=list(range(len(results),self.G)),no_full_sum_formed=True)
                if event_sink:event_sink(dict(event='marginal_pool_failure',call=call,variance=float(v),
                    error_type=type(exc).__name__,message=str(exc),collected_worker_count=len(results)))
                raise  # Preserve infrastructure exception class and abort semantics.
            finally:
                _WORK_CONTEXT=None
            kernel.require(tuple(tuple(vars(r).items()) for r in self.rows)==parent_rows
                and np.array_equal(readonly,p) and not readonly.flags.writeable,'Fork parent input changed')
            kernel.require(original_callables==(private.conditional_terms,private.bb_terms,private.Model.predictors,private.Model.evaluate),
                'Private original callable changed')
            self.last_parallel_worker_records=results
            self.last_parallel_capacity_estimate=dict(workers=8,per_gene=estimates,
                aggregate_numerical_workspaces=aggregate,interpreter_library_parent_overhead_estimate=overhead,
                combined_declared_estimate=aggregate+overhead,not_measured_or_certified=True)
            failed=[r for r in results if not r['success']]
            if failed:
                first=failed[0]
                completed=[r['report'] for r in results if r['gene']<first['gene'] and r['present'] and r['success']]
                self.last_integration_reports=completed
                self.last_integration_failure=dict(call=call,gene=first['gene'],variance=float(v),parameters=np.asarray(p).copy(),
                    completed_gene_reports=completed,refused_gene_ladders=first['refused_gene_ladders'],
                    error_type=first['error_type'],message=first['message'],all_worker_records=results,
                    no_full_sum_formed=True,first_failure_is_ascending_gene=True)
                if gene_terms is not None:
                    for r in results:
                        if r['gene']>=first['gene']:break
                        gene_terms[r['gene']]=r.get('logmass',0.)
                if event_sink:event_sink(dict(event='marginal_refusal',call=call,gene=first['gene'],variance=float(v),
                    with_hessian=with_hessian,seconds=time.perf_counter()-began,error_type=first['error_type'],message=first['message']))
                if first['linalg_error']:
                    raise np.linalg.LinAlgError(first['message'])
                if first['value_error']:
                    exc=kernel.IntegrationFailure(first['message']);exc.integration_reports=first['refused_gene_ladders']
                    raise exc
                raise RuntimeError('Fork worker non-numerical exception: '+first['error_type']+': '+first['message'])
            total,gradient=0.,np.zeros(self.size+1)
            hessian=np.zeros((self.size+1,self.size+1)) if with_hessian else None
            reports=[]
            # Identical ascending gene order and per-array += operations as frozen
            # sequential marginal. No tree reduction or completion-order sum.
            for r in results:
                if not r['present']:
                    if gene_terms is not None:gene_terms[r['gene']]=0.
                    continue
                total+=r['logmass'];gradient+=r['gradient']
                if with_hessian:hessian+=r['hessian']
                if gene_terms is not None:gene_terms[r['gene']]=r['logmass']
                reports.append(r['report'])
            self.last_integration_reports=reports
            errors=dict(logmass=sum(r['refined']['logmass_error_allowance']+r['logmass_refinement_difference'] for r in reports),
                gradient=sum((r['refined']['gradient_error_allowance']+r['gradient_refinement_difference'] for r in reports),np.zeros(self.size+1)),
                hessian=(sum((r['refined']['hessian_error_allowance']+r['hessian_refinement_difference'] for r in reports),np.zeros_like(hessian))
                    if with_hessian else None),empirical_refinement_and_tail_estimates_not_certified=True)
            self.last_total_error_allowance=errors
            kernel.require(errors['logmass']<=1e-8 and np.max(errors['gradient'])<=1e-7,
                'Accumulated gene numerical uncertainty exceeds unchanged full-model allowance')
            if event_sink:event_sink(dict(event='marginal_evaluation',call=call,variance=float(v),with_hessian=with_hessian,
                seconds=time.perf_counter()-began,loglik=float(total),max_abs_score=float(np.max(np.abs(gradient))),
                integration_genes=len(reports),logmass_error_estimate=float(errors['logmass']),
                maximum_gradient_error_estimate=float(np.max(errors['gradient'])),parallel_workers=8))
            return (float(total),gradient,hessian) if with_hessian else (float(total),gradient)
    for name in sequential_module.INHERITED:
        if getattr(ParallelGeneModel,name) is not getattr(private.Model,name):
            raise ValueError('Original fitting/design method changed: '+name)
    return ParallelGeneModel
