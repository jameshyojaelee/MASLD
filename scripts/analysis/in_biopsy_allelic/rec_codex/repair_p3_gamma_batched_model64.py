"""Prepared G8 synthetic64-node GL with eight native-law workers.

The parent keeps the original BB, conditional, null-fit and v=0 callables.
Positive-variance workers use checked batched Gamma conditional mixtures
and gene-local moments, embedding their error arrays before summation. All fitter methods, bounds, starts, ties and restarts are inherited.
This source alone does not admit an integral, fit or calibration run.
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

ORIGINAL_SHA = '41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
GL_SHA = '254a1efce7d8a873f8caeda1f0630cb48e35286e15bd952352246ab26b01e4de'
CONDITIONAL_SHA = '007e46bc796810e1adaf2c7cec0c273706e9e47aa9770acbe636b334600923a1'
INHERITED = ('fit','fit_single','constrained_polish','fit_marginal','fit_marginal_single',
             'bounds','initial','cold_start','informative','predictors','evaluate','sandwich')
_WORK_CONTEXT = None
NUMERICAL_MESSAGES={
    'repair_p3_integral.py':{'Nonfinite derivative moments','Invalid GL integral'},
    'repair_p3_local_gl.py':{'Composite GL evaluation cap exhausted',
        'GL order refinement exceeds unchanged precision allowance',
        'Composite GL relative-tail domain cap exhausted','Composite GL expanded domain exceeds cap',
        'Composite GL tail allowances remain unresolved','Composite GL final domain exceeds cap',
        'Composite GL expanded domain/grid refinement disagrees'},
    'repair_p3_batched_conditional.py':{'Nonfinite batch conditional result'},
    'repair_p3_bb_gamma_batched.py':{'Nonfinite rising-product batch','Nonfinite BB batch'},
    'repair_p3_gamma_batched_model64.py':{'Conditional arithmetic refusal'},
}


class ModelIntegrityError(RuntimeError):
    """Invalid evaluator state; do not count as a scientific refit failure."""


def require_integrity(ok,message):
    if not ok:
        raise ModelIntegrityError(message)


def row_signatures(rows,np):
    """Detached value snapshots, including any stored array-valued fields."""
    def freeze(value):
        if isinstance(value,np.ndarray):
            if value.dtype.hasobject:
                raise ModelIntegrityError('Object arrays are outside the prepared synthetic scope')
            return ('array',value.dtype.str,value.shape,value.tobytes(order='C'))
        if isinstance(value,np.generic):
            return ('numpy_scalar',value.dtype.str,value.tobytes())
        if isinstance(value,float):
            return ('float',value.hex())
        if value is None or isinstance(value,(str,int,bool)):
            return (type(value).__name__,value)
        if isinstance(value,(tuple,list)):
            return (type(value).__name__,tuple(freeze(x) for x in value))
        if isinstance(value,dict):
            return ('dict',tuple((k,freeze(v)) for k,v in sorted(value.items())))
        raise ModelIntegrityError('Unsupported stored row value: '+type(value).__name__)
    return tuple(tuple((key,freeze(value)) for key,value in sorted(vars(row).items())) for row in rows)


def numerical_worker_refusal(exc):
    """Only explicit arithmetic refusals from the pinned helper sources qualify.

    Those helpers predate typed state/arithmetic exceptions. Unknown messages,
    wrong shapes and state checks remain integrity errors. Partial ladders are
    retained when available; a refusal before a report has no inferred count.
    """
    current=exc
    while current is not None:
        tb=current.__traceback__
        while tb is not None:
            origin=Path(tb.tb_frame.f_code.co_filename)
            if (origin.resolve().parent==Path(__file__).resolve().parent
                    and str(current) in NUMERICAL_MESSAGES.get(origin.name,set())):
                return True
            tb=tb.tb_next
        current=current.__cause__
    return False


def first_failed_worker(results):
    """An integrity failure in any gene takes priority over numeric refusals."""
    failed=[r for r in results if not r['success']]
    integrity=[r for r in failed if not r['numerical_integration_refusal'] and not r['linalg_error']]
    return (integrity or failed or [None])[0]


def embed_report_errors(report,model,gene,np,with_hessian):
    """Map local nonnegative allowances into the original parameter coordinates."""
    active=sorted([gene,model.G+gene,2*model.G,
        *range(2*model.G+1,2*model.G+1+model.Z),2*model.G+1+model.Z,
        *range(model.size-3,model.size)])
    indices=active+[model.size]
    require_integrity(len(indices)==15 and len(set(indices))==15
        and report['active_parameter_indices']==active
        and report['augmented_global_indices']==indices
        and report['full_parameter_dimension']==model.size
        and report['local_augmented_dimension']==len(indices),
        'Local/global derivative coordinate metadata differs')
    raw=report['refined']
    lparts=(raw['logmass_error_allowance'],report['logmass_refinement_difference'])
    require_integrity(all(np.isfinite(x) and x>=0 for x in lparts),'Invalid likelihood allowances')
    def components(first,second,shape):
        a,b=np.asarray(first),np.asarray(second)
        require_integrity(a.shape==b.shape==shape and np.isfinite(a).all()
            and np.isfinite(b).all() and np.all(a>=0) and np.all(b>=0),
            'Invalid local error allowance components')
        return a+b
    g=components(raw['gradient_error_allowance'],report['gradient_refinement_difference'],(15,))
    full_g=np.zeros(model.size+1);full_g[indices]=g
    full_h=None
    if with_hessian:
        h=components(raw['hessian_error_allowance'],report['hessian_refinement_difference'],(15,15))
        full_h=np.zeros((model.size+1,model.size+1));full_h[np.ix_(indices,indices)]=h
    return dict(logmass=sum(lparts),gradient=full_g,hessian=full_h,
                augmented_global_indices=indices,local_dimension=15,full_dimension=model.size+1)


def _worker_gene(gene):
    model,p,v,with_hessian,private,kernel,gl,conditional,kernel_path,settings = _WORK_CONTEXT
    rows = tuple(row for row in model.rows if row.gene == gene)
    if not rows:
        return dict(gene=gene,present=False,success=True)
    before = row_signatures(rows,kernel.np)
    before_parameters = p.tobytes()
    original_bb = private.bb_terms
    entered_nodes = 0
    began = time.perf_counter()
    try:
        require_integrity(not p.flags.writeable,'Worker parameters must be read-only')
        Context = conditional.context_class(Path(conditional.__file__).parent)
        context = Context(model,p,rows)
        def checked_callback(row,etas,linear_r,need_hessian=True):
            nonlocal entered_nodes
            # GL increments by the complete batch before calling its first row.
            if row is rows[0]:entered_nodes+=len(etas)
            values,scores,hessians=context.callback(row,etas,linear_r,need_hessian)
            count=len(etas);np=kernel.np
            require_integrity(values.shape==(count,) and scores.shape==(count,2)
                and (not need_hessian or hessians.shape==(count,2,2)),
                'Invalid conditional output shapes')
            if not (np.isfinite(values).all() and np.all(values<=1e-8)
                    and np.isfinite(scores).all()
                    and (not need_hessian or np.isfinite(hessians).all())):
                exc=kernel.IntegrationFailure('Conditional arithmetic refusal')
                exc.integration_reports=[]
                exc.likelihood_evaluations=None
                raise exc
            return values,scores,hessians
        ll,g,h,report = gl.integrate_gene(model,p,v,rows,checked_callback,
            with_hessian=with_hessian,expected_kernel_sha=KERNEL_SHA,
            kernel_path=kernel_path,settings=settings,batch_size=64)
        errors = embed_report_errors(report,model,gene,kernel.np,with_hessian)
        if row_signatures(rows,kernel.np) != before or p.tobytes() != before_parameters or p.flags.writeable or private.bb_terms is not original_bb:
            raise ModelIntegrityError('Worker fixed inputs changed')
        return dict(gene=gene,present=True,success=True,logmass=ll,gradient=g,hessian=h,
            report=dict(gene=gene,**report),error_allowances=errors,
            conditional_BB_calls=context.gamma_BB_calls,support_calls=context.support_calls,
            conditional_chunk_ledger=context.chunk_ledger(),worker_pid=os.getpid(),
            elapsed_seconds=time.perf_counter()-began,fixed_inputs_unchanged=True)
    except Exception as exc:
        state_unchanged=(row_signatures(rows,kernel.np)==before
            and p.tobytes()==before_parameters and not p.flags.writeable
            and private.bb_terms is original_bb)
        raw_count=getattr(exc,'likelihood_evaluations',None)
        count=raw_count if isinstance(raw_count,int) and raw_count>=0 else entered_nodes
        raw_reports=getattr(exc,'integration_reports',None)
        return dict(gene=gene,present=True,success=False,error_type=type(exc).__name__,
            message=str(exc),value_error=isinstance(exc,ValueError),
            numerical_integration_refusal=state_unchanged and numerical_worker_refusal(exc),
            linalg_error=state_unchanged and isinstance(exc,kernel.np.linalg.LinAlgError),
            fixed_inputs_unchanged=state_unchanged,
            refused_gene_ladders=raw_reports if isinstance(raw_reports,list) else [],
            refused_gene_ladders_available=isinstance(raw_reports,list),
            refused_gene_evaluations=count,raw_exception_likelihood_evaluations=raw_count,
            evaluation_count_basis='exception' if count is raw_count else 'first-row GL batches entered',
            traceback=traceback.format_exc(),worker_pid=os.getpid(),
            elapsed_seconds=time.perf_counter()-began)
    finally:
        private.bb_terms = original_bb


def gamma_batched_model_class(private,kernel,gl,conditional,*,kernel_path,event_sink=None):
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise ModelIntegrityError('Individual compute allocation required before numerical imports')
    import numpy as np
    for module,expected in ((private,ORIGINAL_SHA),(kernel,KERNEL_SHA),(gl,GL_SHA),(conditional,CONDITIONAL_SHA)):
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != expected:
            raise ModelIntegrityError('Reviewed native-law numerical source differs')
    kernel_path = Path(kernel_path)
    require_integrity(Path(conditional.__file__).parent.resolve() == kernel_path.parent.resolve()
        == Path(gl.__file__).parent.resolve()==Path(__file__).parent.resolve(),
        'All numerical helpers and wrapper must share the checked archive')
    require_integrity(hashlib.sha256(kernel_path.read_bytes()).hexdigest() == KERNEL_SHA,
                   'Reviewed kernel path differs')
    import ast
    for filename,messages in NUMERICAL_MESSAGES.items():
        # kernel.require raises the messages supplied by the pinned GL caller.
        defining_file='repair_p3_local_gl.py' if filename=='repair_p3_integral.py' else filename
        literals={node.value for node in ast.walk(ast.parse((kernel_path.parent/defining_file).read_text()))
                  if isinstance(node,ast.Constant) and isinstance(node.value,str)}
        require_integrity(messages<=literals,'Pinned numerical refusal messages differ: '+filename)
    require_integrity(private.bb_terms.__name__ == 'bb_terms'
                   and Path(private.bb_terms.__code__.co_filename).resolve() == Path(private.__file__).resolve(),
                   'Parent must retain the original BB function')
    original_callables = (private.bb_terms,private.conditional_terms,private.Model.predictors,private.Model.evaluate)
    settings = kernel.Settings()

    class GammaParallelModel(private.Model):
        def marginal(self,p,v,nodes=20,with_hessian=False,gene_terms=None):
            global _WORK_CONTEXT
            # nodes is compatibility-only for positive GL; inherited fit
            # methods do not compare node counts as a convergence test.
            self.last_integration_failure = None
            self.last_integration_reports = []
            self.last_total_error_allowance = None
            self.last_evaluated_point = None
            require_integrity(original_callables == (private.bb_terms,private.conditional_terms,
                private.Model.predictors,private.Model.evaluate),'Original parent callables changed')
            if v <= 0:
                # Exact original call path, including negative-v refusal and
                # unavailable v-Hessian entries at the zero boundary.
                return super().marginal(p,v,nodes,with_hessian,gene_terms)
            require_integrity(self.G == 8 and self.Z == 7 and self.size == 28
                and 0 < len(self.rows) <= 256
                and all(row.gene in range(8) and 0 <= row.n <= 5000 for row in self.rows),
                'Prepared small-synthetic scope exceeded; do not drop rows to fit this helper')
            require_integrity(all(gl.GL_workspace_estimate_bytes(self,
                tuple(row for row in self.rows if row.gene==gene),64)<=256*1024**2
                for gene in range(8) if any(row.gene==gene for row in self.rows)),
                'Requested64-node GL workspace exceeds unchanged256MiB limit')
            require_integrity(sys.platform.startswith('linux') and threading.active_count() == 1
                and int(os.environ.get('SLURM_CPUS_PER_TASK','0')) >= 8
                and len(os.sched_getaffinity(0)) >= 8,'Single-thread parent and eight allocated CPUs required')
            require_integrity(all(os.environ.get(key) == '1' for key in
                ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
                'Numerical subprocess threads must be fixed to one')
            memory_mb = int(os.environ.get('SLURM_MEM_PER_NODE','0'))
            if not memory_mb:
                memory_mb = int(os.environ.get('SLURM_MEM_PER_CPU','0'))*int(os.environ.get('SLURM_CPUS_PER_TASK','0'))
            require_integrity(memory_mb >= 8192,'Eight GiB declared allocation required')
            require_integrity(_WORK_CONTEXT is None,'Concurrent/reentrant fork evaluation is unsupported')
            require_integrity(np.asarray(p).shape == (28,) and np.asarray(p).dtype == np.dtype('float64')
                and np.isfinite(p).all() and np.isfinite(v) and v > 0,'Finite float64 point required')
            present = sum(any(row.gene == gene for row in self.rows) for gene in range(8))
            entries = [sum(row.n+1 for row in self.rows if row.gene == gene) for gene in range(8)]
            for count in entries:
                require_integrity(count <= 100000,
                    'Prepared per-gene count-axis scope exceeded')
            estimated_overhead = 9*256*1024**2
            estimated_numerical = 512*sum(entries)
            require_integrity(estimated_overhead+estimated_numerical <= 6*1024**3,
                           'Conservative small-synthetic workspace estimate exceeded')
            self.last_parallel_capacity_estimate = dict(workers=8,count_axis_entries=entries,
                interpreter_library_parent_overhead_estimate=estimated_overhead,
                numerical_workspace_estimate=estimated_numerical,
                combined_declared_estimate=estimated_overhead+estimated_numerical,
                not_measured_or_certified=True)
            gene_settings = replace(settings,
                logmass_tolerance=min(settings.logmass_tolerance,1e-8/max(1,present)),
                gradient_tolerance=min(settings.gradient_tolerance,1e-7/max(1,present)),
                hessian_tolerance=min(settings.hessian_tolerance,1e-5/max(1,present)))
            readonly = np.asarray(p).copy()
            readonly.setflags(write=False)
            signatures = row_signatures(self.rows,np)
            began = time.perf_counter()
            call = getattr(self,'marginal_call_count',0)+1
            self.marginal_call_count = call
            self.last_integration_failure = None
            self.last_integration_reports = []
            self.last_total_error_allowance = None
            self.last_evaluated_point = dict(call=call,parameters=readonly.copy(),
                                            variance=float(v),with_hessian=with_hessian)
            _WORK_CONTEXT = (self,readonly,float(v),with_hessian,private,kernel,gl,conditional,kernel_path,gene_settings)
            results = []
            try:
                with ProcessPoolExecutor(max_workers=8,mp_context=multiprocessing.get_context('fork')) as pool:
                    for record in pool.map(_worker_gene,range(8)):
                        results.append(record)
            except Exception as exc:
                self.last_parallel_worker_records = results
                self.last_integration_failure = dict(call=call,parameters=readonly.copy(),variance=float(v),
                    failure_phase='fork_pool_infrastructure',error_type=type(exc).__name__,message=str(exc),
                    collected_workers=results,unavailable_worker_genes=list(range(len(results),8)),
                    traceback=traceback.format_exc(),no_full_sum_formed=True)
                raise
            finally:
                _WORK_CONTEXT = None
            try:
                if not (row_signatures(self.rows,np) == signatures
                    and readonly.tobytes() == np.asarray(p).tobytes() and not readonly.flags.writeable
                    and original_callables == (private.bb_terms,private.conditional_terms,
                        private.Model.predictors,private.Model.evaluate)):
                    self.last_integration_failure = dict(call=call,variance=float(v),
                        parameters=readonly.copy(),failure_phase='fork_parent_integrity',
                        all_worker_records=results,no_full_sum_formed=True)
                    raise ModelIntegrityError('Fork parent state changed')
                self.last_parallel_worker_records = results
                first = first_failed_worker(results)
                if first is not None:
                    prefix_stop=min(r['gene'] for r in results if not r['success'])
                    self.last_integration_reports = [r['report'] for r in results
                        if r['gene'] < prefix_stop and r['present'] and r['success']]
                    self.last_integration_failure = dict(call=call,gene=first['gene'],
                        failure_phase='positive_worker',error_type=first['error_type'],
                        message=first['message'],
                        parameters=readonly.copy(),variance=float(v),all_worker_records=results,
                        refused_gene_ladders=first['refused_gene_ladders'],no_full_sum_formed=True)
                    if gene_terms is not None:
                        for record in results:
                            if record['gene'] >= prefix_stop:
                                break
                            gene_terms[record['gene']] = record.get('logmass',0.)
                    if first['linalg_error']:
                        raise np.linalg.LinAlgError(first['message'])
                    if first['numerical_integration_refusal']:
                        exc = kernel.IntegrationFailure(first['message'])
                        exc.integration_reports = first['refused_gene_ladders']
                        exc.likelihood_evaluations = first['refused_gene_evaluations']
                        raise exc
                    raise ModelIntegrityError('Worker failure: '+first['error_type']+': '+first['message'])
                total,gradient = 0.,np.zeros(29)
                hessian = np.zeros((29,29)) if with_hessian else None
                pending_gene_terms = {}
                for record in results:
                    gene = record['gene']
                    if not record['present']:
                        pending_gene_terms[gene] = 0.
                        continue
                    total += record['logmass']
                    gradient += record['gradient']
                    if with_hessian:
                        hessian += record['hessian']
                    pending_gene_terms[gene] = record['logmass']
                    self.last_integration_reports.append(record['report'])
                reports = self.last_integration_reports
                accepted = [r for r in results if r['present']]
                errors = dict(logmass=sum(r['error_allowances']['logmass'] for r in accepted),
                    gradient=sum((r['error_allowances']['gradient'] for r in accepted),np.zeros(29)),
                    hessian=(sum((r['error_allowances']['hessian'] for r in accepted),np.zeros((29,29)))
                        if with_hessian else None),empirical_estimates_not_certified=True,
                    local_error_arrays_embedded_before_sum=True)
                self.last_total_error_allowance = errors
                if not (np.isfinite(errors['logmass']) and errors['logmass'] >= 0
                        and np.isfinite(errors['gradient']).all() and np.all(errors['gradient'] >= 0)
                        and (not with_hessian or (np.isfinite(errors['hessian']).all()
                                                  and np.all(errors['hessian'] >= 0)))):
                    self.last_integration_failure = dict(call=call,variance=float(v),parameters=readonly.copy(),
                        failure_phase='full_model_error_report_integrity',error_report=errors,
                        all_worker_records=results,full_sum_formed_but_not_returned=True)
                    raise ModelIntegrityError('Invalid nonnegative accumulated numerical error report')
                if not (errors['logmass'] <= 1e-8 and np.max(errors['gradient']) <= 1e-7
                        and (not with_hessian or np.max(errors['hessian']) <= 1e-5)):
                    self.last_integration_failure = dict(call=call,variance=float(v),parameters=readonly.copy(),
                        failure_phase='full_model_numerical_allowance',error_report=errors,
                        all_worker_records=results,full_sum_formed_but_not_returned=True)
                    exc = kernel.IntegrationFailure('Accumulated full-model numerical allowance exceeded')
                    exc.integration_reports = reports
                    exc.likelihood_evaluations = sum(r['likelihood_evaluations'] for r in reports)
                    raise exc
                if event_sink:
                    event_sink(dict(event='marginal_evaluation',call=call,variance=float(v),
                        with_hessian=with_hessian,seconds=time.perf_counter()-began,loglik=float(total),
                        workers=8,integration_genes=len(reports),parent_original_BB_unchanged=True))
                if gene_terms is not None:
                    for gene,value in pending_gene_terms.items():
                        gene_terms[gene] = value
                return (float(total),gradient,hessian) if with_hessian else (float(total),gradient)

            except Exception as exc:
                if self.last_integration_failure is None:
                    self.last_integration_failure = dict(call=call,variance=float(v),parameters=readonly.copy(),
                        failure_phase='parent_post_worker_integrity_or_aggregation',
                        error_type=type(exc).__name__,message=str(exc),all_worker_records=results,
                        completed_reports=self.last_integration_reports,no_successful_return=True)
                raise

    for name in INHERITED:
        require_integrity(getattr(GammaParallelModel,name) is getattr(private.Model,name),
                       'Original inherited method changed: '+name)
    return GammaParallelModel
