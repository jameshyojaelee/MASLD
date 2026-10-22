"""Prepared G8 synthetic positive-variance GL with eight native-law workers.

The parent keeps the original BB, conditional, null-fit and v=0 callables.
Only forked positive-variance workers install the separately checked gamma
arithmetic. All fitter methods, bounds, starts, ties and restarts are inherited.
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
GL_SHA = '45b4423811fded58712bd5e85666f9ffc2b6972abbd9bd71964a584ec6d3c59c'
BB_SHA = '546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'
INHERITED = ('fit','fit_single','constrained_polish','fit_marginal','fit_marginal_single',
             'bounds','initial','cold_start','informative','predictors','evaluate','sandwich')
_WORK_CONTEXT = None


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


def _worker_gene(gene):
    model,p,v,with_hessian,private,kernel,gl,bb,kernel_path,settings = _WORK_CONTEXT
    rows = tuple(row for row in model.rows if row.gene == gene)
    if not rows:
        return dict(gene=gene,present=False,success=True)
    before = row_signatures(rows,kernel.np)
    before_parameters = p.tobytes()
    original_bb = private.bb_terms
    began = time.perf_counter()
    try:
        require_integrity(not p.flags.writeable,'Worker parameters must be read-only')
        private.bb_terms = bb.stable_bb_terms
        ll,g,h,report = gl.integrate_gene(model,p,v,rows,private.conditional_terms,
            with_hessian=with_hessian,expected_kernel_sha=KERNEL_SHA,
            kernel_path=kernel_path,settings=settings)
        if row_signatures(rows,kernel.np) != before or p.tobytes() != before_parameters or p.flags.writeable:
            raise ModelIntegrityError('Worker fixed inputs changed')
        return dict(gene=gene,present=True,success=True,logmass=ll,gradient=g,hessian=h,
            report=dict(gene=gene,**report),worker_pid=os.getpid(),
            elapsed_seconds=time.perf_counter()-began,fixed_inputs_unchanged=True)
    except Exception as exc:
        return dict(gene=gene,present=True,success=False,error_type=type(exc).__name__,
            message=str(exc),value_error=isinstance(exc,ValueError),
            numerical_integration_refusal=(isinstance(exc,ValueError)
                and hasattr(exc,'integration_reports') and hasattr(exc,'likelihood_evaluations')),
            linalg_error=isinstance(exc,kernel.np.linalg.LinAlgError),
            refused_gene_ladders=getattr(exc,'integration_reports',None),
            refused_gene_evaluations=getattr(exc,'likelihood_evaluations',None),
            traceback=traceback.format_exc(),worker_pid=os.getpid(),
            elapsed_seconds=time.perf_counter()-began)
    finally:
        private.bb_terms = original_bb


def gamma_parallel_model_class(private,kernel,gl,bb,*,kernel_path,event_sink=None):
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise ValueError('Individual compute allocation required before numerical imports')
    import numpy as np
    for module,expected in ((private,ORIGINAL_SHA),(kernel,KERNEL_SHA),(gl,GL_SHA),(bb,BB_SHA)):
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed native-law numerical source differs')
    kernel_path = Path(kernel_path)
    require_integrity(hashlib.sha256(kernel_path.read_bytes()).hexdigest() == KERNEL_SHA,
                   'Reviewed kernel path differs')
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
            _WORK_CONTEXT = (self,readonly,float(v),with_hessian,private,kernel,gl,bb,kernel_path,gene_settings)
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
                failed = [r for r in results if not r['success']]
                if failed:
                    first = failed[0]
                    self.last_integration_reports = [r['report'] for r in results
                        if r['gene'] < first['gene'] and r['present'] and r['success']]
                    self.last_integration_failure = dict(call=call,gene=first['gene'],
                        failure_phase='positive_worker',error_type=first['error_type'],
                        message=first['message'],
                        parameters=readonly.copy(),variance=float(v),all_worker_records=results,
                        refused_gene_ladders=first['refused_gene_ladders'],no_full_sum_formed=True)
                    if gene_terms is not None:
                        for record in results:
                            if record['gene'] >= first['gene']:
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
                errors = dict(logmass=sum(r['refined']['logmass_error_allowance']+
                    r['logmass_refinement_difference'] for r in reports),
                    gradient=sum((r['refined']['gradient_error_allowance']+
                        r['gradient_refinement_difference'] for r in reports),np.zeros(29)),
                    hessian=(sum((r['refined']['hessian_error_allowance']+
                        r['hessian_refinement_difference'] for r in reports),np.zeros((29,29)))
                        if with_hessian else None),empirical_estimates_not_certified=True)
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
