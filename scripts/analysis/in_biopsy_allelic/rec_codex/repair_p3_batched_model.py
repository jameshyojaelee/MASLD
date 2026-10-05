"""Isolated small-synthetic marginal implementation with inherited v1.7 fits.

Only the new class's marginal method differs. Positive-variance numerical GL
integration is an explicit revision of the frozen 20-node approximation, not
an objective/estimand change. Fit rules and zero-boundary law are inherited.
"""
from dataclasses import replace
import hashlib
import os
from pathlib import Path
import time

KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
GL_SHA = '8bffaf8057ace5a2593861f879068f7e4a3e3f1ef72d2fecc4199b94c8372c52'
CONDITIONAL_SHA = '184b3ec43d205290e0c4e2d8df0237a9df770af38116d62762a8b1f9811f3d76'
INHERITED = ('fit','fit_single','constrained_polish','fit_marginal','fit_marginal_single',
             'bounds','initial','cold_start','informative','predictors','evaluate','sandwich')


def batched_model_class(private_original, kernel, batched_gl, conditional_module, *,
                        kernel_path, batch_size=32, event_sink=None):
    if not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise ValueError('Individual compute allocation required')
    import numpy as np
    for module, expected in ((kernel,KERNEL_SHA),(batched_gl,GL_SHA),(conditional_module,CONDITIONAL_SHA)):
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()!=expected:
            raise ValueError('Reviewed numerical dependency source differs')
    kernel_path = Path(kernel_path)
    if hashlib.sha256(kernel_path.read_bytes()).hexdigest()!=KERNEL_SHA or batch_size not in (32,64):
        raise ValueError('Require reviewed kernel and fixed32/64 node batch')
    settings = kernel.Settings()
    IntegralModel = kernel.candidate_class(private_original,settings)

    class BatchedModel(IntegralModel):
        def marginal(self,p,v,nodes=20,with_hessian=False,gene_terms=None):
            kernel.require(len(self.rows)<=256 and 0<self.G<=8,'Small synthetic model scope exceeded')
            began = time.perf_counter()
            self.marginal_call_count = getattr(self,'marginal_call_count',0)+1
            call = self.marginal_call_count
            self.last_integration_failure = None
            self.last_evaluated_point = dict(call=call,parameters=np.asarray(p).copy(),
                                            variance=float(v),with_hessian=with_hessian)
            if v<=0:
                answer = super().marginal(p,v,nodes,with_hessian,gene_terms)
                self.last_integration_reports = []
                self.last_total_error_allowance = None
                if event_sink is not None:
                    event_sink(dict(event='marginal_evaluation',call=call,variance=float(v),
                        with_hessian=with_hessian,seconds=time.perf_counter()-began,loglik=float(answer[0]),
                        max_abs_score=float(np.max(np.abs(answer[1]))),integration_genes=0,
                        inherited_exact_zero_boundary=True))
                return answer
            kernel.require(np.asarray(p).shape==(self.size,) and np.isfinite(p).all()
                           and np.isfinite(v) and v>0,'Finite model parameters and variance required')
            total,gradient = 0.,np.zeros(self.size+1)
            hessian = np.zeros((self.size+1,self.size+1)) if with_hessian else None
            reports = []
            self.last_integration_reports = reports
            self.last_total_error_allowance = None
            present = sum(any(row.gene==gene for row in self.rows) for gene in range(self.G))
            # Exactly the reviewed kernel's allocation of full-model LL/score budgets.
            gene_settings = replace(settings,logmass_tolerance=min(settings.logmass_tolerance,1e-8/max(1,present)),
                                    gradient_tolerance=min(settings.gradient_tolerance,1e-7/max(1,present)))
            for gene in range(self.G):
                rows = [row for row in self.rows if row.gene==gene]
                if not rows:
                    if gene_terms is not None:
                        gene_terms[gene] = 0.
                    continue
                kernel.require(sum(row.n+1 for row in rows)<=100000,'Per-gene prepared support scope exceeded')
                setup_start = time.perf_counter()
                context = conditional_module.BatchedConditional(self,p,rows)
                setup_seconds = time.perf_counter()-setup_start
                kernel.require(context.original is private_original
                    and context.original.conditional_terms is private_original.conditional_terms,
                    'Prepared conditional escaped private original')
                combined_estimate = (batched_gl.GL_workspace_estimate_bytes(self,rows,batch_size)
                    +context.memory_estimate['estimated_static_bytes']
                    +context.memory_estimate['estimated_maximum_batch_working_bytes'])
                kernel.require(combined_estimate<=3*1024**3,'Existing declared3GiB workspace allowance exceeded')
                def conditional(row,eta_array,linear_r,need_hessian=True):
                    return context.terms_batch(row,eta_array,linear_r,p,need_hessian)
                try:
                    ll,g,h,report = batched_gl.integrate_gene(self,p,v,rows,conditional,
                        with_hessian=with_hessian,expected_kernel_sha=KERNEL_SHA,kernel_path=kernel_path,
                        settings=gene_settings,batch_size=batch_size)
                except ValueError as exc:
                    self.last_integration_failure = dict(call=call,gene=gene,variance=float(v),
                        parameters=np.asarray(p).copy(),completed_gene_reports=list(reports),
                        refused_gene_ladders=getattr(exc,'integration_reports',None),
                        error_type=type(exc).__name__,message=str(exc))
                    if event_sink is not None:
                        event_sink(dict(event='marginal_refusal',call=call,gene=gene,variance=float(v),
                            with_hessian=with_hessian,seconds=time.perf_counter()-began,
                            error_type=type(exc).__name__,message=str(exc)))
                    raise  # Preserve the inherited numerical failure disposition.
                report['batched_callback'] = dict(setup_seconds=setup_seconds,
                    callback_memory_estimate=context.memory_estimate,
                    combined_declared_allocation_estimate_bytes=combined_estimate,
                    estimate_not_measured_process_peak=True,adaptive_chunk_ledger=context.chunk_ledger())
                reports.append(dict(gene=gene,**report))
                total += ll;gradient += g
                if with_hessian:
                    hessian += h
                if gene_terms is not None:
                    gene_terms[gene] = ll
            errors = dict(logmass=sum(r['refined']['logmass_error_allowance']
                                     +r['logmass_refinement_difference'] for r in reports),
                gradient=sum((r['refined']['gradient_error_allowance']+r['gradient_refinement_difference']
                              for r in reports),np.zeros(self.size+1)),
                hessian=(sum((r['refined']['hessian_error_allowance']+r['hessian_refinement_difference']
                              for r in reports),np.zeros_like(hessian)) if with_hessian else None),
                empirical_refinement_and_tail_estimates_not_certified=True)
            self.last_total_error_allowance = errors
            kernel.require(errors['logmass']<=1e-8 and np.max(errors['gradient'])<=1e-7,
                           'Accumulated gene numerical uncertainty exceeds unchanged full-model allowance')
            # No new aggregate Hessian threshold; retain all per-gene H estimates.
            if event_sink is not None:
                event_sink(dict(event='marginal_evaluation',call=call,variance=float(v),
                    with_hessian=with_hessian,seconds=time.perf_counter()-began,loglik=float(total),
                    max_abs_score=float(np.max(np.abs(gradient))),integration_genes=len(reports),
                    logmass_error_estimate=float(errors['logmass']),
                    maximum_gradient_error_estimate=float(np.max(errors['gradient']))))
            return (float(total),gradient,hessian) if with_hessian else (float(total),gradient)

    for name in INHERITED:
        if getattr(BatchedModel,name) is not getattr(private_original.Model,name):
            raise ValueError('Inherited fitter/design method changed: '+name)
    return BatchedModel
