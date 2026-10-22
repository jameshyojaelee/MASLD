"""Separate synthetic model: fixed-normal Hermite positive integrals only.

Original likelihood at v=0 and every original fitting method are inherited.
No model selection, optimizer changes or alternate backend on refusal.
"""
import multiprocessing
import os
import threading
import time

WORK = None


class HermiteRefusal(ValueError):
    pass


def gene_worker(gene):
    model, p, v, with_hessian, local_module, path_module, conditional_module, integrator = WORK
    rows = tuple(r for r in model.rows if r.gene == gene)
    local = local_module.LocalGeneMoments(model, p, rows, None)
    context = conditional_module.LongDoubleConditional(model, p, rows)
    moments = path_module.PathwiseGeneMoments(local)
    return integrator.integrate(local, context, moments, v, with_hessian=with_hessian)


def model_class(original, local_module, path_module, conditional_module, integrator, event_sink=None):
    import numpy as np
    if not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise RuntimeError('Individual compute required')
    original_bb = original.bb_terms

    class HermiteModel(original.Model):
        def marginal(self, p, v, nodes=20, with_hessian=False, gene_terms=None):
            global WORK
            self.last_integration_failure = None
            if original.bb_terms is not original_bb:
                raise RuntimeError('Original parent arithmetic changed')
            if v <= 0:
                return super().marginal(p, v, nodes, with_hessian, gene_terms)
            if ((self.G, self.Z, self.size) != (8,7,28)
                    or not 0 < len(self.rows) <= 256
                    or not all(r.gene in range(8) and 0 < r.n <= 5000 for r in self.rows)):
                raise RuntimeError('Small synthetic scope exceeded')
            if (threading.active_count() != 1 or int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) < 8
                    or len(os.sched_getaffinity(0)) < 8):
                raise RuntimeError('Eight CPUs and single-thread parent required')
            q = np.asarray(p, dtype=float).copy()
            q.setflags(write=False)
            if q.shape != (28,) or not np.isfinite(q).all() or not np.isfinite(v):
                raise RuntimeError('Invalid marginal point')
            before, row_before = q.tobytes(), repr(self.rows)
            WORK = self, q, v, with_hessian, local_module, path_module, conditional_module, integrator
            genes = sorted({r.gene for r in self.rows})
            started = time.perf_counter()
            with multiprocessing.get_context('fork').Pool(min(8, len(genes))) as pool:
                records = pool.map(gene_worker, genes)
            if q.tobytes() != before or repr(self.rows) != row_before or original.bb_terms is not original_bb:
                raise RuntimeError('Marginal fixed inputs changed')
            self.last_hermite_records = records
            if not all(r['returned'] for r in records):
                self.last_integration_failure = dict(parameters=q.tolist(), variance=float(v), records=records)
                raise HermiteRefusal('Standard Hermite order refinement exceeds unchanged limits')
            ll = sum(r['logmass'] for r in records)
            gradient = np.sum([r['gradient'] for r in records], axis=0)
            hessian = np.sum([r['hessian'] for r in records], axis=0) if with_hessian else None
            errors = [sum(abs(r['orders'][-1]['logmass']-r['orders'][-2]['logmass']) for r in records),
                float(np.max(np.sum([np.abs(r['orders'][-1]['gradient']-r['orders'][-2]['gradient']) for r in records], axis=0))),
                float(np.max(np.sum([np.abs(r['orders'][-1]['hessian']-r['orders'][-2]['hessian']) for r in records], axis=0))) if with_hessian else 0.]
            if not all(a <= b for a,b in zip(errors,(1e-8,1e-7,1e-5))):
                self.last_integration_failure = dict(parameters=q.tolist(), variance=float(v), aggregate_order_differences=errors)
                raise HermiteRefusal('Standard Hermite aggregate order refinement exceeds unchanged limits')
            self.marginal_call_count = getattr(self, 'marginal_call_count', 0)+1
            if gene_terms is not None:
                gene_terms.update({r['gene']:r['logmass'] for r in records})
            if event_sink:
                event_sink(dict(event='marginal_evaluation', call=self.marginal_call_count, variance=float(v),
                    with_hessian=with_hessian, loglik=float(ll), seconds=time.perf_counter()-started,
                    selected_orders=[r['orders'][-1]['order'] for r in records],
                    node_evaluations=sum(r['node_evaluations'] for r in records),
                    aggregate_order_differences=errors))
            return (float(ll), gradient, hessian) if with_hessian else (float(ll), gradient)

    for name in ('fit','fit_single','fit_marginal','fit_marginal_single','constrained_polish','bounds',
                 'initial','cold_start','informative','predictors','evaluate','sandwich'):
        if getattr(HermiteModel,name) is not getattr(original.Model,name):
            raise RuntimeError('Original fitting method changed: '+name)
    return HermiteModel
