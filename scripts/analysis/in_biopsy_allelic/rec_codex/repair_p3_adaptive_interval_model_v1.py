"""Separate small synthetic model using checked adaptive pathwise integrals.

Original zero boundary and fitting methods remain. Each positive call computes
full moments on both domains, including score-only calls. No backend fallback.
"""
import multiprocessing
import os
import threading
import time

WORK = None


class IntervalRefusal(ValueError):
    pass


def gene_worker(gene):
    model, p, variance, modules, kernel, peak_module = WORK
    local = modules['repair_p3_local_moments'].LocalGeneMoments(
        model, p, tuple(r for r in model.rows if r.gene == gene), None)
    conditional = modules['repair_p3_longdouble_conditional_v1'].LongDoubleConditional(model, p, local.rows)
    try:
        scan = peak_module.scan(local, conditional, variance)
    except peak_module.PeakRefusal as exc:
        return dict(gene=gene, returned=False, refusal=str(exc), peak_scan_refused=True)
    kernel.WORK = model, p, variance, {gene: scan['peaks']}, modules
    result = kernel.worker(gene)
    result['peak_scan'] = scan
    return result


def model_class(original, modules, kernel, peak_module, event_sink=None):
    import numpy as np
    if not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise RuntimeError('Individual compute required')
    original_bb = original.bb_terms

    class AdaptiveIntervalModel(original.Model):
        def marginal(self, p, v, nodes=20, with_hessian=False, gene_terms=None):
            global WORK
            self.last_integration_failure = None
            if original.bb_terms is not original_bb:
                raise RuntimeError('Original parent arithmetic changed')
            if v <= 0:
                return super().marginal(p, v, nodes, with_hessian, gene_terms)
            if ((self.G, self.Z, self.size) != (8, 7, 28) or not 0 < len(self.rows) <= 256
                    or not all(0 < r.n <= 5000 for r in self.rows)):
                raise RuntimeError('Small synthetic scope exceeded')
            if (threading.active_count() != 1 or int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) < 8
                    or len(os.sched_getaffinity(0)) < 8):
                raise RuntimeError('Eight CPUs and single-thread parent required')
            q = np.asarray(p, dtype=float).copy()
            q.setflags(write=False)
            if q.shape != (28,) or not np.isfinite(q).all() or not np.isfinite(v):
                raise RuntimeError('Invalid marginal point')
            before, row_before = q.tobytes(), repr(self.rows)
            WORK = self, q, v, modules, kernel, peak_module
            started = time.perf_counter()
            genes = sorted({r.gene for r in self.rows})
            with multiprocessing.get_context('fork').Pool(min(8, len(genes))) as pool:
                records = pool.map(gene_worker, genes)
            if q.tobytes() != before or repr(self.rows) != row_before or original.bb_terms is not original_bb:
                raise RuntimeError('Fixed marginal inputs changed')
            self.last_interval_records = records
            if not all(r['returned'] for r in records):
                self.last_integration_failure = dict(parameters=q.tolist(), variance=float(v), records=records)
                raise IntervalRefusal('Adaptive interval comparison exceeds unchanged limits')
            errors = [sum(r['differences_plus_allowances'][i] for r in records) for i in range(3)]
            if not all(a <= b for a, b in zip(errors, (1e-8, 1e-7, 1e-5))):
                self.last_integration_failure = dict(parameters=q.tolist(), variance=float(v), records=records,
                    aggregate_domain_differences_plus_allowances=errors)
                raise IntervalRefusal('Adaptive aggregate comparison exceeds unchanged limits')
            ll = sum(r['logmass'] for r in records)
            gradient = np.sum([r['gradient'] for r in records], axis=0)
            hessian = np.sum([r['hessian'] for r in records], axis=0)
            self.marginal_call_count = getattr(self, 'marginal_call_count', 0)+1
            if gene_terms is not None:
                gene_terms.update({r['gene']: r['logmass'] for r in records})
            if event_sink:
                event_sink(dict(event='marginal_evaluation', call=self.marginal_call_count,
                    variance=float(v), with_hessian=with_hessian, full_moments_computed=True,
                    loglik=float(ll), seconds=time.perf_counter()-started,
                    node_evaluations=sum(x['evaluations'] for r in records for x in r['records']),
                    peak_counts=[len(r['peak_scan']['peaks']) for r in records],
                    aggregate_domain_differences_plus_allowances=errors))
            return (float(ll), gradient, hessian) if with_hessian else (float(ll), gradient)

    for name in ('fit', 'fit_single', 'fit_marginal', 'fit_marginal_single', 'constrained_polish',
                 'bounds', 'initial', 'cold_start', 'informative', 'predictors', 'evaluate', 'sandwich'):
        if getattr(AdaptiveIntervalModel, name) is not getattr(original.Model, name):
            raise RuntimeError('Original fitting method changed: '+name)
    return AdaptiveIntervalModel
