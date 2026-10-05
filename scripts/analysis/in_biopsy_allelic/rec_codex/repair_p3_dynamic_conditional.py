"""Small-synthetic scalar selected BB law with an exact latest-r cache.

The immutable v2 conditional arithmetic is reused directly. Its stored p is
only an initialization seed here: this explicit four-argument interface owns
dispersion invalidation, and the conditional law depends on row, eta and r.
Eta prefixes remain local to each call. No fit or integration is performed.
"""
import hashlib
import importlib.util
import os
from pathlib import Path
import struct
import sys
import time
import traceback

V2_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
MAX_ROWS = 256
MAX_GENES = 8
MAX_GENE_ENTRIES = 100000
MEMORY_CAP = 1024**3
CHECK_ALLOWANCE = 1e-10
INHERITED = ('fit','fit_single','constrained_polish','fit_marginal','fit_marginal_single',
             'bounds','initial','cold_start','informative','predictors','evaluate','sandwich')


def require(ok, message):
    if not ok:
        raise ValueError(message)


class DynamicConditional:
    def __init__(self, model, seed_p, rows):
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required')
        import numpy as np
        rows = tuple(rows)
        require(0<len(rows)<=MAX_ROWS and len({id(row) for row in rows})==len(rows)
                and 0<model.G<=MAX_GENES,'Small synthetic registered-row scope required')
        per_gene = [0]*model.G
        for row in rows:
            require(all(int(v)==v for v in (row.n,row.a,row.depth,row.minor,row.percent))
                    and row.n>=0 and row.depth>=0 and row.minor>=0 and 0<=row.percent<=50
                    and row.orientation in (-1,1) and row.gene in range(model.G),
                    'Original integer support/orientation required')
            per_gene[row.gene] += int(row.n)+1
        require(max(per_gene)<=MAX_GENE_ENTRIES,'Per-gene count-axis allowance exceeded')
        total = sum(per_gene)
        static_bytes = 64*total*8+1048576
        working_bytes = 64*max(row.n+1 for row in rows)*8+1048576
        require(static_bytes+working_bytes<MEMORY_CAP,'Conservative cache/work estimate exceeds1GiB')
        source = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
        require(hashlib.sha256(source.read_bytes()).hexdigest()==V2_SHA,'Immutable v2 source differs')
        spec = importlib.util.spec_from_file_location('_repair_p3_dynamic_private_v2',source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.prepared = module.PreparedConditional(model,seed_p,rows)
        self.np = np
        self.original = self.prepared.original
        self.model = model
        self._r_keys = {key:struct.pack('!d',float(d['r'])) for key,d in self.prepared.rows.items()}
        self._counts = {key:dict(calls=0,r_rebuilds=0,r_hits=0) for key in self.prepared.rows}
        self.initial_layout = self.layout()
        self.memory_estimate = dict(registered_rows=len(rows),per_gene_count_axis_entries=per_gene,
            static_count_axis_entries=total,estimated_static_bytes=static_bytes,
            estimated_scalar_working_bytes=working_bytes,estimated_cache_plus_work_bytes=static_bytes+working_bytes,
            allocation_estimate_cap_bytes=MEMORY_CAP,estimate_not_measured_process_peak=True,
            persistent_eta_cache=False,persistent_r_entries=len(rows),
            bounds='rows<=256; genes<=8; per-gene sum(n+1)<=100000; estimate<1GiB')

    def layout(self):
        return (tuple(self.prepared.rows),tuple(tuple(d['supports']) for d in self.prepared.rows.values()))

    def terms(self, row, eta, linear_r, need_hessian=True):
        np = self.np
        d = self.prepared.rows.get(id(row))
        require(d is not None and d['row'] is row and d['signature']==tuple(vars(row).items()),
                'Registered row changed or foreign row')
        require(np.asarray(eta).shape==() and np.isfinite(eta)
                and np.asarray(linear_r).shape==() and np.isfinite(linear_r),
                'Finite scalar eta and dispersion required')
        eta,linear_r = float(eta),float(linear_r)
        key = struct.pack('!d',linear_r)
        stats = self._counts[id(row)]
        if key!=self._r_keys[id(row)]:
            expit = self.prepared.expit
            vr,qr = expit(linear_r),expit(-linear_r)
            rho = 1e-6+(1-2e-6)*vr
            rp = (1-2e-6)*vr*qr
            rpp = rp*(qr-vr)
            t = (1-rho)/rho
            tp = -rp/rho**2
            tpp = -rpp/rho**2+2*rp**2/rho**3
            logt = np.log(t)
            rt,dt,qt = self.prepared.rising(d['n'],logt,d['logj'])
            # Prepare all fields before replacing the row's single latest-r entry.
            replacement = dict(r=linear_r,logt=logt,rt=rt[d['n']],dt=dt[d['n']],qt=qt[d['n']],
                               u=tp/t,w=tpp/t-(tp/t)**2)
            require(all(np.isfinite(v) for v in replacement.values()),'Nonfinite rebuilt dispersion state')
            d.update(replacement)
            self._r_keys[id(row)] = key
            stats['r_rebuilds'] += 1
        else:
            stats['r_hits'] += 1
        stats['calls'] += 1
        # v2 creates eta prefixes once and reuses them across all exact supports.
        answer = self.prepared.terms(row,eta,linear_r,self.prepared.p,need_hessian)
        require(all(np.isfinite(v).all() for v in answer if v is not None),'Nonfinite conditional output')
        require(self.layout()==self.initial_layout,'Static row/support layout changed')
        return answer

    def ledger(self):
        return [dict(row_ordinal=i,**self._counts[key]) for i,key in enumerate(self.prepared.rows)]


def install_private(context):
    """Install only on this isolated model's module; preserve method identities."""
    module = sys.modules[context.model.evaluate.__func__.__module__]
    require(module is context.original,'Private conditional module differs')
    identities = {name:getattr(type(context.model),name) for name in INHERITED}
    previous = module.conditional_terms
    module.conditional_terms = context.terms
    require(all(getattr(type(context.model),name) is fn for name,fn in identities.items()),
            'Inherited model method identity changed')
    return previous


def compare_cached_predictors(model, p, rows):
    """Bounded immutable reconstruction hook; retain every point and refusal."""
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required')
    import numpy as np
    from dataclasses import replace
    rows = tuple(rows)
    began = time.perf_counter()
    context = DynamicConditional(model,p,rows)
    setup_seconds = time.perf_counter()-began
    original = context.original.conditional_terms
    original_bb = context.original.bb_terms
    methods = {name:getattr(type(model),name) for name in INHERITED}
    scenarios = (('required_three',np.array([-.5,0.,.5])),
                 ('fixed_grid',np.linspace(-2.,2.,65)),
                 ('prespecified_tails',np.array([-20.,-10.,-5.,5.,10.,20.])))
    offsets = (0.,-.2,0.,.2,0.)
    records,timings = [],[]
    exception_index,exception_traces = {},[]
    maxima = np.zeros(3)
    def retained(value):
        if isinstance(value,np.ndarray): return retained(value.tolist())
        if isinstance(value,(tuple,list)): return [retained(v) for v in value]
        if isinstance(value,(float,np.floating)) and not np.isfinite(value): return str(value)
        return value
    for need_hessian in (False,True):
        for label,us in scenarios:
            ref_seconds,dynamic_seconds = 0.,0.
            for i,row in enumerate(rows):
                for u in us:
                    eta,r = model.predictors(row,p,float(u),need_hessian=False)[:2]
                    for ordinal,offset in enumerate(offsets):
                        rr = float(r+offset)
                        ref,actual,errors,exact,exception = None,None,None,None,None
                        before_key = context._r_keys[id(row)]
                        before_rebuilds = context._counts[id(row)]['r_rebuilds']
                        try:
                            start = time.perf_counter()
                            ref = original(row,float(eta),rr,need_hessian)
                            ref_seconds += time.perf_counter()-start
                            require(all(np.isfinite(a).all() for a in ref if a is not None),
                                    'Nonfinite original conditional output')
                            start = time.perf_counter()
                            actual = context.terms(row,float(eta),rr,need_hessian)
                            dynamic_seconds += time.perf_counter()-start
                            errors = [None if a is None else float(np.max(np.abs(np.asarray(a)-b)))
                                      for a,b in zip(ref,actual)]
                            exact = [a is b if a is None else bool(np.array_equal(a,b)) for a,b in zip(ref,actual)]
                            for j,error in enumerate(errors):
                                if error is not None: maxima[j] = max(maxima[j],error)
                        except Exception:
                            trace = traceback.format_exc()
                            if trace not in exception_index:
                                exception_index[trace] = len(exception_traces)
                                exception_traces.append(trace)
                            exception = exception_index[trace]
                        current_key = struct.pack('!d',rr)
                        transition = (context._r_keys[id(row)]==current_key and
                            context._counts[id(row)]['r_rebuilds']-before_rebuilds==int(before_key!=current_key))
                        within = exception is None and all(e is None or e<=CHECK_ALLOWANCE for e in errors)
                        records.append([label,i,float(u),float(eta),rr,ordinal,need_hessian,
                                        retained(ref),retained(actual),errors,exact,within,exception,transition])
            timings.append(dict(scenario=label,need_hessian=need_hessian,
                original_seconds=ref_seconds,warm_dynamic_seconds=dynamic_seconds,
                setup_seconds=setup_seconds,rows=len(rows),u_nodes=len(us),r_sequence=list(offsets),
                single_run_not_end_to_end_speed=True))
    row = rows[0]
    eta,r = model.predictors(row,p,0.,need_hessian=False)[:2]
    refusals = []
    def refusal(label, candidate_row, candidate_eta, candidate_r):
        rejected,error = False,None
        try: context.terms(candidate_row,candidate_eta,candidate_r)
        except ValueError as exc: rejected,error = True,str(exc)
        except Exception: error = traceback.format_exc()
        refusals.append(dict(case=label,rejected=rejected,error=error))
    refusal('foreign_equal_row',replace(row),eta,r)
    old_omega = row.omega
    try:
        object.__setattr__(row,'omega',old_omega+.125)
        refusal('same_row_field_mutation',row,eta,r)
    finally:
        object.__setattr__(row,'omega',old_omega)
    for value in (np.nan,np.inf,-np.inf):
        refusal('nonfinite_eta_'+str(value),row,value,r)
        refusal('nonfinite_r_'+str(value),row,eta,value)
    refusal('nonscalar_eta',row,np.array([eta]),r)
    refusal('nonscalar_r',row,eta,np.array([r]))
    # Verify changed full parameters produce new eta/r through unchanged predictors.
    # Only the registered gene rows are evaluated; no other cases are introduced.
    chain_records = []
    changed = np.asarray(p).copy()
    changed[model.G+row.gene] += .2
    changed[row.gene] += .1
    for label,q in (('seed',p),('changed_eta_and_r',changed),('returned_seed',p)):
        ref,actual,differences,error = None,None,None,None
        try:
            ref = model.evaluate(q,rows=rows)
            previous = install_private(context)
            try:
                actual = model.evaluate(q,rows=rows)
            finally:
                context.original.conditional_terms = previous
            require(all(np.isfinite(a).all() for a in (*ref,*actual)),'Nonfinite model evaluate comparison')
            differences = [float(np.max(np.abs(np.asarray(a)-b))) for a,b in zip(ref,actual)]
        except Exception:
            error = traceback.format_exc()
        chain_records.append(dict(case=label,parameters=np.asarray(q).copy(),original=retained(ref),
                                  dynamic=retained(actual),max_absolute_differences=differences,exception=error,
                                  original_finite_masks=None if ref is None else
                                      [np.isfinite(a) for a in ref],
                                  dynamic_finite_masks=None if actual is None else
                                      [np.isfinite(a) for a in actual],
                                  no_new_full_model_hessian_threshold=True))
    valid_after = context.terms(row,eta,r)
    unchanged = (context.original.conditional_terms is original and context.original.bb_terms is original_bb
                 and all(getattr(type(model),name) is fn for name,fn in methods.items()))
    return dict(comparison_columns=['scenario','row_ordinal','u','eta','linear_r','r_sequence_ordinal',
        'need_hessian','original','dynamic','max_absolute_differences','exact_equal','within_allowance','exception_trace_index',
        'exact_r_cache_transition_matches'],
        conditional_comparisons=records,timings=timings,setup_seconds=setup_seconds,
        exception_traces=exception_traces,exception_retention='All point failures reference exact interned trace; no truncation',
        nonfinite_value_encoding='Explicit nan/inf strings plus exception; model chain also retains finite masks',
        max_absolute_differences=maxima.tolist(),comparison_allowance=CHECK_ALLOWANCE,
        comparison_records_completed=True,outside_allowance_count=sum(not r[11] for r in records),
        exception_count=sum(r[12] is not None for r in records),
        all_comparisons_within_allowance=all(r[11] for r in records),
        all_exact_r_cache_transitions_match=all(r[13] for r in records),
        guard_checks=refusals,all_required_refusals_observed=all(r['rejected'] for r in refusals),
        valid_context_after_refusal=all(np.isfinite(a).all() for a in valid_after if a is not None),
        original_callables_and_inherited_methods_unchanged=unchanged,
        model_evaluate_chain_comparisons=chain_records,
        memory_estimate=context.memory_estimate,dispersion_cache_ledger=context.ledger(),
        registered_row_metadata=[dict(row_ordinal=i,n=r.n,a=r.a,depth=r.depth,minor=r.minor,percent=r.percent,
            orientation=r.orientation,raw_oriented_genotype_prior=r.prior) for i,r in enumerate(rows)],
        observed_count_boundaries=dict(min_n=min(r.n for r in rows),max_n=max(r.n for r in rows),
            all_registered_rows_compared=True,no_invented_boundary_counts=True),
        r_offset_sequence=list(offsets),fixed_u_points=[-.5,0.,.5],
        deterministic_grid=dict(start=-2.,stop=2.,nodes=65),prespecified_tail_u=[-20.,-10.,-5.,5.,10.,20.],
        scores=2,hessian_shape=[2,2],persistent_eta_cache=False,
        completion_means='All attempted comparison records retained; not automatically accuracy accepted',
        no_integration=True,no_fits=True,no_calibration=True,
        accuracy_limit='Dynamic scalar arithmetic only; integral/fit accuracy budgets unchanged')
