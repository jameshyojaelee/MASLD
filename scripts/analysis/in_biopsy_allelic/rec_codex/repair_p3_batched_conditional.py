"""Fixed-p selected BB conditional rows at multiple nodes, no integration.

Architectural reason: one NumPy prefix pass across nodes reduces repeated Python
per-node setup while preserving the selected probability law and derivatives.
Numerical agreement and speed must be measured on the fixed synthetic cases.
"""
import hashlib
import importlib.util
import os
from pathlib import Path
import sys
import time

PREPARED_SHA = '289e93b651411042093a2b987717d87cdb5ed156ab7c736e1ac5cf820670fd4b'
MAX_NODES = 64
MAX_ENTRIES = 100000
MAX_ROWS = 256
CHECK_ALLOWANCE = 1e-10


def require(ok, message):
    if not ok:
        raise ValueError(message)


class BatchedConditional:
    def __init__(self, model, p, rows):
        require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required')
        import numpy as np
        from scipy.special import expit, logsumexp
        rows = tuple(rows)
        require(0<len(rows)<=MAX_ROWS and len({id(row) for row in rows})==len(rows),
                'Bounded distinct prepared rows required')
        # Check integer support and storage bounds before prefix/context allocation.
        for row in rows:
            require(all(int(value)==value for value in (row.n,row.a,row.depth,row.minor,row.percent))
                    and row.n>=0 and row.depth>=0 and row.minor>=0 and 0<=row.percent<=50
                    and row.orientation in (-1,1),'Original integer support/orientation required')
        static_entries = sum(int(row.n)+1 for row in rows)
        require(static_entries<=MAX_ENTRIES,'Prepared row count-axis storage allowance exceeded')
        source = Path(__file__).with_name('repair_p3_prepared_conditional_v2.py')
        require(hashlib.sha256(source.read_bytes()).hexdigest()==PREPARED_SHA,'Immutable v2 preparation source changed')
        name = '_repair_p3_batched_private_v2'
        spec = importlib.util.spec_from_file_location(name,source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.np,self.expit,self.logsumexp = np,expit,logsumexp
        self.prepared = module.PreparedConditional(model,p,rows)
        self.original = self.prepared.original
        self.memory_estimate = dict(prepared_count_axis_entries=static_entries,
            estimated_static_bytes=64*static_entries*8+1048576,
            maximum_batch_entries=MAX_ENTRIES,maximum_batch_nodes=MAX_NODES,
            estimated_maximum_batch_working_bytes=64*MAX_ENTRIES*8+1048576,
            estimate_not_measured_or_certified=True,
            allocation_guard='effective_nodes*(n+1)<=100000; requested at most64nodes; prepared sum(n+1)<=100000',
            effective_chunk_rule='min(requested_nodes,100000//(n+1)); split internally without dropping nodes')
        self.initial_layout = self.layout()
        self.row_ordinals = {id(row):index for index,row in enumerate(rows)}
        self._chunk_ledger = {}

    def layout(self):
        return (tuple(self.prepared.rows),
                tuple(tuple(d['supports']) for d in self.prepared.rows.values()))

    def rising(self, n, logshape, logj):
        np,expit = self.np,self.expit
        nodes = len(logshape)
        if n==0:
            return tuple(np.zeros((nodes,1)) for _ in range(3))
        # Nodes are rows; cumulative products run only along the integer count axis.
        value = np.empty((nodes,n))
        score = np.empty_like(value)
        curvature = np.empty_like(value)
        value[:,0],score[:,0],curvature[:,0] = logshape,1.,0.
        value[:,1:] = np.logaddexp(logshape[:,None],logj[None,:])
        score[:,1:] = expit(logshape[:,None]-logj[None,:])
        curvature[:,1:] = expit(logshape[:,None]-logj[None,:])*expit(logj[None,:]-logshape[:,None])
        return tuple(np.concatenate((np.zeros((nodes,1)),np.cumsum(a,axis=1)),axis=1)
                     for a in (value,score,curvature))

    def prefixes(self, eta, d):
        np,expit = self.np,self.expit
        m,q = expit(eta),expit(-eta)
        ra,da,qa = self.rising(d['n'],d['logt']-np.logaddexp(0.,-eta),d['logj'])
        rb,db,qb = self.rising(d['n'],d['logt']-np.logaddexp(0.,eta),d['logj'])
        return m,q,ra,da,qa,rb,db,qb

    def bb(self, k, d, prefixes, need_hessian):
        np = self.np
        m,q,ra,da,qa,rb,db,qb = prefixes
        k,choose = d['supports'][tuple(k)]
        nk = d['n']-k
        da,db,qa,qb = da[:,k],db[:,nk],qa[:,k],qb[:,nk]
        mass = choose[None,:]+ra[:,k]+rb[:,nk]-d['rt']
        m,q = m[:,None],q[:,None]
        u,w = d['u'],d['w']
        score = np.stack((da*q-db*m,(da+db-d['dt'])*u),axis=2)
        if not need_hessian:
            return mass,score,None
        hess = np.empty((*mass.shape,2,2))
        hess[:,:,0,0] = qa*q*q+qb*m*m-(da+db)*m*q
        hess[:,:,0,1] = hess[:,:,1,0] = u*(qa*q-qb*m)
        hess[:,:,1,1] = u*u*(qa+qb-d['qt'])+w*(da+db-d['dt'])
        return mass,score,hess

    def sum_terms(self, lp, score, hess):
        np = self.np
        logden = self.logsumexp(lp,axis=1)
        probabilities = np.exp(lp-logden[:,None])
        mean_score = np.einsum('bk,bki->bi',probabilities,score)
        if hess is None:
            return logden,mean_score,None
        den_hess = (np.einsum('bk,bkij->bij',probabilities,hess)
                    +np.einsum('bk,bki,bkj->bij',probabilities,score,score)
                    -np.einsum('bi,bj->bij',mean_score,mean_score))
        return logden,mean_score,den_hess

    def terms_batch(self, row, eta_nodes, linear_r, p, need_hessian=True):
        np = self.np
        self.prepared.validate(p)
        d = self.prepared.rows.get(id(row))
        require(d is not None and d['row'] is row and d['signature']==tuple(vars(row).items()),
                'Prepared row changed')
        eta = np.asarray(eta_nodes,dtype=float)
        require(eta.ndim==1 and 0<len(eta)<=MAX_NODES and np.isfinite(eta).all(),
                'Require finite one-dimensional bounded node batch')
        require(np.asarray(linear_r).shape==() and np.isfinite(linear_r) and linear_r==d['r'],
                'Prepared dispersion changed')
        effective = min(len(eta),MAX_ENTRIES//(d['n']+1))
        require(effective>=1,'Even one node exceeds count-axis allowance')
        chunks = [self._terms_chunk(row,eta[start:start+effective],linear_r,p,need_hessian)
                  for start in range(0,len(eta),effective)]
        answer = tuple(None if chunks[0][j] is None else
                       np.concatenate([chunk[j] for chunk in chunks],axis=0) for j in range(3))
        lengths = tuple(min(effective,len(eta)-start) for start in range(0,len(eta),effective))
        key = (self.row_ordinals[id(row)],len(eta),effective,lengths,bool(need_hessian))
        self._chunk_ledger[key] = self._chunk_ledger.get(key,0)+1
        require(len(self._chunk_ledger)<=MAX_ROWS*MAX_NODES*2,'Bounded chunk ledger exceeded')
        return answer

    terms = terms_batch

    def chunk_ledger(self):
        return [dict(row_ordinal=key[0],requested_nodes=key[1],effective_chunk_limit=key[2],
                     chunk_lengths=list(key[3]),need_hessian=key[4],calls=count,
                     all_requested_nodes_returned=sum(key[3])==key[1])
                for key,count in sorted(self._chunk_ledger.items())]

    def _terms_chunk(self, row, eta_nodes, linear_r, p, need_hessian=True):
        np = self.np
        self.prepared.validate(p)
        d = self.prepared.rows.get(id(row))
        require(d is not None and d['row'] is row and d['signature']==tuple(vars(row).items()),
                'Prepared row changed')
        eta = np.asarray(eta_nodes,dtype=float)
        require(eta.ndim==1 and 0<len(eta)<=MAX_NODES and np.isfinite(eta).all(),
                'Require finite one-dimensional bounded node batch')
        require(len(eta)*(d['n']+1)<=MAX_ENTRIES,'Node-by-count prefix allowance exceeded')
        require(np.asarray(linear_r).shape==() and np.isfinite(linear_r) and linear_r==d['r'],
                'Prepared dispersion changed')
        eta = eta-row.orientation*row.omega
        prefixes = self.prefixes(eta,d)
        nodes = len(eta)
        if not len(d['tails']):
            logc,score = np.zeros(nodes),np.zeros((nodes,2))
            hess = np.zeros((nodes,2,2)) if need_hessian else None
        else:
            logtail,tail_score,tail_hess = self.sum_terms(*self.bb(d['tails'],d,prefixes,need_hessian))
            complement = logtail<np.log(.5)  # The unchanged exact branch, independently at each node.
            logc,score = np.empty(nodes),np.empty((nodes,2))
            hess = np.empty((nodes,2,2)) if need_hessian else None
            if np.any(complement):
                logc[complement] = np.log1p(-np.exp(logtail[complement]))
                ratio = np.exp(logtail[complement]-logc[complement])
                cscore = -ratio[:,None]*tail_score[complement]
                score[complement] = cscore
                if need_hessian:
                    ts = tail_score[complement]
                    hess[complement] = (-ratio[:,None,None]*(tail_hess[complement]+np.einsum('bi,bj->bij',ts,ts))
                                        -np.einsum('bi,bj->bij',cscore,cscore))
            if np.any(~complement):
                selected_prefixes = tuple(array[~complement] for array in prefixes)
                direct = self.sum_terms(*self.bb(np.arange(d['lo'],d['n']-d['lo']+1),
                                                 d,selected_prefixes,need_hessian))
                logc[~complement],score[~complement] = direct[:2]
                if need_hessian:
                    hess[~complement] = direct[2]
        # Keep all three raw oriented genotype priors, with no renormalization.
        logs = d['logprior'][None,:]+np.column_stack((logc,np.full(nodes,d['hom']),np.full(nodes,d['hom'])))
        ds = np.zeros((nodes,3,2));ds[:,0,:] = score
        dh = np.zeros((nodes,3,2,2)) if need_hessian else None
        if need_hessian:
            dh[:,0,:,:] = hess
        logden,ds,dh = self.sum_terms(logs,ds,dh)
        lp,grad,hess = self.bb(np.array([row.a]),d,prefixes,need_hessian)
        parts = np.column_stack((lp[:,0]+d['logprior'][0],
                                 np.full(nodes,d['hom1'][0]+d['logprior'][1]),
                                 np.full(nodes,d['hom2'][0]+d['logprior'][2])))
        mass = self.logsumexp(parts,axis=1)
        wh = np.exp(parts[:,0]-mass)
        ns = wh[:,None]*grad[:,0,:]
        nh = (wh[:,None,None]*hess[:,0,:,:]
              +(wh*(1-wh))[:,None,None]*np.einsum('bi,bj->bij',grad[:,0,:],grad[:,0,:])
              if need_hessian else None)
        answer = mass-logden,ns-ds,nh-dh if need_hessian else None
        require(all(np.isfinite(a).all() for a in answer if a is not None),'Nonfinite batch conditional result')
        require(self.layout()==self.initial_layout,'Persistent row/support cache changed')
        return answer


def compare_cached_predictors(model, p, rows):
    """Immutable synthetic reconstruction hook; no method/global replacement."""
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Compute allocation required')
    import numpy as np
    rows = tuple(rows)
    started = time.perf_counter()
    context = BatchedConditional(model,p,rows)
    preparation_seconds = time.perf_counter()-started
    original = context.original.conditional_terms
    original_bb = context.original.bb_terms
    predictor = type(model).predictors
    scenarios = (('required_three',np.array([-.5,0.,.5])),
                 ('fixed_grid',np.linspace(-2.,2.,65)),
                 ('prespecified_tails',np.array([-20.,-10.,-5.,5.,10.,20.])))
    comparisons,timings = [],[]
    maxima = np.zeros(3)
    for need_hessian in (False,True):
        for label,u_values in scenarios:
            inputs = [(row,np.array([model.predictors(row,p,float(u),need_hessian=False)[0] for u in u_values]),
                       model.predictors(row,p,0.,need_hessian=False)[1]) for row in rows]
            started = time.perf_counter()
            references = [[original(row,float(eta),r,need_hessian) for eta in etas] for row,etas,r in inputs]
            reference_seconds = time.perf_counter()-started
            for batch_size in (32,64):
                started = time.perf_counter()
                actuals = []
                for row,etas,r in inputs:
                    chunks = [context.terms(row,etas[start:start+batch_size],r,p,need_hessian)
                              for start in range(0,len(etas),batch_size)]
                    actuals.append(tuple(None if chunks[0][j] is None else
                                         np.concatenate([chunk[j] for chunk in chunks],axis=0) for j in range(3)))
                batch_seconds = time.perf_counter()-started
                for i,(row,etas,r) in enumerate(inputs):
                    for j,(u,eta,ref) in enumerate(zip(u_values,etas,references[i])):
                        actual = tuple(None if value is None else value[j] for value in actuals[i])
                        errors = [None if a is None else float(np.max(np.abs(np.asarray(a)-b)))
                                  for a,b in zip(ref,actual)]
                        exact = [a is b if a is None else bool(np.array_equal(a,b)) for a,b in zip(ref,actual)]
                        for index,error in enumerate(errors):
                            if error is not None:
                                maxima[index] = max(maxima[index],error)
                        within = all(e is None or e<=CHECK_ALLOWANCE for e in errors)
                        comparisons.append(dict(row_ordinal=i,u=float(u),eta=float(eta),linear_r=float(r),
                            scenario=label,batch_size=batch_size,need_hessian=need_hessian,
                            original=list(ref),batched=list(actual),max_absolute_differences=errors,exact_equal=exact,
                            orientation=row.orientation,raw_oriented_genotype_prior=row.prior,
                            within_allowance=within))
                timings.append(dict(scenario=label,batch_size=batch_size,need_hessian=need_hessian,
                    original_seconds=reference_seconds,warm_batch_seconds=batch_seconds,rows=len(rows),nodes=len(u_values),
                    preparation_seconds=preparation_seconds,single_run_not_end_to_end_speed=True))
    row,etas,r = inputs[0]
    changed = np.asarray(p).copy();changed[0] += .01
    parameter_refused,row_refused,large_batch_refused = False,False,False
    try:
        context.terms(row,etas[:3],r,changed)
    except ValueError:
        parameter_refused = True
    old_omega = row.omega
    try:
        object.__setattr__(row,'omega',old_omega+.125)
        try:
            context.terms(row,etas[:3],r,p)
        except ValueError:
            row_refused = True
    finally:
        object.__setattr__(row,'omega',old_omega)
    try:
        context.terms(row,np.zeros(MAX_NODES+1),r,p)
    except ValueError:
        large_batch_refused = True
    require(parameter_refused and row_refused and large_batch_refused,'Mutation/batch refusal missing')
    context.terms(row,etas[:3],r,p)
    require(context.original.conditional_terms is original and context.original.bb_terms is original_bb
            and type(model).predictors is predictor,'Original callables changed')
    return dict(conditional_comparisons=comparisons,timings=timings,preparation_seconds=preparation_seconds,
        max_absolute_differences=maxima.tolist(),comparison_allowance=CHECK_ALLOWANCE,
        comparison_records_completed=True,
        outside_allowance_count=sum(not record['within_allowance'] for record in comparisons),
        all_comparisons_within_allowance=all(record['within_allowance'] for record in comparisons),
        completion_means='All finite comparison records retained, not automatically accuracy accepted',
        changed_parameters_rejected=parameter_refused,same_row_field_mutation_rejected=row_refused,
        too_many_nodes_rejected=large_batch_refused,valid_context_after_refusal=True,
        original_callables_unchanged=True,persistent_cache_layout_unchanged=context.layout()==context.initial_layout,
        memory_estimate=context.memory_estimate,adaptive_chunk_ledger=context.chunk_ledger(),
        fixed_u_points=[-.5,0.,.5],prespecified_tail_u=[-20.,-10.,-5.,5.,10.,20.],
        deterministic_grid=dict(start=-2.,stop=2.,nodes=65,predeclared_not_outcome_selected=True),
        batch_sizes=[32,64],rows=len(rows),scores=2,hessian_shape=[2,2],
        no_integration=True,no_fits=True,no_calibration=True,
        accuracy_limit='Fixed conditional arithmetic only; integral/fit accuracy budgets unchanged')
