"""Gene-local derivative moments; no fitting/integration or production changes."""
import os
import sys


def compute_guard():
    if not os.environ.get('SLURM_JOB_ID','').isdigit():
        raise RuntimeError('Synthetic local moments require compute')


class LocalGeneMoments:
    def __init__(self, model, p, rows, kernel):
        compute_guard()
        import numpy as np
        self.np,self.model,self.kernel = np,model,kernel
        self.p = np.asarray(p).copy()
        self.rows = tuple(rows)
        if not self.rows or len({r.gene for r in self.rows}) != 1:
            raise ValueError('Exactly one nonempty gene required')
        if self.p.shape != (model.size,) or not np.isfinite(self.p).all():
            raise ValueError('Finite original coordinate vector required')
        self.g,self.G,self.Z = self.rows[0].gene,model.G,model.Z
        g,G,Z = self.g,self.G,self.Z
        self.active = np.array(sorted([g,G+g,2*G,*range(2*G+1,2*G+1+Z),
                                      2*G+1+Z,*range(model.size-3,model.size)]),dtype=int)
        if len(self.active) != Z+7 or len(set(self.active)) != Z+7:
            raise ValueError('Original active-coordinate count differs')
        self.index = {int(global_index):local_index for local_index,global_index in enumerate(self.active)}
        self.augmented = np.r_[self.active,model.size]
        self.signatures = tuple(tuple(vars(row).items()) for row in self.rows)

    def validate(self, p):
        if not self.np.array_equal(p,self.p) or tuple(tuple(vars(r).items()) for r in self.rows) != self.signatures:
            raise ValueError('Local fixed parameters or rows changed')

    def predictors(self, row, p, u=0., need_hessian=True):
        """Archived original Model.predictors formula, active indices only."""
        self.validate(p)
        np,g,G,Z,k = self.np,self.g,self.G,self.Z,2*self.G
        if row.gene != g or not any(row is r for r in self.rows):
            raise ValueError('Row outside fixed gene context')
        factor = 1+p[k]*row.stage+np.dot(p[k+1:k+1+Z],row.z)
        eta = p[g]*factor-row.orientation*p[k+1+Z]*row.stage+u*row.stage
        r = p[G+g]+np.dot(p[-3:],(row.stage,row.background,row.duplication))
        width = len(self.active)
        jac = np.zeros((2,width))
        index = self.index
        jac[0,index[g]] = factor
        jac[0,index[k]] = p[g]*row.stage
        delta = [index[j] for j in range(k+1,k+1+Z)]
        jac[0,delta] = p[g]*np.asarray(row.z)
        jac[0,index[k+1+Z]] = -row.orientation*row.stage
        jac[1,index[G+g]] = 1
        jac[1,[index[j] for j in range(self.model.size-3,self.model.size)]] = [row.stage,row.background,row.duplication]
        eta_hess = None
        if need_hessian:
            eta_hess = np.zeros((width,width))
            eta_hess[index[g],index[k]] = eta_hess[index[k],index[g]] = row.stage
            eta_hess[index[g],delta] = row.z
            eta_hess[delta,index[g]] = row.z
        return eta,r,jac,eta_hess

    def bounds(self):
        """Same row order, cp/hp arithmetic and max values as dense kernel."""
        np = self.np
        cp = np.zeros(len(self.active))
        hp = np.zeros((len(self.active),len(self.active)))
        for row in self.rows:
            _,_,jac,eta_hess = self.predictors(row,self.p,0.,True)
            j = np.abs(jac).sum(axis=0)
            cp += 2*row.n*j
            hp += (16*row.n+8*row.n**2)*np.outer(j,j)+2*row.n*np.abs(eta_hess)
        return cp,hp

    def point(self, p, u, variance, conditional, with_hessian=True):
        self.validate(p)
        np = self.np
        if not np.isfinite(u) or not np.isfinite(variance) or variance <= 0:
            raise ValueError('Finite positive variance point required')
        size = len(self.active)
        ll,gp = 0.,np.zeros(size)
        hh = np.zeros((size+1,size+1)) if with_hessian else None
        for row in self.rows:
            eta,r,jac,eta_hess = self.predictors(row,p,u,with_hessian)
            value,score,row_hess = conditional(row,eta,r,need_hessian=with_hessian)
            ll += value
            gp += jac.T @ score
            if with_hessian:
                hh[:-1,:-1] += jac.T @ row_hess @ jac+score[0]*eta_hess
        z = u/np.sqrt(variance)
        gradient = np.r_[gp,(z*z-1)/(2*variance)]
        if with_hessian:
            hh[-1,-1] = (.5-z*z)/variance**2
        second = hh+np.outer(gradient,gradient) if with_hessian else None
        return ll-z*z/2-.5*np.log(2*np.pi),gradient,hh,second

    def embed_vector(self, value, augmented=True):
        result = self.np.zeros(self.model.size+int(augmented))
        result[self.augmented if augmented else self.active] = value
        return result

    def embed_matrix(self, value, augmented=True):
        index = self.augmented if augmented else self.active
        result = self.np.zeros((self.model.size+int(augmented),)*2)
        result[self.np.ix_(index,index)] = value
        return result


def compare_cached_predictors(model, p, rows):
    """Immutable reconstruction hook; compare local vs original dense points."""
    compute_guard()
    import numpy as np
    import time
    kernel = sys.modules['repair_p3_private_kernel']
    original = sys.modules[model.evaluate.__func__.__module__]
    batch_module = sys.modules['repair_p3_local_validated_batch']
    rows = tuple(rows)
    if len(rows)!=24 or model.G!=8 or model.Z!=7 or model.size!=28 or {r.gene for r in rows}!={2}:
        raise ValueError('Expected exact synthetic G8/P28/Z7/gene2/24-row context')
    # Reconstruction saves the actual immutable variance in its output protocol;
    # the checker passes it via a module-local hook, never as an optimizer start.
    variance = SAVED_VARIANCE
    lower_mass = SAVED_LOGMASS_ESTIMATE
    shapes = [(model,np.asarray(p).copy(),'original_G8')]
    scaled = original.Model(rows,gene_count=379)
    q = np.zeros(scaled.size)
    q[:model.G] = p[:model.G]
    q[scaled.G:scaled.G+model.G] = p[model.G:2*model.G]
    q[2*scaled.G:] = p[2*model.G:]
    if scaled.size!=770:
        raise ValueError('Dimension-only P770 mapping differs')
    shapes.append((scaled,q,'dimension_only_G379'))
    records,all_pass = [],True
    for current,parameters,label in shapes:
        # Dense arrays are allocated sequentially at one point, not per-row cached.
        estimate = 64*8*(current.size+1)**2+1024*64*3001+8*24*2*current.size
        if estimate >= 1024**3 or max(r.n for r in rows)>3000:
            raise ValueError('Declared less-than1GiB synthetic workspace exceeded')
        local = LocalGeneMoments(current,parameters,rows,kernel)
        batch = batch_module.BatchedConditional(current,parameters,rows)
        def callback(row,eta,r,need_hessian=True):
            values,scores,hess = batch.terms_batch(row,np.array([eta]),r,parameters,need_hessian)
            return float(values[0]),scores[0],hess[0] if need_hessian else None
        cp,hp = local.bounds()
        dense_cp,dense_hp = kernel.derivative_bounds(current,parameters,rows)
        bounds_pass = np.array_equal(local.embed_vector(cp,False),dense_cp) and np.array_equal(local.embed_matrix(hp,False),dense_hp)
        all_pass &= bounds_pass
        for with_hessian in (False,True):
            for u in (-.5,0.,.5):
                started = time.perf_counter()
                dense_ll,dense_g,dense_h = current.evaluate(parameters,rows=rows,u=u,need_hessian=with_hessian)
                z = u/np.sqrt(variance)
                dense_gradient = np.r_[dense_g,(z*z-1)/(2*variance)]
                dense_hessian = np.zeros((current.size+1,current.size+1)) if with_hessian else None
                if with_hessian:
                    dense_hessian[:-1,:-1] = dense_h
                    dense_hessian[-1,-1] = (.5-z*z)/variance**2
                dense_second = dense_hessian+np.outer(dense_gradient,dense_gradient) if with_hessian else None
                dense_log = dense_ll-z*z/2-.5*np.log(2*np.pi)
                dense_seconds = time.perf_counter()-started
                started = time.perf_counter()
                actual = local.point(parameters,u,variance,callback,with_hessian)
                local_seconds = time.perf_counter()-started
                gradient_delta = local.embed_vector(actual[1])-dense_gradient
                hessian_delta = local.embed_matrix(actual[2])-dense_hessian if with_hessian else None
                second_delta = local.embed_matrix(actual[3])-dense_second if with_hessian else None
                errors = [abs(actual[0]-dense_log),float(np.max(np.abs(gradient_delta))),
                          float(np.max(np.abs(hessian_delta))) if with_hessian else None,
                          float(np.max(np.abs(second_delta))) if with_hessian else None]
                passed = bounds_pass and errors[0]<=1e-9 and errors[1]<=1e-8 and (not with_hessian or max(errors[2:])<=1e-5)
                tail_records = []
                for bound in (24.,26.):
                    dense_tail = kernel.tail_moments(bound,variance,lower_mass,dense_gradient,dense_second,dense_cp,dense_hp,with_hessian)
                    small_tail = kernel.tail_moments(bound,variance,lower_mass,actual[1],actual[3],cp,hp,with_hessian)
                    gd = local.embed_vector(small_tail['relative_gradient_tail_allowance'])-dense_tail['relative_gradient_tail_allowance']
                    hd = local.embed_matrix(small_tail['relative_hessian_tail_allowance'])-dense_tail['relative_hessian_tail_allowance'] if with_hessian else None
                    tail_pass = np.max(np.abs(gd))<=1e-8 and (not with_hessian or np.max(np.abs(hd))<=1e-5)
                    passed &= bool(tail_pass)
                    tail_records.append(dict(bound=bound,gradient_max_difference=float(np.max(np.abs(gd))),
                        hessian_max_difference=float(np.max(np.abs(hd))) if with_hessian else None,
                        mass_allowance_exact=small_tail['relative_mass_tail_allowance']==dense_tail['relative_mass_tail_allowance'],
                        error_allowances_conditional_on_estimated_mass=True))
                all_pass &= passed
                def sparse_deltas(matrix):
                    if matrix is None:return None
                    i,j = np.nonzero(matrix)
                    return [[int(a),int(b),float(matrix[a,b])] for a,b in zip(i,j)]
                records.append(dict(shape=label,G=current.G,P=current.size,Z=current.Z,u=u,variance=variance,
                    with_hessian=with_hessian,active_indices=local.active.tolist(),augmented_indices=local.augmented.tolist(),
                    errors=errors,within_unchanged_moment_allowances=bool(passed),dense_seconds=dense_seconds,
                    local_seconds=local_seconds,estimated_workspace_bytes=estimate,estimate_not_measured=True,
                    gradient_difference=gradient_delta.tolist(),all_nonzero_hessian_differences=sparse_deltas(hessian_delta),
                    all_nonzero_second_moment_differences=sparse_deltas(second_delta),tail_comparisons=tail_records,
                    cp_hp_embeddings_exact=bounds_pass,maxcp_exact=float(np.max(cp))==float(np.max(dense_cp)),
                    maxhp_exact=float(np.max(hp))==float(np.max(dense_hp)),
                    no_supplied_rows_depths_or_active_parameters_dropped=True))
        changed = parameters.copy();changed[0]+=.01
        rejected = False
        try:local.point(changed,0.,variance,callback)
        except ValueError:rejected=True
        if not rejected:raise ValueError('Local parameter mutation not rejected')
    return dict(comparison_records_completed=True,all_comparisons_within_unchanged_allowances=bool(all_pass),
        all_comparison_failures_retained=True,records=records,record_count=len(records),fixed_u=[-.5,0.,.5],
        actual24_synthetic_rows_unchanged=True,dimension_only_G379_not_fullscale_validation=True,
        no_integration=True,no_fits=True,no_calibration=True,private_sealed_data_not_used=True)
