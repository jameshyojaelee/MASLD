"""Matched-domain GL order diagnostic, preserving the original point arithmetic."""
import os


def diagnostic_ladder(model,p,variance,rows,conditional,kernel,bound):
    """Compare orders at one domain; never accept an integral or fit a model."""
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise ValueError('Individual compute allocation required')
    import numpy as np
    kernel.require(variance > 0 and np.isfinite(variance) and np.isfinite(p).all() and rows,
                   'Positive finite variance/parameters required')
    kernel.require(bound in (15.,16.),'Only the three retained diagnostic domains are authorized')
    cp,hp = kernel.derivative_bounds(model,p,rows)
    width = model.size+1
    sqrtv = np.sqrt(variance)
    evaluations = 0

    def stats(z):
        nonlocal evaluations
        evaluations += 1
        kernel.require(evaluations <= 100000,'Diagnostic evaluation cap exhausted')
        gp = np.zeros(model.size)
        hh = np.zeros((width,width))
        ll = 0.
        for row in rows:
            eta,dispersion,jac,eta_hess = model.predictors(row,p,sqrtv*z,need_hessian=True)
            value,score,row_hess = conditional(row,eta,dispersion,need_hessian=True)
            kernel.require(np.isfinite(value) and value <= 1e-8 and np.isfinite(score).all() and
                           np.isfinite(row_hess).all(),'Invalid conditional value/derivatives')
            ll += value
            gp += jac.T@score
            hh[:-1,:-1] += jac.T@row_hess@jac+score[0]*eta_hess
        g = np.r_[gp,(z*z-1)/(2*variance)]
        hh[-1,-1] = (.5-z*z)/variance**2
        vector = np.r_[1.,g,(hh+np.outer(g,g)).ravel()]
        kernel.require(np.isfinite(vector).all(),'Nonfinite diagnostic moments')
        return ll-z*z/2-.5*np.log(2*np.pi),vector

    reports = []
    spacing = .25
    panels = int(np.ceil(2*bound/spacing))
    edges = np.linspace(-bound,bound,panels+1)
    for order in (16,32,64,128):
        nodes,weights = np.polynomial.legendre.leggauss(order)
        total,correction = np.zeros(1+width+width*width),np.zeros(1+width+width*width)
        shift = -np.inf
        for left,right in zip(edges[:-1],edges[1:]):
            half,middle = (right-left)/2,(right+left)/2
            for node,weight in zip(nodes,weights):
                logdensity,vector = stats(float(middle+half*node))
                exponent = logdensity+np.log(half*weight)
                if exponent > shift:
                    factor = np.exp(shift-exponent)
                    total *= factor
                    correction *= factor
                    shift = exponent
                increment = np.exp(exponent-shift)*vector-correction
                updated = total+increment
                correction = (updated-total)-increment
                total = updated
        kernel.require(total[0] > 0 and np.isfinite(total).all(),'Invalid diagnostic mass/moments')
        ll = shift+np.log(total[0])
        g = total[1:1+width]/total[0]
        second = total[1+width:].reshape(width,width)/total[0]
        h = second-np.outer(g,g)
        tail = kernel.tail_moments(bound,variance,ll,g,second,cp,hp,True)
        reports.append(dict(order=order,z_bound=bound,grid_spacing=spacing,panels=panels,
            node_evaluations=panels*order,cumulative_evaluations=evaluations,
            logmass=float(ll),gradient=g,hessian=h,
            logmass_error_allowance=float(np.log1p(tail['relative_mass_tail_allowance'])),
            gradient_error_allowance=tail['relative_gradient_tail_allowance'],
            hessian_error_allowance=tail['relative_hessian_tail_allowance'],**tail))
    comparisons = []
    for first,second in zip(reports[:-1],reports[1:]):
        gradient = second['gradient']-first['gradient']
        hessian = second['hessian']-first['hessian']
        comparisons.append(dict(first_order=first['order'],second_order=second['order'],
            logmass_difference=second['logmass']-first['logmass'],gradient_difference=gradient,
            hessian_difference=hessian,hessian_vv_difference=float(hessian[-1,-1]),
            gradient_max_absolute_difference=float(np.max(np.abs(gradient))),
            hessian_max_absolute_difference=float(np.max(np.abs(hessian))),
            empirical_logmass_plus_tail_allowance=float(abs(second['logmass']-first['logmass'])+second['logmass_error_allowance']),
            empirical_gradient_plus_tail_allowance=np.abs(gradient)+second['gradient_error_allowance'],
            empirical_hessian_plus_tail_allowance=np.abs(hessian)+second['hessian_error_allowance']))
    return dict(diagnostic_only=True,ladder=reports,adjacent_order_comparisons=comparisons,
        likelihood_evaluations=evaluations,matched_domain=True,grid_spacing=.25,
        no_additional_domain_or_grid_refinement=True,no_accuracy_acceptance=True,
        independent_derivative_accuracy_verified=False,no_model_fit=True,no_calibration=True,
        empirical_error_allowances_not_certified=True)
