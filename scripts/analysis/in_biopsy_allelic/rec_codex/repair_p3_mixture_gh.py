#!/usr/bin/env python3
"""Experimental fixed-point mixture GH revision; no fitting or adoption.

All located maxima and a prior proposal component contribute through a mixture
density denominator. Finite-node/search refinement is empirical, not a certified
error bound or an unconditional guarantee of finding all modes.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


def metadata():
    return dict(method='normalized-prior mixture adaptive Gauss-Hermite',
                orders=[40,80,160,320],search_spacings=[.25,.125],
                initial_z_bound=8.,maximum_z_bound=96.,maximum_evaluations=100000,
                logmass_tolerance=1e-9,gradient_tolerance=1e-8,hessian_tolerance=1e-5,
                prior_policies=['equal','scale_sqrt','scale_precision'],default_prior_policy='equal',
                proposal_weights='positive prior weight; remaining mass equal across located maxima',
                fixed_point_only=True,all_modes_guaranteed=False,certified_error_bound=False)


def load_kernel(path, expected_sha):
    path=Path(path)
    if not expected_sha or hashlib.sha256(path.read_bytes()).hexdigest()!=expected_sha:
        raise ValueError('Reviewed integration kernel source pin differs')
    spec=importlib.util.spec_from_file_location('mixture_gh_reviewed_kernel',path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    return module


def integrate_gene(model,p,variance,rows,conditional,with_hessian=False,*,
                   expected_kernel_sha,kernel_path=None,settings=None,prior_policy='equal'):
    """Compatible four-value fixed-point API: LL, gradient, Hessian, report.

    The Gaussian GH proposal is in z=u/sqrt(v); its importance ratio uses the
    complete mixture q(z), never the single component that generated a node.
    Proposal/node derivatives are not differentiated. Prior scores are fixed-u.
    Report error_allowance fields are empirical node differences plus envelopes,
    not certified integration errors. No original model methods are patched.
    """
    if not os.environ.get('SLURM_JOB_ID','').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
        raise ValueError('Individual compute allocation required before numerical imports')
    kernel_path=Path(__file__).with_name('repair_p3_integral.py') if kernel_path is None else kernel_path
    kernel=load_kernel(kernel_path,expected_kernel_sha)
    import numpy as np
    from scipy.optimize import brentq,minimize_scalar
    from scipy.special import logsumexp
    from scipy.stats import norm
    settings=kernel.Settings() if settings is None else settings
    kernel.require(prior_policy in ('equal','scale_sqrt','scale_precision'),'Unknown prespecified prior proposal policy')
    kernel.require(0<settings.logmass_tolerance<=1e-9 and 0<settings.gradient_tolerance<=1e-8 and
                   0<settings.hessian_tolerance<=1e-5 and 0<settings.maximum_evaluations<=100000 and
                   0<settings.maximum_z_bound<=96 and settings.initial_z_bound==8 and settings.grid_spacing==.25,
                   'Mixture GH settings loosen the fixed precision/domain/evaluation budgets')
    kernel.require(np.isfinite(variance) and variance>0 and np.isfinite(p).all() and rows,
                   'Positive finite variance/parameters required')
    cp,hp=kernel.derivative_bounds(model,p,rows)
    sqrtv=np.sqrt(variance);width=model.size+1
    evaluations=0;reports=[];geometries=[];geometry_cache={};rules={}

    def refuse(message,geometry_failure=None):
        if geometry_failure is not None:
            # Put geometry failures in the existing partial-report channel too;
            # callers already preserve integration_reports on numerical refusal.
            reports.append(dict(phase='proposal_geometry_refusal',**geometry_failure))
        error=kernel.IntegrationFailure(message)
        error.integration_reports=reports
        error.geometry_reports=geometries
        error.likelihood_evaluations=evaluations
        raise error

    def stats(z,full):
        nonlocal evaluations
        z=float(z)
        if abs(z)>settings.maximum_z_bound:
            refuse('Mixture GH node exceeds unchanged standardized domain cap')
        evaluations+=1
        if evaluations>settings.maximum_evaluations:
            refuse('Mixture GH evaluation cap exhausted')
        ll,su,hu=0.,-z,-1.
        gp=np.zeros(model.size) if full else None
        hh=np.zeros((width,width)) if full and with_hessian else None
        for row in rows:
            eta,rr,jac,eta_hess=model.predictors(row,p,sqrtv*z,need_hessian=full and with_hessian)
            value,s,h=conditional(row,eta,rr,need_hessian=True)
            kernel.require(np.isfinite(value) and value<=1e-8 and np.isfinite(s).all() and np.isfinite(h).all(),
                           'Invalid conditional likelihood/derivative')
            ll+=value;su+=sqrtv*row.stage*s[0];hu+=variance*row.stage**2*h[0,0]
            if full:
                gp+=jac.T@s
                if with_hessian:hh[:-1,:-1]+=jac.T@h@jac+s[0]*eta_hess
        logdensity=ll-z*z/2-.5*np.log(2*np.pi)
        if not full:return logdensity,su,hu
        g=np.r_[gp,(z*z-1)/(2*variance)]
        if with_hessian:hh[-1,-1]=(.5-z*z)/variance**2
        kernel.require(np.isfinite(g).all() and (hh is None or np.isfinite(hh).all()),'Nonfinite moments')
        return logdensity,g,hh

    def geometry_stats(z):
        z=float(z)
        if z not in geometry_cache:geometry_cache[z]=stats(z,False)
        return geometry_cache[z]

    def locate(bound,spacing):
        grid=np.linspace(-bound,bound,int(np.ceil(2*bound/spacing))+1)
        values=[geometry_stats(z) for z in grid]
        maxima=[];stationary=[];sampled_refinements=[]
        def details(z,left,right,origin):
            value,score,curvature=geometry_stats(z)
            return dict(z=z,logdensity=value,z_score=score,z_curvature=curvature,
                        stationary_score_threshold=1e-7*(1+abs(curvature)),
                        bracket_left_z=float(left),bracket_right_z=float(right),
                        bracket_left_score=geometry_stats(left)[1],
                        bracket_right_score=geometry_stats(right)[1],
                        origin=origin,z_bound=bound,grid_spacing=spacing)
        for left,right,a,b in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
            if a[1]==0:stationary.append((float(left),float(left),float(right)))
            if a[1]*b[1]<0:
                stationary.append((float(brentq(lambda z:geometry_stats(z)[1],left,right,xtol=1e-11)),
                                   float(left),float(right)))
        if values[-1][1]==0:stationary.append((float(grid[-1]),float(grid[-2]),float(grid[-1])))
        for z,left,right in stationary:
            value,score,curvature=geometry_stats(z)
            if curvature<0:maxima.append((z,value,score,curvature,details(z,left,right,'grid_score_root')))
        for i in range(1,len(grid)-1):
            if values[i][0]>=max(values[i-1][0],values[i+1][0]):
                peak=minimize_scalar(lambda z:-geometry_stats(z)[0],bounds=(grid[i-1],grid[i+1]),
                                     method='bounded',options={'xatol':1e-10})
                left,right=float(grid[i-1]),float(grid[i+1])
                if not peak.success:
                    refuse('Bounded mixture proposal geometry search failed',
                           details(float(peak.x),left,right,'sampled_logdensity_search'))
                raw_z=float(peak.x)
                raw=details(raw_z,left,right,'sampled_logdensity_search')
                # Value optimization locates a peak, but its finite stopping
                # precision does not establish a stationary analytic score.
                # Refine in the SAME sampled window; do not loosen its threshold
                # or discard a failed proposal in favor of another peak.
                if raw['bracket_left_score']<0 or raw['bracket_right_score']>0:
                    refuse('Sampled mixture maximum lacks a positive-to-negative score bracket',raw)
                if raw['bracket_left_score']==0:z=left
                elif raw['bracket_right_score']==0:z=right
                else:z=float(brentq(lambda x:geometry_stats(x)[1],left,right,xtol=1e-11))
                refined=details(z,left,right,'sampled_peak_analytic_score_refinement')
                sampled_refinements.append(dict(raw=raw,refined=refined))
                value,score,curvature=geometry_stats(z)
                if curvature>=0:
                    refuse('Refined sampled mixture maximum has nonnegative curvature',refined)
                maxima.append((z,value,score,curvature,refined))
        maxima.sort(key=lambda x:x[0]);unique=[]
        for peak in maxima:
            z,value,score,curvature,location=peak
            # Numerical geometry criterion only; no model-fit convergence claim.
            if abs(score)>1e-7*(1+abs(curvature)):
                refuse('Located mixture peak has unresolved stationary score',location)
            # Only exact duplicate locations are merged. Nearby genuine modes
            # must not be discarded; overlapping proposals are normalized by q.
            if unique and z==unique[-1][0]:
                if abs(score)<abs(unique[-1][2]):unique[-1]=peak
            else:unique.append(peak)
        components=[dict(kind='prior',mean_z=0.,sd_z=1.)]
        for z,value,score,curvature,location in unique:
            sd=(-curvature)**-.5
            if not np.isfinite(sd) or sd<=0:refuse('Invalid finite mixture proposal width')
            components.append(dict(kind='located_maximum',mean_z=z,sd_z=float(sd),
                logdensity=value,z_score=score,z_curvature=curvature,
                location=location,
                stationary_score_threshold=1e-7*(1+abs(curvature))))
        record=dict(z_bound=bound,grid_spacing=spacing,components=components,
                    sampled_peak_analytic_score_refinements=sampled_refinements,
                    all_grid_points_searched=True,stationary_points_exhaustive=False,
                    nearby_peak_proposals_not_merged=True,
                    mode_search_not_exhaustive=True,positive_prior_proposal_component=True)
        geometries.append(record)
        return record

    def run(geometry,order):
        if order not in rules:
            nodes,weights=np.polynomial.hermite.hermgauss(order)
            kernel.require(np.isfinite(nodes).all() and np.isfinite(weights).all() and np.all(weights>0),
                           'Hermite rule nonfinite or nonpositive')
            rules[order]=(np.sqrt(2)*nodes,np.log(weights)-.5*np.log(np.pi))
        nodes,logweights=rules[order]
        components=geometry['components'];peak_count=len(components)-1
        precision_epsilon=None
        if not peak_count:
            alpha=np.array([1.])
        elif prior_policy=='equal':
            alpha=np.full(len(components),1/len(components))
        else:
            scale=min(1.,min(c['sd_z'] for c in components[1:]))
            if prior_policy=='scale_sqrt':
                a0=min(1/(peak_count+1),scale*np.sqrt(settings.logmass_tolerance))
            else:
                budgets=[settings.logmass_tolerance,
                         settings.gradient_tolerance/(1+np.max(cp))]
                if with_hessian:
                    budgets.append(settings.hessian_tolerance/(1+np.max(hp)+np.max(cp)**2))
                precision_epsilon=.01*min(budgets)
                a0=min(1/(peak_count+1),scale*precision_epsilon)
            # Invalid or underflowed proposal weights are refused, never clipped
            # to zero or replaced by another policy.
            if not np.isfinite(a0) or not 0<a0<1:
                refuse('Prespecified prior proposal weight is invalid or underflowed')
            alpha=np.r_[a0,np.full(peak_count,(1-a0)/peak_count)]
        if not np.isfinite(alpha).all() or not np.all(alpha>0) or abs(alpha.sum()-1)>1e-14:
            refuse('Mixture proposal weights are invalid or underflowed')
        logalpha=np.log(alpha)
        means=np.array([a['mean_z'] for a in components]);sds=np.array([a['sd_z'] for a in components])
        entries=[]
        for component_index,component in enumerate(components):
            for node,lw in zip(nodes,logweights):
                z=component['mean_z']+component['sd_z']*node
                ld,g,h=stats(z,True)
                lq=logsumexp(logalpha-np.log(sds)-.5*np.log(2*np.pi)-.5*((z-means)/sds)**2)
                entries.append((float(logalpha[component_index]+lw+ld-lq),g,h,component_index))
        logs=np.array([x[0] for x in entries]);shift=float(np.max(logs))
        weights=np.exp(logs-shift)
        # Center moments at the largest-weight node before compensated sums.
        anchor=entries[int(np.argmax(logs))][1]
        length=1+width+width*width if with_hessian else 1+width
        total=np.zeros(length);correction=np.zeros(length)
        component_length=1+width+2*width*width if with_hessian else 1+width
        component_totals=np.zeros((len(components),component_length))
        component_corrections=np.zeros_like(component_totals)
        for weight,(_,g,h,index) in zip(weights,entries):
            centered=g-anchor
            vector=np.r_[1.,centered,(h+np.outer(centered,centered)).ravel()] if with_hessian else np.r_[1.,centered]
            increment=weight*vector-correction;updated=total+increment
            correction=(updated-total)-increment;total=updated
            # Diagnostics retain conditional H and centered score outer products
            # separately. They do not change the global summation above.
            vector=np.r_[1.,centered,h.ravel(),np.outer(centered,centered).ravel()] if with_hessian else np.r_[1.,centered]
            increment=weight*vector-component_corrections[index]
            updated=component_totals[index]+increment
            component_corrections[index]=(updated-component_totals[index])-increment
            component_totals[index]=updated
        kernel.require(total[0]>0 and np.isfinite(total).all(),'Nonfinite mixture moment sum')
        if not np.isfinite(component_totals).all() or np.any(component_totals[:,0]<=0):
            refuse('Component importance contribution is nonfinite or underflowed')
        ll=shift+np.log(total[0]);centered_mean=total[1:1+width]/total[0];g=anchor+centered_mean
        h=(total[1+width:].reshape(width,width)/total[0]-np.outer(centered_mean,centered_mean)
           if with_hessian else None)
        second=h+np.outer(g,g) if with_hessian else None
        tail=kernel.tail_moments(geometry['z_bound'],variance,ll,g,second,cp,hp,with_hessian)
        component_reports=[]
        for index,component in enumerate(components):
            moment=component_totals[index];mass_fraction=moment[0]/total[0]
            centered_score=moment[1:1+width]/total[0]
            if not np.isfinite(mass_fraction) or mass_fraction<=0:
                refuse('Normalized component mass is invalid or underflowed')
            component_report=dict(index=index,proposal=component,proposal_weight=float(alpha[index]),
                log_scaled_mass_contribution=float(np.log(moment[0])),global_log_mass_scale=shift,
                log_mass_contribution=float(shift+np.log(moment[0])),
                posterior_normalized_mass_contribution=float(mass_fraction),
                posterior_weighted_score_contribution=centered_score+anchor*mass_fraction,
                posterior_weighted_centered_score_contribution=centered_score)
            if with_hessian:
                component_report.update(
                    posterior_weighted_conditional_H_contribution=moment[1+width:1+width+width*width].reshape(width,width)/total[0],
                    posterior_weighted_centered_score_outer_contribution=moment[1+width+width*width:].reshape(width,width)/total[0])
            component_reports.append(component_report)
        report=dict(z_bound=geometry['z_bound'],grid_spacing=geometry['grid_spacing'],order=order,
                    proposal_geometry=geometry,prior_policy=prior_policy,
                    prior_policy_precision_epsilon=precision_epsilon,
                    proposal_weights=alpha,component_contributions=component_reports,
                    centered_moment_anchor=anchor,posterior_centered_score_mean=centered_mean,
                    component_contributions_are_numerical_diagnostics_not_error_certificates=True,
                    component_count=len(components),quadrature_node_evaluations=len(entries),
                    cumulative_likelihood_evaluations=evaluations,logmass=float(ll),gradient=g,hessian=h,
                    complete_mixture_denominator_used=True,prior_normalization_included_exactly_once=True,
                    centered_compensated_moments=True,importance_quadrature_full_real_line=True,
                    tail_envelope_for_mass_outside_search_domain=True,
                    tail_allowances_use_estimated_mass_not_certified_lower_bound=True,**tail)
        reports.append(report)
        return float(ll),g,h,report

    def allowance(first,second):
        return (abs(second[0]-first[0])+np.log1p(second[3]['relative_mass_tail_allowance']),
                np.abs(second[1]-first[1])+second[3]['relative_gradient_tail_allowance'],
                np.abs(second[2]-first[2])+second[3]['relative_hessian_tail_allowance'] if with_hessian else None)

    def within(errors,fraction):
        return (errors[0]<=settings.logmass_tolerance*fraction and
                np.max(errors[1])<=settings.gradient_tolerance*fraction and
                (not with_hessian or np.max(errors[2])<=settings.hessian_tolerance*fraction))

    bound=settings.initial_z_bound;accepted=None
    for attempt in range(5):
        geometry=locate(bound,settings.grid_spacing)
        ladder=[run(geometry,n) for n in (40,80,160,320)]
        for first,second in zip(ladder[:-1],ladder[1:]):
            second[3]['previous_order_refinement_plus_tail_allowance']=allowance(first,second)
        errors=allowance(ladder[-2],ladder[-1]);r=ladder[-1][3]
        tail_ok=(r['relative_mass_tail_allowance']<=settings.logmass_tolerance/8 and
                 np.max(r['relative_gradient_tail_allowance'])<=settings.gradient_tolerance/8 and
                 (not with_hessian or np.max(r['relative_hessian_tail_allowance'])<=settings.hessian_tolerance/8))
        if tail_ok:
            if not within(errors,.25):refuse('Mixture GH order refinement exceeds unchanged budgets')
            accepted=ladder[-1];break
        requested=ladder[-1][0]+np.log(min(settings.logmass_tolerance,settings.gradient_tolerance,settings.hessian_tolerance))
        requested-=np.log(16*(1+np.max(cp)+np.max(hp)+np.max(cp)**2))
        if np.log(2.)+norm.logsf(settings.maximum_z_bound)>=requested:
            refuse('Mixture GH relative-tail domain cap exhausted')
        desired=brentq(lambda z:np.log(2.)+norm.logsf(z)-requested,0.,settings.maximum_z_bound)
        bound=max(bound+2,float(np.ceil(desired)))
        if bound+2>settings.maximum_z_bound:refuse('Mixture GH expanded domain exceeds cap')
    if accepted is None:refuse('Mixture GH tail allowances unresolved')
    if bound+2>settings.maximum_z_bound:refuse('Mixture GH final domain exceeds cap')
    geometry=locate(bound+2,settings.grid_spacing/2)
    final_ladder=[run(geometry,n) for n in (40,80,160,320)]
    final=final_ladder[-1];order_errors=allowance(final_ladder[-2],final)
    for first,second in zip(final_ladder[:-1],final_ladder[1:]):
        second[3]['previous_order_refinement_plus_tail_allowance']=allowance(first,second)
    errors=(abs(final[0]-accepted[0])+order_errors[0],np.abs(final[1]-accepted[1])+order_errors[1],
            np.abs(final[2]-accepted[2])+order_errors[2] if with_hessian else None)
    if not within(errors,1.):refuse('Mixture GH node/geometry refinement disagrees')
    final[3].update(logmass_error_allowance=float(order_errors[0]),gradient_error_allowance=order_errors[1],
                    hessian_error_allowance=order_errors[2])
    return final[0],final[1],final[2],dict(
        method='normalized-prior mixture adaptive Gauss-Hermite',first=accepted[3],refined=final[3],
        ladder=reports,geometry_ladder=geometries,likelihood_evaluations=evaluations,
        logmass_refinement_difference=float(abs(final[0]-accepted[0])),
        gradient_refinement_difference=np.abs(final[1]-accepted[1]),
        hessian_refinement_difference=np.abs(final[2]-accepted[2]) if with_hessian else None,
        empirical_refinement_plus_tail_allowances=errors,fixed_point_only=True,prior_policy=prior_policy,
        same_probability_estimand=True,numerical_revision_of_prescribed_20_node_approximation=True,
        refinement_differences_not_certified_error_bounds=True,floating_arithmetic_error_not_rigorously_bounded=True,
        unconditional_all_mode_coverage_guaranteed=False,no_model_fits=True,no_calibration=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true',required=True)
    parser.parse_args()
    print(json.dumps(metadata(),sort_keys=True))
