#!/usr/bin/env python3
"""Independent Decimal positive-variance scalar/Hvv diagnostic, synthetic only."""
from decimal import Decimal, localcontext
import math
import os
import time

ORDERS = (20, 40, 80)
DIGITS = (80, 110)
MAXIMUM_EVALUATIONS = 100000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def prepare_probability(rows, p, G, Z, bb_vector):
    """Native selected genotype-mixture probability; no derivative helpers."""
    zero, one, two = Decimal(0), Decimal(1), Decimal(2)
    exact = lambda x: Decimal.from_float(float(x))
    q = [exact(x) for x in p]
    eps, scale = exact(1e-6), exact(1-2e-6)
    prepared = []
    for row in rows:
        n, observed = int(row.n),int(row.a)
        lo = max(int(row.minor),(int(row.percent)*n+99)//100)
        support = tuple(range(lo,n-lo+1))
        if n < int(row.depth) or observed not in support or len(row.z) != Z:
            raise ValueError('Selected support/covariates differ')
        if row.orientation not in (-1,1) or not 0 <= int(row.gene) < G:
            raise ValueError('Gene/orientation differs')
        stage = exact(row.stage)
        eta0 = q[row.gene]*(one+q[2*G]*stage+sum((q[2*G+1+j]*exact(row.z[j]) for j in range(Z)),zero))
        eta0 -= Decimal(row.orientation)*(q[2*G+1+Z]*stage+exact(row.omega))
        rr = q[G+row.gene]+sum((q[-3+j]*exact(x) for j,x in enumerate((row.stage,row.background,row.duplication))),zero)
        rho = eps+scale/(one+(-rr).exp())
        concentration = (one-rho)/rho
        pref,palt,px = (exact(x) for x in (row.p_ref,row.p_alt,row.p_x))
        prior = [two*pref*palt,pref*pref+two*pref*px,palt*palt+two*palt*px]
        if row.orientation == -1:
            prior = [prior[0],prior[2],prior[1]]
        if any(x < zero for x in prior) or sum(prior,zero) <= zero:
            raise ValueError('Invalid genotype mixture')
        logchoose = [Decimal(math.comb(n,k)).ln() for k in range(n+1)]
        error,phi = exact(row.e),exact(row.phi)
        if not zero < error < one or not zero <= phi < one:
            raise ValueError('Invalid homozygote component')

        def homo(mean):
            if phi == zero:
                return [Decimal(math.comb(n,k))*mean**k*(one-mean)**(n-k) for k in support]
            t = (one-phi)/phi
            return [lp.exp() for lp in bb_vector(support,n,mean*t,(one-mean)*t,logchoose)]

        prepared.append((n,support,observed-lo,eta0,stage,concentration,prior,homo(error),homo(one-error),logchoose))

    def loglik(u):
        total = zero
        for n,support,at,eta0,stage,t,prior,homref,homalt,logchoose in prepared:
            eta = eta0+u*stage
            # Separate shapes avoid subtracting a nearly unit mean.
            a,b = t/(one+(-eta).exp()),t/(one+eta.exp())
            het = bb_vector(support,n,a,b,logchoose)
            masses = [prior[0]*lp.exp()+prior[1]*homref[j]+prior[2]*homalt[j]
                      for j,lp in enumerate(het)]
            total += masses[at].ln()-sum(masses,zero).ln()
        return total


    return loglik


def hermite_pair(n, z):
    """Probabilists He_n and He_(n-1), normalized standard-normal rule."""
    previous, current = Decimal(1), z
    for k in range(1, n):
        previous, current = current, z*current-Decimal(k)*previous
    return current, previous


def decimal_rule(n, digits, seed_nodes):
    require(n in ORDERS and len(seed_nodes) == n, 'Unexpected GH order/seeds')
    tolerance = Decimal(10)**(-(digits-15))
    positives = []
    corrections = []
    for seed in seed_nodes:
        if seed <= 0:
            continue
        z = Decimal.from_float(float(seed))
        for iteration in range(100):
            value, previous = hermite_pair(n, z)
            delta = value/(Decimal(n)*previous)
            z -= delta
            if abs(delta) <= tolerance:
                break
        else:
            raise ValueError('Decimal Hermite Newton iteration cap exhausted')
        value, previous = hermite_pair(n, z)
        correction = abs(value/(Decimal(n)*previous))
        require(correction <= tolerance, 'Decimal root correction unresolved')
        positives.append(z)
        corrections.append(correction)
    require(len(positives) == n//2 and all(z > 0 for z in positives),
            'Positive root family differs')
    positives.sort()
    roots = [-z for z in reversed(positives)]+positives
    require(all(a < b for a,b in zip(roots[:-1],roots[1:])), 'Repeated or unordered roots')
    factorial = Decimal(math.factorial(n-1))
    weights = [factorial/(Decimal(n)*hermite_pair(n,z)[1]**2) for z in roots]
    require(all(w.is_finite() and w > 0 for w in weights), 'Invalid Decimal GH weight')
    moment_checks = []
    for k in range(9):
        expected = Decimal(0) if k%2 else Decimal(math.prod(range(1,k,2)))
        actual = sum((w*z**k for z,w in zip(roots,weights)), Decimal(0))
        error = abs(actual-expected)
        require(error <= Decimal(10)**(-(digits-20))*(1+abs(expected)),
                'Normalized Gaussian rule moment check failed')
        moment_checks.append(dict(power=k,actual=str(actual),expected=str(expected),absolute_error=str(error)))
    require(all(a == -b for a,b in zip(roots,reversed(roots))), 'Root symmetry differs')
    require(all(a == b for a,b in zip(weights,reversed(weights))), 'Weight symmetry differs')
    return roots,weights,dict(order=n,decimal_digits=digits,max_root_correction=str(max(corrections)),
        normalized_gaussian_moments=moment_checks,roots=[str(z) for z in roots],
        weights=[str(w) for w in weights],double_used_for_initial_guesses_only=True,
        probabilists_weight_formula='(n-1)! / (n * He_(n-1)(z)^2)',
        root_and_moment_checks_are_not_quadrature_accuracy_checks=True)


def tail_diagnostic(bound, v, logmass, score):
    """Conservative polynomial envelopes; divide by an uncertified mass estimate.

    For B>0, use phi(B)<exp(-B²/2), Mills M0 <=2exp(-B²/2)/B.
    Integration-by-parts bounds give M2 and M4. These bound probability-law
    tails with L<=1; GH rule error inside the node range is not bounded.
    """
    one, two = Decimal(1), Decimal(2)
    exponential = (-bound*bound/two).exp()
    m0 = two*exponential/bound
    m2 = two*bound*exponential+m0
    m4 = two*bound**3*exponential+Decimal(3)*m2
    mass = logmass.exp()
    prior = m0/mass
    derivative = (m2+m0)/(two*v*mass)
    rawvv = (m4+Decimal(6)*m2+Decimal(3)*m0)/(Decimal(4)*v*v*mass)
    return dict(node_range_bound=str(bound),relative_mass_envelope=str(prior),
        absolute_score_tail_envelope=str(derivative+abs(score)*prior),
        absolute_raw_vv_tail_envelope=str(rawvv),
        conditional_on_estimated_mass_not_certified_lower_bound=True,
        GH_is_full_real_line_not_a_truncated_integral=True,
        interior_quadrature_error_not_bounded=True,tail_diagnostic_not_accuracy_acceptance=True)


def check(rows, p, variance, G, Z, bb_vector, seed_factory, save):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required')
    require(len(p) == 2*G+Z+5 and rows and variance > 0, 'Fixed synthetic shape/variance differs')
    evaluations, records, began = 0, [], time.monotonic()
    for digits in DIGITS:
        with localcontext() as context:
            context.prec = digits
            v = Decimal.from_float(float(variance))
            rootv = v.sqrt()
            probability = prepare_probability(rows,p,G,Z,bb_vector)
            for order in ORDERS:
                phase_began = time.monotonic()
                nodes,weights,rule = decimal_rule(order,digits,seed_factory(order))
                save('rule_%d_%d' % (digits,order),rule)
                values = []
                for z in nodes:
                    evaluations += 1
                    require(evaluations <= MAXIMUM_EVALUATIONS, 'Node evaluation cap exhausted')
                    value = probability(rootv*z)
                    require(value.is_finite() and value <= Decimal('1e-25'), 'Invalid selected log probability')
                    values.append(value)
                shift = max(values)
                masses = [w*(ell-shift).exp() for w,ell in zip(weights,values)]
                mass = sum(masses,Decimal(0))
                require(mass.is_finite() and mass > 0, 'Invalid scaled Decimal mass')
                score = sum((w*(z*z-1)/(2*v) for w,z in zip(masses,nodes)),Decimal(0))/mass
                raw = sum((w*(z**4-6*z*z+3)/(4*v*v) for w,z in zip(masses,nodes)),Decimal(0))/mass
                hessian = raw-score*score
                logmass = shift+mass.ln()
                record = dict(decimal_digits=digits,order=order,logmass=str(logmass),
                    variance_score=str(score),raw_variance_second_moment=str(raw),Hvv=str(hessian),
                    fixed_u_prior_score='H2(z)/(2v)',fixed_u_raw_vv='H4(z)/(4v^2)',
                    tail=tail_diagnostic(max(abs(z) for z in nodes),v,logmass,score),
                    elapsed_seconds=time.monotonic()-phase_began,
                    cumulative_probability_evaluations=evaluations,independent_Hvv_accuracy_accepted=False)
                records.append(record)
                save('point_%d_%d' % (digits,order),record)
    with localcontext() as context:
        context.prec = max(DIGITS)
        fields = ('logmass','variance_score','Hvv')
        refinements = []
        for digits in DIGITS:
            series = [r for r in records if r['decimal_digits']==digits]
            for a,b in zip(series[:-1],series[1:]):
                refinements.append(dict(kind='order',digits=digits,first_order=a['order'],second_order=b['order'],
                    differences={k:str(Decimal(b[k])-Decimal(a[k])) for k in fields}))
        for order in ORDERS:
            a,b = [r for r in records if r['order']==order]
            refinements.append(dict(kind='digits',order=order,first_digits=a['decimal_digits'],second_digits=b['decimal_digits'],
                differences={k:str(Decimal(b[k])-Decimal(a[k])) for k in fields}))
    save('refinements',refinements)
    return dict(completed=True,records=records,refinements=refinements,
        probability_evaluations=evaluations,maximum_probability_evaluations=MAXIMUM_EVALUATIONS,
        elapsed_seconds=time.monotonic()-began,arbitrary_u_predictors_in_Decimal=True,
        binary_float_inputs_converted_exactly=True,no_original_derivative_helpers=True,
        GH20_accuracy_assumed=False,independent_Hvv_accuracy_accepted=False,
        full_Hessian_accuracy_verified=False,no_model_fit=True,no_calibration=True,
        arithmetic_and_quadrature_error_not_certified=True)

