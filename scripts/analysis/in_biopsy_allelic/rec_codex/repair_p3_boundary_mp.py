#!/usr/bin/env python3
"""Independent high-precision synthetic boundary derivative; no data loading."""
import argparse
import hashlib
import json
import os
from pathlib import Path


def boundary_check(rows, p, G, Z):
    """Right variance score from a literal conditional probability calculation.

    All inputs must be the caller's already reconstructed synthetic fixed point.
    No original likelihood, score, Hessian, or quadrature helpers are imported.
    Existing declared finite-difference records remain necessary and unchanged.
    """
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    import mpmath as mp
    if mp.__version__ != '1.3.0':
        raise RuntimeError('Reviewed mpmath version differs')
    rows = tuple(rows)
    if not rows or len(p) != 2*G+Z+5:
        raise ValueError('Fixed-point shape differs')
    records = []
    for digits in (80, 110):
        with mp.workdps(digits):
            # mp.mpf(float) preserves the exact input binary value, not its
            # shortened decimal display. Integer counts retain integer values.
            q = [mp.mpf(float(x)) for x in p]
            eps = mp.mpf(float(1e-6))
            rho_scale = mp.mpf(float(1-2e-6))

            def rising(shape, count):
                return mp.fsum(mp.log(shape+j) for j in range(count))

            def bb(k, n, mean, concentration):
                return (mp.log(mp.binomial(n,k))
                        + rising(mean*concentration,k)
                        + rising((1-mean)*concentration,n-k)
                        - rising(concentration,n))

            def sigmoid(x):
                return 1/(1+mp.exp(-x))

            prepared = []
            for row in rows:
                n, observed = int(row.n), int(row.a)
                lo = max(int(row.minor),(int(row.percent)*n+99)//100)
                support = tuple(range(lo,n-lo+1))
                if n < int(row.depth) or observed not in support or len(row.z) != Z:
                    raise ValueError('Selected synthetic support/covariates differ')
                if row.orientation not in (-1,1) or not 0 <= int(row.gene) < G:
                    raise ValueError('Gene/orientation differs')
                stage = mp.mpf(float(row.stage))
                z = [mp.mpf(float(x)) for x in row.z]
                eta0 = q[row.gene]*(1+q[2*G]*stage+mp.fsum(q[2*G+1+j]*z[j] for j in range(Z)))
                eta0 -= row.orientation*q[2*G+1+Z]*stage+row.orientation*mp.mpf(float(row.omega))
                rr = q[G+row.gene]+mp.fsum(q[-3+j]*mp.mpf(float(x)) for j,x in enumerate((row.stage,row.background,row.duplication)))
                rho = eps+rho_scale*sigmoid(rr)
                concentration = (1-rho)/rho
                pref, palt, px = (mp.mpf(float(x)) for x in (row.p_ref,row.p_alt,row.p_x))
                prior = [2*pref*palt,pref**2+2*pref*px,palt**2+2*palt*px]
                if row.orientation == -1:
                    prior = [prior[0],prior[2],prior[1]]
                if any(x < 0 for x in prior) or mp.fsum(prior) <= 0:
                    raise ValueError('Invalid genotype mixture')
                error, phi = mp.mpf(float(row.e)), mp.mpf(float(row.phi))
                if not 0 < error < 1 or not 0 <= phi < 1:
                    raise ValueError('Invalid homozygote component')

                def homo(mean):
                    if phi == 0:
                        return [mp.binomial(n,k)*mean**k*(1-mean)**(n-k) for k in support]
                    return [mp.exp(bb(k,n,mean,(1-phi)/phi)) for k in support]

                prepared.append((n,support,observed-lo,eta0,stage,concentration,prior,homo(error),homo(1-error)))

            def loglik(u):
                contributions = []
                for n,support,at,eta0,stage,concentration,prior,homref,homalt in prepared:
                    mean = sigmoid(eta0+u*stage)
                    masses = [prior[0]*mp.exp(bb(k,n,mean,concentration))
                              +prior[1]*homref[j]+prior[2]*homalt[j]
                              for j,k in enumerate(support)]
                    contributions.append(mp.log(masses[at])-mp.log(mp.fsum(masses)))
                return mp.fsum(contributions)

            ell0 = loglik(mp.mpf(0))
            second = mp.diff(lambda u: mp.exp(loglik(u)-ell0),mp.mpf(0),2)
            records.append(dict(decimal_digits=digits,
                fixed_gene_logmass=mp.nstr(ell0,digits),
                normalized_probability_second_u_derivative=mp.nstr(second,digits),
                right_variance_score=mp.nstr(second/2,digits)))
    with mp.workdps(110):
        agreement = abs(mp.mpf(records[0]['right_variance_score'])-mp.mpf(records[1]['right_variance_score']))
        difference = mp.nstr(agreement,110)
    return dict(method='Gaussian heat-semigroup: right variance score = half normalized probability second u derivative',
        records=records,absolute_precision_difference=difference,mpmath_version=mp.__version__,
        input_float_conversion='exact binary floats via mp.mpf(float)',
        fixed_synthetic_rows=len(rows),no_model_fit=True,no_statistical_calibration=True,
        declared_finite_difference_failures_preserved=True,
        no_original_derivative_helpers=True,no_quadrature=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        numerical_execution=False,data_reads=False,no_model_fit=True,no_calibration=True,
        intended_decimal_precisions=[80,110])))


if __name__ == '__main__':
    main()
