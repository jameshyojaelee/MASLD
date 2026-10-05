#!/usr/bin/env python3
"""Independent Decimal boundary diagnostic with shared integer rising products."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path


def bb_vector(support, n, a, b, logchoose):
    """Full integer rising prefixes; floating Decimal addition order changes."""
    zero = Decimal(0)
    def prefix(shape):
        values = [zero]
        for j in range(n):
            values.append(values[-1]+(shape+Decimal(j)).ln())
        return values
    ra, rb, rt = prefix(a), prefix(b), prefix(a+b)
    return [logchoose[k]+ra[k]+rb[n-k]-rt[n] for k in support]


def boundary_check(rows, p, G, Z):
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    rows = tuple(rows)
    if not rows or len(p) != 2*G+Z+5:
        raise ValueError('Fixed-point shape differs')
    records = []
    for digits in (80,110):
        with localcontext() as context:
            context.prec = digits
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

            ell0 = loglik(zero)
            ladder = []
            for t_text in ('1e-15','1e-20','1e-25'):
                t = Decimal(t_text)
                plus,minus = loglik(t)-ell0,loglik(-t)-ell0
                score = (plus.exp()+minus.exp()-two)/(two*t*t)
                ladder.append(dict(u_step=t_text,right_variance_score=str(score),
                                   plus_log_probability_ratio=str(plus),minus_log_probability_ratio=str(minus)))
            records.append(dict(decimal_digits=digits,fixed_gene_logmass=str(ell0),
                                right_variance_score=ladder[-1]['right_variance_score'],ladder=ladder))
    with localcontext() as context:
        context.prec = 110
        differences = [dict(u_step=a['u_step'],absolute_precision_difference=str(abs(Decimal(a['right_variance_score'])-Decimal(b['right_variance_score']))))
                       for a,b in zip(records[0]['ladder'],records[1]['ladder'])]
    return dict(method='Gaussian heat-semigroup via symmetric normalized probability curvature',
        records=records,precision_comparisons=differences,
        absolute_precision_difference=differences[-1]['absolute_precision_difference'],
        input_float_conversion='exact binary floats via Decimal.from_float(float)',
        fixed_synthetic_rows=len(rows),no_model_fit=True,no_statistical_calibration=True,
        declared_finite_difference_failures_preserved=True,
        no_original_derivative_helpers=True,no_quadrature=True,stdlib_only=True,
        shared_integer_rising_prefixes=True,arithmetic_addition_order_changed=True,
        arithmetic_roundoff_and_finite_difference_truncation_not_certified=True)


def compare_small():
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    records = []
    for digits in (80,110):
        with localcontext() as context:
            context.prec = digits
            zero,one = Decimal(0),Decimal(1)
            for n in (5,20,100):
                logchoose = [Decimal(math.comb(n,k)).ln() for k in range(n+1)]
                for a,b in ((Decimal('0.000001'),Decimal('3.25')),
                            (Decimal('2.5'),Decimal('2.5')),
                            (Decimal('10000'),Decimal('0.000001'))):
                    support = tuple(range(n+1))
                    def rise(shape,count):
                        return sum(((shape+Decimal(j)).ln() for j in range(count)),zero)
                    literal = [logchoose[k]+rise(a,k)+rise(b,n-k)-rise(a+b,n) for k in support]
                    fast = bb_vector(support,n,a,b,logchoose)
                    selected = slice(n//5,n-n//7+1)
                    refden = sum((x.exp() for x in literal[selected]),zero).ln()
                    fastden = sum((x.exp() for x in fast[selected]),zero).ln()
                    conditional = max(abs((x-refden)-(y-fastden)) for x,y in zip(literal[selected],fast[selected]))
                    records.append(dict(decimal_digits=digits,n=n,a=str(a),b=str(b),
                        max_absolute_logpmf_difference=str(max(abs(x-y) for x,y in zip(literal,fast))),
                        max_absolute_conditional_logpmf_difference=str(conditional)))
    return dict(fixed_synthetic_arithmetic=True,no_fit=True,no_calibration=True,
                arithmetic_addition_order_changed=True,cases=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true')
    parser.add_argument('--compare',action='store_true')
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.compare:
        if args.inspect or args.output is None or args.output.exists():
            raise ValueError('Comparison requires a new output file')
        result = compare_small()
        result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        with args.output.open('x') as handle:
            json.dump(result,handle,indent=2,allow_nan=False)
            handle.write('\n')
        return
    if not args.inspect:
        raise ValueError('Use metadata inspection or compute comparison')
    print(json.dumps(dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        numerical_execution=False,data_reads=False,no_model_fit=True,no_calibration=True,
        intended_decimal_precisions=[80,110],u_steps=['1e-15','1e-20','1e-25'],stdlib_only=True)))


if __name__ == '__main__':
    main()

