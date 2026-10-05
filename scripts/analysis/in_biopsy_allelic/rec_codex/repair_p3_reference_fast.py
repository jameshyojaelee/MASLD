#!/usr/bin/env python3
"""Isolated full-support BB reference arithmetic comparison; no fits or data reads."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time


def fast_bb(k, n, eta, rr):
    """Literal integer rising products, shared prefix sums over complete support.

    Imports occur only on numerical execution. No likelihood/score evaluator is
    reused. This preserves the reference rho, shapes, and combinatorial factor;
    cumulative summation changes floating-point addition order.
    """
    import numpy as np
    from scipy.special import expit, gammaln
    rho = 1e-6 + (1-2e-6)*expit(rr)
    logt = np.log((1-rho)/rho)
    la, lb = logt-np.logaddexp(0.,-eta), logt-np.logaddexp(0.,eta)

    def prefixes(shape):
        terms = np.empty(n, dtype=float)
        if n:
            terms[0] = shape
            terms[1:] = np.logaddexp(shape, np.log(np.arange(1, n, dtype=float)))
        return np.concatenate(([0.], np.cumsum(terms)))

    k = np.asarray(k, dtype=int)
    if n < 0 or np.any(k < 0) or np.any(k > n):
        raise ValueError('Invalid integer support')
    return (gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
            + prefixes(la)[k]+prefixes(lb)[n-k]-prefixes(logt)[n])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.inspect:
        print(json.dumps(dict(source_sha256=digest, numerical_execution=False,
                              data_reads=False, fitting=False, calibration=False)))
        return
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    if args.output is None or args.output.exists():
        raise ValueError('New output file required')
    import numpy as np
    import scipy
    from scipy.special import expit, gammaln, logsumexp

    def rise(shape, count):
        return 0. if count == 0 else shape + np.logaddexp(shape, np.log(np.arange(1,count,dtype=float))).sum()

    def literal(k, n, eta, rr):
        rho = 1e-6 + (1-2e-6)*expit(rr)
        logt = np.log((1-rho)/rho)
        la, lb = logt-np.logaddexp(0.,-eta), logt-np.logaddexp(0.,eta)
        return np.array([gammaln(n+1)-gammaln(a+1)-gammaln(n-a+1)+rise(la,int(a))+rise(lb,int(n-a))-rise(logt,n) for a in k])

    cases = []
    # Full support and asymmetric truncation together detect reversed indexing.
    for n in (0, 5, 20, 100, 500, 2000):
        for eta in (-1000., -30., -6., 0., 6., 30., 1000.):
            for rr in (-1000., -15., -3., 5., 1000.):
                k = np.arange(n+1)
                start = time.perf_counter()
                reference = literal(k,n,eta,rr)
                literal_seconds = time.perf_counter()-start
                start = time.perf_counter()
                candidate = fast_bb(k,n,eta,rr)
                fast_seconds = time.perf_counter()-start
                lo, hi = n//5, n-n//7
                sel = slice(lo,hi+1)
                ref_cond = reference[sel]-logsumexp(reference[sel])
                cand_cond = candidate[sel]-logsumexp(candidate[sel])
                cases.append(dict(n=n,eta=eta,linear_r=rr,
                    max_absolute_logpmf_difference=float(np.max(np.abs(reference-candidate))),
                    max_absolute_conditional_logpmf_difference=float(np.max(np.abs(ref_cond-cand_cond))),
                    literal_full_mass=float(np.exp(logsumexp(reference))),
                    candidate_full_mass=float(np.exp(logsumexp(candidate))),
                    finite=bool(np.isfinite(reference).all() and np.isfinite(candidate).all()),
                    literal_seconds=literal_seconds,fast_seconds=fast_seconds))
    result = dict(source_sha256=digest, environment=dict(python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__), fixed_synthetic_arithmetic=True,
        stochastic_seed=None, no_fits=True, no_calibration=True,
        no_integration_accuracy_claim=True, addition_order_changed=True, cases=cases)
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')


if __name__ == '__main__':
    main()
