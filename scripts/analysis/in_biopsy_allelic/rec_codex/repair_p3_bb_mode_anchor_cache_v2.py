"""One-point full-support reuse for the reviewed native BB arithmetic.

The original conditional law requests tails, sometimes the central support,
then the observed count at identical n/eta/r. Retain one normalized support
for those consecutive requests. No quadrature, probability or fit rule changes.
Each instance belongs to one sequential evaluator or one process worker.
"""
import hashlib
import inspect
import os
from pathlib import Path
import threading

SCALAR_SHA = 'a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266'
MAXIMUM_N = 100000
MAXIMUM_RETAINED_BYTES = 8*1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)


class OnePointBBTerms:
    def __init__(self, scalar):
        require(os.environ.get('SLURM_JOB_ID', '').isdigit()
                and not os.environ.get('SLURM_ARRAY_TASK_ID'),
                'Individual compute allocation required before numerical imports')
        source = Path(inspect.getsourcefile(scalar))
        require(scalar.__name__ == 'stable_bb_terms'
                and hashlib.sha256(source.read_bytes()).hexdigest() == SCALAR_SHA,
                'Reviewed scalar BB source differs')
        self.scalar = scalar
        self.owner_process = os.getpid()
        self.owner_thread = threading.get_ident()
        self.key = None
        self.values = None
        self.calls = self.hits = self.misses = self.maximum_retained_bytes = 0

    def __call__(self, k, n, eta, linear_r, need_hessian=True):
        require(os.environ.get('SLURM_JOB_ID', '').isdigit()
                and not os.environ.get('SLURM_ARRAY_TASK_ID'),
                'Individual compute allocation required')
        require(os.getpid() == self.owner_process and threading.get_ident() == self.owner_thread,
                'Support reuse requires its constructing process and thread')
        import numpy as np
        def exact_float64_input(value):
            return (isinstance(value, (float, np.float64))
                    or (type(value) is int and abs(value) <= 2**53))
        require(exact_float64_input(eta) and exact_float64_input(linear_r),
                'Predictors require float64 or exactly represented bounded Python integers')
        k = np.atleast_1d(np.asarray(k))
        require(isinstance(n, (int, float, np.integer, np.floating))
                and np.isfinite(n) and int(n) == n and 0 <= n <= MAXIMUM_N
                and k.ndim == 1 and np.isfinite(k).all()
                and np.all(k == np.floor(k)) and np.all((k >= 0) & (k <= n))
                and np.isfinite(eta) and np.isfinite(linear_r),
                'Invalid bounded integer BB input')
        n, k = int(n), k.astype(int)
        # Exact float encodings include signed zero. Hessian requests cannot
        # reuse score-only values. No rounding or approximate key is used.
        key = (n, float(eta).hex(), float(linear_r).hex(), bool(need_hessian))
        retained_estimate = (56 if need_hessian else 24)*(n+1)
        require(retained_estimate <= MAXIMUM_RETAINED_BYTES,
                'Single retained full-support capacity exceeded')
        self.calls += 1
        if key != self.key:
            self.misses += 1
            # Discard the preceding support before allocating its replacement.
            # Failed computation leaves an empty cache and preserves refusal.
            self.key, self.values = None, None
            values = self.scalar(np.arange(n+1), n, eta, linear_r, need_hessian)
            lp, score, hessian = values
            require(lp.shape == (n+1,) and score.shape == (n+1, 2)
                    and ((hessian is None and not need_hessian)
                         or (need_hessian and hessian.shape == (n+1, 2, 2)))
                    and np.isfinite(lp).all() and np.isfinite(score).all()
                    and (hessian is None or np.isfinite(hessian).all()),
                    'Invalid reviewed full-support result')
            retained_bytes = sum(v.nbytes for v in values if v is not None)
            require(retained_bytes == retained_estimate,
                    'Retained array size differs from declared capacity')
            for value in values:
                if value is not None:
                    value.setflags(write=False)
            self.values, self.key = values, key
            self.maximum_retained_bytes = max(self.maximum_retained_bytes, retained_bytes)
        else:
            self.hits += 1
        lp, score, hessian = self.values
        # Advanced indexing returns independent arrays. A caller cannot mutate
        # the retained support or poison a later request at the same point.
        return lp[k], score[k], hessian[k] if need_hessian else None

    def statistics(self):
        return dict(calls=self.calls, hits=self.hits, misses=self.misses,
                    maximum_retained_bytes=self.maximum_retained_bytes,
                    retained_supports=0 if self.values is None else 1,
                    maximum_retained_bytes_allowed=MAXIMUM_RETAINED_BYTES,
                    not_a_total_process_memory_bound=True)
