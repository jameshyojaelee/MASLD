"""Synthetic check of the N2a cross-fit statistic (n2a_lib.site_record). Simulated reads only; no study data.
Usage: python -I test_crossfit.py
Scenarios: unchanged skew (0.8, 0.8), balanced (0.5, 0.5), weakened (0.8, 0.6), skew gained from a balanced biopsy 1
(0.5, 0.8; the known blind spot), flipped dominant X (0.8, 0.2). 'naive' phases by the whole biopsy-1 majority and
measures biopsy 1 on the same reads."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from n2a_lib import site_record

rng = np.random.default_rng(20260923)
SITES, REPS = 60, 2000


def sim_counts(n, p):
    n0 = rng.binomial(n, 0.5)
    r0 = rng.binomial(n0, p); r1 = rng.binomial(n - n0, p)
    return (r0, n0 - r0, r1, (n - n0) - r1, 0, 0)


def one_woman(th1, th2):
    num_cf = [0, 0, 0, 0]; naive = [0, 0, 0, 0]
    for _ in range(SITES):
        z = 1 if rng.random() < 0.5 else -1
        p1 = th1 if z == 1 else 1 - th1
        p2 = th2 if z == 1 else 1 - th2
        n1, n2 = 20 + rng.poisson(40), 20 + rng.poisson(40)
        a, b = sim_counts(n1, p1), sim_counts(n2, p2)
        rec = site_record(a, b)
        if rec is None:
            continue
        _, x1, m1, x2, m2 = rec
        for i, v in enumerate((x1, m1, x2, m2)):
            num_cf[i] += v
        ref1, alt1 = a[0] + a[2], a[1] + a[3]
        if ref1 != alt1:
            phase_ref = ref1 > alt1
            naive[0] += ref1 if phase_ref else alt1; naive[1] += ref1 + alt1
            naive[2] += (b[0] + b[2]) if phase_ref else (b[1] + b[3]); naive[3] += n2
    cf = num_cf[2] / num_cf[3] - num_cf[0] / num_cf[1] if num_cf[1] and num_cf[3] else np.nan
    nv = naive[2] / naive[3] - naive[0] / naive[1] if naive[1] and naive[3] else np.nan
    return cf, nv


print("scenario            true_change  mean_cross_fit  sd_cross_fit  mean_naive")
bad = []
for name, th1, th2 in [("unchanged_0.8", 0.8, 0.8), ("balanced_0.5", 0.5, 0.5), ("weakened_0.8_0.6", 0.8, 0.6),
                       ("gain_from_balanced", 0.5, 0.8), ("flipped_0.8_0.2", 0.8, 0.2)]:
    res = np.array([one_woman(th1, th2) for _ in range(REPS)])
    cf, nv = res[:, 0], res[:, 1]
    print("%-19s %+8.2f %14.4f %13.4f %11.4f" % (name, abs(th2 - 0.5) - abs(th1 - 0.5), np.nanmean(cf), np.nanstd(cf), np.nanmean(nv)))
    if name in ("unchanged_0.8", "balanced_0.5", "gain_from_balanced") and abs(np.nanmean(cf)) > 0.01:
        bad.append(name)
    if name == "flipped_0.8_0.2" and np.nanmean(cf) > -0.5:
        bad.append(name)
print("FAILED: " + ",".join(bad) if bad else "all checks passed")
sys.exit(1 if bad else 0)
