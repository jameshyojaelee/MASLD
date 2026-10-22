"""Pre-launch checks for REC v1.8.4 (no fits). Writes <out>/check.json; exits 1 on any failure.

1. select_k at K = 2 equals the v1.8 §2 item 3 table on a grid of cases, including differences of exactly eps,
   with differences computed as doubles (addendum item 1). likelihood_v18.select_start compares l_a >= l_b - eps
   in its one-converged branch, which differs from the difference form by one rounding at the band edge
   (first check run, job 22289597); v1.8.4 uses select_k for both section 3 and section 4.
2. Probe start: 2G + 5 uniforms consumed in the item 4 order; values equal an independent draw; holds and
   clip applied (a model with gene 0's rows removed); delta, gamma keep the §7 values.
3. Observed starts 0 and 1 are bitwise equal to REC v1.8.3's (observed_null.json start0_point, start1_point).
"""
import argparse
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
for d in (HERE, PARENT / 'v183', PARENT):
    sys.path.insert(0, str(d))
V183_OBS = Path('/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/'
                'executions/codex-rec-20260929T142434Z/rec-v1.8.3-20261016T234346Z/observed/observed_null.json')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    import numpy as np
    from bfix_v1 import load
    from likelihood_v18 import clip
    from likelihood_v184 import ModelV184, select_k, LOGIT02
    res = {}
    fails = []

    # 1  reference: v1.8 §2 item 3 table written with differences d = l_a - l_b computed as doubles (addendum item 1)
    def table_k2(fits, eps):
        (lw, cw), (lc, cc) = [(f['loglik'], f['converged']) for f in fits]
        if cw and cc:
            d = lc - lw
            return 1 if d >= eps else (0 if -d >= eps else 0)
        if cw or cc:
            k = 0 if cw else 1
            other = fits[1 - k]['loglik']
            if np.isfinite(other) and other - fits[k]['loglik'] > eps:
                return None
            return k
        return None
    eps = 1e-6
    vals = [-1.0, -1.0 + eps, -1.0 - eps, -1.0 + 2 * eps, -1.0 + .5 * eps, float('-inf')]
    n = 0
    for lw, lc, cw, cc in itertools.product(vals, vals, (True, False), (True, False)):
        fits = [{'loglik': lw, 'converged': cw and np.isfinite(lw)}, {'loglik': lc, 'converged': cc and np.isfinite(lc)}]
        k, _ = select_k(fits, eps)
        ref = table_k2(fits, eps)
        n += 1
        if k != ref:
            fails.append({'check': 1, 'lw': lw, 'lc': lc, 'cw': cw, 'cc': cc, 'select_k': k, 'v18': ref})
    res['select_k_K2_cases'] = n

    # 2
    design, settings = load()
    model = ModelV184.from_model(design.model())
    G, Z, size = model.G, model.Z, model.size
    E = (20260923, 13, 3, 7, 1)
    x = model.probe_start(E)
    rng = np.random.default_rng(np.random.SeedSequence(list(E)))
    a = rng.uniform(-1.5, 1.5, G); r = rng.uniform(-6.0, -2.0, G); k = rng.uniform(-0.5, 0.5); q = rng.uniform(-1.0, 1.0, 4)
    after = rng.random()
    rng2 = np.random.default_rng(np.random.SeedSequence(list(E)))
    rng2.random(2 * G + 5)
    consumed_ok = rng2.random() == after
    cold = model.cold_start()
    expect = cold.copy()
    expect[:G] = a; expect[G:2 * G] = r; expect[2 * G] = k
    expect[size - 3], expect[size - 2], expect[size - 1] = q[0], q[1], q[2]
    expect[2 * G + 1 + Z] = q[3]
    expect = clip(model.apply_holds(expect), model.bounds())
    res['probe'] = {'G': G, 'Z': Z, 'size': size, 'uniforms_consumed_2G_plus_5': bool(consumed_ok),
                    'equal_independent_draw': bool(np.array_equal(x, expect)),
                    'delta_equal_cold': bool(np.array_equal(x[2 * G + 1:2 * G + 1 + Z], cold[2 * G + 1:2 * G + 1 + Z])),
                    'probe_start': x.tolist()}
    if not (consumed_ok and np.array_equal(x, expect)):
        fails.append({'check': 2, 'what': 'probe draw'})
    sub = ModelV184([row for row in model.rows if row.gene != 0], G)
    xs = sub.probe_start(E)
    held = xs[0] == 0. and xs[G] == LOGIT02 and np.array_equal(xs[1:G], x[1:G])
    res['probe_hold_gene0_removed'] = {'alpha0': xs[0], 'r0': xs[G], 'others_equal': bool(np.array_equal(xs[1:G], x[1:G]))}
    if not held:
        fails.append({'check': 2, 'what': 'holds'})

    # 3
    prev = json.loads(V183_OBS.read_text())
    x0 = clip(model.apply_holds(clip(model.cold_start(), model.bounds())), model.bounds())
    b0 = model._section3_bounds(x0)
    free = np.array([not (lo is not None and hi is not None and lo == hi) for lo, hi in b0])
    x1 = model._jitter(x0, model._scales(model._hdiag3(x0), free), free, (20260923, 7, 0), b0)
    s0 = np.array_equal(x0, np.asarray(prev['start0_point'], float))
    s1 = np.array_equal(x1, np.asarray(prev['start1_point'], float))
    res['observed_starts_bitwise_equal_v183'] = {'start0': bool(s0), 'start1': bool(s1)}
    if not (s0 and s1):
        fails.append({'check': 3})
    res['failures'] = fails
    (out / 'check.json').write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps({k: v for k, v in res.items() if k != 'probe'}, indent=1))
    if fails:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
