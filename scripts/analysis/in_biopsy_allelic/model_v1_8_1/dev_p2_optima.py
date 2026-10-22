"""Development check: P2 resample refits that differ between optimizer paths. For chosen draws,
refit section 3 from several starts and print log-likelihood, kappa, bound coordinates and the
smallest Hessian eigenvalues."""
import sys

import numpy as np

import fitting as ft
import fixture_io as fio
import model_a as ma

FIX = ("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/"
       "executions/model-a-bfix-20260930T115521Z")

prob, glob = fio.load_problem(FIX)
L = prob.L
draws = fio.load_draws(FIX, prob)
D = prob.data()
fit3 = ft.fit_section3(D)
for b in [int(x) for x in sys.argv[1:]]:
    Db = prob.data(mult=draws["mult"][b])
    starts = {"warm": fit3.theta, "cold": None}
    rng = np.random.default_rng(b)
    for k in range(4):
        s = fit3.theta.copy()
        s[L.kappa] = rng.uniform(-0.6, 0.6)
        starts[f"kappa{k}"] = s
    for name, s in starts.items():
        f = ft.fit_section3(Db, start=s)
        ll, g, H = ma.loglik3(f.theta, Db, 2)
        idx = np.flatnonzero(Db.free3 & (f.theta > L.lb) & (f.theta < L.ub))
        ev = np.linalg.eigvalsh(-H[np.ix_(idx, idx)])
        print(b, name, repr(f.loglik), repr(f.theta[L.kappa]), f.converged, f.proj_max_grad,
              [L.names[i] for i in f.at_bound], ev[:3], flush=True)
        print("   theta", {n: round(float(x), 4) for n, x in zip(L.names, f.theta)
                           if n.startswith(("alpha", "r[", "delta", "rho", "omega"))}, flush=True)
