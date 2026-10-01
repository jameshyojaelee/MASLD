"""Smoke test for fit.py / b3_tests.py on simulated data (small SVI and MC sizes)."""
import numpy as np
import pandas as pd
import b1_simulate as b1
import fit
import b3_tests as b3

rng = np.random.default_rng(1)
d, w, lam = b1.simulate(rng)
cohorts = [c for c, _, _ in b1.DESIGN]
frame = pd.DataFrame(np.asarray(d["x"]), columns=[f"hotspot_{i}" for i in range(d["x"].shape[1])])
frame["cohort"] = [cohorts[i] for i in d["cohort"]]
for k, v in d["labels"].items():
    frame[k] = np.where(v >= 0, v, np.nan)
frame["nas"] = np.where(d["nas"] >= 0, d["nas"], np.nan)
frame["fib_lo"] = d["fib_lo"]; frame["fib_hi"] = d["fib_hi"]
prog = [c for c in frame.columns if c.startswith("hotspot_")]
a = fit.arrays(frame, prog, cohorts)
idx = {c: i for i, c in enumerate(cohorts)}
post, loss = fit.fit_svi(a, 2, steps=300, draws=20)
print("svi loss", loss)
te = b3.score_as(fit.subset(a, np.asarray(a["cohort"]) == 0), {"fibrosis", "nas"})
print("lpd K2 on cohort 0 (fib+nas)", fit.cohort_logpred(post, 2, te, rng, n_eps=8, n_shift=2, draws=5))
te2 = b3.score_as(fit.subset(a, np.asarray(a["cohort"]) == 1), {"components"})
print("lpd K2 on cohort 1 (components)", fit.cohort_logpred(post, 2, te2, rng, n_eps=8, n_shift=2, draws=5))
b = fit.simulate_labels(fit.posterior_mean(post), 2, a, rng)
print("simulated label counts", {k: np.bincount(np.asarray(v)[np.asarray(v) >= 0]).tolist() for k, v in b["labels"].items()})
print("smoke ok")
