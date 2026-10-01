#!/usr/bin/env python3
"""B3: PRESPEC_B tests 1-2 with their parametric-bootstrap nulls (Amendment 2).

observed: statistic on real labels; also fits K = 1, 2, 3 on all development data and
          saves posterior means used to simulate null datasets.
null:     statistics on datasets simulated from the fitted K-1 model (same x, same
          label-missingness pattern), refit end to end with the same code.
summary:  decision per test: gain >= 0.02 nats per participant AND above the null's
          95th percentile.

Test 1: K = 2 vs K = 1, leave-one-cohort-out over the six development cohorts; the
        held-out cohort is scored on fibrosis + NAS only (component cohorts collapsed
        to their NAS sum for scoring).
Test 2: K = 3 vs K = 2 and K = 4 vs K = 3, leave-one-cohort-out between GSE130970 and
        GSE267145; the held-out cohort is scored on the three component labels only.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import jax
import numpy as np

import fit

DEV = ["GSE130970", "GSE267145", "GSE135251", "GSE162694", "GSE193066", "GSE213621", "GSE240729"]
COMPONENT_COHORTS = ["GSE130970", "GSE267145"]
MARGIN = 0.02


def check_prespec():
    here = Path(__file__).resolve().parent
    want = (here / "PRESPEC.sha256").read_text().split()[0]
    got = hashlib.sha256((here / "PRESPEC_B.md").read_bytes()).hexdigest()
    if want != got:
        raise SystemExit("PRESPEC_B.md does not match PRESPEC.sha256")


def score_as(a, keep):
    """Copy of a with only the label types in `keep` retained for scoring."""
    b = dict(a)
    comp = {c: np.asarray(v) for c, v in a["labels"].items()}
    has_comp = (comp["steatosis"] >= 0) & (comp["ballooning"] >= 0) & (comp["inflammation"] >= 0)
    nas = np.asarray(a["nas"]).copy()
    if "nas" in keep:
        nas = np.where(has_comp, comp["steatosis"] + comp["ballooning"] + comp["inflammation"], nas)
    else:
        nas[:] = -1
    labels = {c: (v if "components" in keep else -np.ones_like(v)) for c, v in comp.items()}
    b["labels"] = {c: fit.jnp.asarray(v) for c, v in labels.items()}
    b["nas"] = fit.jnp.asarray(nas)
    if "fibrosis" not in keep:
        b["fib_lo"] = fit.jnp.asarray(-np.ones_like(np.asarray(a["fib_lo"])))
        b["fib_hi"] = fit.jnp.asarray(-np.ones_like(np.asarray(a["fib_hi"])))
    return b


def loco(a, cohorts, idx, k_pair, keep, rng):
    out, num, den = {}, 0.0, 0
    for c in cohorts:
        test = np.asarray(a["cohort"]) == idx[c]
        tr = fit.subset(a, ~test)
        te = score_as(fit.subset(a, test), keep)
        lp = {}
        for k in k_pair:
            post, _ = fit.fit_svi(tr, k)
            lp[k] = fit.cohort_logpred(post, k, te, rng)
        jax.clear_caches()  # each subset shape compiles anew; without this, memory grows until the task dies
        n = int(test.sum())
        g = (lp[k_pair[1]] - lp[k_pair[0]]) / n
        out[c] = {"n": n, "gain_per_participant": g, "lpd": {str(k): v for k, v in lp.items()}}
        num += g * n; den += n
    return num / den, out


def statistics(a, idx, rng):
    t1, d1 = loco(a, DEV, idx, (1, 2), {"fibrosis", "nas"}, rng)
    t2a, d2a = loco(a, COMPONENT_COHORTS, idx, (2, 3), {"components"}, rng)
    t2b, d2b = loco(a, COMPONENT_COHORTS, idx, (3, 4), {"components"}, rng)
    return {"test1_K2_vs_K1": t1, "test2_K3_vs_K2": t2a, "test2_K4_vs_K3": t2b,
            "detail": {"test1": d1, "test2_K3_vs_K2": d2a, "test2_K4_vs_K3": d2b}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["observed", "null", "summary"])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--n", type=int, default=1)
    a_ = ap.parse_args()
    check_prespec()
    out = Path(a_.out); out.mkdir(parents=True, exist_ok=True)
    d, prog = fit.load_inputs(a_.inputs)
    d = d[d["cohort"].isin(DEV)].reset_index(drop=True)
    idx = {c: i for i, c in enumerate(DEV)}
    a = fit.arrays(d, prog, DEV)

    if a_.mode == "observed":
        rng = np.random.default_rng(fit.SEED)
        res = statistics(a, idx, rng)
        res["n_participants"] = {c: int((d["cohort"] == c).sum()) for c in DEV}
        (out / "observed.json").write_text(json.dumps(res, indent=2))
        for k in (1, 2, 3):
            post, loss = fit.fit_svi(a, k)
            np.savez(out / f"fit_full_K{k}.npz", **fit.posterior_mean(post))
        print(json.dumps({k: v for k, v in res.items() if k != "detail"}, indent=2))
    elif a_.mode == "null" and a_.n > 1:
        # one process per draw: XLA's compiled code accumulates within a process until
        # mmap fails ("LLVM compilation error: Cannot allocate memory") at ~12 GB RSS
        for i in range(a_.start, a_.start + a_.n):
            if not (out / "null" / f"sim_{i:04d}.json").exists():
                subprocess.run([sys.executable, __file__, "null", "--inputs", a_.inputs, "--out", str(out),
                                "--start", str(i), "--n", "1"], check=True)
    elif a_.mode == "null":
        means = {k: dict(np.load(out / f"fit_full_K{k}.npz")) for k in (1, 2, 3)}
        rows = []
        for i in range(a_.start, a_.start + a_.n):
            if (out / "null" / f"sim_{i:04d}.json").exists():  # a rerun after a failure keeps finished draws
                continue
            rng = np.random.default_rng(fit.SEED + 1000 + i)
            # test 1 null from K = 1; test 2 nulls from K = 2 (for 3 vs 2) and K = 3 (for 4 vs 3)
            s1 = loco(fit.simulate_labels(means[1], 1, a, rng), DEV, idx, (1, 2), {"fibrosis", "nas"}, rng)[0]
            s2a = loco(fit.simulate_labels(means[2], 2, a, rng), COMPONENT_COHORTS, idx, (2, 3), {"components"}, rng)[0]
            s2b = loco(fit.simulate_labels(means[3], 3, a, rng), COMPONENT_COHORTS, idx, (3, 4), {"components"}, rng)[0]
            rows.append({"sim": i, "test1_K2_vs_K1": s1, "test2_K3_vs_K2": s2a, "test2_K4_vs_K3": s2b})
            (out / "null").mkdir(exist_ok=True)
            (out / "null" / f"sim_{i:04d}.json").write_text(json.dumps(rows[-1]))
    else:
        obs = json.loads((out / "observed.json").read_text())
        sims = [json.loads(p.read_text()) for p in sorted((out / "null").glob("sim_*.json"))]
        summ = {"n_null": len(sims)}
        if len(sims) != 200:
            raise SystemExit(f"null incomplete: {len(sims)} of 200 simulated datasets")
        for t in ("test1_K2_vs_K1", "test2_K3_vs_K2", "test2_K4_vs_K3"):
            null = np.array([s[t] for s in sims])
            q95 = float(np.quantile(null, 0.95)) if len(null) else None
            summ[t] = {"observed": obs[t], "null_q95": q95, "null_mean": float(null.mean()) if len(null) else None,
                       "p_upper": float((1 + (null >= obs[t]).sum()) / (1 + len(null))) if len(null) else None,
                       "keep_larger_K": bool(len(null) and obs[t] >= MARGIN and obs[t] > q95)}
        (out / "summary.json").write_text(json.dumps(summ, indent=2))
        print(json.dumps(summ, indent=2))


if __name__ == "__main__":
    main()
