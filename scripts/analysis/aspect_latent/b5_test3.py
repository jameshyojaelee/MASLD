#!/usr/bin/env python3
"""B5: PRESPEC_B test 3, component-specific signal, with its parametric-bootstrap null.

Leave-one-cohort-out between the two component cohorts (GSE130970, GSE267145): fit the
kept K (2) by SVI on all development participants outside the held-out cohort; for each
label j, predicted_j = x @ W @ lambda_j averaged over the 200 posterior draws (a
monotone function of the expected ordinal label, so its ranks suffice). Statistic per
label: Spearman partial correlation between predicted_j and observed j given the other
three observed labels (ranks residualized on the other labels' ranks by OLS), then the
participant-weighted mean over the two held-out cohorts.

Null: 1,000 datasets simulated from the fitted K - 1 model (fit_full_K1.npz from B3),
same x and label-missingness pattern, statistic recomputed from the refit. One-sided
p = (1 + #null >= observed) / 1001; BH across the four labels.

modes: observed | null (--start/--n, one process per draw) | summary
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import jax
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit  # noqa: E402
from b3_tests import COMPONENT_COHORTS, DEV, check_prespec  # noqa: E402
from model import LABELS  # noqa: E402

K_KEPT, N_NULL = 2, 1000


def observed_labels(a):
    lab = {c: np.asarray(v) for c, v in a["labels"].items()}
    lab["fibrosis"] = np.where(np.asarray(a["fib_lo"]) == np.asarray(a["fib_hi"]), np.asarray(a["fib_lo"]), -1)
    return lab


def partial_spearman(pred, y, others):
    ok = (y >= 0) & np.all(np.stack([o >= 0 for o in others]), 0)
    if ok.sum() < 10:
        return np.nan, int(ok.sum())
    r = lambda v: stats.rankdata(v[ok])
    Z = np.column_stack([np.ones(ok.sum())] + [r(o) for o in others])
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    return float(np.corrcoef(res(r(pred)), res(r(y)))[0, 1]), int(ok.sum())


def statistic(a, idx):
    per = {name: [] for name in LABELS}
    for c in COMPONENT_COHORTS:
        test = np.asarray(a["cohort"]) == idx[c]
        post, _ = fit.fit_svi(fit.subset(a, ~test), K_KEPT)
        jax.clear_caches()
        te = fit.subset(a, test)
        lam = np.asarray(fit._lam(post, K_KEPT))  # draws x labels x K
        W = np.asarray(post["W"])  # draws x P x K
        x = np.asarray(te["x"])
        lab = observed_labels(te)
        for j, name in enumerate(LABELS):
            pred = np.einsum("np,dpk,dk->n", x, W, lam[:, j, :]) / W.shape[0]
            others = [lab[o] for o in LABELS if o != name]
            rho, n = partial_spearman(pred, lab[name], others)
            per[name].append((c, rho, n))
    out = {}
    for name, rows in per.items():
        num = sum(r * n for _, r, n in rows if np.isfinite(r)); den = sum(n for _, r, n in rows if np.isfinite(r))
        out[name] = {"stat": num / den if den else np.nan, "by_cohort": {c: {"rho": r, "n": n} for c, r, n in rows}}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["observed", "null", "summary"])
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--b3-dir", required=True, help="B3 directory with fit_full_K1.npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--n", type=int, default=1)
    a_ = ap.parse_args()
    check_prespec()
    out = Path(a_.out); out.mkdir(parents=True, exist_ok=True)
    if a_.mode == "null" and a_.n > 1:  # one process per draw (XLA compiled code accumulates in a process)
        for i in range(a_.start, a_.start + a_.n):
            if not (out / "null" / f"sim_{i:04d}.json").exists():
                subprocess.run([sys.executable, __file__, "null", "--inputs", a_.inputs, "--b3-dir", a_.b3_dir,
                                "--out", str(out), "--start", str(i), "--n", "1"], check=True)
        return
    if a_.mode == "summary":
        obs = json.loads((out / "observed.json").read_text())
        sims = [json.loads(p.read_text()) for p in sorted((out / "null").glob("sim_*.json"))]
        if len(sims) != N_NULL:
            raise SystemExit(f"null incomplete: {len(sims)} of {N_NULL}")
        rows = {}
        for name in LABELS:
            null = np.array([s[name] for s in sims], float)
            p = (1 + np.sum(null >= obs[name]["stat"])) / (1 + len(null))
            rows[name] = {"observed": obs[name]["stat"], "by_cohort": obs[name]["by_cohort"],
                          "null_mean": float(np.nanmean(null)), "null_q95": float(np.nanquantile(null, 0.95)), "p": float(p)}
        ps = np.array([rows[n]["p"] for n in LABELS]); o = np.argsort(ps)
        q = np.empty(len(ps)); q[o] = np.minimum.accumulate((ps[o] * len(ps) / np.arange(1, len(ps) + 1))[::-1])[::-1]
        for n, qq in zip(LABELS, q):
            rows[n]["q_bh_four_labels"] = float(min(qq, 1.0)); rows[n]["pass"] = bool(qq < 0.05)
        rows["p_floor"] = 1 / (N_NULL + 1)
        (out / "summary.json").write_text(json.dumps(rows, indent=2)); print(json.dumps(rows, indent=2))
        return
    d, prog = fit.load_inputs(a_.inputs)
    d = d[d["cohort"].isin(DEV)].reset_index(drop=True)
    a = fit.arrays(d, prog, DEV)
    idx = {c: i for i, c in enumerate(DEV)}
    if a_.mode == "observed":
        res = statistic(a, idx)
        (out / "observed.json").write_text(json.dumps(res, indent=2, default=float)); print(json.dumps(res, indent=2, default=float))
        return
    k1 = dict(np.load(Path(a_.b3_dir) / f"fit_full_K{K_KEPT - 1}.npz"))
    i = a_.start
    rng = np.random.default_rng(fit.SEED + 50000 + i)
    res = statistic(fit.simulate_labels(k1, K_KEPT - 1, a, rng), idx)
    (out / "null").mkdir(exist_ok=True)
    (out / "null" / f"sim_{i:04d}.json").write_text(json.dumps({"sim": i, **{n: res[n]["stat"] for n in LABELS}}))


if __name__ == "__main__":
    main()
