"""Time one section-3 and one section-4 fit on a synthetic data set of the planned real size
(558 individuals in five cohorts, 400 genes, about 5x10^4 included rows). Synthetic only: depths,
genotypes and counts are drawn here; no MASLD data are read.

Usage: scale_timing.py --out <dir> [--genes 400] [--seed 20260930]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.special import expit, logit

import fitting as ft
import model_a as ma
from het_rule import het

COHORTS = [("A", 208, (10, 1, 10)), ("B", 160, (20, 2, 10)), ("C", 83, (10, 2, 15)),
           ("D", 54, (30, 3, 5)), ("E", 53, (10, 1, 5))]


def simulate(G, seed):
    rng = np.random.default_rng(seed)
    coh = np.concatenate([[c] * m for c, m, _ in COHORTS])
    rule_of = {c: r for c, _, r in COHORTS}
    n_ind = len(coh)
    S = rng.integers(0, 4, n_ind).astype(float)
    zraw = np.column_stack([rng.normal(-6, 0.4, n_ind), rng.normal(0.25, 0.08, n_ind),
                            rng.normal(0, 0.003, n_ind), rng.integers(0, 2, n_ind),
                            rng.normal(5, 1, n_ind), rng.normal(0.003, 0.002, n_ind)])
    St = ma.centre_stage(coh, S, np.ones(n_ind))
    alpha = np.clip(rng.normal(0, 0.6, G), -3, 3)
    u = rng.normal(0, 0.05, G)
    p_alt = rng.uniform(0.05, 0.5, G)
    orient_g = rng.choice([-1, 1], G)
    ii, gg = np.meshgrid(np.arange(n_ind), np.arange(G), indexing="ij")
    ii, gg = ii.ravel(), gg.ravel()
    n = np.minimum(6821, np.round(np.exp(rng.normal(np.log(30), 1.2, len(ii))))).astype(int)
    pr, pa = 1 - p_alt[gg], p_alt[gg]
    w = np.column_stack([2 * pr * pa, pr ** 2, pa ** 2])
    cls = (rng.random(len(ii))[:, None] > np.cumsum(w, 1)).sum(1)
    e = rng.uniform(5e-4, 5e-3, n_ind)[ii]
    s_t = -orient_g[gg]
    x1 = alpha[gg] * (1 - 0.08 * St[ii]) + u[gg] * St[ii] + s_t * 0.1
    p_het = rng.beta(expit(x1) * 32, expit(-x1) * 32)
    alt_frac = np.where(cls == 0, np.where(orient_g[gg] == 1, p_het, 1 - p_het),
                        np.where(cls == 1, e, 1 - e))
    alt = rng.binomial(n, rng.beta(alt_frac * 22, (1 - alt_frac) * 22))
    ref = n - alt
    rule = np.array([rule_of[c] for c in coh[ii]])
    keep = het(ref, alt, rule[:, 0], rule[:, 1], rule[:, 2])
    ii, gg, n, ref, alt, rule, e = ii[keep], gg[keep], n[keep], ref[keep], alt[keep], rule[keep], e[keep]
    orient = orient_g[gg]
    a = np.where(orient == 1, alt, ref)
    R = ma.Rows(n, a, rule[:, 0], rule[:, 1], rule[:, 2], orient, 1 - p_alt[gg], p_alt[gg],
                np.zeros(len(n)), e, 0.0439)
    ids = [f"I{i:04d}" for i in range(n_ind)]
    prob = ma.Problem(ids, coh, S, zraw, ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female",
                                          "z_age10", "z_c"], np.array(ids), R, ii, gg,
                      [f"S{g:04d}" for g in range(G)], np.full(len(n), 0.1))
    return prob


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20260930)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    t = time.perf_counter()
    prob = simulate(args.genes, args.seed)
    D = prob.data()
    L = prob.L
    rec = dict(genes=args.genes, rows=int(D.nr), individuals=len(prob.ind_ids), coordinates=L.size,
               max_n=int(D.rows.n.max()), tail_terms=int(len(D.rows.tail_a)),
               simulate_seconds=time.perf_counter() - t)
    th = ft.start_values(D)
    for order in (0, 1, 2):
        t = time.perf_counter()
        ma.loglik3(th, D, order)
        rec[f"loglik3_order{order}_seconds"] = time.perf_counter() - t
    print(rec, flush=True)
    fit3 = ft.fit_section3(D)
    rec.update(section3_seconds=fit3.seconds, section3_converged=fit3.converged,
               section3_attempts=fit3.attempts, section3_loglik=fit3.loglik,
               section3_kappa=float(fit3.theta[L.kappa]), section3_proj_max_grad=fit3.proj_max_grad)
    print(rec, flush=True)
    th4 = fit3.theta.copy()
    th4[L.v] = 0.05 ** 2
    for order in (0, 1, 2):
        t = time.perf_counter()
        ma.loglik4(th4, D, order)
        rec[f"loglik4_order{order}_seconds"] = time.perf_counter() - t
    print(rec, flush=True)
    fit4 = ft.fit_section4(D, fit3.theta, 0.05 ** 2, 0.05 ** 2)
    rec.update(section4_seconds=fit4.seconds, section4_converged=fit4.converged,
               section4_attempts=fit4.attempts, section4_v=float(fit4.theta[L.v]),
               section4_loglik=fit4.loglik, section4_proj_max_grad=fit4.proj_max_grad,
               tau_planted_sd=0.05, kappa_planted=-0.08, logit_rho_planted=float(logit(1 / 33)))
    print(rec, flush=True)
    (out / "scale_timing.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
