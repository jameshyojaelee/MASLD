"""Time Model A fits on synthetic data at the v1.8 section 4 planning bounds (Gate B cost input).
Synthetic only: depths, genotypes and counts are drawn here; no MASLD data are read.

Design at the bound: 558 individuals in five cohorts (208, 160, 83, 54, 53), G genes, every
individual included for every gene (558 rows per gene, the v1.8 rows-per-gene bound), 7 base z
columns each partly missing in some cohort (design width 14 = 7 + 7 missing indicators), depths
log-normal around 30 with a tail capped at the T4 limit n = 10,000 (one row per gene at 10,000).

Timed: one section-3 fit; section-4 evaluations at v = tau_max^2 (fresh state, so the v1.8 global
check and the reference comparison run for every gene; then warmed; then every gene on the
reference integrator); per-gene global-check and reference costs on 10 genes at v = tau_max^2,
0.25 and 1 (extrapolated to all genes, an estimate); one full section-4 fit from v = tau_max^2.

Usage: scale_timing.py --out <dir> [--genes 400] [--seed 20261015] [--tau-max 0.064373]
"""
import argparse
import json
import platform
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import expit, logit

import fitting as ft
import model_a as ma
from het_rule import truncation_bounds
from betabinom import bb_logpmf

COHORTS = [("A", 208, (10, 1, 10)), ("B", 160, (20, 2, 10)), ("C", 83, (10, 2, 15)),
           ("D", 54, (30, 3, 5)), ("E", 53, (10, 1, 5))]
BASE_Z = ["z_b_log_e", "z_d_dup", "z_ffpe", "z_sex_female", "z_age10", "z_anc1", "z_anc2"]
N_MAX = 10_000


def maxrss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2


def simulate(G, seed, tau_sd=0.05, kappa=-0.08):
    rng = np.random.default_rng(np.random.SeedSequence([seed, G]))
    coh = np.concatenate([[c] * m for c, m, _ in COHORTS])
    rule_of = {c: r for c, _, r in COHORTS}
    n_ind = len(coh)
    S = rng.integers(0, 3, n_ind).astype(float)          # F0-1 / F2 / F3-4
    zraw = np.column_stack([rng.normal(-6, 0.4, n_ind), rng.normal(0.25, 0.08, n_ind),
                            rng.integers(0, 2, n_ind).astype(float),
                            rng.integers(0, 2, n_ind).astype(float), rng.normal(5, 1, n_ind),
                            rng.normal(0, 1, n_ind), rng.normal(0, 1, n_ind)])
    # each base column partly missing in one cohort -> one missing indicator per column
    for j in range(zraw.shape[1]):
        c = COHORTS[j % len(COHORTS)][0]
        sel = np.flatnonzero(coh == c)
        zraw[rng.choice(sel, max(2, len(sel) // 5), replace=False), j] = np.nan
    alpha = np.clip(rng.normal(0, 0.6, G), -3, 3)
    u = rng.normal(0, tau_sd, G)
    p_alt = rng.uniform(0.05, 0.5, G)
    orient_g = rng.choice([-1, 1], G)
    ii, gg = np.meshgrid(np.arange(n_ind), np.arange(G), indexing="ij")
    ii, gg = ii.ravel(), gg.ravel()
    order = np.lexsort((gg, ii))
    ii, gg = ii[order], gg[order]
    rule = np.array([rule_of[c] for c in coh[ii]])
    n = np.round(np.exp(rng.normal(np.log(30), 1.2, len(ii)))).astype(np.int64)
    n = np.clip(n, rule[:, 0].max(), N_MAX)
    first = np.flatnonzero(ii == 0)
    n[first] = N_MAX                                      # one row per gene at the depth limit
    St = ma.centre_stage(coh, S, np.ones(n_ind))
    s_t = -orient_g[gg]
    x1 = alpha[gg] * (1 + kappa * St[ii]) + u[gg] * St[ii] + s_t * 0.1
    lo, hi, _ = truncation_bounds(n, rule[:, 0], rule[:, 1], rule[:, 2])
    mu, rho = expit(x1), 1 / 33.0
    a = np.empty(len(n), np.int64)
    U01 = rng.random(len(n))
    for i in range(len(n)):                               # every row is a het call (H class)
        k = np.arange(lo[i], hi[i] + 1)
        lf = bb_logpmf(k, n[i], mu[i], rho)
        pm = np.exp(lf - lf.max())
        cdf = np.cumsum(pm) / pm.sum()
        a[i] = k[min(int(np.searchsorted(cdf, U01[i])), len(k) - 1)]
    R = ma.Rows(n, a, rule[:, 0], rule[:, 1], rule[:, 2], orient_g[gg], 1 - p_alt[gg], p_alt[gg],
                np.zeros(len(n)), np.full(len(n), 0.002), 0.0439)
    ids = [f"I{i:04d}" for i in range(n_ind)]
    prob = ma.Problem(ids, coh, S, zraw, BASE_Z, np.array(ids), R, ii, gg,
                      [f"S{g:04d}" for g in range(G)], np.full(len(n), 0.1))
    return prob


def timed(rec, key, f):
    t = time.perf_counter()
    try:
        out = f()
        rec[key + "_seconds"] = time.perf_counter() - t
        return out
    except Exception as err:            # noqa: BLE001 - record and continue
        rec[key + "_seconds"] = time.perf_counter() - t
        rec[key + "_error"] = repr(err)
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--genes", type=int, default=400)
    ap.add_argument("--seed", type=int, default=20261015)
    ap.add_argument("--tau-max", type=float, default=0.064373)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "scale_timing.json"
    if path.exists():
        raise SystemExit(f"{path} exists; write to a new directory")
    rec = dict(argv=sys.argv, host=platform.node(), seed=args.seed, tau_max=args.tau_max)

    def save():
        rec["max_rss_gb"] = maxrss_gb()
        path.write_text(json.dumps(rec, indent=1, default=float))
        print(json.dumps(rec, default=float), flush=True)

    prob = timed(rec, "simulate", lambda: simulate(args.genes, args.seed))
    D = prob.data()
    L = prob.L
    gsum = np.bincount(D.gene, D.rows.n + 1, L.G)
    rec.update(genes=args.genes, rows=int(D.nr), individuals=len(prob.ind_ids),
               design_width=L.p, coordinates=L.size, free_section3=int(D.free3.sum()),
               max_n=int(D.rows.n.max()), mean_n=float(D.rows.n.mean()),
               max_gene_sum_n_plus_1=int(gsum.max()), tail_terms=int(len(D.rows.tail_a)))
    save()
    th = ft.start_values(D)
    for order in (0, 1, 2):
        timed(rec, f"loglik3_order{order}", lambda: ma.loglik3(th, D, order))
    fit3 = timed(rec, "section3_fit", lambda: ft.fit_section3(D))
    rec.update(section3_converged=fit3.converged, section3_status=fit3.status,
               section3_attempts=fit3.attempts, section3_loglik=fit3.loglik,
               section3_kappa=float(fit3.theta[L.kappa]),
               section3_proj_max_grad=fit3.proj_max_grad)
    save()
    v_max = args.tau_max ** 2
    # full evaluations at v = tau_max^2: fresh state (first evaluation: global check and reference
    # comparison for every gene), then the same state warmed, then every gene on the reference
    th4 = fit3.theta.copy()
    th4[L.v] = v_max
    st = ma.Sec4State(L.G)
    res = timed(rec, "loglik4_vtaumax_fresh_order0",
                lambda: ma.loglik4(th4, D, 0, detail=True, state=st))
    if res is not None:
        det = res[1]
        rec["loglik4_vtaumax_reference_genes"] = int(np.sum(det["integrator"] == ma.REF))
        rec["loglik4_vtaumax_events_first"] = dict(st.events)
        for order in (0, 1, 2):
            timed(rec, f"loglik4_vtaumax_warm_order{order}",
                  lambda: ma.loglik4(th4, D, order, state=st))
        rec["loglik4_vtaumax_events_after_warm"] = dict(st.events)
    save()
    timed(rec, "loglik4_vtaumax_all_reference_order0",
          lambda: ma.loglik4(th4, D, 0, all_reference=True))
    save()
    # per-gene cost of the global check and the reference integrator at larger v, on the first
    # SUBSET genes (a full evaluation at v = 0.25 or 1 is estimated at many hours)
    SUBSET = 10
    for label, v in (("vtaumax", v_max), ("v0.25", 0.25), ("v1", 1.0)):
        th4 = fit3.theta.copy()
        th4[L.v] = v
        x1b, x2, _, _ = ma.predictors(th4, D)
        ctx = ma._Ctx(D, x1b, x2, v)
        stg = ma.Sec4State(L.G)
        per = []
        for g in range(SUBSET):
            r = dict(gene=g, rows=int(len(D.gene_idx(g))))
            t = time.perf_counter()
            h0, h10, h20 = ctx.h_gene(g, 0.0, 2)
            U = float(np.sqrt(2 * v * -h0[0]))
            hU = ctx.h_gene(g, np.array([-U, U]), 0)[0]
            uu, hh, hh2, ok, _ = ma._local_search(lambda x: ctx.h_gene(g, x, 2), [-U], [0.0], [U],
                                                  [hU[0]], h0, [hU[1]], h10, h20)
            cands = [(float(uu[0]), float(hh[0]), float(hh2[0]))] if ok[0] else []
            r["local_search_seconds"] = time.perf_counter() - t
            t = time.perf_counter()
            cands, cap = ma._global_check(ctx, g, U, cands, stg)
            r.update(global_check_seconds=time.perf_counter() - t, grid_cap=bool(cap), U=U,
                     grid_points=stg.events["max_grid_points"], n_modes=len(cands))
            best = ma._pick(cands)
            md = dict(modes={g: cands}, sigma={g: ma._sigma(best[2])}, h_hat={g: best[1]})
            r["sigma_hat"] = md["sigma"][g]
            t = time.perf_counter()
            try:
                rr = ma._reference(ctx, g, md, stg)
                r.update(reference_ell=rr["ell"], reference_panels=rr["panels"])
            except RuntimeError as err:
                r["reference_error"] = repr(err)
            r["reference_seconds"] = time.perf_counter() - t
            per.append(r)
            rec[f"per_gene_{label}"] = per
            save()
        n_ok = [x for x in per if "reference_error" not in x]
        rec[f"per_gene_{label}_summary"] = dict(
            genes=SUBSET, reference_failed=SUBSET - len(n_ok),
            mean_global_check_seconds=float(np.mean([x["global_check_seconds"] for x in per])),
            mean_reference_seconds=float(np.mean([x["reference_seconds"] for x in per])),
            extrapolated_first_evaluation_hours_all_genes=float(
                L.G * np.mean([x["global_check_seconds"] + x["reference_seconds"] for x in per])
                / 3600))
        save()
    fit4 = timed(rec, "section4_fit", lambda: ft.fit_section4(D, fit3.theta, v_max, v_max))
    if fit4 is not None:
        rec.update(section4_converged=fit4.converged, section4_status=fit4.status,
                   section4_attempts=fit4.attempts, section4_v=float(fit4.theta[L.v]),
                   section4_loglik=fit4.loglik, section4_loglik_ref=fit4.loglik_ref,
                   section4_proj_max_grad=fit4.proj_max_grad, section4_events=fit4.events,
                   section4_repolished=fit4.repolished,
                   section4_n_reference_genes=fit4.n_reference_genes)
    rec.update(tau_planted_sd=0.05, kappa_planted=-0.08, rho_planted=1 / 33.0)
    save()


if __name__ == "__main__":
    main()
