"""B-MODEL on the B-FIX shared fixture (spec v1 section 11): log-likelihood and gradient at
theta1..4, the section-3 and section-4 fits, and P1, P2, P3 on the supplied draw files.

Usage: run_fixture.py --fixture <model-a-bfix dir> --out <new dir> [--quick]
"""
import argparse
import csv
import json
import time
from pathlib import Path

import numpy as np

import drivers as dv
import fitting as ft
import fixture_io as fio
import model_a as ma


def fit_record(fit, L):
    return dict(loglik=fit.loglik, theta=L.to_dict(fit.theta), converged=fit.converged,
                proj_max_grad=fit.proj_max_grad, attempts=fit.attempts, seconds=fit.seconds,
                at_bound=[L.names[i] for i in fit.at_bound], lbfgs_message=fit.lbfgs_message)


def write_tsv(path, recs, cols):
    if path.exists():
        raise FileExistsError(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in recs:
            w.writerow(["NA" if r.get(c) is None else (repr(r[c]) if isinstance(r[c], float)
                                                        else r[c]) for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--quick", action="store_true", help="stop after the observed fits")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = "quick_" if args.quick else ""
    if (out / f"{stem}results.json").exists():
        raise SystemExit(f"{out} already holds results; write to a new directory")
    t_all = time.perf_counter()

    prob, glob = fio.load_problem(args.fixture)
    params = fio.load_params(args.fixture, prob.L)
    L = prob.L
    D = prob.data()
    tau_max, Delta = glob["tau_max_fixture"], glob["Delta"]
    v_max = tau_max ** 2
    res = dict(fixture=str(Path(args.fixture).resolve()), design_columns=prob.design_columns,
               coordinates=L.names, n_rows=int(D.nr), n_individuals=len(prob.ind_ids),
               n_clusters=int(len(np.unique(prob.cluster))), tau_max=tau_max, Delta=Delta)

    # 1. log-likelihood and gradient at the parameter files
    res["at_params"] = {}
    for name, th in params.items():
        if th[L.v] > 0:
            ll, g = ma.loglik4(th, D, 1)
            ll40 = ma.loglik4(th, D, 0, nodes=40)
            rec = dict(form="section 4 (adaptive GH, 20 nodes)", v=float(th[L.v]), loglik=ll,
                       loglik_gh40=ll40, gradient=L.to_dict(g))
        else:
            ll, g = ma.loglik3(th, D, 1)
            g[L.v] = ma.boundary_score(th, D)
            rec = dict(form="section 3 (v = 0)", v=0.0, loglik=ll, gradient=L.to_dict(g),
                       gradient_note="the v entry is the boundary score dl/dv at v = 0")
        res["at_params"][name] = rec
        print(name, rec["form"], repr(ll), flush=True)

    # 2. section-3 MLE from the section 7 starts (timed)
    start = ft.start_values(D)
    start_no_omega = ft.start_values(D, include_omega_c=False)
    fit3 = ft.fit_section3(D)
    z, se, cond = dv.kappa_z(fit3.theta, D)
    alpha_hat = fit3.theta[L.alpha]
    res["section3_mle"] = dict(fit_record(fit3, L), kappa_hat=float(fit3.theta[L.kappa]),
                               se_sandwich=se, z=z, info_condition_number=cond,
                               start_alpha=start[L.alpha].tolist(),
                               start_alpha_if_omega_c_zero=start_no_omega[L.alpha].tolist(),
                               tau_max_rule_value=float(Delta * np.median(np.abs(alpha_hat))))
    print("section3", repr(fit3.loglik), fit3.theta[L.kappa], fit3.converged, fit3.seconds, flush=True)

    # 3. section-4 fits from v = 0 and v = tau_max^2 (timed), LR
    best, f0, f1 = ft.fit_section4_two_starts(D, fit3.theta, v_max)
    LR = ft.lr_statistic(fit3, best)
    res["section4"] = dict(v_hat=float(best.theta[L.v]), tau_hat=float(np.sqrt(best.theta[L.v])),
                           loglik=best.loglik, LR=LR,
                           kappa_hat_section4=float(best.theta[L.kappa]),
                           boundary_score_at_section3_mle=ma.boundary_score(fit3.theta, D),
                           start_kept="v0" if best is f0 else "vmax",
                           fit_start_v0=fit_record(f0, L), fit_start_vmax=fit_record(f1, L))
    if best.theta[L.v] > 0:
        res["section4"]["loglik_gh40_at_mle"] = ma.loglik4(best.theta, D, 0, nodes=40)
    print("section4", best.theta[L.v], repr(best.loglik), LR, f0.seconds, f1.seconds, flush=True)
    res["timings_seconds"] = dict(section3_fit_cold_start=fit3.seconds,
                                  section4_fit_start_v0=f0.seconds,
                                  section4_fit_start_vmax=f1.seconds)
    if args.quick:
        (out / f"{stem}results.json").write_text(json.dumps(res, indent=1))
        return

    draws = fio.load_draws(args.fixture, prob)
    res["draw_files"] = dict(p2_strata_match_rule=draws["strata_match_rule"],
                             B_p1=len(draws["Sstar"]), B_p2=len(draws["mult"]),
                             B_p3=int(draws["U"].shape[0]))
    if not fit3.converged:
        raise SystemExit("observed section-3 fit failed: P1 not computed (spec 5.1)")
    # the observed section-3 fit from 20 further random starts (multimodality check)
    rng = np.random.default_rng(np.random.SeedSequence([20260923, 91]))
    probe = []
    for k in range(20):
        s = ft.start_values(D)
        s[L.alpha] = rng.uniform(-1.5, 1.5, L.G)
        s[L.r] = rng.uniform(-6.0, -2.0, L.G)
        s[L.kappa] = rng.uniform(-0.5, 0.5)
        s[[L.rho_S, L.rho_b, L.rho_d, L.omega_S]] = rng.uniform(-1.0, 1.0, 4)
        f = ft.fit_section3(D, start=s)
        probe.append(dict(start=k, loglik=f.loglik, kappa=float(f.theta[L.kappa]),
                          converged=f.converged, at_bound=[L.names[i] for i in f.at_bound]))
    res["section3_mle"]["random_start_probe"] = dict(
        best_loglik=max(p["loglik"] for p in probe if p["converged"]),
        n_higher_than_section7_start=sum(p["converged"] and p["loglik"] > fit3.loglik + 1e-6
                                         for p in probe), fits=probe)

    t = time.perf_counter()
    p1 = dv.run_p1(prob, fit3, draws["Sstar"])
    p1["seconds"] = time.perf_counter() - t
    print("P1", p1["z"], p1["p_warm"], p1["p_best"], p1["seconds"], flush=True)
    t = time.perf_counter()
    p2 = dv.run_p2(prob, fit3, best, draws["mult"], draws["lambdas"], tau_max, Delta, v_max)
    p2["seconds"] = time.perf_counter() - t
    print("P2", p2["p_warm"], p2["p_best"], p2["seconds"], flush=True)
    t = time.perf_counter()
    p3 = dv.run_p3(prob, fit3, LR, draws["U"], v_max)
    p3["seconds"] = time.perf_counter() - t
    print("P3", p3["p_warm"], p3["p_best"], p3["seconds"], flush=True)
    res["P1"], res["P2"], res["P3"] = p1, p2, p3
    res["total_seconds"] = time.perf_counter() - t_all
    (out / "results.json").write_text(json.dumps(res, indent=1))

    write_tsv(out / "p1_draws.tsv", p1["draws"],
              ["draw", "kappa_warm", "se_warm", "z_warm", "exceed_warm", "in_band_warm",
               "converged_warm", "loglik_warm", "proj_max_grad_warm", "loglik_cold", "kappa_cold",
               "cold_better", "kappa_best", "z_best", "exceed_best"])
    srcs = sorted(p2["sources"])
    p2rows = []
    for r0 in p2["draws"]:
        r = dict(r0)
        b = r["draw"]
        for s in srcs:
            for pol in ("warm", "best"):
                comp = p2["sources"][s][pol]
                r[f"ratio_{pol}_{s}"] = comp["ratio_draws"][b]
                r[f"up_{pol}_{s}"] = comp["up_indicators"][b]
                r[f"lo_{pol}_{s}"] = comp["lo_indicators"][b]
        p2rows.append(r)
    write_tsv(out / "p2_draws.tsv", p2rows,
              ["draw", "kappa_warm", "loglik_warm", "converged_warm", "kappa_cold", "loglik_cold",
               "cold_better", "kappa_best"]
              + [f"{x}_{pol}_{s}" for pol in ("warm", "best") for s in srcs
                 for x in ("ratio", "up", "lo")])
    write_tsv(out / "p3_draws.tsv", p3["draws"],
              ["draw", "LR", "exceed", "in_band", "loglik3", "loglik4", "v_hat", "kappa3",
               "start_kept", "converged", "n_class_H", "n_class_R", "n_class_A", "sum_a",
               "loglik3_cold", "kappa3_cold", "cold_better", "LR_best", "v_hat_best",
               "exceed_best"])


if __name__ == "__main__":
    main()
