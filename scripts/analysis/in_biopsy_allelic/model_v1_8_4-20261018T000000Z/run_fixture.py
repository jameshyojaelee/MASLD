"""B-MODEL on the B-FIX shared fixture (spec v1 section 11, rerun under v1.8 as a new labelled
family): log-likelihood and gradient at theta1..4, the section-3 and section-4 fits, and P1, P2,
P3 on the supplied draw files.

v1.8.4: observed section-3 fit from K_obs = 3 starts (probe on stream 10, [20260923, 13, 0, 0, 1]); every refit at
K = 3; the observed and P3 section-4 pairs of item 6.1; LR from l3* (items 6.3-6.4). Every fit end point (each start
of each refit, each pair member, each observed start, the P2 profile fits and the 20-start probe) is written to
fits.tsv with its full theta (repr doubles, which round-trip) and its l values.

Usage: run_fixture.py --fixture <model-a-bfix dir> --out <new dir> [--quick] [--only p1,p2,p3]
"""
import argparse
import csv
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np

import drivers as dv
import fitting as ft
import fixture_io as fio
import model_a as ma


def fit_record(fit, L):
    return dict(loglik=fit.loglik, loglik_ref=fit.loglik_ref, theta=L.to_dict(fit.theta),
                converged=fit.converged, status=fit.status, hessian_ok=fit.hessian_ok,
                proj_max_grad=fit.proj_max_grad, attempts=fit.attempts, seconds=fit.seconds,
                at_bound=[L.names[i] for i in fit.at_bound], lbfgs_message=fit.lbfgs_message,
                repolished=fit.repolished, n_reference_genes=fit.n_reference_genes,
                v_at_upper=fit.v_at_upper, events=fit.events)


def run_header():
    """v1.8.3a (H-12 m6): host, CPU model, glibc, thread variables, backend and versions in results.json."""
    import ctypes
    import scipy
    model = next((ln.split(":", 1)[1].strip() for ln in open("/proc/cpuinfo") if ln.startswith("model name")), None)
    try:
        libc = ctypes.CDLL("libc.so.6")
        libc.gnu_get_libc_version.restype = ctypes.c_char_p
        glibc = libc.gnu_get_libc_version().decode()
    except Exception as err:  # noqa: BLE001
        glibc = repr(err)
    head = dict(argv=sys.argv, host=platform.node(), cpu_model=model, glibc=glibc, python=sys.version.split()[0],
                numpy=np.__version__, scipy=scipy.__version__, python_executable=sys.executable,
                backend=getattr(ma, "BACKEND", "cpu"),
                env={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                                     "SLURM_JOB_ID", "MODEL_A_TORCH_DEVICE", "CUDA_VISIBLE_DEVICES")},
                started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    if head["backend"] == "torch":
        import backend_torch as bt
        head["torch_env"] = bt.environment()
    return head


FIT_COLS = ["test", "phase", "draw", "role", "kind", "v_start", "loglik", "loglik_ref", "l3_at_v0", "converged",
            "status", "hessian_ok", "proj_max_grad", "attempts", "repolished", "v_at_upper", "n_reference_genes",
            "at_bound", "attempt_log", "seconds"]


def fit_row(test, phase, draw, role, fit, L, extra=None):
    extra = extra or {}
    return dict(test=test, phase=phase, draw=draw, role=role, kind=extra.get("kind", ""),
                v_start=extra.get("v_start"), loglik=float(fit.loglik), loglik_ref=float(fit.loglik_ref),
                l3_at_v0=extra.get("l3_at_v0"), converged=bool(fit.converged), status=str(fit.status),
                hessian_ok=bool(fit.hessian_ok), proj_max_grad=float(fit.proj_max_grad), attempts=int(fit.attempts),
                repolished=bool(fit.repolished), v_at_upper=bool(fit.v_at_upper),
                n_reference_genes=int(fit.n_reference_genes), at_bound=";".join(L.names[i] for i in fit.at_bound),
                attempt_log=json.dumps([[float(a), bool(c)] for a, c in fit.attempt_log]),
                seconds=float(fit.seconds), theta=np.array(fit.theta, float).copy())


def write_fits(path, rows, L):
    if path.exists():
        raise FileExistsError(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(FIT_COLS + [f"theta:{n}" for n in L.names])
        for r in rows:
            w.writerow(["NA" if r[c] is None else (repr(float(r[c])) if isinstance(r[c], (float, np.floating))
                                                   else r[c]) for c in FIT_COLS]
                       + [repr(float(x)) for x in r["theta"]])


def write_tsv(path, recs, cols):
    if path.exists():
        raise FileExistsError(path)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in recs:
            w.writerow(["NA" if r.get(c) is None else (repr(r[c]) if isinstance(r[c], float)
                                                        else r[c]) for c in cols])


START_COLS = [c for n in ("warm", "cold", "probe1")
              for c in (f"loglik_{n}", f"dl_{n}_vs_warm", f"converged_{n}", f"attempts_{n}", f"seconds_{n}",
                        f"status_{n}")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--quick", action="store_true", help="stop after the observed fits")
    ap.add_argument("--only", default="p1,p2,p3", help="comma list of tests to run")
    args = ap.parse_args()
    only = set(args.only.split(","))
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
    res = dict(header=run_header(), fixture=str(Path(args.fixture).resolve()), design_columns=prob.design_columns,
               coordinates=L.names, n_rows=int(D.nr), n_individuals=len(prob.ind_ids),
               n_clusters=int(len(np.unique(prob.cluster))), tau_max=tau_max, Delta=Delta)

    # 1. log-likelihood and gradient at the parameter files
    res["at_params"] = {}
    for name, th in params.items():
        if th[L.v] > 0:
            ll, g, det = ma.loglik4(th, D, 1, detail=True)
            ll40 = ma.loglik4(th, D, 0, nodes=40)
            llref = ma.loglik4(th, D, 0, all_reference=True)
            rec = dict(form="section 4 (adaptive GH, 20 nodes, v1.8 mode search)", v=float(th[L.v]),
                       loglik=ll, loglik_gh40=ll40, loglik_reference_all_genes=llref,
                       integrator=[int(x) for x in det["integrator"]],
                       u_hat=[float(x) for x in det["u_hat"]], gradient=L.to_dict(g))
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
    test = ",".join(sorted(only))
    FR = []   # every fit end point
    fit3 = ft.fit_section3_observed(D, ft.FIXTURE_OBS_START_SEED,
                                    ft.probe_entropy_fixed_observed(1))  # v1.8.4 item 5: K_obs = 3
    for nm, f in zip(("start0_section7", "start1_jitter", "start2_probe1"), fit3.multistart_fits):
        FR.append(fit_row(test, "observed", "", f"observed_s3/{nm}", f, L, dict(kind="section3", v_start=0.0)))
    z, se, cond = dv.kappa_z(fit3.theta, D)
    alpha_hat = fit3.theta[L.alpha]
    res["section3_mle"] = dict(fit_record(fit3, L), kappa_hat=float(fit3.theta[L.kappa]),
                               se_sandwich=se, z=z, info_condition_number=cond,
                               start_alpha=start[L.alpha].tolist(),
                               start_alpha_if_omega_c_zero=start_no_omega[L.alpha].tolist(),
                               tau_max_rule_value=float(Delta * np.median(np.abs(alpha_hat))),
                               observed_multistart=fit3.multistart)
    print("section3", repr(fit3.loglik), fit3.theta[L.kappa], fit3.converged, fit3.seconds, flush=True)

    # 3. v1.8.4 item 6: the observed section-4 pair from the reported section-3 theta, l3*, LR
    pr = dv.section4_pair_lr(D, fit3.theta, fit3.loglik, v_max)
    best, f0, f1, rule4 = pr["kept4"], pr["member0"], pr["member1"], pr["rule4"]
    for i, (nm, m) in enumerate((("v0", f0), ("vmax", f1))):
        FR.append(fit_row(test, "observed", "", f"observed_pair/{nm}", m, L,
                          dict(kind="section4", v_start=(0.0, v_max)[i], l3_at_v0=pr["l3star_contrib"][i][1])))
    if best is None:
        write_fits(out / "fits.tsv", FR, L)
        (out / "results_partial.json").write_text(json.dumps(res, indent=1, default=float))
        raise SystemExit(f"observed section-4 pair failed ({rule4}): P2/P3 not computed (v1.8.4 item 6.2)")
    LR = pr["LR"]
    res["section4"] = dict(v_hat=float(best.theta[L.v]), tau_hat=float(np.sqrt(best.theta[L.v])),
                           loglik=best.loglik, loglik_ref=best.loglik_ref, LR=LR,
                           l3=fit3.loglik, l3star=pr["l3star"], l3star_minus_l3=pr["l3star"] - fit3.loglik,
                           missed_v0=pr["missed_v0"], l3star_contrib=pr["l3star_contrib"],
                           LR_gh=pr["LR_gh"], converged=best.converged,
                           kappa_hat_section4=float(best.theta[L.kappa]),
                           boundary_score_at_section3_mle=ma.boundary_score(fit3.theta, D),
                           start_kept="v0" if best is f0 else "vmax", rule=rule4,
                           v_at_upper=best.v_at_upper,
                           fit_start_v0=fit_record(f0, L), fit_start_vmax=fit_record(f1, L))
    if best.theta[L.v] > 0:
        res["section4"]["loglik_gh40_at_mle"] = ma.loglik4(best.theta, D, 0, nodes=40)
    print("section4", best.theta[L.v], repr(best.loglik), LR, f0.seconds, f1.seconds, flush=True)
    # v1.8.3a (H-12 m4): section3_fit_cold_start is the reported start's own time (Fit.seconds);
    # section3_observed_fit_both_starts is the whole two-start observed fit.
    res["timings_seconds"] = dict(section3_fit_cold_start=fit3.seconds,
                                  section3_observed_fit_both_starts=fit3.multistart["seconds_total"],
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
        write_fits(out / "fits.tsv", FR, L)
        (out / "results_partial.json").write_text(json.dumps(res, indent=1, default=float))
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
        FR.append(fit_row(test, "observed_probe", k, "random_start_probe", f, L, dict(kind="section3", v_start=0.0)))
        probe.append(dict(start=k, loglik=f.loglik, kappa=float(f.theta[L.kappa]),
                          converged=f.converged, at_bound=[L.names[i] for i in f.at_bound]))
    res["section3_mle"]["random_start_probe"] = dict(
        best_loglik=max(p["loglik"] for p in probe if p["converged"]),
        n_higher_than_section7_start=sum(p["converged"] and p["loglik"] > fit3.loglik + 1e-6
                                         for p in probe), fits=probe)

    if "p1" in only:
        t = time.perf_counter()
        p1 = dv.run_p1(prob, fit3, draws["Sstar"])
        p1["seconds"] = time.perf_counter() - t
        FR += [fit_row(test, "P1", b, role, f, L, ex) for b, role, f, ex in p1.pop("_fits")]
        res["P1"] = p1
        print("P1", p1["z"], p1["p"], p1["failed_draws"], p1["seconds"], flush=True)
        write_tsv(out / "p1_draws.tsv", p1["draws"],
                  ["draw", "kappa", "se", "z", "exceed", "in_band", "failed", "kept", "kept_index", "rule",
                   "probe_entropy"] + START_COLS)
    if "p2" in only:
        t = time.perf_counter()
        p2 = dv.run_p2(prob, fit3, best, draws["mult"], draws["lambdas"], tau_max, Delta, v_max)
        p2["seconds"] = time.perf_counter() - t
        FR += [fit_row(test, "P2", b, role, f, L, ex) for b, role, f, ex in p2.pop("_fits")]
        res["P2"] = p2
        print("P2", p2["p"], p2["failed_draws"], p2["seconds"], flush=True)
        srcs = sorted(p2["sources"])
        p2rows = []
        for r0 in p2["draws"]:
            r = dict(r0)
            b = r["draw"]
            for s_ in srcs:
                comp = p2["sources"][s_]
                r[f"ratio_{s_}"] = comp["ratio_draws"][b]
                r[f"up_{s_}"] = comp["up_indicators"][b]
                r[f"lo_{s_}"] = comp["lo_indicators"][b]
            p2rows.append(r)
        write_tsv(out / "p2_draws.tsv", p2rows,
                  ["draw", "kappa", "failed", "kept", "kept_index", "rule", "probe_entropy"]
                  + [f"kappa_{n}" for n in ft.START_NAMES] + START_COLS
                  + [f"{x}_{s_}" for s_ in srcs for x in ("ratio", "up", "lo")])
    if "p3" in only:
        t = time.perf_counter()
        p3 = dv.run_p3(prob, fit3, LR, draws["U"], v_max)
        p3["seconds"] = time.perf_counter() - t
        FR += [fit_row(test, "P3", b, role, f, L, ex) for b, role, f, ex in p3.pop("_fits")]
        counts = p3.pop("_counts")
        res["P3"] = p3
        print("P3", p3["p"], p3["failed_draws"], p3["seconds"], flush=True)
        write_tsv(out / "p3_draws.tsv", p3["draws"],
                  ["draw", "LR", "LR_gh", "exceed", "in_band", "failed", "probe_entropy", "rule3", "kept3",
                   "kept3_index", "loglik3", "kappa3", "rule4", "kept4", "loglik4_ref", "loglik4_gh", "v_hat",
                   "v_at_upper", "l3star", "l3star_minus_l3", "missed_v0", "n_reference_genes4", "repolished4"]
                  + ["s3_" + c for c in START_COLS]
                  + [f"m{i}_{c}" for i in (0, 1) for c in ("v_hat", "loglik_ref", "loglik_gh", "converged",
                                                           "attempts", "seconds", "status", "l3_at_v0")]
                  + ["n_class_H", "n_class_R", "n_class_A", "sum_a"])
        R = prob.rows
        with open(out / "p3_counts.tsv", "x", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["draw", "row", "individual_id", "gene", "n", "a_observed", "a_star", "class_HRA"])
            for b, (a_b, c_b) in enumerate(counts):
                for j in range(len(a_b)):
                    w.writerow([b, j, prob.ind_ids[prob.row_ind[j]], L.gene_names[prob.row_gene[j]], int(R.n[j]),
                                int(R.a[j]), int(a_b[j]), "HRA"[int(c_b[j])]])
    # v1.8 section 1 item 10: fits with v_hat within 1e-8 of 1, with the gene count
    G_fit = int(D.n_genes_with_rows)
    vb = []
    if best.v_at_upper:
        vb.append(dict(fit="observed section 4", genes=G_fit))
    for d in res.get("P3", {}).get("draws", []):
        if d.get("v_at_upper"):
            vb.append(dict(fit=f"P3 draw {d['draw']} section 4", genes=G_fit, LR=d["LR"]))
    res["v_upper_bound_fits"] = vb
    write_fits(out / "fits.tsv", FR, L)
    res["n_fit_rows"] = len(FR)
    res["total_seconds"] = time.perf_counter() - t_all
    (out / "results.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
