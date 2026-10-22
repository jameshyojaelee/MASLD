"""Spec section 4 accuracy test read with v1.8 section 1 item 1 (1e-8 relative in exp(ell_g) =
1e-8 absolute in ell_g). Five synthetic genes of 5, 20, 50, 200 and 500 rows (the B-FIX fixture
genes have 3-30 rows) and the 8 fixture genes, at tau in {0.25, 1, 4} x tau_max (fixture
tau_max 0.05): adaptive GH20 with the v1.8 mode, GH40, the v1.8 reference integrator, and
scipy.integrate.quad (epsabs 1e-12). Synthetic inputs and the B-FIX fixture only.

Usage: accuracy_v18.py --fixture <model-a-bfix dir> --out <new dir>
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

import fixture_io as fio
import model_a as ma

TAU_MAX = 0.05
TOL = 1e-8
HERE = Path(__file__).resolve().parent


def _synthetic():
    spec = importlib.util.spec_from_file_location("tmod", HERE / "tests" / "test_model_a.py")
    tmod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tmod)
    return tmod._synthetic_problem([5, 20, 50, 200, 500])


def table(D, th, names, label):
    L = D.L
    out = []
    for mult in (0.25, 1.0, 4.0):
        v = (mult * TAU_MAX) ** 2
        t = th.copy()
        t[L.v] = v
        _, d20 = ma.loglik4(t, D, 0, detail=True)
        _, d40 = ma.loglik4(t, D, 0, nodes=40, detail=True)
        _, dref = ma.loglik4(t, D, 0, detail=True, all_reference=True)
        for g in range(L.G):
            q, qerr = ma.gh_gene_integral_quad(t, D, g, v)
            r = dict(set=label, gene=names[g], rows=int(len(D.gene_idx(g))), tau_mult=mult, v=v,
                     u_hat=float(d20["u_hat"][g]), sigma_hat=float(d20["sigma"][g]),
                     integrator_v18=int(d20["integrator"][g]), gh20=float(d20["ell_gh"][g]),
                     gh40=float(d40["ell_gh"][g]), reference=float(dref["ell_g"][g]),
                     reference_panels=dref["reference_panels"].get(g), quad=float(q),
                     quad_abserr=float(qerr))
            r.update(gh20_minus_quad=r["gh20"] - q, reference_minus_quad=r["reference"] - q,
                     gh20_minus_gh40=r["gh20"] - r["gh40"])
            r["pass"] = bool(abs(r["gh20_minus_quad"]) <= TOL and abs(r["reference_minus_quad"]) <= TOL
                             and abs(r["gh20_minus_gh40"]) <= TOL)
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / "accuracy_v18.json").exists():
        raise SystemExit("accuracy_v18.json exists; write to a new directory")
    prob, th = _synthetic()
    rows = table(prob.data(), th, prob.L.gene_names, "synthetic_5_to_500_rows")
    fprob, _ = fio.load_problem(args.fixture)
    params = fio.load_params(args.fixture, fprob.L)
    for name in ("theta3_truth_vpos", "theta4_perturbed"):
        rows += table(fprob.data(), params[name], fprob.L.gene_names, f"fixture_{name}")
    res = dict(tolerance_abs_ell=TOL, tau_max=TAU_MAX, n=len(rows),
               n_pass=sum(r["pass"] for r in rows),
               max_abs_gh20_minus_quad=max(abs(r["gh20_minus_quad"]) for r in rows),
               max_abs_reference_minus_quad=max(abs(r["reference_minus_quad"]) for r in rows),
               max_abs_gh20_minus_gh40=max(abs(r["gh20_minus_gh40"]) for r in rows), rows=rows)
    (out / "accuracy_v18.json").write_text(json.dumps(res, indent=1))
    print({k: v for k, v in res.items() if k != "rows"}, flush=True)


if __name__ == "__main__":
    main()
