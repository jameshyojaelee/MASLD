#!/usr/bin/env python3
"""H3 step 5: independent re-derivation of selected cohort estimates (numpy, not the R code path).
Usage: 05_verify.py CALLS_DIR SCORES_DIR TESTS_DIR"""
import csv, sys
import numpy as np
import pandas as pd

calls, scores, tests = sys.argv[1:4]
ph = pd.read_csv(scores + "/person_phenotypes.tsv", sep="\t")
gt = pd.read_csv(calls + "/person_genotypes.tsv", sep="\t", na_values=["NA"], keep_default_na=False)
d0 = ph.merge(gt.drop(columns=["cohort"]), on="person_id")
d0["sex2"] = d0["sex"].where(d0["sex"].isin(["F", "M"]))


def merged_stage(stage):
    s = stage.astype(float).copy()
    while True:
        vc = s.value_counts()
        if len(vc) < 2 or (vc >= 3).all():
            return s
        small = vc.idxmin() if (vc == vc.min()).sum() == 1 else min(vc[vc == vc.min()].index)
        others = sorted(set(vc.index) - {small})
        dist = [abs(o - small) for o in others]
        s[s == small] = others[int(np.argmin(dist))]


def hc3(X, y, j):
    n, p = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    e = y - X @ beta
    w = (e / (1 - np.diag(H))) ** 2
    XtXi = np.linalg.inv(X.T @ X)
    V = XtXi @ (X.T * w) @ X @ XtXi
    return beta[j], np.sqrt(V[j, j])


def cohort_fit(co, dose, yname, extra=None):
    d = d0[(d0.cohort == co) & d0[dose].notna() & d0.stage.notna() & d0.sex2.notna()].copy()
    d["st"] = merged_stage(d.stage)
    dum = pd.get_dummies(d["st"], drop_first=True).astype(float).values
    X = np.column_stack([np.ones(len(d)), d[dose].values, dum, (d.sex2 == "F").astype(float).values])
    return hc3(X, d[yname].values, 1), len(d)


def dose_col(rs):
    return {"rs72613567": lambda c: c.map({"CARRIER": 1.0, "NC": 0.0})}.get(rs, lambda c: c.map({"RR": 0.0, "RA": 1.0, "AA": 2.0}))


for rs in ["rs72613567", "rs738409", "rs58542926", "rs2642438", "rs1260326"]:
    d0["dose_" + rs] = dose_col(rs)(d0[rs + "_call"])
d0["dose_pi"] = d0["pi_status"].map({"MZ": 1.0, "MM": 0.0})

est = pd.read_csv(tests + "/h3_cohort_estimates.tsv", sep="\t")
rows = []
for rs, co, pg in [("rs738409", "GSE213621", "ecm"), ("rs738409", "GSE162694", "lipid"), ("rs58542926", "GSE135251", "ecm"),
                   ("rs2642438", "GSE130970", "lipid"), ("rs1260326", "GSE193066", "ecm"), ("rs72613567", "GSE162694", "lipid")]:
    (b, s), n = cohort_fit(co, "dose_" + rs, "score_" + pg)
    r = est[(est.rsid == rs) & (est.cohort == co) & (est.program == pg)].iloc[0]
    rows.append((rs, co, pg, n, b, r.beta, s, r.se))
h31 = pd.read_csv(tests + "/h31_cohort_estimates.tsv", sep="\t")
for co in ["GSE213621", "GSE130970"]:
    (b, s), n = cohort_fit(co, "dose_pi", "score_upr")
    r = h31[h31.cohort == co].iloc[0]
    rows.append(("Pi*MZ", co, "upr", n, b, r.beta, s, r.se))
print("%-11s %-10s %-5s %4s %10s %10s %10s %10s" % ("variant", "cohort", "prog", "n", "beta_py", "beta_R", "se_py", "se_R"))
worst = 0.0
for r in rows:
    print("%-11s %-10s %-5s %4d %10.5f %10.5f %10.5f %10.5f" % r)
    worst = max(worst, abs(r[4] - r[5]), abs(r[6] - r[7]))
print("largest absolute difference:", worst)
assert worst < 1e-8
