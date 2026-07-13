#!/usr/bin/env python3
"""
76b_snatac_permnull.py — Phase-6 WS-4: permutation null for the Decima argmax-vs-snATAC
corroboration (red-team critique L3: "19/29 has no interpretable null").

The Decima cell-type argmax is corroborated when a variant's top compartment is among
the compartments whose snATAC peaks overlap the variant (atac_accessible_celltypes).
Observed statistic = #(argmax in accessible set) among variants with >=1 accessible peak.

NULL (preserves peak-overlap MULTIPLICITY + compartment availability): for each testable
variant, draw a random argmax compartment and count whether it lands in that variant's
accessible set. The chance for variant i is |accessible_i| / n_compartments, so the null
naturally preserves each variant's accessible-set SIZE. Two nulls reported:
  (1) uniform over the 5 scored compartments;
  (2) drawn from the EMPIRICAL argmax frequency (preserves Decima's compartment bias, e.g.
      hepatocyte picked more often) — the more conservative, appropriate null.

This is a semi-independent SANITY CHECK, not disjoint validation (same Human-Multiome cohort;
peak-presence vs disease-DA). It does NOT rehabilitate the aggregate hep/non-parenchymal %.
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SF = os.path.join(ROOT, "GWAS/finemapping/results/seqfunc")
TSV = os.path.join(SF, "decima_celltype.tsv")
SUMMARY = os.path.join(SF, "decima_celltype_summary.json")
COMPARTMENTS = ["hepatocyte", "LSEC_endothelial", "kupffer_macrophage",
                "lymphoid_T_NK_B", "stellate_fibroblast"]
NPERM = 100_000
SEED = 42


def parse_set(x):
    if not isinstance(x, str) or not x.strip() or x.strip() in ("[]", "set()", "nan"):
        return set()
    return set(t.strip().strip("'\"") for t in
               x.replace("[", "").replace("]", "").replace("{", "").replace("}", "").split(",")
               if t.strip())


def main():
    df = pd.read_csv(TSV, sep="\t")
    # accessible set restricted to the 5 SCORED argmax compartments (argmax can never be
    # cholangiocyte — no eff_ column — so a cholangiocyte-only peak can never match).
    df["acc"] = df["atac_accessible_celltypes"].map(parse_set).map(
        lambda s: {c for c in s if c in COMPARTMENTS})
    # testable = any snATAC peak overlap at all (preserves the reported denominator of 29)
    test = df[df["atac_any_peak"] == True].copy()      # noqa: E712
    n_test = len(test)
    obs = int((test["argmax_matches_atac"] == True).sum())   # noqa: E712  (== reported 19)
    acc_sizes = np.array([len(s) for s in test["acc"]])

    rng = np.random.default_rng(SEED)
    # empirical argmax frequency over the 5 compartments (from ALL scored variants)
    freq = df["argmax_celltype"].value_counts()
    p_emp = np.array([freq.get(c, 0) for c in COMPARTMENTS], float)
    p_emp = p_emp / p_emp.sum()

    def perm_p(prob):
        acc_masks = np.array([[1 if c in s else 0 for c in COMPARTMENTS] for s in test["acc"]])
        hits = np.empty(NPERM, int)
        for b in range(NPERM):
            draw = rng.choice(len(COMPARTMENTS), size=n_test, p=prob)
            hits[b] = int(acc_masks[np.arange(n_test), draw].sum())
        p = (1 + int(np.sum(hits >= obs))) / (1 + NPERM)
        return p, float(hits.mean()), float(np.percentile(hits, 97.5))

    p_uni, mean_uni, hi_uni = perm_p(np.full(len(COMPARTMENTS), 1 / len(COMPARTMENTS)))
    p_empd, mean_emp, hi_emp = perm_p(p_emp)

    expected_uniform = float((acc_sizes / len(COMPARTMENTS)).sum())

    res = dict(
        n_testable=n_test, observed_matches=obs,
        observed_rate=round(obs / n_test, 4),
        expected_matches_uniform_analytic=round(expected_uniform, 3),
        null_uniform=dict(mean_matches=round(mean_uni, 3), p95=hi_uni,
                          perm_p=p_uni, nperm=NPERM),
        null_empirical_argmax_freq=dict(mean_matches=round(mean_emp, 3), p95=hi_emp,
                                        perm_p=p_empd, nperm=NPERM,
                                        note="preserves Decima compartment bias; conservative"),
        interpretation=("Semi-independent sanity check (same Human-Multiome cohort; peak-presence "
                        "vs disease-DA, NOT disjoint validation). Does NOT rehabilitate the "
                        "aggregate hepatocyte/non-parenchymal proportion."))
    print(json.dumps(res, indent=2))

    # merge into the summary json (additive; do not disturb other keys)
    with open(SUMMARY) as fh:
        summ = json.load(fh)
    prev = summ.get("snATAC_corroboration")
    summ["snATAC_corroboration"] = {"legacy_value": prev, "permutation_test": res}
    with open(SUMMARY, "w") as fh:
        json.dump(summ, fh, indent=2)
    print(f"\nupdated {SUMMARY}")


if __name__ == "__main__":
    main()
