#!/usr/bin/env python
"""CGB2: risk-allele share of transcript in heterozygotes vs stage (P1) and within-stage continuum (P2).
Follows PRESPEC_CGB2.md."""
import argparse
import json
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

H3 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H3-20261009T180331Z/results"
COUNTS = f"{H3}/geno-20261009T180507Z/person_variant_counts.tsv"
GENO = f"{H3}/calls-20261009T181453Z/person_genotypes.tsv"
VARIANTS = {"rs738409": "PNPLA3", "rs58542926": "TM6SF2", "rs1260326": "GCKR", "rs2642438": "MTARC1"}
P1_COHORTS = ["GSE130970", "GSE135251", "GSE162694", "GSE213621"]
P2_COHORTS = ["GSE162694", "GSE213621"]


def load(units_path):
    u = pd.read_csv(units_path, sep="\t")
    g = pd.read_csv(GENO, sep="\t")
    lib = g[["person_id", "libraries"]].assign(lib=g.libraries.str.split(";|,")).explode("lib")
    lib["lib"] = lib.lib.str.strip()
    pu = lib.merge(u, left_on="lib", right_on="sample_id")
    pu = pu.drop_duplicates("person_id")
    c = pd.read_csv(COUNTS, sep="\t")
    c = c[c.rsid.isin(VARIANTS) & (c.depth >= 20)].copy()
    c["f"] = c.alt_reads / c.depth
    c = c[(c.f >= 0.2) & (c.f <= 0.8)]
    d = c.merge(pu, on="person_id")
    d = d[d.dataset.isin(P1_COHORTS) & d.fibrosis_stage.notna()]
    d = d[~((d.dataset == "GSE213621") & (d.fibrosis_stage == 0))].copy()
    d["sex"] = (d.inferred_sex == "M").astype(int)
    d["stage"] = d.fibrosis_stage.astype(int)
    for ds, idx in d.groupby("dataset").groups.items():
        v = u[(u.dataset == ds)].continuum_z
        d.loc[idx, "cz"] = (d.loc[idx, "continuum_z"] - v.mean()) / v.std(ddof=1)
    return d


def qb(formula, df, term):
    # quasibinomial on the alt fraction with depth as variance weight; a two-column count endog
    # makes statsmodels' Pearson scale ignore depth and understates the SE
    f = smf.glm(formula, df, family=sm.families.Binomial(), var_weights=df.depth).fit(scale="X2")
    return dict(beta=f.params[term], se=f.bse[term], p=f.pvalues[term], dispersion=f.scale, n=len(df))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    d = load(a.units)
    d.to_csv(f"{a.out}/heterozygote_counts.tsv", sep="\t", index=False)
    rows = []
    for rs, gene in VARIANTS.items():
        s = d[d.rsid == rs]
        per_stage = s.groupby(["dataset", "stage"]).size()
        row = dict(variant=rs, gene=gene, n_het=len(s), mean_alt_fraction=s.f.mean(),
                   min_hets_per_stage_group=int(per_stage.min()))
        row["p1_testable"] = bool(per_stage.min() >= 10)
        p2_groups = s[s.dataset.isin(P2_COHORTS)].groupby(["dataset", "stage"]).size()
        row["p2_testable"] = bool(p2_groups.min() >= 10)
        if per_stage.min() < 10:
            row["p1_note"] = "some stage group has <10 heterozygotes; groups listed in heterozygote_counts.tsv"
        p1 = qb("f ~ stage + C(dataset) + sex", s, "stage")
        row.update({f"p1_{k}": v for k, v in p1.items()})
        signs = {ds: np.sign(qb("f ~ stage + sex", t, "stage")["beta"]) for ds, t in s.groupby("dataset")
                 if t.stage.nunique() > 1 and len(t) >= 10}
        row["p1_cohort_signs"] = signs
        row["p1_same_sign_cohorts"] = int(sum(v == np.sign(p1["beta"]) for v in signs.values()))
        hep = qb("f ~ stage + C(dataset) + sex + bp_Hepatocytes", s, "stage")
        row["p1_hep_adjusted_beta"], row["p1_hep_adjusted_p"] = hep["beta"], hep["p"]
        s2 = s[s.dataset.isin(P2_COHORTS)]
        p2 = qb("f ~ cz + C(stage) + C(dataset) + sex", s2, "cz")
        row.update({f"p2_{k}": v for k, v in p2.items()})
        rows.append(row)
    t = pd.DataFrame(rows)
    t["p1_holm"] = multipletests(t.p1_p, method="holm")[1]
    t["p2_holm"] = multipletests(t.p2_p, method="holm")[1]
    t["stage_dependent"] = t.p1_testable & (t.p1_holm < 0.05) & (t.p1_same_sign_cohorts >= 3)
    t.to_csv(f"{a.out}/allelic_results.tsv", sep="\t", index=False)
    json.dump(t.to_dict(orient="records"), open(f"{a.out}/results.json", "w"), indent=1, default=str)
    print(t.drop(columns=["p1_cohort_signs"]).to_string())


if __name__ == "__main__":
    main()
