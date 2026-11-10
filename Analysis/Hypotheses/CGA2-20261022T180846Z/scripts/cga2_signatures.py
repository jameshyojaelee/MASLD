#!/usr/bin/env python
"""CGA2: published signatures vs matched random sets; within-stage continuum tracking (E1) and
stage information beyond the continuum (E2). Follows PRESPEC_CGA2.md."""
import argparse
import json
import os

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PANELS = {"govaere_2020": ["GSE162694", "GSE213621"], "hoang_2019": ["GSE162694", "GSE213621"],
          "moylan_2014": ["GSE162694", "GSE213621"], "arendt_2015": ["GSE162694", "GSE213621"],
          "suppli_2019": ["GSE162694", "GSE213621"], "pantano_2021": ["GSE213621"]}
SEED, NRAND = 20261022, 1000


def read_panel(name, sym2g):
    p = pd.read_csv(f"{ROOT}/data/published_gene_panels/{name}_panel.tsv", sep="\t", comment="#")
    p = p[p.direction.isin(["up", "down"])]
    p["g"] = p.gene_symbol.map(sym2g)
    return p


def score(z, up, dn):
    s = z.loc[up].mean(0) if len(up) else 0
    if len(dn):
        s = s - z.loc[dn].mean(0)
    return (s - s.mean()) / s.std(ddof=1)


def e1(df, s):
    d = df.assign(s=s.reindex(df.sample_id).values)
    f = smf.ols("continuum_z ~ s + C(stage) + sex", d).fit()
    return f.params["s"], f.pvalues["s"]


def e2(df, s, ll0):
    d = df.assign(s=s.reindex(df.sample_id).values)
    m = OrderedModel(d.stage, d[["continuum_z", "sex", "s"]], distr="logit").fit(method="bfgs", disp=False)
    return 2 * (m.llf - ll0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--units", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--logcpm", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=False)
    rng = np.random.default_rng(SEED)
    genes = pd.read_csv(a.genes, sep="\t")
    sig = set(genes.g[genes.signature_gene.fillna(False).astype(bool)])
    sym2g = genes.dropna(subset=["gene_name"]).drop_duplicates("gene_name").set_index("gene_name").g
    units = pd.read_csv(a.units, sep="\t")
    rows, nulls = [], []
    for ds in ["GSE162694", "GSE213621"]:
        x = pd.read_csv(f"{a.logcpm}/{ds}_logcpm.tsv.gz", sep="\t").set_index("g")
        u = units[(units.dataset == ds) & units.fibrosis_stage.notna()]
        u = u[~((ds == "GSE213621") & (u.fibrosis_stage == 0))].copy()
        u["stage"] = u.fibrosis_stage.astype(int)
        u["sex"] = (u.inferred_sex == "M").astype(int)
        u["continuum_z"] = (u.continuum_z - u.continuum_z.mean()) / u.continuum_z.std(ddof=1)
        x = x[u.sample_id]
        z = x.sub(x.mean(1), axis=0).div(x.std(1, ddof=1), axis=0).dropna()
        pool = z.index.difference(list(sig))
        dec = pd.qcut(x.loc[pool].mean(1), 10, labels=False)
        bins = {k: dec.index[dec == k].to_numpy() for k in range(10)}
        ll0 = OrderedModel(u.stage, u[["continuum_z", "sex"]], distr="logit").fit(method="bfgs", disp=False).llf
        for name, cohorts in PANELS.items():
            if ds not in cohorts:
                continue
            p = read_panel(name, sym2g)
            n_listed, n_sig = len(p), int(p.g.isin(sig).sum())
            p = p[p.g.isin(pool)]
            up, dn = list(p.g[p.direction == "up"]), list(p.g[p.direction == "down"])
            s = score(z, up, dn)
            b1, p1 = e1(u, s)
            chi = e2(u, s, ll0)
            rb, rc = np.empty(NRAND), np.empty(NRAND)
            for i in range(NRAND):
                def pick(gs):
                    out = []
                    for k, n in dec.loc[gs].value_counts().items():
                        out.extend(rng.choice(bins[k], n, replace=False))
                    return out
                rs = score(z, pick(up), pick(dn))
                rb[i] = e1(u, rs)[0]
                rc[i] = e2(u, rs, ll0)
            nulls.append(pd.DataFrame({"dataset": ds, "signature": name, "rand_e1": rb, "rand_e2_chi2": rc}))
            rows.append(dict(dataset=ds, signature=name, n_units=len(u), n_listed=n_listed,
                             n_signature_genes_removed=n_sig, n_unmapped_or_unexpressed=n_listed - n_sig - len(p), n_up=len(up), n_down=len(dn),
                             e1_beta=b1, e1_p=p1, e1_emp_p=(1 + (rb >= b1).sum()) / (NRAND + 1),
                             e1_rand_median=np.median(rb), e2_chi2=chi, e2_lr_p=stats.chi2.sf(chi, 1),
                             e2_emp_p=(1 + (rc >= chi).sum()) / (NRAND + 1), e2_rand_median=np.median(rc),
                             e2_rand_frac_lr_p_lt_0_05=float((stats.chi2.sf(rc, 1) < 0.05).mean())))
            print(rows[-1], flush=True)
    t = pd.DataFrame(rows)
    t["e1_bh"] = multipletests(t.e1_emp_p, method="fdr_bh")[1]
    t["e2_bh"] = multipletests(t.e2_emp_p, method="fdr_bh")[1]
    adds = {}
    for name, cohorts in PANELS.items():
        r = t[t.signature == name]
        adds[name] = bool(len(r) == len(cohorts) and (r.e2_bh < 0.05).all())
    t.to_csv(f"{a.out}/signature_results.tsv", sep="\t", index=False)
    pd.concat(nulls).to_csv(f"{a.out}/random_set_null.tsv.gz", sep="\t", index=False)
    json.dump({"adds_stage_info_beyond_continuum": adds, "seed": SEED, "n_random": NRAND},
              open(f"{a.out}/results.json", "w"), indent=1)


if __name__ == "__main__":
    main()
