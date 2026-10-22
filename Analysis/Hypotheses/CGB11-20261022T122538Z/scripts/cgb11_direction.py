#!/usr/bin/env python
"""CGB11: colocalized risk-haplotype expression direction vs disease direction.
Follows PRESPEC_CGB11.md (sha256 in PRESPEC_CGB11.md.sha256)."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
DIR = "Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z/coloc_direction/full-20261015T225130Z/direction_by_trait.tsv"
SEED = 20261022
NPERM = 10000


def gene_signs(d):
    d = d[d.release == "adopted_r2_2026-10-14"].copy()
    n = len(d)
    vc = d.direction_sign.value_counts(dropna=False).to_dict()
    d["direction_sign"] = pd.to_numeric(d.direction_sign, errors="coerce")
    dd = d[d.direction_sign.isin([1, -1])].copy()
    dd["g"] = dd.ensembl.astype(str).str.split(".").str[0]
    dd["chrom"] = dd.hit1.astype(str).str.split(":").str[0]
    dd["pos"] = pd.to_numeric(dd.hit1.astype(str).str.split(":").str[1], errors="coerce")
    agg = dd.groupby("g").agg(gene=("gene", "first"), n_rows=("direction_sign", "size"),
                              s_min=("direction_sign", "min"), s_max=("direction_sign", "max"),
                              pp4=("pair_pp4", "max"), chrom=("chrom", "first"), pos=("pos", "median"),
                              traits=("trait", lambda x: ";".join(sorted(set(x)))))
    agg["mixed"] = agg.s_min != agg.s_max
    agg["gen_sign"] = np.where(agg.mixed, 0, agg.s_max)
    return n, vc, agg.reset_index()


def assign_loci(df):
    df = df.sort_values(["chrom", "pos"]).copy()
    loc, cur, last_c, last_p = [], -1, None, None
    for c, p in zip(df.chrom, df.pos):
        if c != last_c or p - last_p > 1_000_000:
            cur += 1
        loc.append(cur)
        last_c, last_p = c, p
    df["locus"] = loc
    return df


def prune(df):
    return (df.sort_values(["locus", "pp4", "g"], ascending=[True, False, True])
              .drop_duplicates("locus"))


def mh_test(df, dcol):
    df = df.copy()
    df["q"] = pd.qcut(df.AveExpr.rank(method="first"), 5, labels=False)
    tabs = []
    for _, s in df.groupby("q"):
        a = int(((s.gen_sign > 0) & (s[dcol] > 0)).sum())
        b = int(((s.gen_sign > 0) & (s[dcol] < 0)).sum())
        c = int(((s.gen_sign < 0) & (s[dcol] > 0)).sum())
        d = int(((s.gen_sign < 0) & (s[dcol] < 0)).sum())
        if (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
            tabs.append(np.array([[a, b], [c, d]]) + 0.0)
    if len(tabs) == 0:
        return dict(n=len(df), or_mh=np.nan, lo=np.nan, hi=np.nan, p=np.nan)
    st = StratifiedTable([t + 0.5 if (t == 0).any() else t for t in tabs])
    lo, hi = st.oddsratio_pooled_confint()
    conc = float((np.sign(df.gen_sign) == np.sign(df[dcol])).mean())
    return dict(n=int(len(df)), n_strata=len(tabs), concordant_fraction=conc,
                or_mh=float(st.oddsratio_pooled), lo=float(lo), hi=float(hi),
                p=float(st.test_null_odds().pvalue))


def perm_corr(df, dcol, rng):
    x = df.gen_sign.values.astype(float)
    y = df[dcol].values.astype(float)
    obs = stats.pearsonr(x, y)[0]
    q = pd.qcut(df.AveExpr.rank(method="first"), 5, labels=False).values
    idx = [np.where(q == k)[0] for k in np.unique(q)]
    null = np.empty(NPERM)
    for i in range(NPERM):
        yp = y.copy()
        for ix in idx:
            yp[ix] = y[rng.permutation(ix)]
        null[i] = np.corrcoef(x, yp)[0, 1]
    p = (1 + np.sum(np.abs(null) >= abs(obs))) / (NPERM + 1)
    return dict(n=int(len(df)), r=float(obs), perm_p=float(p))


def mde(df, dcol, rng, nsim=2000):
    """Smallest OR with 80% power at alpha 0.025 (Holm first step), observed margins."""
    p_up_g = (df.gen_sign > 0).mean()
    p_up_d = (df[dcol] > 0).mean()
    n = len(df)
    for orv in [1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0]:
        # solve cell probabilities for target OR with fixed margins (Plackett)
        a, b = p_up_g, p_up_d
        if orv == 1:
            p11 = a * b
        else:
            s = 1 + (orv - 1) * (a + b)
            p11 = (s - np.sqrt(s * s - 4 * orv * (orv - 1) * a * b)) / (2 * (orv - 1))
        probs = np.array([p11, a - p11, b - p11, 1 - a - b + p11])
        hits = 0
        for _ in range(nsim):
            c = rng.multinomial(n, probs)
            t = np.array([[c[0], c[1]], [c[2], c[3]]]) + 0.5
            lor = np.log(t[0, 0] * t[1, 1] / (t[0, 1] * t[1, 0]))
            se = np.sqrt((1 / t).sum())
            hits += abs(lor / se) > stats.norm.ppf(1 - 0.0125)
        if hits / nsim >= 0.8:
            return orv
    return np.inf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.chdir(ROOT)
    os.makedirs(a.out, exist_ok=False)
    rng = np.random.default_rng(SEED)
    raw = pd.read_csv(DIR, sep="\t", low_memory=False)
    n, vc, gs = gene_signs(raw)
    cur = raw[raw.release == "adopted_r2_2026-10-14"]
    counts = (len(cur), cur.ensembl.nunique(), int((cur.direction_sign == 1).sum()),
              int((cur.direction_sign == -1).sum()))
    res = {"input_rows_genes_up_down": counts}
    if counts[:2] != (589, 436):
        print("S2 stop: input counts differ", counts)
        sys.exit(2)
    genes = pd.read_csv(a.genes, sep="\t")
    gs = gs.merge(genes[["g", "AveExpr", "disease_t", "disease_padj", "cont_z", "cont_q", "deg"]],
                  on="g", how="left")
    res["genes_directional_any"] = int(len(gs))
    res["genes_mixed_sign_excluded"] = int(gs.mixed.sum())
    gs = gs[~gs.mixed & gs.AveExpr.notna()].copy()
    gs = assign_loci(gs)
    gs.to_csv(f"{a.out}/gene_directions.tsv", sep="\t", index=False)
    pr = prune(gs)
    res["genes_usable"] = int(len(gs))
    res["loci"] = int(pr.locus.nunique())

    out = {}
    for name, dcol, qcol in [("P1_disease", "disease_t", "disease_padj"), ("P2_continuum", "cont_z", "cont_q")]:
        sub = pr[pr[qcol] < 0.05]
        if len(sub) < 20:
            out[name] = {"status": "untestable_S1", "n": int(len(sub))}
            continue
        r = mh_test(sub, dcol)
        r["mde_or"] = mde(sub, dcol, rng)
        r["unpruned"] = mh_test(gs[gs[qcol] < 0.05], dcol)
        r["unthresholded_perm"] = perm_corr(pr[pr[dcol].notna()], dcol, rng)
        out[name] = r
    ps = [out[k].get("p", np.nan) for k in out]
    ok = [i for i, p in enumerate(ps) if np.isfinite(p)]
    if ok:
        adj = multipletests([ps[i] for i in ok], method="holm")[1]
        for i, q in zip(ok, adj):
            out[list(out)[i]]["holm_p"] = float(q)
    res["tests"] = out
    conc26 = gs[gs.deg.fillna(False).astype(bool)][["gene", "gen_sign", "disease_t", "cont_z", "traits"]]
    conc26.to_csv(f"{a.out}/deg_and_colocalized_directions.tsv", sep="\t", index=False)
    res["deg_and_colocalized_directional_n"] = int(len(conc26))
    with open(f"{a.out}/results.json", "w") as f:
        json.dump(res, f, indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
