#!/usr/bin/env python3
"""Model A T6: predictive allelic model, development cross-fitting (PRESPEC_A section 9).

Pairs: genes with main-effect q < 0.05 (T5a) and >= 20 development participants
with composition (guard-passing). Target: each held-out participant's oriented
alt-haplotype read count a out of n at the gene's chosen het tag.

Models (beta-binomial, overdispersion fitted on the training fold):
  GTEx  p = 2^aFC / (1 + 2^aFC) from GTEx v8 liver log2 aFC (no fitted mean).
  L0    one alt-haplotype fraction per gene.
  L1    p_d = sum_k w_dk * pi_k: w_dk = theta_dk * e_gk / sum_j theta_dj * e_gj with
        theta_dk the participant's InstaPrism mRNA fraction and e_gk the gene's CPM in
        the atlas pseudobulk of lineage k. Components: hepatocytes plus the gene's two
        other lineages with the largest median w; remaining weight joins hepatocytes.
        Fitted only for genes with >= 20% median non-hepatocyte weight; other genes
        keep L0. logit(pi_k) has a ridge penalty (lambda = 1) toward the L0 logit.
  L2    fitted only if G3 passes (not in this script).

Cross-fitting: 5 participant folds by sha256("MASLD-A-T6" + individual) mod 5,
first-degree and same-individual relatives in one fold. Metric: held-out
deviance D = -2 * sum log BB(a | n, p, rho) with the binomial coefficient.
L1 beats L0 if (D_L0 - D_L1) / D_L0 >= 0.005 and the gain exceeds the 95th
percentile of 200 refits with composition rows shuffled across participants
within cohort (seed 20260923).
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t5a_allelic_table import bb_fit  # noqa: E402

SEED, N_SHUF, K_FOLD, MIN_PAIR, MIN_NONHEP, LAMBDA, MARGIN = 20260923, 200, 5, 20, 0.20, 1.0, 0.005


def bb_loglik(k, n, p, rho):
    s = (1 - rho) / rho
    a_, b_ = p * s, (1 - p) * s
    return (special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1)
            + special.betaln(k + a_, n - k + b_) - special.betaln(a_, b_))


def fit_rho_given_p(k, n, p):
    f = lambda g: -np.sum(bb_loglik(k, n, p, special.expit(g[0])))
    return float(special.expit(optimize.minimize(f, [-3.0], method="L-BFGS-B", bounds=[(-12, 6)]).x[0]))


def fit_l0(k, n):
    par, _ = bb_fit(k, n)
    return float(special.expit(par[0])), float(special.expit(par[-1]))


def fit_l1(k, n, W, p0):
    K = W.shape[1]
    t0 = special.logit(np.clip(p0, 1e-4, 1 - 1e-4))

    def f(x):
        pi, rho = special.expit(x[:K]), special.expit(x[K])
        p = np.clip(W @ pi, 1e-6, 1 - 1e-6)
        return -np.sum(bb_loglik(k, n, p, rho)) + LAMBDA * np.sum((x[:K] - t0) ** 2)
    r = optimize.minimize(f, np.r_[np.full(K, t0), -3.0], method="L-BFGS-B", bounds=[(-10, 10)] * K + [(-12, 6)])
    return special.expit(r.x[:K]), float(special.expit(r.x[K])), r


def lineage_weights(theta, e_g, comps):
    """theta: participants x lineages (mRNA fraction); e_g: lineage CPM of the gene; comps: component lineage names."""
    raw = theta.to_numpy() * e_g[theta.columns].to_numpy()[None, :]
    raw = raw / raw.sum(1, keepdims=True)
    W = pd.DataFrame(raw, columns=theta.columns)
    out = W[comps].copy()
    out[comps[0]] += W.drop(columns=comps).sum(1)
    return out.to_numpy()


def folds(ids, kin, xwalk):
    h = {i: hashlib.sha256(("MASLD-A-T6" + i).encode()).hexdigest() for i in ids}
    run_ind = xwalk.set_index("run")["individual_id"]
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    rel = kin[kin["relation"].isin(["first_degree", "same_individual"])]
    for a_, b_ in zip(rel["#IID1"].map(run_ind), rel["IID2"].map(run_ind)):
        if a_ in parent and b_ in parent:
            ra, rb = find(a_), find(b_)
            if ra != rb:
                parent[max(ra, rb, key=lambda z: h[z])] = min(ra, rb, key=lambda z: h[z])
    return {i: int(h[find(i)], 16) % K_FOLD for i in ids}


def crossfit_gene(g, W, fold, afc):
    """Held-out log-likelihood per participant for GTEx, L0 and (if W is not None) L1."""
    k, n = g["a_o"].to_numpy(float), g["n"].to_numpy(float)
    ll = {"gtex": np.full(len(g), np.nan), "l0": np.full(len(g), np.nan), "l1": np.full(len(g), np.nan)}
    for f in range(K_FOLD):
        te, tr = fold == f, fold != f
        if te.sum() == 0 or tr.sum() < 10:
            continue
        p0, rho0 = fit_l0(k[tr], n[tr])
        ll["l0"][te] = bb_loglik(k[te], n[te], p0, rho0)
        if not np.isnan(afc):
            pg = 2 ** afc / (1 + 2 ** afc)
            ll["gtex"][te] = bb_loglik(k[te], n[te], pg, fit_rho_given_p(k[tr], n[tr], pg))
        if W is not None:
            pi, rho1, _ = fit_l1(k[tr], n[tr], W[tr], p0)
            ll["l1"][te] = bb_loglik(k[te], n[te], np.clip(W[te] @ pi, 1e-6, 1 - 1e-6), rho1)
    return ll


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--participant-gene", required=True)
    ap.add_argument("--t5a-dir", required=True)
    ap.add_argument("--identity-dir", required=True)
    ap.add_argument("--pseudobulk", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    me = pd.read_csv(Path(a.t5a_dir) / "main_allelic_effects.tsv", sep="\t")
    pg = pd.read_csv(a.participant_gene, sep="\t")
    pb = pd.read_csv(a.pseudobulk, sep="\t", index_col=0)
    cpm = pb / pb.sum(0) * 1e6
    lineages = list(cpm.columns)
    pg = pg[pg["passes_guard"]].dropna(subset=lineages)
    pairs = me[me["q"] < 0.05]
    counts = pg.groupby("gene_id").size()
    genes = sorted(g for g in pairs["gene_id"] if counts.get(g, 0) >= MIN_PAIR)
    afc = me.set_index("gene_id")["gtex_log2_aFC"]
    xwalk = pd.read_csv(Path(a.identity_dir) / "sample_to_individual.tsv", sep="\t")
    kin = pd.read_csv(Path(a.identity_dir) / "kinship_related_pairs.tsv", sep="\t")
    fold_of = folds(sorted(pg["individual_id"].unique()), kin, xwalk)

    per_gene, W_of, comps_of = [], {}, {}
    ll_rows = []
    for gid in genes:
        g = pg[pg["gene_id"] == gid].reset_index(drop=True)
        name = g["gene_name"].iloc[0]
        W = None
        nonhep = np.nan
        if name in cpm.index and cpm.loc[name].sum() > 0:
            e_g = cpm.loc[name]
            full = lineage_weights(g[lineages], e_g, lineages)
            med = pd.Series(np.median(full, 0), index=lineages)
            nonhep = float(1 - med["Hepatocytes"])
            if nonhep >= MIN_NONHEP:
                comps = ["Hepatocytes"] + list(med.drop("Hepatocytes").sort_values(ascending=False).index[:2])
                W = lineage_weights(g[lineages], e_g, comps)
                W_of[gid], comps_of[gid] = W, comps
        fold = g["individual_id"].map(fold_of).to_numpy()
        ll = crossfit_gene(g, W, fold, afc.get(gid, np.nan))
        ll_rows.append(pd.DataFrame({"gene_id": gid, "cohort": g["cohort"], **{f"ll_{m}": v for m, v in ll.items()}}))
        per_gene.append({"gene_id": gid, "gene_name": name, "n": len(g), "median_nonhep_weight": nonhep,
                         "l1_fitted": W is not None, "l1_components": ",".join(comps_of.get(gid, [])),
                         **{f"D_{m}": float(-2 * np.nansum(v)) if not np.all(np.isnan(v)) else np.nan for m, v in ll.items()}})
    pgdf = pd.DataFrame(per_gene)
    llp = pd.concat(ll_rows, ignore_index=True)
    llp.to_csv(out / "heldout_loglik_per_participant.tsv.gz", sep="\t", index=False, compression="gzip")

    def rel_gain(sub, base, new):
        s = sub.dropna(subset=[f"ll_{base}", f"ll_{new}"])
        d0, d1 = -2 * s[f"ll_{base}"].sum(), -2 * s[f"ll_{new}"].sum()
        return float((d0 - d1) / d0) if d0 else np.nan, int(s["gene_id"].nunique())

    l1g = [g for g in genes if g in W_of]
    sub = llp[llp["gene_id"].isin(l1g)]
    gain_l1, n_l1 = rel_gain(sub, "l0", "l1")
    gain_l0_vs_gtex, n_g = rel_gain(llp, "gtex", "l0")

    # shuffle null for L1: composition rows permuted across participants within cohort
    rng = np.random.default_rng(SEED)
    null = []
    for b in range(N_SHUF):
        rows = []
        for gid in l1g:
            g = pg[pg["gene_id"] == gid].reset_index(drop=True)
            perm = np.arange(len(g))
            for c, ii in g.groupby("cohort").indices.items():
                perm[ii] = ii[rng.permutation(len(ii))]
            e_g = cpm.loc[g["gene_name"].iloc[0]]
            Wb = lineage_weights(g[lineages].iloc[perm].reset_index(drop=True), e_g, comps_of[gid])
            ll = crossfit_gene(g, Wb, g["individual_id"].map(fold_of).to_numpy(), np.nan)
            rows.append(pd.DataFrame({"gene_id": gid, "ll_l0": ll["l0"], "ll_l1": ll["l1"]}))
        null.append(rel_gain(pd.concat(rows), "l0", "l1")[0] if rows else np.nan)
    null = np.array(null)
    summ = {"pairs": len(genes), "l1_genes": n_l1,
            "gain_L0_over_GTEx_relative_deviance": gain_l0_vs_gtex, "genes_with_gtex_aFC": n_g,
            "gain_L1_over_L0_relative_deviance": gain_l1,
            "shuffle_null_95th": float(np.nanpercentile(null, 95)) if len(null) else None,
            "shuffle_null_median": float(np.nanmedian(null)) if len(null) else None,
            "l1_beats_l0": bool(len(null) and gain_l1 >= MARGIN and gain_l1 > np.nanpercentile(null, 95)),
            "gain_L1_over_L0_by_cohort": {c: rel_gain(d, "l0", "l1")[0] for c, d in sub.groupby("cohort")},
            "fold_sizes": pd.Series(fold_of).value_counts().sort_index().to_dict(),
            "settings": {"seed": SEED, "n_shuffle": N_SHUF, "folds": K_FOLD, "min_participants": MIN_PAIR,
                         "min_nonhep_weight": MIN_NONHEP, "ridge_lambda": LAMBDA, "margin": MARGIN}}
    pgdf.to_csv(out / "t6_per_gene.tsv", sep="\t", index=False)
    np.savetxt(out / "l1_shuffle_null_gains.txt", null)
    (out / "t6_summary.json").write_text(json.dumps(summ, indent=2, default=str))
    print(json.dumps(summ, indent=2, default=str))


if __name__ == "__main__":
    main()
