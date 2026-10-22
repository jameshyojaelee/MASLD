#!/usr/bin/env python3
"""H48 step 3: score knockouts against human signatures; N1 gene-label permutation null,
N2 NC pseudo-knockout null, calibration, BH. Writes OUT/step3/."""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata
from statsmodels.stats.multitest import multipletests

from h48_stats import empirical_two_sided, fraction_below, read_null_array

SEED = 48
NPERM = int(os.environ.get("H48_NPERM", 20000))
NPERM_CAL = int(os.environ.get("H48_NPERM_CAL", 2000))

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
a = ap.parse_args()
s1, s2, s3 = (os.path.join(a.out, x) for x in ("step1", "step2", "step3"))
os.makedirs(s3, exist_ok=False)

univ = pd.read_csv(os.path.join(s1, "gene_universe.tsv"), sep="\t")
dims = json.load(open(os.path.join(s2, "step2_dims.json")))
kos = dims["kos"]; G = dims["n_genes"]; B = dims["B"]; BF = dims["B_full"]; nK = len(kos)
eff = pd.read_csv(os.path.join(s2, "effects_real_log2fc.tsv.gz"), sep="\t").set_index("gene_id")
effnf = pd.read_csv(os.path.join(s2, "effects_real_dropA_log2fc.tsv.gz"), sep="\t").set_index("gene_id")
fitset = pd.read_csv(os.path.join(s2, "fit_set.tsv"), sep="\t").set_index("target_gene")
assert list(eff.index) == list(univ.gene_id) and list(eff.columns) == kos
E = eff.values  # G x K
E_nf = effnf.values
tested_cols = [c for c in kos if fitset.loc[c, "analysis"] == "tested"]
E_cm = E - eff[tested_cols].values.mean(1, keepdims=True)
sig_union = [l.strip() for l in open(os.path.join(s2, "sig_union_genes.txt"))]
nullfull_path = os.path.join(s2, "null_full_f32.bin")
nullsig_path = os.path.join(s2, "null_sig_f32.bin")
elig = pd.read_csv(os.path.join(s1, "knockout_eligibility.tsv"), sep="\t").set_index("target_gene")
gid_pos = {g: i for i, g in enumerate(univ.gene_id)}
sym2gid = dict(zip(univ.gene_name, univ.gene_id))
bins = univ.expr_bin.values


def spearman_rows(M, w_rank):
    """Spearman of each row of M with a fixed vector whose ranks are w_rank."""
    R = rankdata(M, axis=1)
    R = R - R.mean(1, keepdims=True)
    wr = w_rank - w_rank.mean()
    return (R @ wr) / (np.sqrt((R ** 2).sum(1)) * np.sqrt((wr ** 2).sum()))


def cosine(x, w):
    return float(x @ w / (np.linalg.norm(x) * np.linalg.norm(w)))


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); m = ~np.isnan(p)
    if m.any():
        q[m] = multipletests(p[m], method="fdr_bh")[1]
    return q


summary = {}
for sigfile in sorted(f for f in os.listdir(s1) if f.startswith("signature_")):
    label = sigfile.replace("signature_", "").replace(".tsv", "")
    sig = pd.read_csv(os.path.join(s1, sigfile), sep="\t")
    sg = sig.mouse_gene_id.values
    w = sig.human_t.values.astype(float)
    pos = np.array([gid_pos[g] for g in sg])
    upos = np.array([sig_union.index(g) for g in sg])
    # N1 permutation index matrix: expression-bin-matched, without replacement within a permutation
    rng = np.random.default_rng(SEED)
    P = np.empty((NPERM, len(sg)), dtype=np.int64)
    bin_members = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    sb = bins[pos]
    for b in np.unique(sb):
        cols = np.where(sb == b)[0]
        mem = bin_members[b]
        for r in range(NPERM):
            P[r, cols] = rng.choice(mem, size=len(cols), replace=False)
    rows = []
    cal_p = []
    for k, g in enumerate(kos):
        tg = sym2gid.get(g)
        keep = np.array([x != tg for x in sg])
        wk = w[keep]; wr = rankdata(wk)
        x = E[pos[keep], k]
        obs = spearman_rows(x[None, :], wr)[0]
        nul = spearman_rows(E[:, k][P[:, keep]], wr)
        p1 = empirical_two_sided(obs, nul)
        # common-mode-removed and no-fast scores, same N1 perms
        x_cm = E_cm[pos[keep], k]
        obs_cm = spearman_rows(x_cm[None, :], wr)[0]
        nul_cm = spearman_rows(E_cm[:, k][P[:, keep]], wr)
        p1_cm = empirical_two_sided(obs_cm, nul_cm)
        obs_nf = spearman_rows(E_nf[pos[keep], k][None, :], wr)[0] if not np.isnan(E_nf[0, k]) else np.nan
        tested = fitset.loc[g, "analysis"] == "tested"
        # N2: NC pseudo-knockouts (tested knockouts only)
        mu = sd = z2 = p2 = p2_emp = np.nan; loo = [np.nan]
        if tested:
            nullsig = read_null_array(nullsig_path, (nK, B, len(sig_union)), index=k)
            n2 = spearman_rows(nullsig[:, upos[keep]].astype(float), wr)
            mu, sd = n2.mean(), n2.std(ddof=1)
            if np.isfinite(sd) and sd > 0:
                z2 = (obs - mu) / sd
                p2 = 2 * norm.sf(abs(z2))
            p2_emp = empirical_two_sided(obs - mu, n2 - mu)
            # N2 leave-one-draw-out calibration
            loo = []
            for b in range(B):
                o = np.delete(n2, b)
                loo.append(2 * norm.sf(abs((n2[b] - o.mean()) / o.std(ddof=1))))
            # N1 calibration on NC pseudo-knockout draws (full-universe draws)
            nullfull = read_null_array(nullfull_path, (nK, BF, G), index=k)
            for b in range(BF):
                xb = nullfull[b].astype(float)
                ob = spearman_rows(xb[pos[keep]][None, :], wr)[0]
                nb = spearman_rows(xb[P[:NPERM_CAL][:, keep]], wr)
                cal_p.append(dict(ko=g, draw=b, score=ob, p_N1=empirical_two_sided(ob, nb)))
        rows.append(dict(knockout=g, analysis=fitset.loc[g, "analysis"], n_cells=int(elig.loc[g, "n_cells_total"]),
                         n_guides_ge10=int(elig.loc[g, "n_guides_ge10"]), n_ko_samples=int(fitset.loc[g, "n_ko_samples"]),
                         letters=elig.loc[g, "letters"], cells_per_guide_letter=elig.loc[g, "cells_per_guide_letter"],
                         n_sig_genes=int(keep.sum()), spearman=obs, cosine=cosine(x, wk), p_N1=p1,
                         N2_null_mean=mu, N2_null_sd=sd, z_N2=z2, p_N2=p2, p_N2_empirical=p2_emp,
                         N2_loo_frac_p_lt_005=fraction_below(loo, 0.05) if tested else np.nan,
                         spearman_commonmode_removed=obs_cm, p_N1_commonmode_removed=p1_cm,
                         spearman_drop_letter_A=obs_nf))
    df = pd.DataFrame(rows)
    t = df.analysis == "tested"
    for c in ["q_N1", "q_N2", "q_N1_commonmode_removed"]:
        df[c] = np.nan
    df.loc[t, "q_N1"] = bh(df.loc[t, "p_N1"]); df.loc[t, "q_N2"] = bh(df.loc[t, "p_N2"])
    df.loc[t, "q_N1_commonmode_removed"] = bh(df.loc[t, "p_N1_commonmode_removed"])
    df["direction"] = np.where(df.spearman > 0, "toward_MASLD", "away_from_MASLD")
    df["pass_N1"] = (df.q_N1 < 0.05) & t
    df["robust_N1_and_N2"] = df.pass_N1 & (df.q_N2 < 0.05) & (np.sign(df.z_N2) == np.sign(df.spearman))
    df = df.sort_values("spearman", ascending=False)
    df["rank_by_spearman"] = np.arange(1, len(df) + 1)
    df.to_csv(os.path.join(s3, f"knockout_scores_{label}.tsv"), sep="\t", index=False)
    cal = pd.DataFrame(cal_p)
    cal.to_csv(os.path.join(s3, f"calibration_N1_on_NC_draws_{label}.tsv.gz"), sep="\t", index=False)
    df_all = df
    df = df_all[df_all.analysis == "tested"]
    nf = df.dropna(subset=["spearman_drop_letter_A"])
    summary[label] = dict(
        n_tested=int(len(df)), n_descriptive_only=int((df_all.analysis != "tested").sum()), n_sig_genes_max=int(df.n_sig_genes.max()),
        n_pass_N1_toward=int((df.pass_N1 & (df.spearman > 0)).sum()),
        n_pass_N1_away=int((df.pass_N1 & (df.spearman < 0)).sum()),
        n_robust_toward=int((df.robust_N1_and_N2 & (df.spearman > 0)).sum()),
        n_robust_away=int((df.robust_N1_and_N2 & (df.spearman < 0)).sum()),
        n_pass_N1_commonmode_removed=int((df.q_N1_commonmode_removed < 0.05).sum()),
        median_spearman_real=float(df.spearman.median()),
        median_N2_null_mean=float(df.N2_null_mean.median()),
        median_N2_null_sd=float(df.N2_null_sd.median()),
        N1_calibration_frac_p_lt_005=fraction_below(cal.p_N1, 0.05),
        N1_calibration_frac_p_lt_001=fraction_below(cal.p_N1, 0.01),
        N1_calibration_n=int(len(cal)),
        N2_loo_frac_p_lt_005_mean=float(df.N2_loo_frac_p_lt_005.mean()),
        frac_real_p_N2_lt_005=fraction_below(df.p_N2, 0.05),
        spearman_scores_vs_dropA_fit=float(pd.Series(nf.spearman).corr(pd.Series(nf.spearman_drop_letter_A), method="spearman")),
        spearman_scores_vs_commonmode_removed=float(df.spearman.corr(df.spearman_commonmode_removed, method="spearman")),
    )
json.dump(summary, open(os.path.join(s3, "summary.json"), "w"), indent=1)
print(json.dumps(summary, indent=1))
