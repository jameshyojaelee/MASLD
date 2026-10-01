#!/usr/bin/env python3
"""Model A T5b: stage dependence and lineage dependence of oriented allelic ratios.

Run only after T5a has written stage_mde.tsv. Development participants only.

Per gene, y = log2((a + 0.5) / (r + 0.5)) with a / r the reads on the lead eQTL's
alt / ref haplotype at the participant's chosen het tag; v = delta-method variance;
weights 1 / (v + tau^2), tau^2 by Paule-Mandel from the model without the tested
term. Covariates: hepatocyte and two other lineage mRNA fractions (from T5a), log
reads at the tag, one indicator per cohort.

Stage (PRESPEC_A 8): slope per step of S (0 / 1 / 2) in staged first-biopsy
participants; genes with >= 30. Null: Freedman-Lane, residuals of the stage-free
weighted fit permuted within cohort, 1,000 permutations, seed 20260923; p =
(1 + #|t_perm| >= |t_obs|) / 1001; BH across genes; Storey pi1 at lambda 0.5.
Two-indicator sensitivity (F2 and F3-4 vs F0-1) reported without a null.

Lineage: F-statistic of the three lineage fractions over log reads + cohort, all
development participants with composition; genes with >= 30; Freedman-Lane as above.

G3 (PRESPEC_A 7): two fixed participant halves by hash parity with relatives kept
together; per gene stage t in each half; slope of half-2 t on half-1 t across genes;
chromosome block bootstrap (1,000 draws, seed 20260923). The same statistic for
negative contexts that carry no stage: log reads at the tag, and stage shuffled
within cohort x hepatocyte-fraction tertile (seed 20260923). G3 passes if the
bootstrap 2.5% bound of the stage slope is above 0 and above the point value of
every negative context.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from t5a_allelic_table import bh, design, log2_ratio, paule_mandel  # noqa: E402

SEED, N_PERM, N_BOOT, MIN_N, MIN_HALF = 20260923, 1000, 1000, 30, 15


def whitened(g, lineages, term):
    """Whitened y, full X, reduced X, cohort codes for one gene; term in {'S', 'lineage', ctx}."""
    y, v = log2_ratio(g["a_o"].to_numpy(float), g["r_o"].to_numpy(float))
    if term == "lineage":
        Xr = np.column_stack([np.log(g["n"].to_numpy(float)), pd.get_dummies(g["cohort"]).to_numpy(float)])
        Xf = np.column_stack([g[l].to_numpy(float) for l in lineages] + [Xr])
    else:
        Xr = design(g, lineages, with_stage=False)
        Xf = np.column_stack([g[term].to_numpy(float), Xr])
    tau2 = paule_mandel(y, v, Xr)
    sw = 1 / np.sqrt(v + tau2)
    return y * sw, Xf * sw[:, None], Xr * sw[:, None], g["cohort"].to_numpy(), tau2


def ols_stats(X, Y, k_test):
    """t of column 0 (k_test == 1) or F of the first k_test columns, for each column of Y."""
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    s2 = (R ** 2).sum(0) / (n - p)
    if k_test == 1:
        return B[0] / np.sqrt(s2 * XtX_inv[0, 0])
    V = XtX_inv[:k_test, :k_test]
    Vi = np.linalg.pinv(V)
    b = B[:k_test]
    return np.einsum("ij,ik,kj->j", b, Vi, b) / k_test / s2


def freedman_lane(yw, Xf, Xr, cohort, k_test, rng):
    obs = ols_stats(Xf, yw[:, None], k_test)[0]
    br = np.linalg.lstsq(Xr, yw, rcond=None)[0]
    fit, res = Xr @ br, yw - Xr @ br
    idx = {c: np.where(cohort == c)[0] for c in np.unique(cohort)}
    Y = np.empty((len(yw), N_PERM))
    for b in range(N_PERM):
        perm = np.arange(len(yw))
        for ii in idx.values():
            perm[ii] = ii[rng.permutation(len(ii))]
        Y[:, b] = fit + res[perm]
    null = ols_stats(Xf, Y, k_test)
    stat = np.abs if k_test == 1 else (lambda z: z)
    p = (1 + np.sum(stat(null) >= stat(obs))) / (N_PERM + 1)
    return float(obs), float(p)


def storey_pi1(p, lam=0.5):
    p = np.asarray(p)
    return float(1 - min(1.0, np.mean(p > lam) / (1 - lam))) if len(p) else None


def halves(ids, kin, xwalk):
    """Hash-parity halves; first-degree and same-individual relatives share the lowest-hash member's half."""
    h = {i: hashlib.sha256(("MASLD-A-G3" + i).encode()).hexdigest() for i in ids}
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
    return {i: int(h[find(i)][-1], 16) % 2 for i in ids}


def split_slope(t1, t2):
    t1, t2 = np.asarray(t1), np.asarray(t2)
    return float(np.polyfit(t1, t2, 1)[0]) if len(t1) > 2 and np.std(t1) > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--participant-gene", required=True)
    ap.add_argument("--t5a-dir", required=True)
    ap.add_argument("--identity-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    t5a = json.loads((Path(a.t5a_dir) / "t5a_summary.json").read_text())
    lineages = t5a["lineages"]
    mde = pd.read_csv(Path(a.t5a_dir) / "stage_mde.tsv", sep="\t")
    main_eff = pd.read_csv(Path(a.t5a_dir) / "main_allelic_effects.tsv", sep="\t")
    pg = pd.read_csv(a.participant_gene, sep="\t")
    pg = pg[pg["passes_guard"]].dropna(subset=lineages)
    ss = np.random.SeedSequence(SEED)

    # stage
    staged = pg[pg["stage_eligible"].astype(bool)].copy()
    genes = sorted(set(mde["gene_id"]))
    rngs = dict(zip(genes, [np.random.default_rng(s) for s in ss.spawn(len(genes))]))
    rows = []
    for gid in genes:
        g = staged[staged["gene_id"] == gid]
        if len(g) < MIN_N:
            continue
        yw, Xf, Xr, coh, tau2 = whitened(g, lineages, "S")
        t_obs, p = freedman_lane(yw, Xf, Xr, coh, 1, rngs[gid])
        b = np.linalg.lstsq(Xf, yw, rcond=None)[0][0]
        g2 = g.assign(S1=(g["S"] == 1).astype(float), S2=(g["S"] == 2).astype(float))
        yv, vv = log2_ratio(g2["a_o"].to_numpy(float), g2["r_o"].to_numpy(float))
        X2 = np.column_stack([g2["S1"], g2["S2"], design(g2, lineages, with_stage=False)]) / np.sqrt(vv + tau2)[:, None]
        b2 = np.linalg.lstsq(X2, yv / np.sqrt(vv + tau2), rcond=None)[0]
        rows.append({"gene_id": gid, "gene_name": g["gene_name"].iloc[0], "chrom": g["chrom"].iloc[0], "n": len(g),
                     "n_by_S": "/".join(str(int((g["S"] == s_).sum())) for s_ in (0, 1, 2)),
                     "tau2": tau2, "stage_slope_log2_per_step": float(b), "t": t_obs, "p_perm": p,
                     "F2_vs_F01_log2": float(b2[0]), "F34_vs_F01_log2": float(b2[1])})
    st = pd.DataFrame(rows)
    st["q"] = bh(st["p_perm"])
    st = st.merge(main_eff[["gene_id", "log2_allelic_ratio", "q", "mde80_log2_per_step"]].rename(
        columns={"q": "main_q", "log2_allelic_ratio": "main_log2_allelic_ratio"}), on="gene_id", how="left")
    # a slope toward larger |ratio| with stage is positive in this column
    st["slope_toward_larger_effect"] = st["stage_slope_log2_per_step"] * np.sign(st["main_log2_allelic_ratio"])
    st.sort_values("p_perm").to_csv(out / "stage_tests.tsv", sep="\t", index=False)
    m = len(st)
    summ = {"stage": {"genes": m, "q_lt_0.05": int((st["q"] < 0.05).sum()), "q_lt_0.10": int((st["q"] < 0.10).sum()),
                      "pi1": storey_pi1(st["p_perm"]), "p_floor": 1 / (N_PERM + 1),
                      "bh_k_min_at_q_0.05": math.ceil((1 / (N_PERM + 1)) * m / 0.05) if m else None,
                      "median_abs_slope": float(st["stage_slope_log2_per_step"].abs().median()) if m else None,
                      "median_slope_toward_larger_effect_main_q_lt_0.05": float(st.loc[st["main_q"] < 0.05, "slope_toward_larger_effect"].median()) if m else None}}

    # lineage
    lin_all = pg.copy()
    lgenes = sorted(g for g, d in lin_all.groupby("gene_id") if len(d) >= MIN_N)
    rngs = dict(zip(lgenes, [np.random.default_rng(s) for s in np.random.SeedSequence([SEED, 1]).spawn(len(lgenes))]))
    rows = []
    for gid in lgenes:
        g = lin_all[lin_all["gene_id"] == gid]
        yw, Xf, Xr, coh, tau2 = whitened(g, lineages, "lineage")
        F, p = freedman_lane(yw, Xf, Xr, coh, len(lineages), rngs[gid])
        b = np.linalg.lstsq(Xf, yw, rcond=None)[0][:len(lineages)]
        rows.append({"gene_id": gid, "gene_name": g["gene_name"].iloc[0], "n": len(g), "tau2": tau2, "F": F, "p_perm": p,
                     **{f"coef_{l}": float(c) for l, c in zip(lineages, b)}})
    ln = pd.DataFrame(rows)
    ln["q"] = bh(ln["p_perm"])
    ln.sort_values("p_perm").to_csv(out / "lineage_tests.tsv", sep="\t", index=False)
    summ["lineage"] = {"genes": len(ln), "q_lt_0.05": int((ln["q"] < 0.05).sum()), "pi1": storey_pi1(ln["p_perm"])}

    # G3
    xwalk = pd.read_csv(Path(a.identity_dir) / "sample_to_individual.tsv", sep="\t")
    kin = pd.read_csv(Path(a.identity_dir) / "kinship_related_pairs.tsv", sep="\t")
    half = halves(sorted(staged["individual_id"].unique()), kin, xwalk)
    staged["half"] = staged["individual_id"].map(half)
    rng = np.random.default_rng([SEED, 2])
    staged["logn"] = np.log(staged["n"])
    fake = staged[["individual_id", "cohort", "S", "Hepatocytes"]].drop_duplicates("individual_id").sort_values("individual_id")
    fake["hep_tertile"] = fake.groupby("cohort")["Hepatocytes"].transform(
        lambda z: pd.qcut(z.rank(method="first"), 3, labels=False))
    fake["S_fake"] = fake.groupby(["cohort", "hep_tertile"])["S"].transform(lambda z: rng.permutation(z.to_numpy()))
    staged = staged.merge(fake[["individual_id", "S_fake"]], on="individual_id")
    contexts = {"stage": "S", "log_reads": "logn", "stage_shuffled_within_hepatocyte_tertile": "S_fake"}
    tt = {k: [] for k in contexts}
    chrom_of = {}
    for gid, g in staged.groupby("gene_id"):
        g0, g1 = g[g["half"] == 0], g[g["half"] == 1]
        if len(g0) < MIN_HALF or len(g1) < MIN_HALF:
            continue
        chrom_of[gid] = g["chrom"].iloc[0]
        for k, col in contexts.items():
            ts = []
            for gh in (g0, g1):
                yw, Xf, Xr, coh, _ = whitened_ctx(gh, lineages) if col == "logn" else whitened(gh, lineages, col)
                ts.append(ols_stats(Xf, yw[:, None], 1)[0])
            tt[k].append((gid, ts[0], ts[1]))
    g3 = {"genes": len(chrom_of), "min_per_half": MIN_HALF}
    chroms = sorted(set(chrom_of.values()))
    boot_idx = [rng.choice(len(chroms), len(chroms), replace=True) for _ in range(N_BOOT)]
    for k in contexts:
        d = pd.DataFrame(tt[k], columns=["gene_id", "t_half0", "t_half1"])
        d["chrom"] = d["gene_id"].map(chrom_of)
        d.to_csv(out / f"g3_split_t_{k}.tsv", sep="\t", index=False)
        by = {c: d[d["chrom"] == c] for c in chroms}
        est = split_slope(d["t_half0"], d["t_half1"])
        boots = []
        for bi in boot_idx:
            dd = pd.concat([by[chroms[i]] for i in bi])
            boots.append(split_slope(dd["t_half0"], dd["t_half1"]))
        g3[k] = {"slope": est, "boot_2.5": float(np.nanpercentile(boots, 2.5)), "boot_97.5": float(np.nanpercentile(boots, 97.5))}
    g3["pass"] = bool(g3["stage"]["boot_2.5"] > 0 and all(g3["stage"]["boot_2.5"] > g3[k]["slope"] for k in contexts if k != "stage"))
    summ["g3"] = g3
    (out / "t5b_summary.json").write_text(json.dumps(summ, indent=2, default=str))
    print(json.dumps(summ, indent=2, default=str))


def whitened_ctx(g, lineages):
    """Negative context 'log reads': log reads is the tested column; the model drops it from the covariates."""
    y, v = log2_ratio(g["a_o"].to_numpy(float), g["r_o"].to_numpy(float))
    Xr = np.column_stack([g[l].to_numpy(float) for l in lineages] + [pd.get_dummies(g["cohort"]).to_numpy(float)])
    Xf = np.column_stack([np.log(g["n"].to_numpy(float)), Xr])
    tau2 = paule_mandel(y, v, Xr)
    sw = 1 / np.sqrt(v + tau2)
    return y * sw, Xf * sw[:, None], Xr * sw[:, None], g["cohort"].to_numpy(), tau2


if __name__ == "__main__":
    main()
