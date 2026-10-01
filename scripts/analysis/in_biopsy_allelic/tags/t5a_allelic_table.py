#!/usr/bin/env python3
"""Model A T5a: per-participant allelic table, main allelic effects, GTEx replication.

Development participants only; sealed individuals are dropped before any count is
read. Reads no statistic that relates allelic ratio to stage.

Steps (PRESPEC_A sections 6-8):
1. One library per genetic individual: its home-cohort library, first biopsy in
   GSE193066, then most identity sites called. Possibly mixed libraries excluded.
2. Stage S (0 = F0-1, 1 = F2, 2 = F3-4) from Kleiner stage; GSE213621 from its
   condition label (its `fibrosis_stage` column is an ordinal code with Control = 0).
   Source controls (C = 1) are not staged; GSE130970 controls are the 6 from
   Hoang 2019 Table S1.
3. Tags whose position lies in exons of more than one gene, tags at REDIportal v3
   RNA-editing sites, and tags in HLA genes (dense polymorphism, mapping bias) are dropped.
4. Het calls at tags with the T2 rule. FFPE gate G1-FFPE for C>T / G>A tags.
5. Per participant x gene, the het tag with the most reads (ties: lowest position);
   reads oriented to the lead eQTL's alt haplotype with the T1 LD sign.
6. Detection guard: per gene, the depth n_min at which the T2 rule calls a het with
   probability >= 0.9 under binomial sampling at the gene's median imbalance.
7. Main effect per gene (>= 10 participants): beta-binomial MLE of the oriented
   alt-haplotype fraction, likelihood-ratio test against 0.5; BH across genes.
   Cohort heterogeneity: likelihood-ratio test of cohort-specific means.
8. G2: sign agreement with GTEx v8 liver slope; Spearman with GTEx log2 aFC.
9. Stage MDE per gene (design only, before T5b): WLS SE of the stage slope with
   hepatocyte + the two other lineages with the largest median mRNA fraction in
   staged participants + log reads at the tag + cohort; tau^2 from the stage-free fit.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

FFPE = ["GSE213621", "GSE240729"]
STAGED = ["GSE213621", "GSE135251", "GSE162694", "GSE130970", "GSE193066", "GSE240729", "GSE174478"]
KLEINER = ["GSE135251", "GSE162694", "GSE130970", "GSE193066", "GSE240729", "GSE174478"]
DEAMINATION = {("C", "T"), ("G", "A")}
MIN_PARTICIPANTS_MAIN, MIN_STAGED = 10, 30
LN2 = np.log(2.0)


# ---------- participants ----------
def choose_libraries(xwalk, mixed, seal, placement):
    x = xwalk[xwalk["individual_id"].notna() & ~xwalk["run"].isin(mixed)].copy()
    x = x.merge(seal, on="individual_id", how="inner")
    x = x[(x["cohort"] == x["home_cohort"]) & ~x["sealed"]]
    first = placement.set_index("sample_id")["is_first_biopsy"].astype(bool)
    x["first_biopsy"] = x["run"].map(first).where(x["cohort"] == "GSE193066", True).fillna(False).astype(bool)
    x = x.sort_values(["individual_id", "first_biopsy", "autosomal_called", "run"],
                      ascending=[True, False, False, True]).drop_duplicates("individual_id")
    return x


def stage_coding(meta, ctrl130970):
    m = meta.rename(columns={"sample_id": "run", "dataset": "cohort"})[["run", "cohort", "condition", "fibrosis_stage"]].copy()
    m["S"], m["C"] = np.nan, 0
    kle = m["cohort"].isin(KLEINER)
    fs = m["fibrosis_stage"]
    m.loc[kle & fs.isin([0, 1]), "S"] = 0
    m.loc[kle & (fs == 2), "S"] = 1
    m.loc[kle & fs.isin([3, 4]), "S"] = 2
    g = m["cohort"] == "GSE213621"
    m.loc[g, "S"] = m.loc[g, "condition"].map({"Fibrosis_F0F1": 0, "Fibrosis_F2": 1, "Fibrosis_F3F4": 2})
    ctrl = (m["condition"] == "Control") & (m["cohort"] != "GSE130970")
    ctrl |= m["run"].isin(ctrl130970)
    m.loc[ctrl, "C"] = 1
    m.loc[ctrl, "S"] = np.nan
    return m[["run", "S", "C"]]


# ---------- tags ----------
def multi_gene_positions(tags, gtf):
    ex = []
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] == "exon":
            ex.append((f[0], int(f[3]), int(f[4]), f[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".")[0]))
    ex = pd.DataFrame(ex, columns=["chrom", "start", "end", "gene"])
    n_genes = {}
    for c, t in tags.drop_duplicates(["chrom", "pos"]).groupby("chrom"):
        e = ex[ex["chrom"] == c]
        s, en, g = e["start"].to_numpy(), e["end"].to_numpy(), e["gene"].to_numpy()
        for p in t["pos"]:
            n_genes[(c, p)] = len(set(g[(s <= p) & (en >= p)]))
    return n_genes


def het_call(ref, alt, rule):
    dp = ref + alt
    minor = np.minimum(ref, alt)
    return (dp >= rule["min_dp"]) & (minor >= rule["min_minor_reads"]) & (minor >= rule["min_minor_frac"] * dp)


def call_probability(n, pi, rule):
    x = np.arange(n + 1)
    pmf = stats.binom.pmf(x, n, pi)
    minor = np.minimum(x, n - x)
    ok = (n >= rule["min_dp"]) & (minor >= rule["min_minor_reads"]) & (minor >= rule["min_minor_frac"] * n)
    return float(pmf[ok].sum())


def n_min_for(pi, rule, target=0.9):
    grid = list(range(int(rule["min_dp"]), 200)) + list(range(200, 10001, 50))
    for n in grid:
        if call_probability(n, pi, rule) >= target:
            return n
    return None


# ---------- beta-binomial ----------
def bb_negll(params, k, n, groups):
    b, g = params[:-1], params[-1]
    p = special.expit(b[groups])
    rho = special.expit(g)
    s = (1 - rho) / rho
    a_, b_ = p * s, (1 - p) * s
    return -np.sum(special.betaln(k + a_, n - k + b_) - special.betaln(a_, b_))


def bb_fit(k, n, groups=None, fix_b=None):
    groups = np.zeros(len(k), int) if groups is None else groups
    G = groups.max() + 1
    if fix_b is not None:
        f = lambda gpar: bb_negll(np.r_[np.full(G, fix_b), gpar], k, n, groups)
        r = optimize.minimize(f, x0=[-3.0], method="L-BFGS-B", bounds=[(-12, 6)])
        return np.r_[np.full(G, fix_b), r.x], r.fun
    p0 = np.clip((k.sum() + 0.5) / (n.sum() + 1), 1e-3, 1 - 1e-3)
    x0 = np.r_[np.full(G, special.logit(p0)), -3.0]
    r = optimize.minimize(bb_negll, x0=x0, args=(k, n, groups), method="L-BFGS-B",
                          bounds=[(-12, 12)] * G + [(-12, 6)])
    return r.x, r.fun


def bb_se_intercept(params, k, n):
    h = 1e-4
    f = lambda v: bb_negll(v, k, n, np.zeros(len(k), int))
    H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            e_i, e_j = np.eye(2)[i] * h, np.eye(2)[j] * h
            H[i, j] = (f(params + e_i + e_j) - f(params + e_i - e_j) - f(params - e_i + e_j) + f(params - e_i - e_j)) / (4 * h * h)
    try:
        return float(np.sqrt(np.linalg.inv(H)[0, 0]))
    except np.linalg.LinAlgError:
        return np.nan


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        o = np.argsort(p[ok])
        r = p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty(ok.sum()); qq[o] = np.minimum(r, 1)
        q[ok] = qq
    return q


# ---------- WLS helpers (shared with T5b) ----------
def log2_ratio(a, r):
    y = np.log2((a + 0.5) / (r + 0.5))
    v = (1 / (a + 0.5) + 1 / (r + 0.5)) / LN2 ** 2
    return y, v


def paule_mandel(y, v, X):
    """tau^2 >= 0 at which the weighted residual sum of squares equals n - p."""
    n, p = X.shape

    def q(t):
        w = 1 / (v + t)
        beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        return np.sum(w * (y - X @ beta) ** 2) - (n - p)

    if q(0.0) <= 0:
        return 0.0
    hi = 1.0
    while q(hi) > 0 and hi < 1e4:
        hi *= 4
    return float(optimize.brentq(q, 0.0, hi, xtol=1e-8)) if q(hi) <= 0 else hi


def design(df, lineages, with_stage=True):
    """[S] + lineage mRNA fractions + log reads at the tag + one indicator per cohort (no separate intercept)."""
    cols = []
    if with_stage:
        cols.append(df["S"].to_numpy(float))
    for l in lineages:
        cols.append(df[l].to_numpy(float))
    cols.append(np.log(df["n"].to_numpy(float)))
    coh = pd.get_dummies(df["cohort"]).to_numpy(float)
    return np.column_stack(cols + [coh])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", required=True)
    ap.add_argument("--ad", required=True)
    ap.add_argument("--t2", required=True)
    ap.add_argument("--identity-dir", required=True)
    ap.add_argument("--seal-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--gse130970-controls", required=True)
    ap.add_argument("--gse193066-placement", required=True)
    ap.add_argument("--composition", required=True)
    ap.add_argument("--gtex-egenes", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--rediportal-hits", required=True, help="chrom<TAB>pos of tag positions listed in REDIportal")
    ap.add_argument("--restricted-out", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sex-het-threshold", type=float, default=0.05)
    a = ap.parse_args()
    out, rout = Path(a.out), Path(a.restricted_out)
    out.mkdir(parents=True, exist_ok=True); rout.mkdir(parents=True, exist_ok=True)
    summary = {}

    rule = json.loads(Path(a.t2).read_text())["chosen_rule"]
    if rule is None:
        raise SystemExit("G1' failed: no T2 rule with FDP <= 0.02")
    summary["t2_rule"] = rule
    summary["g1_prime_pass"] = bool(rule["het_sensitivity"] >= 0.85)
    if not summary["g1_prime_pass"]:
        raise SystemExit(f"G1' failed: sensitivity {rule['het_sensitivity']:.3f} < 0.85")

    # participants
    xwalk = pd.read_csv(Path(a.identity_dir) / "sample_to_individual.tsv", sep="\t")
    mixed = set(pd.read_csv(Path(a.seal_dir) / "excluded_possibly_mixed_libraries.tsv", sep="\t")["run"])
    seal = pd.read_csv(Path(a.seal_dir) / "sealed_individuals.tsv", sep="\t")
    placement = pd.read_csv(a.gse193066_placement, sep="\t")
    lib = choose_libraries(xwalk, mixed, seal, placement)
    meta = pd.read_csv(a.metadata, keep_default_na=False, na_values=[""])
    c130 = pd.read_csv(a.gse130970_controls, sep="\t")
    stage = stage_coding(meta, set(c130.loc[c130["source_control_status"] == "Control", "run"]))
    lib = lib.merge(stage, on="run", how="left")
    lib["stage_eligible"] = lib["cohort"].isin(STAGED) & lib["S"].notna() & (lib["C"] == 0) & lib["first_biopsy"]
    summary["development_participants"] = lib.groupby("cohort").size().to_dict()
    summary["development_staged_participants"] = {f"{c}|S{int(s_)}": int(n) for (c, s_), n in
                                                  lib[lib["stage_eligible"]].groupby(["cohort", "S"]).size().items()}

    # sex check (G2 part)
    s = lib[lib["chrX_called"] >= 30].copy()
    s["geno_sex"] = np.where(s["chrX_het_rate"] >= a.sex_het_threshold, "F", "M")
    both = s[s["expr_sex"].isin(["M", "F"])]
    summary["sex_check"] = {"libraries_compared": int(len(both)),
                            "concordant_fraction": float((both["expr_sex"] == both["geno_sex"]).mean()) if len(both) else None,
                            "chrX_het_rate_by_expr_sex_median": both.groupby("expr_sex")["chrX_het_rate"].median().to_dict(),
                            "threshold": a.sex_het_threshold}

    # composition
    comp = pd.read_csv(a.composition, sep="\t")
    ctypes = [c for c in comp.columns if c not in ("cohort", "run")]
    lib = lib.merge(comp.drop(columns="cohort"), on="run", how="left")
    staged_comp = lib.loc[lib["stage_eligible"], ctypes]
    others = staged_comp.drop(columns="Hepatocytes").median().sort_values(ascending=False)
    lineages = ["Hepatocytes"] + list(others.index[:2])
    summary["lineages"] = lineages
    summary["lineage_median_mrna_fraction_staged"] = staged_comp.median()[lineages].to_dict()

    # tags
    tags = pd.read_csv(a.tags, sep="\t")
    ng = multi_gene_positions(tags, a.gtf)
    tags["n_genes_at_pos"] = [ng[(c, p)] for c, p in zip(tags["chrom"], tags["pos"])]
    summary["tag_rows_in_multi_gene_exons_dropped"] = int((tags["n_genes_at_pos"] > 1).sum())
    tags = tags[tags["n_genes_at_pos"] == 1]
    red = pd.read_csv(a.rediportal_hits, sep="\t", header=None, names=["chrom", "pos"])
    at_red = pd.Series(list(zip(tags["chrom"], tags["pos"]))).isin(set(zip(red["chrom"], red["pos"]))).to_numpy()
    hla = tags["gene_name"].astype(str).str.startswith("HLA-").to_numpy()
    summary["tag_rows_at_rediportal_sites_dropped"] = int(at_red.sum())
    summary["tag_rows_in_hla_genes_dropped"] = int((hla & ~at_red).sum())
    tags = tags[~at_red & ~hla]
    summary["tag_snvs_used"] = int(tags.drop_duplicates(["chrom", "pos"]).shape[0])
    summary["egenes_with_tags_used"] = int(tags["gene_id"].nunique())
    tags["deamination"] = [(r, x) in DEAMINATION for r, x in zip(tags["ref"], tags["alt"])]

    ad = pd.read_csv(a.ad, sep="\t")
    ad = ad[ad["run"].isin(set(lib["run"]))]
    d = ad.merge(tags[["gene_id", "gene_name", "chrom", "pos", "ref", "alt", "orientation", "deamination"]],
                 on=["chrom", "pos", "ref", "alt"])
    d = d.merge(lib[["run", "individual_id"]], on="run")
    d["het"] = het_call(d["ref_reads"].to_numpy(), d["alt_reads"].to_numpy(), rule)
    h = d[d["het"]].copy()
    h["alt_frac"] = h["alt_reads"] / (h["ref_reads"] + h["alt_reads"])
    summary["het_observations"] = int(len(h))
    summary["median_ref_fraction_at_called_hets"] = {"all": float(1 - h["alt_frac"].median()),
                                                     **{c: float(1 - v.median()) for c, v in h.groupby("cohort")["alt_frac"]}}

    # G1-FFPE
    g1f = {}
    non = h[~h["cohort"].isin(FFPE) & h["deamination"]].groupby(["chrom", "pos"])["alt_frac"].median()
    for c in FFPE:
        f = h[(h["cohort"] == c) & h["deamination"]].groupby(["chrom", "pos"])["alt_frac"].median()
        common = f.index.intersection(non.index)
        diff = float(f[common].median() - non[common].median()) if len(common) else np.nan
        ctrl_f = h[(h["cohort"] == c) & ~h["deamination"]].groupby(["chrom", "pos"])["alt_frac"].median()
        ctrl_n = h[~h["cohort"].isin(FFPE) & ~h["deamination"]].groupby(["chrom", "pos"])["alt_frac"].median()
        cc = ctrl_f.index.intersection(ctrl_n.index)
        g1f[c] = {"deamination_tags_compared": int(len(common)), "median_alt_frac_diff_vs_nonffpe": diff,
                  "other_tags_compared": int(len(cc)),
                  "other_tags_median_alt_frac_diff": float(ctrl_f[cc].median() - ctrl_n[cc].median()) if len(cc) else None,
                  "pass": bool(len(common) >= 20 and abs(diff) <= 0.02)}
        if not g1f[c]["pass"]:
            h = h[~((h["cohort"] == c) & h["deamination"])]
    summary["g1_ffpe"] = g1f

    # orient and pick one tag per participant x gene
    pos_or = h["orientation"] == 1
    h["a_o"] = np.where(pos_or, h["alt_reads"], h["ref_reads"])
    h["r_o"] = np.where(pos_or, h["ref_reads"], h["alt_reads"])
    h["n"] = h["a_o"] + h["r_o"]
    h = h.sort_values(["individual_id", "gene_id", "n", "pos"], ascending=[True, True, False, True])
    pg = h.drop_duplicates(["individual_id", "gene_id"]).copy()
    pg = pg.merge(lib[["individual_id", "S", "C", "stage_eligible"] + ctypes], on="individual_id", how="left")

    # detection guard
    guard = []
    for gid, g in pg.groupby("gene_id"):
        m = float(np.median(np.abs(g["a_o"] / g["n"] - 0.5)))
        nm = n_min_for(0.5 - m, rule)
        guard.append((gid, m, nm))
    guard = pd.DataFrame(guard, columns=["gene_id", "median_abs_imbalance", "n_min"])
    pg = pg.merge(guard, on="gene_id")
    pg["passes_guard"] = pg["n_min"].notna() & (pg["n"] >= pg["n_min"].fillna(np.inf))
    pg.to_csv(rout / "participant_gene_allelic.tsv.gz", sep="\t", index=False, compression="gzip")
    ex = pg[pg["stage_eligible"]].assign(excluded=~pg["passes_guard"]).groupby(["gene_id", "S"])["excluded"].sum().unstack(fill_value=0)
    ex.columns = [f"guard_excluded_S{int(c)}" for c in ex.columns]
    guard = guard.merge(ex, left_on="gene_id", right_index=True, how="left")

    # main effects
    egenes = pd.read_csv(a.gtex_egenes, sep="\t", usecols=["gene_id", "slope", "log2_aFC", "qval"])
    rows = []
    for gid, g in pg[pg["passes_guard"]].groupby("gene_id"):
        if len(g) < MIN_PARTICIPANTS_MAIN:
            continue
        k, n = g["a_o"].to_numpy(float), g["n"].to_numpy(float)
        par1, nll1 = bb_fit(k, n)
        _, nll0 = bb_fit(k, n, fix_b=0.0)
        codes = g["cohort"].where(g["cohort"].map(g["cohort"].value_counts()) >= 5, "other")
        grp = pd.factorize(codes)[0]
        _, nllc = bb_fit(k, n, grp)
        lrt = max(0.0, 2 * (nll0 - nll1))
        het_lrt = max(0.0, 2 * (nll1 - nllc))
        rows.append({"gene_id": gid, "gene_name": g["gene_name"].iloc[0], "n_participants": len(g),
                     "n_cohorts": g["cohort"].nunique(), "median_reads": float(np.median(n)),
                     "log2_allelic_ratio": par1[0] / LN2, "se": bb_se_intercept(par1, k, n) / LN2,
                     "rho": float(special.expit(par1[-1])), "p": float(stats.chi2.sf(lrt, 1)),
                     "cohort_heterogeneity_p": float(stats.chi2.sf(het_lrt, grp.max())) if grp.max() > 0 else np.nan,
                     "n_staged": int(g["stage_eligible"].sum())})
    me = pd.DataFrame(rows)
    me["q"] = bh(me["p"])
    me = me.merge(guard, on="gene_id", how="left").merge(egenes.rename(columns={"slope": "gtex_slope", "log2_aFC": "gtex_log2_aFC", "qval": "gtex_q"}), on="gene_id", how="left")

    # G2 replication
    same = lambda d_: float((np.sign(d_["log2_allelic_ratio"]) == np.sign(d_["gtex_slope"])).mean()) if len(d_) else None
    g2set = me[(me["n_participants"] >= 30) & (me["gtex_slope"].abs() >= 0.3)]
    sig = me[me["q"] < 0.05]
    summary["g2"] = {"genes_tested": int(len(me)), "genes_q_lt_0.05": int(len(sig)),
                     "g2_genes": int(len(g2set)), "g2_sign_agreement": same(g2set),
                     "g2_pass": bool(len(g2set) >= 20 and same(g2set) >= 0.8),
                     "sign_agreement_q_lt_0.05": same(sig), "sign_agreement_all": same(me),
                     "spearman_vs_gtex_log2_aFC_all": float(stats.spearmanr(me["log2_allelic_ratio"], me["gtex_log2_aFC"], nan_policy="omit")[0]),
                     "spearman_vs_gtex_log2_aFC_q_lt_0.05": float(stats.spearmanr(sig["log2_allelic_ratio"], sig["gtex_log2_aFC"], nan_policy="omit")[0]) if len(sig) > 2 else None,
                     "gtex_aFC_sign_equals_slope_sign": float((np.sign(egenes["log2_aFC"]) == np.sign(egenes["slope"])).mean())}
    summary["g2_sex_pass"] = bool(summary["sex_check"]["concordant_fraction"] is not None and summary["sex_check"]["concordant_fraction"] >= 0.98)

    # stage MDE (design only)
    mde = []
    st = pg[pg["passes_guard"] & pg["stage_eligible"]].dropna(subset=lineages)
    for gid, g in st.groupby("gene_id"):
        if len(g) < MIN_STAGED:
            continue
        y, v = log2_ratio(g["a_o"].to_numpy(float), g["r_o"].to_numpy(float))
        Xr = design(g, lineages, with_stage=False)
        tau2 = paule_mandel(y, v, Xr)
        X = design(g, lineages, with_stage=True)
        w = 1 / (v + tau2)
        try:
            cov = np.linalg.pinv(X.T @ (X * w[:, None]))
            se = float(np.sqrt(cov[0, 0]))
        except np.linalg.LinAlgError:
            se = np.nan
        mde.append({"gene_id": gid, "n_staged": len(g), "n_by_S": "/".join(str(int((g["S"] == s_).sum())) for s_ in (0, 1, 2)),
                    "tau2": tau2, "se_stage_slope_log2": se, "mde80_log2_per_step": 2.8 * se})
    mde = pd.DataFrame(mde)
    me = me.merge(mde[["gene_id", "tau2", "se_stage_slope_log2", "mde80_log2_per_step"]], on="gene_id", how="left")
    summary["stage_mde"] = {"genes": int(len(mde)), "median_mde80_log2_per_step": float(mde["mde80_log2_per_step"].median()) if len(mde) else None}

    me = me.sort_values("p")
    me[me["n_participants"] >= MIN_PARTICIPANTS_MAIN].to_csv(out / "main_allelic_effects.tsv", sep="\t", index=False)
    mde.to_csv(out / "stage_mde.tsv", sep="\t", index=False)
    guard.to_csv(out / "detection_guard.tsv", sep="\t", index=False)
    (out / "t5a_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
