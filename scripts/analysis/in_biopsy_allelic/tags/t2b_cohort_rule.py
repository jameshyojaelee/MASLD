#!/usr/bin/env python3
"""Model A-final, B-T2: per-cohort het rule R_c and homozygote-error overdispersion phi.

Spec: docs/technical/agent_exchange/2026-09-29_model_A_likelihood_spec_v1.md (sections 1.2, 1.3, 2.3).

Geuvadis step (30 LCL RNA libraries against 1kGP DNA; t2_work/chr*.ad_vs_dna.tsv.gz):
  - e_j: pooled other-allele fraction of individual j at its DNA-hom rows with n >= 20.
  - phi: beta-binomial ML of other-allele reads o ~ BB(n, e_j, phi) at DNA-hom rows, 10 <= n <= 1000.
  - (mu_het, rho_het): beta-binomial ML of REF reads at DNA-het rows, 10 <= n <= 1000.
  - For every rule, the empirical Geuvadis FDP and het sensitivity beside the model values
    computed from (e_j, phi, mu_het, rho_het), as a check of the model used for the liver cohorts.

Liver step (depth only; no allele fraction is formed):
  - Unit libraries per cohort and set (development: the spec fit set F in the five staged cohorts,
    all development unit libraries in GSE126848, GSE167523, PRJNA512027; sealed: sealed F).
  - Chosen tag per library x gene: the largest n = REF + ALT reads among the gene's tags kept in the
    spec 1.2 tag table (build_tag_table.py; kept_primary, and drop_in_<cohort> for FFPE transitions),
    ties to the smallest position.
  - Row het prior pi = het_ref_alt of the tag in the library's ancestry call (JPT for EAS-called
    libraries of the Japanese cohorts; nearest superpopulation for intermediates, provisional);
    rows whose tag has no frequency row are dropped and counted.
  - Background: the cohort's 90th-percentile e_i for that set (QC v2 cohort_error_rate.tsv).
  - Per rule: FDP = sum (1-pi) P(call|hom,n) / sum [(1-pi) P(call|hom,n) + pi P(call|het,n)],
    sensitivity = sum pi P(call|het,n) / sum pi, expected true-het calls ETP = sum pi P(call|het,n),
    all over rows with n >= d.
  - Choice: among rules with FDP <= 0.02 and sensitivity >= 0.85 (G1'), the largest ETP; ties go
    to the larger d, then the larger k. No such rule: the cohort fails G1' (untestable).

Gate A inputs (aggregate; for the design simulation, spec 3.4 and 5.2): per cohort x set x gene, rows at
the chosen rule's depth, expected het rows (sum pi), expected included rows (sum pi P(call|het,n)) and
depth quantiles; the same by stage; individuals per stage; and per-cohort inputs (R_c, e quantiles,
omega_c = logit of the null-site median REF fraction, unit_all, depth_20_3_0.10; non-estimable cohorts
take the logit of the median over estimable development cohorts), phi, mu_het, rho_het.

Outputs (aggregate only) in --out; nothing per library is written.
"""
import argparse
import functools
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, gammaln

# d, k, F as in T2 (2026-09-25), extended to k <= 5 and F <= 20 on 2026-09-29 before any allelic data were read:
# at the original ceiling (k 3, F 10) two high-background cohorts had no rule with FDP <= 0.02.
GRID = [(d, k, F) for d, k, F in itertools.product((10, 20, 30), (1, 2, 3, 4, 5), (1, 2, 5, 10, 15, 20))]
FDP_MAX, SENS_MIN = 0.02, 0.85
STAGED5 = ["GSE213621", "GSE135251", "GSE162694", "GSE130970", "GSE240729"]
MAIN_ONLY = ["GSE126848", "GSE167523", "PRJNA512027"]
JAPANESE = {"GSE167523", "GSE174478", "GSE193066"}


def lo_of(n, k, F):
    """Smallest minor count the rule calls het at depth n (spec 1.3, integer arithmetic)."""
    return np.maximum(k, -(-F * n // 100))


def bb_logpmf(a, n, mu, rho):
    al, be = mu * (1 - rho) / rho, (1 - mu) * (1 - rho) / rho
    return (gammaln(n + 1) - gammaln(a + 1) - gammaln(n - a + 1)
            + betaln(a + al, n - a + be) - betaln(al, be))


@functools.lru_cache(maxsize=None)
def call_prob(n, rule, mu, rho):
    """P(het call | n) for a count a ~ BB(n, mu, rho): sum over K(n) = {lo..n-lo}."""
    d, k, F = rule
    if n < d:
        return 0.0
    lo = int(lo_of(n, k, F))
    if lo > n - lo:
        return 0.0
    a = np.arange(lo, n - lo + 1)
    return float(np.exp(bb_logpmf(a, n, mu, rho)).sum())


def geuvadis(work):
    hom_parts, het_parts, emp = [], [], []
    for f in sorted(Path(work).glob("chr*.ad_vs_dna.tsv.gz")):
        x = pd.read_csv(f, sep="\t", usecols=["individual", "ref_reads", "alt_reads", "dna"])
        x["n"] = x["ref_reads"] + x["alt_reads"]
        x = x[x["n"] >= 10]
        hom = x[x["dna"] != 1]
        hom_parts.append(pd.DataFrame({"ind": hom["individual"], "n": hom["n"],
                                       "o": np.where(hom["dna"] == 0, hom["alt_reads"], hom["ref_reads"])}))
        het = x[x["dna"] == 1]
        het_parts.append(pd.DataFrame({"n": het["n"], "a": het["ref_reads"]}))
        emp.append(x[["n", "ref_reads", "alt_reads", "dna"]])
    hom = pd.concat(hom_parts, ignore_index=True)
    het = pd.concat(het_parts, ignore_index=True)
    deep = hom[hom["n"] >= 20].groupby("ind")[["o", "n"]].sum()
    e_j = (deep["o"] / deep["n"]).clip(lower=1e-6)
    h = hom[hom["n"] <= 1000].groupby(["ind", "n", "o"]).size().rename("w").reset_index()
    h["e"] = h["ind"].map(e_j)

    def nll_phi(lphi):
        return -(h["w"] * bb_logpmf(h["o"].to_numpy(), h["n"].to_numpy(), h["e"].to_numpy(), np.exp(lphi))).sum()
    phi = float(np.exp(minimize_scalar(nll_phi, bounds=(np.log(1e-6), np.log(0.5)), method="bounded").x))

    g = het[het["n"] <= 1000].groupby(["n", "a"]).size().rename("w").reset_index()

    def nll_het(p):
        mu, rho = 1 / (1 + np.exp(-p[0])), 1 / (1 + np.exp(-p[1]))
        return -(g["w"] * bb_logpmf(g["a"].to_numpy(), g["n"].to_numpy(), mu, rho)).sum()
    r = minimize(nll_het, x0=[0.0, np.log(0.05 / 0.95)], method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-6})
    mu_het, rho_het = float(1 / (1 + np.exp(-r.x[0]))), float(1 / (1 + np.exp(-r.x[1])))

    # empirical vs model FDP and sensitivity per rule, Geuvadis depth <= 1000
    x = pd.concat(emp, ignore_index=True)
    x = x[x["n"] <= 1000]
    minor = np.minimum(x["ref_reads"], x["alt_reads"]).to_numpy()
    n = x["n"].to_numpy()
    is_het = (x["dna"] == 1).to_numpy()
    homrows = hom[hom["n"] <= 1000]
    hom_n = homrows.groupby(["ind", "n"]).size().rename("w").reset_index()
    het_n = het[het["n"] <= 1000].groupby("n").size()
    rows = []
    for rule in GRID:
        d, k, F = rule
        call = (n >= d) & (minor >= k) & (100 * minor >= F * n)
        at = n >= d
        emp_fdp = float((call & ~is_het).sum() / max(call.sum(), 1))
        emp_sens = float((call & is_het).sum() / max((is_het & at).sum(), 1))
        hn = hom_n[hom_n["n"] >= d]
        fp = sum(w * call_prob(int(nn), rule, e_j[ind], phi) for ind, nn, w in hn.itertuples(index=False))
        hh = het_n[het_n.index >= d]
        tp = sum(w * call_prob(int(nn), rule, mu_het, rho_het) for nn, w in hh.items())
        rows.append({"min_dp": d, "min_minor_reads": k, "min_minor_pct": F,
                     "geuvadis_fdp_empirical": emp_fdp, "geuvadis_fdp_model": fp / max(fp + tp, 1e-300),
                     "geuvadis_sens_empirical": emp_sens, "geuvadis_sens_model": tp / max(hh.sum(), 1)})
    fit = {"phi": phi, "mu_het": mu_het, "rho_het": rho_het, "e_j_median": float(e_j.median()),
           "e_j_p90": float(e_j.quantile(0.9)), "hom_rows_n10_1000": int(len(homrows)),
           "het_rows_n10_1000": int(het_n.sum()), "individuals": int(len(e_j))}
    return fit, pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geuvadis-work", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--ad", required=True, help="T4 v2 long table (restricted); only n = REF + ALT is used")
    ap.add_argument("--tag-table", required=True, help="build_tag_table.py tag_table.tsv (spec 1.2 exclusions as flags)")
    ap.add_argument("--null-bias", required=True, help="QC v2 cohort_null_site_ref_bias.tsv (omega_c for the Gate A table)")
    ap.add_argument("--error-rate", required=True, help="QC v2 cohort_error_rate.tsv")
    ap.add_argument("--e-calibrated", help="B-QC cohort_model_inputs.tsv with e_p90_calibrated (spec 9 male chrX rule); "
                                           "when given, it replaces e_p90 as the rule background")
    ap.add_argument("--tag-freq", required=True, help="B-Ancestry tag_2pq_by_superpop.tsv")
    ap.add_argument("--ancestry", required=True, help="B-Ancestry library_ancestry_pcs.tsv (restricted)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    for d, k, F in GRID:  # the integer predicate against direct enumeration
        for n in range(0, 301):
            lo = int(lo_of(n, k, F))
            direct = [m for m in range(n + 1) if n >= d and min(m, n - m) >= k and 100 * min(m, n - m) >= F * n]
            assert direct == ([m for m in range(lo, n - lo + 1)] if n >= d else []), (d, k, F, n)
    fit, geu = geuvadis(a.geuvadis_work)
    geu.to_csv(out / "geuvadis_rule_check.tsv", sep="\t", index=False)
    print(json.dumps(fit, indent=1))

    xw = pd.read_csv(a.crosswalk, sep="\t", dtype=str)
    xw = xw[xw["unit_library"] == "True"]
    staged = (xw["C"] == "0") & xw["S"].notna() & (xw["S"] != "")
    dev = (xw["run_role"] == "development") & ((xw["cohort"].isin(STAGED5) & staged) | xw["cohort"].isin(MAIN_ONLY))
    sea = (xw["run_role"] == "sealed") & staged & xw["cohort"].isin(STAGED5 + ["GSE174478"])
    units = pd.concat([xw[dev].assign(set="development"), xw[sea].assign(set="sealed")])[["run", "cohort", "set", "S"]]
    print(units.groupby(["set", "cohort"]).size().to_string())

    tags = pd.read_csv(a.tag_table, sep="\t", dtype={"pos": int})
    ffpe_drop = sorted(c[len("drop_in_"):] for c in tags.columns if c.startswith("drop_in_"))
    print("transition tags dropped in", ffpe_drop)
    tags = tags[tags["kept_primary"]]

    ad = pd.read_csv(a.ad, sep="\t", usecols=["run", "chrom", "pos", "ref_reads", "alt_reads"])
    ad["n"] = ad["ref_reads"] + ad["alt_reads"]
    ad = ad.drop(columns=["ref_reads", "alt_reads"]).merge(units, on="run")
    rows = ad.merge(tags[["gene_id", "chrom", "pos", "transition"]], on=["chrom", "pos"])
    rows = rows[~(rows["cohort"].isin(ffpe_drop) & rows["transition"])]
    rows = rows.sort_values(["run", "gene_id", "n", "pos"], ascending=[True, True, False, True])
    chosen = rows.drop_duplicates(["run", "gene_id"])
    chosen = chosen[chosen["n"] >= min(d for d, _, _ in GRID)]

    anc = pd.read_csv(a.ancestry, sep="\t", usecols=["run", "ancestry_call", "nearest_superpop", "intermediate"])
    anc["pop"] = np.where(anc["intermediate"].astype(str) == "True", anc["nearest_superpop"], anc["ancestry_call"])
    chosen = chosen.merge(anc[["run", "pop"]], on="run", how="left")
    chosen.loc[chosen["cohort"].isin(JAPANESE) & (chosen["pop"] == "EAS"), "pop"] = "JPT"
    fr = pd.read_csv(a.tag_freq, sep="\t")
    long = fr.melt(id_vars=["chrom", "pos"], value_vars=[c for c in fr.columns if c.startswith("het_ref_alt_")],
                   var_name="pop", value_name="pi")
    long["pop"] = long["pop"].str.replace("het_ref_alt_", "", regex=False)
    chosen = chosen.merge(long, on=["chrom", "pos", "pop"], how="left")
    dropped = chosen[chosen["pi"].isna()].groupby(["cohort", "set"]).size().rename("rows_no_freq")
    chosen = chosen[chosen["pi"].notna()]

    err = pd.read_csv(a.error_rate, sep="\t")
    err["set"] = err["set"].map({"unit_development": "development", "unit_sealed": "sealed"})
    e90 = err.dropna(subset=["set"]).set_index(["cohort", "set"])["e_p90"]
    if a.e_calibrated:
        cal = pd.read_csv(a.e_calibrated, sep="\t")
        cal = cal[cal["rule_set"].isin(["development", "sealed"])].set_index(["cohort", "rule_set"])["e_p90_calibrated"]
        cal.index = cal.index.set_names(["cohort", "set"])
        e90 = cal.reindex(e90.index).fillna(e90)
        print("rule background = calibrated e_p90:\n" + e90.to_string())

    grid_rows, chosen_rules = [], []
    for (c, s), g in chosen.groupby(["cohort", "set"]):
        e = float(e90[(c, s)])
        cnt = g.groupby(["n"])["pi"].agg(["sum", "size"])
        cnt["hom"] = cnt["size"] - cnt["sum"]
        for rule in GRID:
            d, k, F = rule
            sub = cnt[cnt.index >= d]
            p_hom = np.array([call_prob(int(nn), rule, e, fit["phi"]) for nn in sub.index])
            p_het = np.array([call_prob(int(nn), rule, fit["mu_het"], fit["rho_het"]) for nn in sub.index])
            fp, tp = float((sub["hom"] * p_hom).sum()), float((sub["sum"] * p_het).sum())
            grid_rows.append({"cohort": c, "set": s, "min_dp": d, "min_minor_reads": k, "min_minor_pct": F,
                              "e_p90": e, "rows_at_depth": int(sub["size"].sum()), "expected_hets_at_depth": float(sub["sum"].sum()),
                              "fdp_model": fp / max(fp + tp, 1e-300), "sens_model": tp / max(float(sub["sum"].sum()), 1e-300),
                              "expected_true_het_calls": tp})
        gr = pd.DataFrame([r for r in grid_rows if r["cohort"] == c and r["set"] == s])
        ok = gr[(gr["fdp_model"] <= FDP_MAX) & (gr["sens_model"] >= SENS_MIN)]
        if len(ok):
            best = ok.sort_values(["expected_true_het_calls", "min_dp", "min_minor_reads"], ascending=False).iloc[0]
            verdict = "pass"
        else:
            best = gr.sort_values(["fdp_model"]).iloc[0]
            verdict = "fail_untestable"
        chosen_rules.append({"cohort": c, "set": s, "g1prime": verdict, "min_dp": int(best["min_dp"]),
                             "min_minor_reads": int(best["min_minor_reads"]), "min_minor_pct": int(best["min_minor_pct"]),
                             "fdp_model": best["fdp_model"], "sens_model": best["sens_model"],
                             "expected_true_het_calls": best["expected_true_het_calls"], "e_p90": e,
                             "transition_tags_dropped": c in ffpe_drop,
                             "rows_no_freq": int(dropped.get((c, s), 0))})
    pd.DataFrame(grid_rows).to_csv(out / "cohort_rule_grid.tsv", sep="\t", index=False)

    # Gate A aggregate inputs
    rule_of = {(r["cohort"], r["set"]): (r["min_dp"], r["min_minor_reads"], r["min_minor_pct"]) for r in chosen_rules}
    q = [0.1, 0.25, 0.5, 0.75, 0.9]
    gene_rows, stage_rows = [], []
    for (c, s_), g in chosen.groupby(["cohort", "set"]):
        rule = rule_of[(c, s_)]
        g = g[g["n"] >= rule[0]].copy()
        g["p_inc"] = [call_prob(int(nn), rule, fit["mu_het"], fit["rho_het"]) for nn in g["n"]]
        g["e_inc"] = g["pi"] * g["p_inc"]
        for gid, h in g.groupby("gene_id"):
            r = {"cohort": c, "set": s_, "gene_id": gid, "rows_at_depth": len(h), "expected_het_rows": h["pi"].sum(),
                 "expected_included_rows": h["e_inc"].sum(), "distinct_chosen_tags": h["pos"].nunique()}
            r.update({f"n_q{int(x * 100)}": float(np.quantile(h["n"], x)) for x in q})
            gene_rows.append(r)
        for (gid, st), h in g.groupby(["gene_id", "S"], dropna=False):
            stage_rows.append({"cohort": c, "set": s_, "gene_id": gid, "S": st, "rows_at_depth": len(h),
                               "median_n": float(h["n"].median()), "expected_included_rows": h["e_inc"].sum()})
    pd.DataFrame(gene_rows).to_csv(out / "gateA_gene.tsv", sep="\t", index=False)
    pd.DataFrame(stage_rows).to_csv(out / "gateA_gene_by_stage.tsv", sep="\t", index=False)
    units.groupby(["cohort", "set", "S"], dropna=False).size().rename("individuals").reset_index() \
        .to_csv(out / "gateA_stage_counts.tsv", sep="\t", index=False)
    nb = pd.read_csv(a.null_bias, sep="\t")
    nb = nb[(nb["set"] == "unit_all") & (nb["het_rule"] == "depth_20_3_0.10") & (nb["site_set"] == "null_sites")].set_index("cohort")
    est = nb["estimable"].astype(str) == "True"
    dev_est = [c for c in STAGED5 + MAIN_ONLY if c in nb.index and est[c]]
    pooled = float(np.median(nb.loc[dev_est, "ref_frac_median"]))
    inp = []
    for r in chosen_rules:
        c, s_ = r["cohort"], r["set"]
        ok = c in nb.index and bool(est[c])
        rf = float(nb.loc[c, "ref_frac_median"]) if ok else pooled
        e = err[(err["cohort"] == c) & (err["set"] == s_)].iloc[0]
        inp.append({"cohort": c, "set": s_, "min_dp": r["min_dp"], "min_minor_reads": r["min_minor_reads"],
                    "min_minor_pct": r["min_minor_pct"], "e_p10": e["e_p10"], "e_p50": e["e_p50"], "e_p90": e["e_p90"],
                    "eps_p50": e["eps_p50"], "c_p50": e["c_p50"], "omega_c": float(np.log(rf / (1 - rf))),
                    "omega_source": "own" if ok else "pooled_estimable_development", "phi": fit["phi"],
                    "mu_het_geuvadis": fit["mu_het"], "rho_het_geuvadis": fit["rho_het"]})
    pd.DataFrame(inp).to_csv(out / "gateA_cohort_inputs.tsv", sep="\t", index=False)
    cr = pd.DataFrame(chosen_rules)
    cr.to_csv(out / "cohort_rules.tsv", sep="\t", index=False)
    summary = {"geuvadis_fit": fit, "grid": [list(r) for r in GRID], "fdp_max": FDP_MAX, "sens_min": SENS_MIN,
               "ffpe_transition_drop_cohorts": ffpe_drop, "tag_table": a.tag_table,
               "tag_table_kept_rows": int(len(tags))}
    (out / "t2b_summary.json").write_text(json.dumps(summary, indent=1))
    print(cr.to_string(index=False))


if __name__ == "__main__":
    main()
