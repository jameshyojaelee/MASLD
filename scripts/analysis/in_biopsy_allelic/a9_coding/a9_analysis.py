#!/usr/bin/env python3
"""A9: named coding / risk variants genotyped from RNA reads (PRESPEC_A section 10).

Development participants only (one library per genetic individual, sealed excluded;
same selection and stage coding as T5a). Genotype from the reads at the site with the
T2 rule: het by the rule; otherwise hom for the majority allele if depth >= the
rule's minimum; otherwise missing. Call rate is reported by stage because depth
depends on expression, which can change with stage.

1. Positive control: PNPLA3 rs738409 G (148M) dosage higher at S = 2 than S = 0.
   Logistic regression of S2-vs-S0 on dosage + cohort indicators, one-sided Wald p.
   The same model for the other five variants (secondary; BH within the six).
2. Program associations: for variants with >= 20 minor-allele carriers among staged
   participants with program scores, OLS of each of the 117 cohort-centred program
   scores (Model B inputs, ALT-corrected counts) on dosage + S indicators + cohort;
   BH within the 117 programs per variant. Null: dosage permuted within cohort x S,
   1,000 permutations (seed 20260923), for the number of programs at q < 0.05.
   Reported as cross-sectional associations within stage.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tags"))
from t5a_allelic_table import bh, choose_libraries, het_call, stage_coding  # noqa: E402

SEED, N_PERM, MIN_CARRIERS = 20260923, 1000, 20


def read_ad(ad_dir, variants):
    rows = []
    for f in sorted(Path(ad_dir).glob("*.ad.tsv.gz")):
        d = pd.read_csv(f, sep="\t", dtype=str)
        libs = [col.split("]", 1)[1].rsplit("/", 1)[-1].replace(".bam:AD", "") for col in d.columns[4:]]
        for rec in d.itertuples(index=False):
            key = (rec[0], int(rec[1]))
            if key not in variants:
                continue
            ref, alt = variants[key]
            if rec[2] != ref:
                raise SystemExit(f"REF mismatch at {key}")
            alts = rec[3].split(",")
            j = 1 + alts.index(alt) if alt in alts else None
            for lib, ad in zip(libs, rec[4:]):
                v = [int(x) if x != "." else 0 for x in ad.split(",")]
                rows.append((key[0], key[1], lib, v[0], v[j] if j is not None else 0))
    return pd.DataFrame(rows, columns=["chrom", "pos", "run", "ref_reads", "alt_reads"])


def genotype(ref, alt, rule):
    dp = ref + alt
    het = het_call(ref, alt, rule)
    g = np.where(het, 1.0, np.where(alt > ref, 2.0, 0.0))
    return np.where(dp >= rule["min_dp"], g, np.nan)


def ols_t(y, X):
    """t of column 0 for each column of y (n x m)."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ y
    R = y - X @ B
    s2 = (R ** 2).sum(0) / (len(X) - X.shape[1])
    return B[0] / np.sqrt(s2 * XtX_inv[0, 0]), B[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ad-dir", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--t2", required=True)
    ap.add_argument("--identity-dir", required=True)
    ap.add_argument("--seal-dir", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--gse130970-controls", required=True)
    ap.add_argument("--gse193066-placement", required=True)
    ap.add_argument("--programs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rule = json.loads(Path(a.t2).read_text())["chosen_rule"]
    var = pd.read_csv(a.variants, sep="\t")
    vmap = {(c, p): (r, x) for c, p, r, x in zip(var["chrom"], var["pos"], var["ref"], var["alt"])}

    xwalk = pd.read_csv(Path(a.identity_dir) / "sample_to_individual.tsv", sep="\t")
    mixed = set(pd.read_csv(Path(a.seal_dir) / "excluded_possibly_mixed_libraries.tsv", sep="\t")["run"])
    seal = pd.read_csv(Path(a.seal_dir) / "sealed_individuals.tsv", sep="\t")
    lib = choose_libraries(xwalk, mixed, seal, pd.read_csv(a.gse193066_placement, sep="\t"))
    c130 = pd.read_csv(a.gse130970_controls, sep="\t")
    stage = stage_coding(pd.read_csv(a.metadata, keep_default_na=False, na_values=[""]),
                         set(c130.loc[c130["source_control_status"] == "Control", "run"]))
    lib = lib.merge(stage, on="run", how="left")
    lib["stage_eligible"] = lib["S"].notna() & (lib["C"] == 0) & lib["first_biopsy"]

    ad = read_ad(a.ad_dir, vmap)
    ad = ad[ad["run"].isin(set(lib["run"]))]
    ad["dosage"] = genotype(ad["ref_reads"].to_numpy(), ad["alt_reads"].to_numpy(), rule)
    ad = ad.merge(var[["variant", "gene", "chrom", "pos"]], on=["chrom", "pos"])
    G = ad.pivot_table(index="run", columns="variant", values="dosage", aggfunc="first")
    lib = lib.merge(G, left_on="run", right_index=True, how="left")

    summary, pc = {"rule": rule, "participants": int(len(lib))}, []
    st = lib[lib["stage_eligible"]]
    for v in var["variant"]:
        d = lib[v] if v in lib else pd.Series(np.nan, index=lib.index)
        callrate = {f"S{int(s)}": float(st.loc[st["S"] == s, v].notna().mean()) if v in st else 0.0 for s in (0, 1, 2)}
        called = lib[lib[v].notna()] if v in lib else lib.iloc[0:0]
        af = float(called[v].mean() / 2) if len(called) else np.nan
        # HWE chi-square on called development participants
        n0, n1, n2 = [(called[v] == k).sum() for k in (0, 1, 2)]
        n = n0 + n1 + n2
        p_ = (n1 + 2 * n2) / (2 * n) if n else np.nan
        exp = np.array([(1 - p_) ** 2, 2 * p_ * (1 - p_), p_ ** 2]) * n if n else np.zeros(3)
        hwe = float(stats.chi2.sf(np.sum((np.array([n0, n1, n2]) - exp) ** 2 / np.where(exp > 0, exp, 1)), 1)) if n else np.nan
        row = {"variant": v, "called": int(n), "alt_allele_freq": af, "hom_ref": int(n0), "het": int(n1), "hom_alt": int(n2),
               "hwe_chi2_p": hwe, **{f"call_rate_{k}": x for k, x in callrate.items()}}
        s02 = st[st["S"].isin([0, 2])].dropna(subset=[v]) if v in st else st.iloc[0:0]
        if len(s02) > 30 and s02[v].std() > 0:
            X = pd.get_dummies(s02["cohort"], drop_first=True).astype(float)
            X.insert(0, "dosage", s02[v].to_numpy())
            X = sm.add_constant(X)
            fit = sm.Logit((s02["S"] == 2).astype(float).to_numpy(), X).fit(disp=0)
            z = fit.params["dosage"] / fit.bse["dosage"]
            row.update({"n_S0": int((s02["S"] == 0).sum()), "n_S2": int((s02["S"] == 2).sum()),
                        "log_or_S2_vs_S0_per_alt": float(fit.params["dosage"]), "se": float(fit.bse["dosage"]),
                        "p_one_sided_alt_higher_at_S2": float(stats.norm.sf(z)), "p_two_sided": float(2 * stats.norm.sf(abs(z)))})
        pc.append(row)
    pc = pd.DataFrame(pc)
    pc["q_two_sided_within_six"] = bh(pc.get("p_two_sided", pd.Series(np.nan, index=pc.index)))
    pc.to_csv(out / "coding_variants_stage.tsv", sep="\t", index=False)
    pn = pc[pc["variant"] == "rs738409"].iloc[0].to_dict()
    summary["pnpla3_positive_control"] = {k: pn.get(k) for k in ("n_S0", "n_S2", "log_or_S2_vs_S0_per_alt", "p_one_sided_alt_higher_at_S2")}
    summary["pnpla3_positive_control_pass"] = bool(pn.get("p_one_sided_alt_higher_at_S2", 1) < 0.05)

    # program associations within stage
    prog = pd.read_csv(a.programs, sep="\t")
    pcols = [c for c in prog.columns if c.startswith("hotspot_")]
    sp = st.merge(prog[["sample"] + pcols], left_on="run", right_on="sample")
    rng = np.random.default_rng(SEED)
    res = []
    for v in var["variant"]:
        if v not in sp:
            continue
        d = sp.dropna(subset=[v] + pcols)
        carriers = int((d[v] > 0).sum()) if d[v].mean() < 1 else int((d[v] < 2).sum())  # carriers of the minor allele
        if carriers < MIN_CARRIERS or d[v].std() == 0:
            continue
        cov = np.column_stack([(d["S"] == 1).astype(float), (d["S"] == 2).astype(float),
                               pd.get_dummies(d["cohort"]).to_numpy(float)])
        X = np.column_stack([d[v].to_numpy(float), cov])
        Y = d[pcols].to_numpy(float)
        t, b = ols_t(Y, X)
        p = 2 * stats.t.sf(np.abs(t), len(d) - X.shape[1])
        q = bh(p)
        strata = d.groupby(["cohort", "S"]).indices
        nulls = []
        for _ in range(N_PERM):
            g = d[v].to_numpy(float).copy()
            for ii in strata.values():
                g[ii] = g[ii][rng.permutation(len(ii))]
            tp, _ = ols_t(Y, np.column_stack([g, cov]))
            nulls.append(int((bh(2 * stats.t.sf(np.abs(tp), len(d) - X.shape[1])) < 0.05).sum()))
        n_hit = int((q < 0.05).sum())
        res.append(pd.DataFrame({"variant": v, "program": pcols, "beta_per_alt": b, "t": t, "p": p, "q": q}))
        summary[f"programs_{v}"] = {"n": int(len(d)), "minor_carriers": carriers, "programs_q_lt_0.05": n_hit,
                                    "null_hits_95th": float(np.percentile(nulls, 95)),
                                    "p_count_vs_null": float((1 + np.sum(np.array(nulls) >= n_hit)) / (N_PERM + 1))}
    if res:
        pd.concat(res).sort_values("p").to_csv(out / "coding_variant_program_associations.tsv", sep="\t", index=False)
    (out / "a9_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
