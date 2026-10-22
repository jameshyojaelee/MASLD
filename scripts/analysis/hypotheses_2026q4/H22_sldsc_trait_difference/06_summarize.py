#!/usr/bin/env python3
"""H22 step 6: tau*, paired-jackknife differences, BH families, seeded-null
calibration, controls and stop rules. Reads ldsc outputs only.

tau*_c = tau_c * sd_c * M / h2g (Gazal 2017), M = total common SNPs (baselineLD
M_5_50 base), sd_c = SD of binary annotation c over common SNPs, h2g = M_5_50 .
coef from the same regression. Delete values per block give the jackknife
variance and the cross-trait covariance (identical blocks for PDFF and cirrhosis).

Usage: 06_summarize.py <OUT> <BASELINE_PREFIX>
"""
import json
import os
import re
import sys

import numpy as np
import pandas as pd

TRAITS = ["pdff", "cirr", "alt", "height"]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def norm_p(z):
    from math import erfc, sqrt
    return np.array([erfc(abs(x) / sqrt(2)) for x in np.atleast_1d(z)])


def jk_cov(a, b):
    B = len(a)
    return (B - 1) / B * np.sum((a - a.mean()) * (b - b.mean()))


class Ctx:
    def __init__(self, out, base_prefix):
        self.out = out
        self.base_names = None
        bm = []
        for c in range(1, 23):
            with open(f"{base_prefix}{c}.l2.M_5_50") as fh:
                bm.append([float(x) for x in fh.readline().split()])
        self.base_m550 = np.sum(bm, axis=0)
        hdr = pd.read_csv(f"{base_prefix}1.l2.ldscore.gz", sep="\t", nrows=0).columns[3:]
        self.base_names = [h[:-2] for h in hdr]
        assert len(self.base_names) == len(self.base_m550)
        self.M_tot = float(self.base_m550[0])
        mm = [pd.read_csv(os.path.join(out, f"ldsc/packed/M.{c}.tsv"), sep="\t") for c in range(1, 23)]
        self.m550 = pd.concat(mm).groupby("annotation", sort=False)["M_5_50"].sum()
        st = pd.concat([pd.read_csv(os.path.join(out, f"ldsc/annot_stats/{c}.tsv"), sep="\t")
                        for c in range(1, 23)])
        g = st.groupby("annotation", sort=False)[["n_common_in", "n_common_total", "n_snp_in"]].sum()
        p = g["n_common_in"] / g["n_common_total"]
        self.sd = np.sqrt(p * (1 - p))
        self.prop_common = p
        self.n_common_total = int(st.groupby("annotation")["n_common_total"].sum().iloc[0])
        self.models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", dtype=str,
                                  keep_default_na=False).set_index("model")
        self.cache = {}

    def run(self, trait, model):
        key = (trait, model)
        if key in self.cache:
            return self.cache[key]
        pre = os.path.join(self.out, "ldsc/h2", trait, model)
        cols = [] if model == "BASE" else self.models.loc[model, "columns"].split(",")
        names = self.base_names + cols
        M = np.concatenate([self.base_m550, [self.m550[c] for c in cols]])
        log = open(pre + ".log").read()
        if os.path.exists(pre + ".results"):
            res = pd.read_csv(pre + ".results", sep="\t")
            coef = res["Coefficient"].to_numpy(float)
            ldsc_se = res["Coefficient_std_error"].to_numpy(float)
            enr = res[["Enrichment", "Enrichment_std_error", "Enrichment_p"]].to_numpy(float)
        else:
            m = re.search(r"Coefficients:(.*?)Coefficient SE:(.*?)Lambda GC", log, re.S)
            coef = np.array(m.group(1).split(), float)
            ldsc_se = np.array(m.group(2).split(), float)
            enr = None
        assert len(coef) == len(names), (trait, model, len(coef), len(names))
        dele = np.loadtxt(pre + ".part_delete", ndmin=2)
        assert dele.shape[1] == len(names)
        se_jk = np.sqrt([jk_cov(dele[:, i], dele[:, i]) for i in range(len(names))])
        h2 = float(M @ coef)
        h2_del = dele @ M
        n_snp = int(re.findall(r"(\d+) SNPs remain", log)[-1])
        chi_rm = re.search(r"Removed (\d+) SNPs with chi\^2 >", log)
        r = dict(names=names, cols=cols, M=M, coef=coef, ldsc_se=ldsc_se, se_jk=se_jk, enr=enr,
                 dele=dele, h2=h2, h2_del=h2_del, n_snp=n_snp, B=dele.shape[0],
                 chisq_removed=int(chi_rm.group(1)) if chi_rm else None)
        tau_star, tau_star_del = {}, {}
        for c in cols:
            i = names.index(c)
            k = self.sd[c] * self.M_tot
            tau_star[c] = coef[i] * k / h2
            tau_star_del[c] = dele[:, i] * k / h2_del
        r["tau_star"], r["tau_star_del"] = tau_star, tau_star_del
        self.cache[key] = r
        return r


def diff_traits(ctx, model, ann):
    a, b = ctx.run("pdff", model), ctx.run("cirr", model)
    assert a["n_snp"] == b["n_snp"] and a["B"] == b["B"], f"blocks differ for {model}"
    da, db = a["tau_star_del"][ann], b["tau_star_del"][ann]
    va, vb, cab = jk_cov(da, da), jk_cov(db, db), jk_cov(da, db)
    d = a["tau_star"][ann] - b["tau_star"][ann]
    se = np.sqrt(va + vb - 2 * cab)
    return dict(tau_star_pdff=a["tau_star"][ann], se_pdff=np.sqrt(va),
                tau_star_cirr=b["tau_star"][ann], se_cirr=np.sqrt(vb),
                cov=cab, corr=cab / np.sqrt(va * vb), D=d, se_D=se, z=d / se,
                se_D_if_independent=np.sqrt(va + vb), n_snp=a["n_snp"])


def diff_gain_lost(ctx, trait, model, g, l):
    r = ctx.run(trait, model)
    dg, dl = r["tau_star_del"][g], r["tau_star_del"][l]
    d = r["tau_star"][g] - r["tau_star"][l]
    se = np.sqrt(jk_cov(dg, dg) + jk_cov(dl, dl) - 2 * jk_cov(dg, dl))
    return dict(tau_star_gain=r["tau_star"][g], tau_star_lost=r["tau_star"][l], D=d, se_D=se, z=d / se)


def main(out, base_prefix):
    ctx = Ctx(out, base_prefix)
    R = os.path.join(out, "results")
    os.makedirs(R, exist_ok=True)
    meta = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t").set_index("annotation")
    models = ctx.models
    lineages = sorted(models["lineage"].dropna().unique())
    summary = {"M_tot_common": ctx.M_tot, "n_common_total_frq": ctx.n_common_total}

    # ---- stop rules and controls ----
    rows = []
    for t in TRAITS:
        r = ctx.run(t, "BASE")
        se = np.sqrt(jk_cov(r["h2_del"], r["h2_del"]))
        rows.append(dict(trait=t, h2=r["h2"], h2_se=se, h2_z=r["h2"] / se, n_snp=r["n_snp"],
                         chisq_removed=r["chisq_removed"]))
    base = pd.DataFrame(rows)
    base.to_csv(os.path.join(R, "baseline_h2.tsv"), sep="\t", index=False)
    stop = bool((base.set_index("trait").loc[["pdff", "cirr"], "h2_z"] < 4).any())

    ctl = []
    for t in TRAITS:
        for L in lineages:
            m = f"CT_{L}"
            if m not in models.index or not os.path.exists(os.path.join(out, "ldsc/h2", t, m + ".log")):
                continue
            r = ctx.run(t, m)
            c = f"SPEC_{L}"
            i = r["names"].index(c)
            ctl.append(dict(trait=t, lineage=L, annotation=c, tau=r["coef"][i], tau_se_ldsc=r["ldsc_se"][i],
                            tau_se_jk=r["se_jk"][i], tau_z=r["coef"][i] / r["ldsc_se"][i],
                            tau_star=r["tau_star"][c],
                            tau_star_se=np.sqrt(jk_cov(r["tau_star_del"][c], r["tau_star_del"][c])),
                            enrichment=r["enr"][i, 0], enrichment_se=r["enr"][i, 1], enrichment_p=r["enr"][i, 2],
                            prop_common_snps=ctx.prop_common[c]))
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(os.path.join(R, "cell_type_spec_tau.tsv"), sep="\t", index=False)
    hep = ctl.set_index(["trait", "lineage"])
    pos_z = float(hep.loc[("alt", "hepatocyte"), "tau_z"])
    neg_z = float(hep.loc[("height", "hepatocyte"), "tau_z"])
    controls = {"positive_alt_hepatocyte_tau_z": pos_z, "positive_pass": pos_z >= 1.96,
                "negative_height_hepatocyte_tau_z": neg_z, "negative_pass": neg_z < 1.96}
    summary.update(controls)
    summary["stop_h2"] = stop
    summary["stop"] = bool(stop or not controls["positive_pass"])
    if os.path.exists(os.path.join(out, "ldsc/STOP.json")):
        summary["stop"] = True
    if summary["stop"]:
        json.dump(summary, open(os.path.join(R, "summary.json"), "w"), indent=2, default=float)
        print("STOP: contrast untestable")
        return

    # ---- checks: jackknife SE reproduces LDSC SE; identical pair SNP sets ----
    chk = []
    for (t, m), r in list(ctx.cache.items()):
        rel = np.max(np.abs(r["se_jk"] - r["ldsc_se"]) / r["ldsc_se"])
        chk.append(dict(trait=t, model=m, max_rel_diff_se=rel, n_snp=r["n_snp"], chisq_removed=r["chisq_removed"]))

    # ---- descriptive tau* for every non-null model ----
    desc = []
    for m, row in models[models["kind"] != "null"].iterrows():
        for t in TRAITS:
            if not os.path.exists(os.path.join(out, "ldsc/h2", t, m + ".log")):
                continue
            r = ctx.run(t, m)
            for c in r["cols"]:
                i = r["names"].index(c)
                se = np.sqrt(jk_cov(r["tau_star_del"][c], r["tau_star_del"][c]))
                desc.append(dict(trait=t, model=m, annotation=c, n_lifted_peaks=int(meta.loc[c, "n_lifted"]),
                                 prop_common_snps=ctx.prop_common[c], tau=r["coef"][i],
                                 tau_z=r["coef"][i] / r["ldsc_se"][i], tau_star=r["tau_star"][c],
                                 tau_star_se=se, tau_star_z=r["tau_star"][c] / se,
                                 enrichment=r["enr"][i, 0], enrichment_p=r["enr"][i, 2]))
    desc = pd.DataFrame(desc)
    desc.to_csv(os.path.join(R, "tau_star_all.tsv"), sep="\t", index=False)

    # ---- families ----
    fam = {}
    rows = []
    for L in lineages:
        if f"STEA_{L}" in models.index:
            for c in models.loc[f"STEA_{L}", "columns"].split(",")[1:]:
                rows.append(dict(lineage=L, annotation=c, direction="gain" if "_GAIN_" in c else "lost",
                                 **diff_traits(ctx, f"STEA_{L}", c)))
    fam["F1"] = pd.DataFrame(rows)
    rows = []
    for L in lineages:
        if f"CT_{L}" in models.index:
            rows.append(dict(lineage=L, annotation=f"SPEC_{L}", **diff_traits(ctx, f"CT_{L}", f"SPEC_{L}")))
    fam["F2"] = pd.DataFrame(rows)
    rows = []
    for L in lineages:
        m = f"STEA_{L}"
        if m in models.index:
            cols = models.loc[m, "columns"].split(",")
            g, l = f"STEA_GAIN_{L}", f"STEA_LOST_{L}"
            if g in cols and l in cols:
                for t in ("pdff", "cirr"):
                    rows.append(dict(lineage=L, trait=t, **diff_gain_lost(ctx, t, m, g, l)))
    fam["F3"] = pd.DataFrame(rows)
    rows = []
    for L in lineages:
        if f"FIB_{L}" in models.index:
            for c in models.loc[f"FIB_{L}", "columns"].split(",")[1:]:
                rows.append(dict(lineage=L, annotation=c, direction="gain" if "_GAIN_" in c else "lost",
                                 **diff_traits(ctx, f"FIB_{L}", c)))
    fam["F4"] = pd.DataFrame(rows)
    for k, df in fam.items():
        df["p"] = norm_p(df["z"].values)
        df["q_BH"] = bh(df["p"].values)
        df["family_size"] = len(df)

    # ---- seeded null ----
    nrows = []
    for m, row in models[models["kind"] == "null"].iterrows():
        L, k = row["lineage"], int(row["draw"])
        cols = row["columns"].split(",")
        g, l = f"NULL_GAIN_{L}_{k}", f"NULL_LOST_{L}_{k}"
        for c, dirn in ((g, "gain"), (l, "lost")):
            if c in cols:
                d = diff_traits(ctx, m, c)
                nrows.append(dict(lineage=L, draw=k, stat="F1", direction=dirn, trait="pdff-cirr", z=d["z"], D=d["D"]))
        if g in cols and l in cols:
            for t in ("pdff", "cirr"):
                d = diff_gain_lost(ctx, t, m, g, l)
                nrows.append(dict(lineage=L, draw=k, stat="F3", direction="gain-lost", trait=t, z=d["z"], D=d["D"]))
        for t in ("pdff", "cirr"):
            r = ctx.run(t, m)
            rel = np.max(np.abs(r["se_jk"] - r["ldsc_se"]) / r["ldsc_se"])
            chk.append(dict(trait=t, model=m, max_rel_diff_se=rel, n_snp=r["n_snp"], chisq_removed=r["chisq_removed"]))
        for t in ("pdff", "cirr"):
            ctx.cache.pop((t, m), None)
    nulls = pd.DataFrame(nrows)
    nulls.to_csv(os.path.join(R, "null_draws.tsv"), sep="\t", index=False)
    pd.DataFrame(chk).to_csv(os.path.join(R, "checks.tsv"), sep="\t", index=False)

    f1n = nulls[nulls["stat"] == "F1"]["z"].abs()
    f3n = nulls[nulls["stat"] == "F3"]["z"].abs()
    calib = {"F1_null_n": int(len(f1n)), "F1_null_frac_abs_z_ge_1.96": float((f1n >= 1.96).mean()),
             "F3_null_n": int(len(f3n)), "F3_null_frac_abs_z_ge_1.96": float((f3n >= 1.96).mean())}
    summary.update(calib)
    anti = calib["F1_null_frac_abs_z_ge_1.96"] > 0.10
    summary["F1_jackknife_anticonservative"] = bool(anti)

    def emp(df, stat, keyfun):
        out_p, out_n = [], []
        for _, r in df.iterrows():
            sel = nulls[(nulls["stat"] == stat) & (nulls["lineage"] == r["lineage"]) & keyfun(r)]
            n = len(sel)
            out_p.append((1 + (sel["z"].abs() >= abs(r["z"])).sum()) / (n + 1) if n else np.nan)
            out_n.append(n)
        df["p_empirical"] = out_p
        df["n_null"] = out_n

    emp(fam["F1"], "F1", lambda r: nulls["direction"] == r["direction"])
    emp(fam["F3"], "F3", lambda r: nulls["trait"] == r["trait"])
    if anti:
        fam["F1"]["q_BH_empirical"] = bh(fam["F1"]["p_empirical"].values)
        fam["F1"]["supported"] = fam["F1"]["q_BH_empirical"] < 0.05
    else:
        fam["F1"]["supported"] = (fam["F1"]["q_BH"] < 0.05) & (fam["F1"]["p_empirical"] < 0.05)
    fam["F3"]["supported"] = (fam["F3"]["q_BH"] < 0.05) & (fam["F3"]["p_empirical"] < 0.05)
    fam["F2"]["jackknife_only_q_lt_0.05"] = fam["F2"]["q_BH"] < 0.05
    fam["F4"]["jackknife_only_q_lt_0.05"] = fam["F4"]["q_BH"] < 0.05
    for k, df in fam.items():
        df.to_csv(os.path.join(R, f"{k}.tsv"), sep="\t", index=False)
        summary[f"{k}_family_size"] = int(len(df))
        summary[f"{k}_n_q_lt_0.05"] = int((df["q_BH"] < 0.05).sum())
        if "supported" in df:
            summary[f"{k}_n_supported"] = int(df["supported"].sum())

    # ---- overlap check ----
    rg_log = os.path.join(out, "ldsc/rg/pdff_cirr.log")
    if os.path.exists(rg_log):
        lines = open(rg_log).read().splitlines()
        i = [j for j, x in enumerate(lines) if x.startswith("Summary of Genetic Correlation Results")][0]
        hdr, val = lines[i + 1].split(), lines[i + 2].split()
        summary["rg_pdff_cirr"] = dict(zip(hdr, val))
    summary["pair_report"] = json.load(open(os.path.join(out, "ldsc/sumstats/pair_report.json")))
    summary["max_rel_diff_jackknife_vs_ldsc_se"] = float(pd.DataFrame(chk)["max_rel_diff_se"].max())
    json.dump(summary, open(os.path.join(R, "summary.json"), "w"), indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
