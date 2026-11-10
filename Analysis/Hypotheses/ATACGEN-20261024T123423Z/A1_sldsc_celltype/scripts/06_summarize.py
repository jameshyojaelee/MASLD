#!/usr/bin/env python3
"""A1 step 6: tau*, standardized coefficients, paired-jackknife differences, BH
families, seeded-null calibration and controls (adapted from H22 06). Reads ldsc
outputs only.

tau*_c = tau_c * sd_c * M / h2g (Gazal 2017); std coef_c = tau_c * M / h2g.
M = baselineLD M_5_50 of the base annotation; sd_c over common SNPs (binary or
continuous); h2g = sum_k tau_k M_5_50,k of the same regression. Block delete
values give the jackknife (co)variances, B = 200.

Usage: 06_summarize.py <OUT> <BASELINE_PREFIX> <TRAITS_TSV> <H22_DIR>
"""
import json
import os
import re
import sys
from math import erfc, sqrt

import numpy as np
import pandas as pd

V3_LINEAGES = ["cholangiocyte", "hepatocyte", "macrophage", "stellate", "t_nk"]
DA_LINEAGES = ["hepatocyte", "stellate"]
H22_LINEAGES = ["b_cell", "cholangiocyte", "endothelial_cell", "fibroblast", "hepatocyte", "macrophage", "t_cell"]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def p_two(z):
    return np.array([erfc(abs(x) / sqrt(2)) for x in np.atleast_1d(z)])


def p_one_pos(z):
    return np.array([0.5 * erfc(x / sqrt(2)) for x in np.atleast_1d(z)])


def jk_cov(a, b):
    B = len(a)
    return (B - 1) / B * np.sum((a - a.mean()) * (b - b.mean()))


class Ctx:
    def __init__(self, out, base_prefix, h22):
        self.out = out
        bm = []
        for c in range(1, 23):
            with open(f"{base_prefix}{c}.l2.M_5_50") as fh:
                bm.append([float(x) for x in fh.readline().split()])
        self.base_m550 = np.sum(bm, axis=0)
        hdr = pd.read_csv(f"{base_prefix}1.l2.ldscore.gz", sep="\t", nrows=0).columns[3:]
        self.base_names = [h[:-2] for h in hdr]
        self.M_tot = float(self.base_m550[0])
        # A1 annotations
        mm = pd.concat([pd.read_csv(os.path.join(out, f"ldsc/packed/M.{c}.tsv"), sep="\t") for c in range(1, 23)])
        self.m550 = mm.groupby("annotation", sort=False)["M_5_50"].sum().to_dict()
        st = pd.concat([pd.read_csv(os.path.join(out, f"ldsc/annot_stats/{c}.tsv"), sep="\t") for c in range(1, 23)])
        g = st.groupby("annotation", sort=False)[["n_common_nonzero", "sum_common", "sumsq_common",
                                                   "n_common_total", "n_snp_nonzero"]].sum()
        n = g["n_common_total"]
        mean = g["sum_common"] / n
        self.sd = np.sqrt(g["sumsq_common"] / n - mean ** 2).to_dict()
        self.prop_common = (g["n_common_nonzero"] / n).to_dict()
        self.n_common_in = g["n_common_nonzero"].to_dict()
        self.n_common_total = int(n.iloc[0])
        # H22 GSE296875 annotations (read-only)
        mm = pd.concat([pd.read_csv(os.path.join(h22, f"ldsc/packed/M.{c}.tsv"), sep="\t") for c in range(1, 23)])
        h22m = mm.groupby("annotation", sort=False)["M_5_50"].sum()
        st = pd.concat([pd.read_csv(os.path.join(h22, f"ldsc/annot_stats/{c}.tsv"), sep="\t") for c in range(1, 23)])
        g = st.groupby("annotation", sort=False)[["n_common_in", "n_common_total"]].sum()
        p = g["n_common_in"] / g["n_common_total"]
        for a in ["UNION"] + [f"SPEC_{L}" for L in H22_LINEAGES]:
            k = "H22:" + a
            self.m550[k] = float(h22m[a])
            self.sd[k] = float(np.sqrt(p[a] * (1 - p[a])))
            self.prop_common[k] = float(p[a])
            self.n_common_in[k] = int(g.loc[a, "n_common_in"])
        self.models = pd.read_csv(os.path.join(out, "ldsc/models.tsv"), sep="\t", dtype=str,
                                  keep_default_na=False).set_index("model")
        self.cache = {}

    def cols(self, model):
        if model == "BASE":
            return []
        if model.startswith("CT296875_"):
            return ["H22:UNION", "H22:SPEC_" + model[len("CT296875_"):]]
        return self.models.loc[model, "columns"].split(",")

    def run(self, trait, model):
        key = (trait, model)
        if key in self.cache:
            return self.cache[key]
        pre = os.path.join(self.out, "ldsc/h2", trait, model)
        if not os.path.exists(pre + ".log"):
            return None
        cols = self.cols(model)
        names = self.base_names + cols
        M = np.concatenate([self.base_m550, [self.m550[c] for c in cols]])
        log = open(pre + ".log").read()
        enr = None
        if os.path.exists(pre + ".results"):
            res = pd.read_csv(pre + ".results", sep="\t")
            coef = res["Coefficient"].to_numpy(float)
            ldsc_se = res["Coefficient_std_error"].to_numpy(float)
            enr = res[["Prop._SNPs", "Prop._h2", "Enrichment", "Enrichment_std_error", "Enrichment_p"]].to_numpy(float)
        else:
            m = re.search(r"Coefficients:(.*?)Coefficient SE:(.*?)Lambda GC", log, re.S)
            coef = np.array(m.group(1).split(), float)
            ldsc_se = np.array(m.group(2).split(), float)
        assert len(coef) == len(names), (trait, model, len(coef), len(names))
        dele = np.loadtxt(pre + ".part_delete", ndmin=2)
        assert dele.shape[1] == len(names)
        se_jk = np.sqrt([jk_cov(dele[:, i], dele[:, i]) for i in range(len(names))])
        h2 = float(M @ coef)
        h2_del = dele @ M
        n_snp = int(re.findall(r"(\d+) SNPs remain", log)[-1])
        chi_rm = re.search(r"Removed (\d+) SNPs with chi\^2 >", log)
        r = dict(names=names, cols=cols, coef=coef, ldsc_se=ldsc_se, se_jk=se_jk, enr=enr, dele=dele,
                 h2=h2, h2_del=h2_del, n_snp=n_snp, B=dele.shape[0],
                 chisq_removed=int(chi_rm.group(1)) if chi_rm else 0)
        r["ts"], r["ts_del"], r["sc"], r["sc_del"] = {}, {}, {}, {}
        for c in cols:
            i = names.index(c)
            r["ts"][c] = coef[i] * self.sd[c] * self.M_tot / h2
            r["ts_del"][c] = dele[:, i] * self.sd[c] * self.M_tot / h2_del
            r["sc"][c] = coef[i] * self.M_tot / h2
            r["sc_del"][c] = dele[:, i] * self.M_tot / h2_del
        self.cache[key] = r
        return r


def annot_row(ctx, r, c):
    i = r["names"].index(c)
    se_ts = np.sqrt(jk_cov(r["ts_del"][c], r["ts_del"][c]))
    d = dict(annotation=c, tau=r["coef"][i], tau_se=r["ldsc_se"][i], tau_z=r["coef"][i] / r["ldsc_se"][i],
             tau_star=r["ts"][c], tau_star_se=se_ts, tau_star_z=r["ts"][c] / se_ts,
             prop_common_snps=ctx.prop_common[c], n_common_snps=ctx.n_common_in[c],
             n_regression_snps=r["n_snp"], n_blocks=r["B"], h2g=r["h2"])
    if r["enr"] is not None:
        d.update(prop_h2=r["enr"][i, 1], enrichment=r["enr"][i, 2], enrichment_se=r["enr"][i, 3],
                 enrichment_p=r["enr"][i, 4])
    return d


def diff(r, a, b, which="ts"):
    da, db = r[which + "_del"][a], r[which + "_del"][b]
    d = r[which][a] - r[which][b]
    se = np.sqrt(jk_cov(da, da) + jk_cov(db, db) - 2 * jk_cov(da, db))
    return d, se


def main(out, base_prefix, traits_path, h22):
    ctx = Ctx(out, base_prefix, h22)
    R = os.path.join(out, "results")
    os.makedirs(R, exist_ok=True)
    tr = pd.read_csv(traits_path, sep="\t")
    testab = pd.read_csv(os.path.join(out, "ldsc/phase1_testability.tsv"), sep="\t")
    tr = tr.merge(testab, on="trait_id")
    ann = pd.read_csv(os.path.join(out, "peaks/annotations.tsv"), sep="\t", keep_default_na=False).set_index("annotation")
    summary = {"M_tot_common": ctx.M_tot, "n_common_total_frq": ctx.n_common_total}

    # ---- baseline h2 ----
    rows = []
    for _, t in tr.iterrows():
        r = ctx.run(t["trait_id"], "BASE")
        se = np.sqrt(jk_cov(r["h2_del"], r["h2_del"]))
        rows.append(dict(trait_id=t["trait_id"], trait=t["trait"], study=t["study"], role=t["role"],
                         primary=t["primary"], ancestry=t["ancestry"], h2=r["h2"], h2_se_jk=se,
                         h2_z_ldsc_log=t["h2_z"], testable=t["testable"], n_regression_snps=r["n_snp"],
                         chisq_removed=r["chisq_removed"], n_blocks=r["B"]))
    pd.DataFrame(rows).to_csv(os.path.join(R, "baseline_h2.tsv"), sep="\t", index=False)
    ok = tr[tr["testable"]]
    test_ids = list(ok.loc[ok["role"] == "test", "trait_id"])
    ctrl_ids = list(ok.loc[ok["role"] == "control", "trait_id"])
    desc_ids = list(ok.loc[ok["role"] == "descriptive", "trait_id"])
    meta = tr.set_index("trait_id")

    def tag(t):
        return dict(trait_id=t, trait=meta.loc[t, "trait"], study=meta.loc[t, "study"],
                    role=meta.loc[t, "role"], primary=bool(meta.loc[t, "primary"]))

    # ---- Q1 families (one-sided, tau > 0) ----
    def q1(ids, prefix, lineages, colname):
        rows = []
        for t in ids:
            for L in lineages:
                m = f"{prefix}_{L}"
                r = ctx.run(t, m)
                if r is None:
                    continue
                c = colname(L)
                d = annot_row(ctx, r, c)
                d["annotation"] = c.replace("H22:", "")
                u = annot_row(ctx, r, r["cols"][0])
                rows.append(dict(**tag(t), lineage=L, model=m, **d,
                                 union_tau_z=u["tau_z"], union_tau_star=u["tau_star"]))
        df = pd.DataFrame(rows)
        if len(df):
            df["p_one_sided"] = p_one_pos(df["tau_z"].values)
            df["p_two_sided"] = p_two(df["tau_z"].values)
        return df

    fam = {}
    fam["F1"] = q1(test_ids, "ONLY", V3_LINEAGES, lambda L: f"ONLY_{L}")
    fam["F1_ALL"] = q1(test_ids, "ALL", V3_LINEAGES, lambda L: f"ALL_{L}")
    fam["F1_REP"] = q1(test_ids, "CT296875", H22_LINEAGES, lambda L: f"H22:SPEC_{L}")
    fam["F1_CTRL"] = q1(ctrl_ids, "ONLY", V3_LINEAGES, lambda L: f"ONLY_{L}")
    for k in ("F1", "F1_ALL", "F1_REP", "F1_CTRL"):
        df = fam[k]
        df["q_BH"] = bh(df["p_one_sided"].values)
        df["family_size"] = len(df)
    desc = {}
    desc["Q1_controls_ALL"] = q1(ctrl_ids, "ALL", V3_LINEAGES, lambda L: f"ALL_{L}")
    desc["Q1_controls_REP"] = q1(ctrl_ids, "CT296875", H22_LINEAGES, lambda L: f"H22:SPEC_{L}")
    desc["Q1_TG_ONLY"] = q1(desc_ids, "ONLY", V3_LINEAGES, lambda L: f"ONLY_{L}")
    desc["Q1_TG_ALL"] = q1(desc_ids, "ALL", V3_LINEAGES, lambda L: f"ALL_{L}")
    desc["Q1_TG_REP"] = q1(desc_ids, "CT296875", H22_LINEAGES, lambda L: f"H22:SPEC_{L}")

    # ---- Q3: promoter vs distal ----
    def q3(ids):
        rows = []
        for t in ids:
            for L in V3_LINEAGES:
                r = ctx.run(t, f"PD_{L}")
                if r is None:
                    continue
                pr, di = f"PROM_{L}", f"DIST_{L}"
                a, b = annot_row(ctx, r, pr), annot_row(ctx, r, di)
                D, se = diff(r, di, pr, "sc")
                se_p = np.sqrt(jk_cov(r["sc_del"][pr], r["sc_del"][pr]))
                se_d = np.sqrt(jk_cov(r["sc_del"][di], r["sc_del"][di]))
                rows.append(dict(**tag(t), lineage=L, std_coef_prom=r["sc"][pr], std_coef_prom_se=se_p,
                                 std_coef_dist=r["sc"][di], std_coef_dist_se=se_d, D_dist_minus_prom=D, se_D=se,
                                 z=D / se, tau_star_prom=a["tau_star"], tau_star_prom_se=a["tau_star_se"],
                                 tau_star_dist=b["tau_star"], tau_star_dist_se=b["tau_star_se"],
                                 tau_z_prom=a["tau_z"], tau_z_dist=b["tau_z"],
                                 enrichment_prom=a.get("enrichment"), enrichment_prom_p=a.get("enrichment_p"),
                                 enrichment_dist=b.get("enrichment"), enrichment_dist_p=b.get("enrichment_p"),
                                 n_common_prom=a["n_common_snps"], n_common_dist=b["n_common_snps"],
                                 n_regression_snps=r["n_snp"]))
        df = pd.DataFrame(rows)
        if len(df):
            df["p"] = p_two(df["z"].values)
        return df

    fam["F3"] = q3(test_ids)
    fam["F3"]["q_BH"] = bh(fam["F3"]["p"].values)
    fam["F3"]["family_size"] = len(fam["F3"])
    desc["Q3_TG"] = q3(desc_ids)

    # ---- Q2 ----
    def q2_obs(ids):
        cont, vs, dd, sign, top = [], [], [], [], []
        for t in ids:
            for L in DA_LINEAGES:
                r = ctx.run(t, f"CONT_{L}")
                if r is not None:
                    d = annot_row(ctx, r, f"RESP_{L}")
                    cont.append(dict(**tag(t), lineage=L, statistic="tau_star_RESP", **d, z=d["tau_star_z"]))
                r = ctx.run(t, f"SIGN_{L}")
                if r is not None:
                    d = annot_row(ctx, r, f"SIGN_{L}")
                    sign.append(dict(**tag(t), lineage=L, statistic="tau_star_SIGN", **d, z=d["tau_star_z"]))
                r = ctx.run(t, f"TERT_{L}")
                if r is not None:
                    for s in ("GAIN", "LOST"):
                        d = annot_row(ctx, r, f"{s}_{L}")
                        vs.append(dict(**tag(t), lineage=L, statistic=f"tau_star_{s}_vs_stable", **d,
                                       z=d["tau_star_z"]))
                    D, se = diff(r, f"GAIN_{L}", f"LOST_{L}")
                    dd.append(dict(**tag(t), lineage=L, statistic="D_gain_minus_lost",
                                   tau_star_gain=r["ts"][f"GAIN_{L}"], tau_star_lost=r["ts"][f"LOST_{L}"],
                                   D=D, se_D=se, z=D / se, n_regression_snps=r["n_snp"]))
                r = ctx.run(t, f"TOP10_{L}")
                if r is not None:
                    for s in ("GAIN10", "LOST10"):
                        d = annot_row(ctx, r, f"{s}_{L}")
                        top.append(dict(**tag(t), lineage=L, statistic=f"tau_star_{s}_vs_rest", **d,
                                        z=d["tau_star_z"]))
        return [pd.DataFrame(x) for x in (cont, vs, dd, sign, top)]

    fam["F2"], fam["F2_VS"], fam["F2_D"], desc["Q2_SIGN"], desc["Q2_TOP10"] = q2_obs(test_ids)
    tgq2 = q2_obs(desc_ids)
    desc["Q2_TG_CONT"], desc["Q2_TG_VS"], desc["Q2_TG_D"] = tgq2[0], tgq2[1], tgq2[2]

    # nulls
    nrows, chk = [], []
    for m, row in ctx.models[ctx.models["kind"].isin(["nullcont", "nulltert"])].iterrows():
        L, k = row["lineage"], int(row["draw"])
        for t in test_ids:
            r = ctx.run(t, m)
            if r is None:
                continue
            if row["kind"] == "nullcont":
                c = f"NULLRESP_{L}_{k}"
                se = np.sqrt(jk_cov(r["ts_del"][c], r["ts_del"][c]))
                nrows.append(dict(trait_id=t, lineage=L, draw=k, statistic="tau_star_RESP", z=r["ts"][c] / se))
            else:
                for s, nm in (("GAIN", "NULLGAIN"), ("LOST", "NULLLOST")):
                    c = f"{nm}_{L}_{k}"
                    se = np.sqrt(jk_cov(r["ts_del"][c], r["ts_del"][c]))
                    nrows.append(dict(trait_id=t, lineage=L, draw=k, statistic=f"tau_star_{s}_vs_stable",
                                      z=r["ts"][c] / se))
                D, se = diff(r, f"NULLGAIN_{L}_{k}", f"NULLLOST_{L}_{k}")
                nrows.append(dict(trait_id=t, lineage=L, draw=k, statistic="D_gain_minus_lost", z=D / se))
            rel = np.max(np.abs(r["se_jk"] - r["ldsc_se"]) / r["ldsc_se"])
            chk.append(dict(trait_id=t, model=m, max_rel_diff_se=rel, n_snp=r["n_snp"]))
            ctx.cache.pop((t, m), None)
    nulls = pd.DataFrame(nrows)
    nulls.to_csv(os.path.join(R, "q2_null_draws.tsv"), sep="\t", index=False)
    calib = {}
    for s in ("tau_star_RESP", "tau_star_GAIN_vs_stable", "tau_star_LOST_vs_stable", "D_gain_minus_lost"):
        z = nulls.loc[nulls["statistic"] == s, "z"].abs()
        calib[s] = dict(n=int(len(z)), frac_abs_z_ge_1_96=float((z >= 1.96).mean()) if len(z) else None)
    summary["q2_null_calibration"] = calib
    anti = calib["tau_star_RESP"]["frac_abs_z_ge_1_96"] is not None and calib["tau_star_RESP"]["frac_abs_z_ge_1_96"] > 0.10
    summary["F2_jackknife_anticonservative"] = bool(anti)

    def emp(df):
        ps, ns = [], []
        for _, r in df.iterrows():
            sel = nulls[(nulls["trait_id"] == r["trait_id"]) & (nulls["lineage"] == r["lineage"]) &
                        (nulls["statistic"] == r["statistic"])]
            ps.append((1 + (sel["z"].abs() >= abs(r["z"])).sum()) / (len(sel) + 1) if len(sel) else np.nan)
            ns.append(len(sel))
        df["p_empirical"] = ps
        df["n_null"] = ns

    for k in ("F2", "F2_VS", "F2_D"):
        df = fam[k]
        df["p"] = p_two(df["z"].values)
        df["q_BH"] = bh(df["p"].values)
        df["family_size"] = len(df)
        emp(df)
        if anti and k == "F2":
            df["q_BH_empirical"] = bh(df["p_empirical"].values)
            df["supported"] = df["q_BH_empirical"] < 0.05
        else:
            df["supported"] = (df["q_BH"] < 0.05) & (df["p_empirical"] < 0.05)

    # ---- controls ----
    f1 = fam["F1"].set_index(["trait_id", "lineage"])
    ctl = fam["F1_CTRL"].set_index(["trait_id", "lineage"])
    pos = float(f1.loc[("alt_panukbb", "hepatocyte"), "tau_z"])
    neg = float(ctl.loc[("height_panukbb", "hepatocyte"), "tau_z"]) if ("height_panukbb", "hepatocyte") in ctl.index else None
    hst = float(ctl.loc[("height_panukbb", "stellate"), "tau_z"]) if ("height_panukbb", "stellate") in ctl.index else None
    summary["controls"] = {"positive_alt_hepatocyte_only_tau_z": pos, "positive_pass": pos >= 1.96,
                           "negative_height_hepatocyte_only_tau_z": neg,
                           "negative_pass": (neg is not None and neg < 1.96),
                           "height_stellate_only_tau_z": hst}

    # ---- checks ----
    for (t, m), r in list(ctx.cache.items()):
        rel = np.max(np.abs(r["se_jk"] - r["ldsc_se"]) / r["ldsc_se"])
        chk.append(dict(trait_id=t, model=m, max_rel_diff_se=rel, n_snp=r["n_snp"]))
    chk = pd.DataFrame(chk)
    chk.to_csv(os.path.join(R, "checks.tsv"), sep="\t", index=False)
    summary["max_rel_diff_jackknife_vs_ldsc_se"] = float(chk["max_rel_diff_se"].max())
    summary["regression_snps_by_trait"] = chk.groupby("trait_id")["n_snp"].agg(["min", "max"]).to_dict("index")

    for k, df in fam.items():
        df.to_csv(os.path.join(R, f"{k}.tsv"), sep="\t", index=False)
        summary[f"{k}_family_size"] = int(len(df))
        summary[f"{k}_n_q_lt_0.05"] = int((df["q_BH"] < 0.05).sum())
        if "supported" in df:
            summary[f"{k}_n_supported"] = int(df["supported"].sum())
    for k, df in desc.items():
        if len(df) and "p" not in df and "z" in df:
            df["p"] = p_two(df["z"].values)
        df.to_csv(os.path.join(R, f"desc_{k}.tsv"), sep="\t", index=False)

    # set sizes
    st = ann.loc[ann["draw"].astype(str) == "", ["kind", "set", "lineage", "n_peaks", "n_lifted", "testable"]].copy()
    st["n_common_snps"] = [ctx.n_common_in.get(a, np.nan) for a in st.index]
    st["prop_common_snps"] = [ctx.prop_common.get(a, np.nan) for a in st.index]
    st.to_csv(os.path.join(R, "annotation_sizes.tsv"), sep="\t")
    json.dump(summary, open(os.path.join(R, "summary.json"), "w"), indent=2, default=float)
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:5])
