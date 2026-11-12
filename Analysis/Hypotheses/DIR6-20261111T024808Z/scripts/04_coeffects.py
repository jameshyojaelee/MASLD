#!/usr/bin/env python3
"""DIR6 amendment 1 (PRESPEC_DIR6_AMENDMENT1_coeffects.md): LDL, T2D and CAD co-effects of the
enzyme-lowering allele at the 250 H26 instruments, plus the six-trait BH family.

Order: (1) md5 checked by the sbatch wrapper; (2) known-locus orientation and build asserts;
(3) instrument match and frequency agreement with Pan-UKBB EUR (direction not read); (4) co-effects.
Usage: python -I 04_coeffects.py <project_root> <dir6_run> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT, RUN, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
DL = RUN / "downloads_untrusted"
S1 = RUN / "out_23190132"
S2 = RUN / "out_23210381"
GHOUSE = ROOT / "GWAS/MR_Data/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz"
FILES = {"ldl_glgc": "GCST90239658", "t2d_mvp": "GCST90475667", "cad_aragam": "GCST90132314"}
# known loci: (rsid, allele, expected sign of its effect)
KNOWN = {"ldl_glgc": [("rs11591147", "T", -1), ("rs6511720", "T", -1), ("rs7412", "T", -1)],
         "cad_aragam": [("rs10455872", "G", 1), ("rs1333049", "C", 1), ("rs6511720", "T", -1)],
         "t2d_mvp": [("rs7903146", "T", 1)]}
# liver-fat-raising alleles of the parent DIR6 controls (PRESPEC_DIR6.md section 7); co-effects are
# reported descriptively, outside the BH family, and do not stop the job
CONTROLS = [("PNPLA3", "rs738409", "G"), ("TM6SF2", "rs58542926", "T"), ("GCKR", "rs1260326", "T")]


def log(*a):
    print(*a, flush=True)


def pair_eq(a1, a2, b1, b2):
    return ((a1 == b1) & (a2 == b2)) | ((a1 == b2) & (a2 == b1))


def scan(path, keys, rsids):
    hdr = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    rscol = "rsid" if "rsid" in hdr else ("hm_rsid" if "hm_rsid" in hdr else None)
    want = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "standard_error"]
    want += [c for c in ["beta", "odds_ratio", "effect_allele_frequency", "ci_upper", "ci_lower",
                         "p_value", "hm_code", rscol] if c and c in hdr]
    log(f"{path.name}: columns used {want}")
    out = []
    n_all = n_se_na = n_p0 = 0
    for ch in pd.read_csv(path, sep="\t", usecols=want, dtype=str, chunksize=2_000_000):
        n_all += len(ch)
        n_se_na += int(pd.to_numeric(ch["standard_error"], errors="coerce").isna().sum())
        if "p_value" in ch:
            n_p0 += int((pd.to_numeric(ch["p_value"], errors="coerce") == 0).sum())
        k = ch["chromosome"] + ":" + ch["base_pair_location"]
        m = k.isin(keys) | (ch[rscol].isin(rsids) if rscol else False)
        out.append(ch[m].assign(key=k[m]))
    d = pd.concat(out, ignore_index=True)
    if rscol and rscol != "rsid":
        d = d.rename(columns={rscol: "rsid"})
    for c in ["effect_allele", "other_allele"]:
        d[c] = d[c].str.upper()
    b = pd.to_numeric(d.get("beta"), errors="coerce") if "beta" in d else pd.Series(np.nan, index=d.index)
    if "odds_ratio" in d:
        lor = np.log(pd.to_numeric(d["odds_ratio"], errors="coerce"))
        d["effect_source"] = np.where(b.notna(), "beta", "ln(odds_ratio)")
        b = b.fillna(lor)
    else:
        d["effect_source"] = "beta"
    d["b"] = b
    d["s"] = pd.to_numeric(d["standard_error"], errors="coerce")
    d["se_source"] = np.where(d["s"].notna(), "standard_error", "missing")
    # Deviation (logged): the MVP T2D file has standard_error = NA; SE of ln OR is taken from the
    # 95% CI, (ln ci_upper - ln ci_lower) / (2 * 1.959964), only where standard_error is missing.
    if "ci_upper" in d and "ci_lower" in d:
        se_ci = (np.log(pd.to_numeric(d["ci_upper"], errors="coerce"))
                 - np.log(pd.to_numeric(d["ci_lower"], errors="coerce"))) / (2 * stats.norm.isf(0.025))
        fill = d["s"].isna() & np.isfinite(se_ci)
        d.loc[fill, "s"] = se_ci[fill]
        d.loc[fill, "se_source"] = "95% CI of odds_ratio"
    d["p_reported"] = pd.to_numeric(d["p_value"], errors="coerce") if "p_value" in d else np.nan
    # second fallback where SE and CI are both missing: SE = |b| / |Z|, Z from the reported p
    zp = pd.Series(stats.norm.isf(pd.to_numeric(d["p_reported"]) / 2), index=d.index)
    fill2 = d["s"].isna() & np.isfinite(zp) & (zp > 0) & d["b"].notna()
    d.loc[fill2, "s"] = np.abs(d.loc[fill2, "b"]) / zp[fill2]
    d.loc[fill2, "se_source"] = "abs(effect) / Z from p_value"
    log(f"{path.name}: full file rows {n_all}, standard_error missing {n_se_na}, p_value == 0 {n_p0}, "
        f"ci columns present {'ci_upper' in hdr and 'ci_lower' in hdr}; "
        f"kept rows se_source {d['se_source'].value_counts().to_dict()}")
    d["eaf"] = pd.to_numeric(d["effect_allele_frequency"], errors="coerce") if "effect_allele_frequency" in d else np.nan
    return d


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    inst = pd.read_csv(S1 / "instruments.tsv", sep="\t", dtype={"chr": str, "vid": str})
    u = inst.drop_duplicates("vid")[["vid", "chr", "pos38", "a_raise", "a_other", "eaf_raise_panukbb"]].copy()
    u["key"] = u["chr"] + ":" + u["pos38"].astype(str)
    keys = set(u["key"])
    known_rs = {r for v in KNOWN.values() for r, _, _ in v} | {r for _, r, _ in CONTROLS}

    # build reference: GRCh38 positions of the known rsIDs in the Ghouse file
    gh = []
    for ch in pd.read_csv(GHOUSE, sep="\t", usecols=["chromosome", "base_pair_location", "rsid"], dtype=str, chunksize=2_000_000):
        gh.append(ch[ch["rsid"].isin(known_rs)])
    gh = pd.concat(gh).drop_duplicates("rsid").set_index("rsid")

    checks, rows, agree, ctrl, matches = [], [], [], [], []
    for tr, acc in FILES.items():
        d = scan(DL / acc / f"{acc}.h.tsv.gz", keys, known_rs)
        # parent controls (descriptive only)
        for gene, rs, allele in CONTROLS:
            r = d[d["rsid"] == rs] if "rsid" in d else d.iloc[:0]
            r = r[(r["effect_allele"] == allele) | (r["other_allele"] == allele)]
            if len(r) != 1:
                ctrl.append(dict(trait=tr, gene=gene, rsid=rs, liver_fat_raising_allele=allele, n_rows=len(r)))
                continue
            r = r.iloc[0]
            e = r["b"] if r["effect_allele"] == allele else -r["b"]
            pos_ok = rs in gh.index and f"{gh.loc[rs, 'chromosome']}:{gh.loc[rs, 'base_pair_location']}" == r["key"]
            ctrl.append(dict(trait=tr, gene=gene, rsid=rs, liver_fat_raising_allele=allele, n_rows=1,
                             position=r["key"], build_pass=bool(pos_ok), effect=e, se=r["s"],
                             se_source=r["se_source"], z=e / r["s"], p=2 * stats.norm.sf(abs(e / r["s"]))))
        # 1-2. known-locus orientation and build
        for rs, allele, sign in KNOWN[tr]:
            r = d[d.get("rsid") == rs] if "rsid" in d else d.iloc[:0]
            r = r[(r["effect_allele"] == allele) | (r["other_allele"] == allele)]
            if len(r) != 1:
                checks.append(dict(trait=tr, rsid=rs, allele=allele, n_rows=len(r), effect=np.nan, orientation_pass=False, build_pass=False))
                continue
            r = r.iloc[0]
            e = r["b"] if r["effect_allele"] == allele else -r["b"]
            bpass = rs in gh.index and f"{gh.loc[rs, 'chromosome']}:{gh.loc[rs, 'base_pair_location']}" == r["key"]
            checks.append(dict(trait=tr, rsid=rs, allele=allele, n_rows=1, effect=e, z=e / r["s"],
                               orientation_pass=bool(np.sign(e) == sign), build_pass=bool(bpass), position=r["key"]))
        # 3. instrument match and frequency agreement (direction not read)
        m = u.merge(d, on="key")
        m = m[pair_eq(m["effect_allele"], m["other_allele"], m["a_raise"], m["a_other"])]
        m = m[~m["vid"].duplicated(keep=False)]
        pal = m.apply(lambda x: {x["a_raise"], x["a_other"]} in ({"A", "T"}, {"C", "G"}), axis=1)
        f_raise = np.where(m["effect_allele"] == m["a_raise"], m["eaf"], 1 - m["eaf"])
        ok = (~pal) & np.isfinite(f_raise) & np.isfinite(m["eaf_raise_panukbb"])
        within = np.abs(f_raise[ok] - m.loc[ok, "eaf_raise_panukbb"]) <= 0.10
        # descriptive only (not a gate): palindromic instruments are kept and aligned on the harmonised
        # effect_allele label, as in the parent run; their frequency agreement and hm_code are reported
        okp = pal & np.isfinite(f_raise) & np.isfinite(m["eaf_raise_panukbb"])
        within_p = np.abs(f_raise[okp] - m.loc[okp, "eaf_raise_panukbb"]) <= 0.10
        lp = np.abs(np.log10(2 * stats.norm.sf(np.abs(m["b"] / m["s"]))) - np.log10(m["p_reported"]))
        agree.append(dict(trait=tr, accession=acc, instruments=len(u), matched=len(m), frac_matched=len(m) / len(u),
                          freq_checkable=int(ok.sum()), freq_within_010=int(within.sum()),
                          frac_freq_agree=float(within.mean()) if ok.sum() else np.nan,
                          effect_source=",".join(sorted(m["effect_source"].unique())),
                          se_source=",".join(sorted(m["se_source"].unique())),
                          median_abs_log10p_diff=float(np.nanmedian(lp)) if np.isfinite(lp).any() else np.nan,
                          palindromic_matched=int(pal.sum()), palindromic_freq_checkable=int(okp.sum()),
                          palindromic_freq_within_010=int(within_p.sum()),
                          palindromic_hm_codes=",".join(sorted(m.loc[pal, "hm_code"].astype(str).unique())) if "hm_code" in m else ""))
        matches.append(m.assign(trait=tr, palindromic=pal.to_numpy(), eaf_raise_file=f_raise))
        for x in m.itertuples(index=False):
            b_raise = x.b if x.effect_allele == x.a_raise else -x.b
            rows.append(dict(id=x.vid, trait=tr, beta_al=b_raise, se=x.s))
    ck = pd.DataFrame(checks); ag = pd.DataFrame(agree)
    ck.to_csv(OUT / "known_locus_checks.tsv", sep="\t", index=False)
    pd.DataFrame(ctrl).to_csv(OUT / "controls_coeffects.tsv", sep="\t", index=False)
    pd.concat(matches, ignore_index=True).to_csv(OUT / "instrument_matches.tsv", sep="\t", index=False)
    log(pd.DataFrame(ctrl).to_string())
    ag.to_csv(OUT / "instrument_match_freq_agreement.tsv", sep="\t", index=False)
    log(ck.to_string()); log(ag.to_string())
    if not (ck["orientation_pass"].all() and ck["build_pass"].all()):
        sys.exit("KNOWN-LOCUS ORIENTATION OR BUILD CHECK FAILED; stopping before co-effects")
    bad = ag[ag["frac_freq_agree"].notna() & (ag["frac_freq_agree"] < 0.90)]
    if len(bad):
        sys.exit(f"FREQUENCY AGREEMENT BELOW 0.90: {bad['trait'].tolist()}; stopping")

    # 4. co-effects, six-trait family
    new = pd.DataFrame(rows)
    old = pd.read_csv(S1 / "rows.tsv", sep="\t", dtype={"id": str})
    old = old[old["trait"].isin(["tg_glgc", "bmi_yengo", "homair_manning"]) & ~old["id"].str.startswith("ctrl_")]
    a_raise = u.set_index("vid")["a_raise"]
    old["beta_al"] = np.where(old["ea"] == old["id"].map(a_raise), old["beta"], -old["beta"])
    co = pd.concat([old[["id", "trait", "beta_al", "se"]], new], ignore_index=True)
    co = co[np.isfinite(co["beta_al"]) & np.isfinite(co["se"]) & (co["se"] > 0)]
    co["beta_lower"] = -co["beta_al"]
    co["z"] = co["beta_lower"] / co["se"]
    co["p"] = 2 * stats.norm.sf(np.abs(co["z"]))
    p = co["p"].to_numpy(); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    co["q6"] = np.minimum(q, 1)
    co["flag"] = np.where(co["q6"] < 0.05, np.where(co["beta_lower"] > 0, "enzyme-lowering allele raises (q < 0.05)",
                                                    "enzyme-lowering allele lowers (q < 0.05)"), "no evidence (q >= 0.05)")
    co.to_csv(OUT / "coeffects_six_trait_family.tsv", sep="\t", index=False)
    log(f"six-trait family: {n} tests; q < 0.05: {(co['q6'] < 0.05).sum()}")

    # power: effect detectable with 80% power at alpha 0.05 / n (two-sided)
    zc = stats.norm.isf(0.05 / n / 2) + stats.norm.isf(0.20)
    pw = co.groupby("trait").agg(tests=("se", "size"), median_se=("se", "median")).reset_index()
    pw["mde_80pct_alpha_bonf"] = zc * pw["median_se"]
    pw["units"] = pw["trait"].map({"ldl_glgc": "SD (inverse-normal LDL) per allele", "t2d_mvp": "log odds per allele",
                                   "cad_aragam": "log odds per allele", "tg_glgc": "SD (inverse-normal log TG) per allele",
                                   "bmi_yengo": "SD BMI per allele", "homair_manning": "log HOMA-IR per allele"})
    pw.to_csv(OUT / "power.tsv", sep="\t", index=False); log(pw.to_string())

    # per-gene and per-locus counts
    g = pd.read_csv(S2 / "dir6_gene_classes.tsv", sep="\t", dtype={"vid": str})[["symbol", "call", "vid", "locus_id", "combined"]]
    wide = co.pivot(index="id", columns="trait", values="flag")
    out = []
    for tr in ["ldl_glgc", "t2d_mvp", "cad_aragam", "tg_glgc", "bmi_yengo", "homair_manning"]:
        g[tr] = g["vid"].map(wide[tr]) if tr in wide else np.nan
        g[tr] = g[tr].fillna("variant not in file")
        for call in ["all", "inhibit", "activate"]:
            s = g if call == "all" else g[g["call"] == call]
            for unit in ["gene", "independent locus"]:
                if unit == "gene":
                    vals = s[tr]
                else:  # a locus has a reading if any gene there does; flagged if any gene there is flagged
                    vals = s.groupby("locus_id")[tr].agg(lambda v: next((x for x in v if "q < 0.05" in x),
                                                                       next((x for x in v if x != "variant not in file"), "variant not in file")))
                vc = vals.value_counts()
                out.append(dict(trait=tr, call=call, unit=unit, n=len(vals),
                                with_reading=int((vals != "variant not in file").sum()),
                                lowering_allele_raises_q05=int(vc.get("enzyme-lowering allele raises (q < 0.05)", 0)),
                                lowering_allele_lowers_q05=int(vc.get("enzyme-lowering allele lowers (q < 0.05)", 0)),
                                no_evidence=int(vc.get("no evidence (q >= 0.05)", 0))))
    pd.DataFrame(out).to_csv(OUT / "coeffect_counts.tsv", sep="\t", index=False)
    g.to_csv(OUT / "gene_coeffects.tsv", sep="\t", index=False)
    log(pd.DataFrame(out)[lambda x: x.unit == "gene"].to_string())
    log("COEFF DONE")


if __name__ == "__main__":
    main()
