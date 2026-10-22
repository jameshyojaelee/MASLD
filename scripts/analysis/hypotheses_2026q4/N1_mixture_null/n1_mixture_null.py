#!/usr/bin/env python3
"""N1 step one: technical null for the GSE193066 "possibly mixed" call.

Implements Analysis/Hypotheses/N1-<timestamp>/PRESPEC.md. Reads existing allele
depths (AD) at the identity-panel SNVs; reads no BAM and reruns no published analysis.
Unit: library (nested in person).
"""
import argparse
import gzip
import json
import multiprocessing as mp
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import betaln, logsumexp

COHORTS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523", "GSE174478",
           "GSE193066", "GSE213621", "GSE240729", "PRJNA512027"]
TARGET = "GSE193066"
MIN_DEPTH = 10
MIN_SITES = 500
FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30]
STEP = {"thin_ref": 1, "pick_ref": 2, "pc": 3, "nc": 4}


# ---------------------------------------------------------------- data
def site_class(ref, alt):
    """0: C>T or G>A (deamination makes REF read as ALT); 1: T>C or A>G; 2: transversion."""
    if (ref, alt) in (("C", "T"), ("G", "A")):
        return 0
    if (ref, alt) in (("T", "C"), ("A", "G")):
        return 1
    return 2


def load_targets(path):
    t = pd.read_csv(path, sep="\t")
    t = t[t["chrom"].str.fullmatch(r"chr[0-9]+")].reset_index(drop=True)
    t["cls"] = [site_class(r, a) for r, a in zip(t["ref"], t["alt"])]
    index = {(c, int(p)): i for i, (c, p) in enumerate(zip(t["chrom"], t["pos"]))}
    return t, index


def load_vcf(path, index, n_sites):
    samples, ref, alt = None, None, None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = f[9:]
                ref = np.zeros((n_sites, len(samples)), dtype=np.int32)
                alt = np.zeros((n_sites, len(samples)), dtype=np.int32)
                continue
            i = index.get((f[0], int(f[1])))
            if i is None:
                continue
            ad_pos = f[8].split(":").index("AD")
            for s, cell in enumerate(f[9:]):
                parts = cell.split(":")
                if len(parts) <= ad_pos:
                    continue
                ad = parts[ad_pos].split(",")
                if ad[0] not in (".", ""):
                    ref[i, s] = int(ad[0])
                if len(ad) > 1 and ad[1] not in (".", ""):
                    alt[i, s] = int(ad[1])
    return samples, ref, alt


def naive_dose(ref, alt, min_depth=MIN_DEPTH):
    n = ref + alt
    frac = np.where(n > 0, alt / np.maximum(n, 1), np.nan)
    dose = np.full(n.shape, np.nan)
    ok = n >= min_depth
    dose[ok & (frac < 0.15)] = 0
    dose[ok & (frac >= 0.15) & (frac <= 0.85)] = 1
    dose[ok & (frac > 0.85)] = 2
    return dose


def alt_frequency(dose):
    n = np.sum(~np.isnan(dose), axis=1)
    return (np.nansum(dose, axis=1) + 1.0) / (2.0 * n + 2.0), n


# ---------------------------------------------------------------- likelihood
def log_bb(a, n, pi, rho):
    pi = np.clip(pi, 1e-6, 1 - 1e-6)
    s = (1.0 - rho) / rho
    al, be = pi * s, (1.0 - pi) * s
    return betaln(a + al, n - a + be) - betaln(al, be)


def loglik(theta, data, model, rho_het_fixed, rho_hom_fixed):
    n, a, p, cls = data
    k = 0
    e_ct, e_tc, e_tv = 10.0 ** theta[0], 10.0 ** theta[1], 10.0 ** theta[2]
    b = theta[3]
    k = 4
    if rho_hom_fixed is None:
        rho_hom = 10.0 ** theta[k]
        k += 1
    else:
        rho_hom = rho_hom_fixed
    if rho_het_fixed is None:
        rho_het = 10.0 ** theta[k]
        k += 1
    else:
        rho_het = rho_het_fixed
    f = theta[k] if model != "M0" else 0.0
    e_ra = np.where(cls == 0, e_ct, np.where(cls == 1, e_tc, e_tv))
    e_ar = np.where(cls == 0, e_tc, np.where(cls == 1, e_ct, e_tv))
    prior = np.stack([2 * np.log1p(-p), np.log(2) + np.log(p) + np.log1p(-p), 2 * np.log(p)])
    m_major = (0.0, b, 1.0)
    terms = []
    for g in range(3):
        rho = rho_het if g == 1 else rho_hom
        if model == "M2":
            for h in range(3):
                m = (1 - f) * m_major[g] + f * h / 2.0
                pi = m * (1 - e_ar) + (1 - m) * e_ra
                terms.append(prior[g] + prior[h] + log_bb(a, n, pi, rho))
        else:
            m = (1 - f) * m_major[g] + (f * p if model == "MP" else 0.0)
            pi = m * (1 - e_ar) + (1 - m) * e_ra
            terms.append(prior[g] + log_bb(a, n, pi, rho))
    return float(np.sum(logsumexp(np.stack(terms), axis=0)))


def fit(data, model, rho_het_fixed, rho_hom_fixed, start=None):
    x0 = [-2.5, -2.5, -2.5, 0.5]
    bounds = [(-5, np.log10(0.25))] * 3 + [(0.3, 0.7)]
    if rho_hom_fixed is None:
        x0.append(-2.0)
        bounds.append((-4, np.log10(0.5)))
    if rho_het_fixed is None:
        x0.append(-1.0)
        bounds.append((-3, np.log10(0.8)))
    if start is not None:
        x0 = list(start)
    best = None
    f_starts = [None] if model == "M0" else [0.03, 0.15]
    for f0 in f_starts:
        xs, bs = list(x0), list(bounds)
        if f0 is not None:
            xs.append(f0)
            bs.append((0.0, 0.5))
        res = minimize(lambda th: -loglik(th, data, model, rho_het_fixed, rho_hom_fixed),
                       np.array(xs), method="L-BFGS-B", bounds=bs)
        if best is None or res.fun < best.fun:
            best = res
    return -best.fun, best.x


def unpack(x, model, rho_het_fixed, rho_hom_fixed, prefix):
    out = {f"{prefix}e_CT": 10 ** x[0], f"{prefix}e_TC": 10 ** x[1], f"{prefix}e_tv": 10 ** x[2],
           f"{prefix}b": x[3]}
    k = 4
    if rho_hom_fixed is None:
        out[f"{prefix}rho_hom"] = 10 ** x[k]
        k += 1
    if rho_het_fixed is None:
        out[f"{prefix}rho_het"] = 10 ** x[k]
        k += 1
    if model != "M0":
        out[f"{prefix}f"] = x[k]
    return out


def fit_free_rho(task):
    """Pass A on reference: M0 with rho_het and rho_hom free."""
    key, data = task
    ll, x = fit(data, "M0", None, None)
    out = {"key": key, "n_sites": int(len(data[0])), "ll_M0_free": ll}
    out.update(unpack(x, "M0", None, None, "free_"))
    return out


def fit_models(task):
    """M0, M2, MP with rho_het fixed (primary); again with rho_hom fixed (sensitivity)."""
    key, data, rho_het, rho_hom_ref = task
    out = {"key": key, "n_sites": int(len(data[0])), "panel_depth_used": int(data[0].sum())}
    for tag, rho_hom_fixed in (("", None), ("strict_", rho_hom_ref)):
        ll0, x0 = fit(data, "M0", rho_het, rho_hom_fixed)
        out.update(unpack(x0, "M0", rho_het, rho_hom_fixed, f"{tag}M0_"))
        res = {}
        for model in ("M2", "MP"):
            ll, x = fit(data, model, rho_het, rho_hom_fixed, start=x0)
            if ll < ll0:  # nested at f = 0
                ll, x = ll0, np.append(x0, 0.0)
            res[model] = (ll, x)
            out[f"{tag}ll_{model}"] = ll
            out.update(unpack(x, model, rho_het, rho_hom_fixed, f"{tag}{model}_"))
        out[f"{tag}ll_M0"] = ll0
        better = "M2" if res["M2"][0] >= res["MP"][0] else "MP"
        out[f"{tag}LLR"] = res[better][0] - ll0
        out[f"{tag}f_hat"] = float(res[better][1][-1])
        out[f"{tag}M2_minus_MP"] = res["M2"][0] - res["MP"][0]
    return out


# ---------------------------------------------------------------- helpers
def usable(ref, alt, p, cls):
    n = ref + alt
    keep = n >= MIN_DEPTH
    return (n[keep].astype(float), alt[keep].astype(float), p[keep], cls[keep])


def thin(rng, ref, alt, rate):
    rate = min(1.0, rate)
    return rng.binomial(ref, rate), rng.binomial(alt, rate)


def model_free(ref, alt, p, cls):
    """Other-allele reads at naively homozygous sites (depth >= 20), by direction, and the
    depth-weighted slope of other-allele fraction on the population frequency of the other allele."""
    n = ref + alt
    frac = np.where(n > 0, alt / np.maximum(n, 1), 0.0)
    hr = (n >= 20) & (frac < 0.15)
    ha = (n >= 20) & (frac > 0.85)
    other = np.concatenate([alt[hr], ref[ha]]).astype(float)
    depth = np.concatenate([n[hr], n[ha]]).astype(float)
    q = np.concatenate([p[hr], 1 - p[ha]])
    # direction of the miscall that would create the other allele
    d_ct = np.concatenate([cls[hr] == 0, cls[ha] == 1])
    d_tc = np.concatenate([cls[hr] == 1, cls[ha] == 0])
    d_tv = np.concatenate([cls[hr] == 2, cls[ha] == 2])
    out = {"hom_sites_dp20": int(len(depth)),
           "other_frac_hom": other.sum() / depth.sum() if depth.sum() else np.nan}
    for name, d in (("CT", d_ct), ("TC", d_tc), ("tv", d_tv)):
        out[f"other_frac_{name}_direction"] = other[d].sum() / depth[d].sum() if depth[d].sum() else np.nan
    if len(depth) >= 100:
        w = depth
        y = other / depth
        xm, ym = np.average(q, weights=w), np.average(y, weights=w)
        sxx = np.sum(w * (q - xm) ** 2)
        slope = np.sum(w * (q - xm) * (y - ym)) / sxx if sxx > 0 else np.nan
        out["slope_on_other_allele_freq"] = slope
        out["intercept_at_freq0"] = ym - slope * xm
    return out


def star_facts(root, cohort, run):
    d = root / cohort / "alignments" / "star" / run
    out = {"run": run, "cohort": cohort}
    log = d / f"{run}.Log.final.out"
    if log.exists():
        txt = log.read_text()
        for name, pat in (("input_reads", r"Number of input reads \|\s+(\d+)"),
                          ("avg_input_read_length", r"Average input read length \|\s+(\d+)"),
                          ("unique_pct", r"Uniquely mapped reads % \|\s+([\d.]+)%"),
                          ("mismatch_rate_pct", r"Mismatch rate per base, % \|\s+([\d.]+)%"),
                          ("multi_pct", r"% of reads mapped to multiple loci \|\s+([\d.]+)%")):
            m = re.search(pat, txt)
            out[name] = float(m.group(1)) if m else np.nan
    rpg = d / f"{run}.ReadsPerGene.out.tab"
    if rpg.exists():
        t = pd.read_csv(rpg, sep="\t", header=None, skiprows=4)
        fwd, rev = t[2].sum(), t[3].sum()
        out["strand_fwd_fraction"] = fwd / (fwd + rev) if fwd + rev else np.nan
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--identity", required=True, help="sample_to_individual.tsv")
    ap.add_argument("--kinship", required=True, help="kinship_related_pairs.tsv")
    ap.add_argument("--results-root", required=True, help="Patient_Cohorts/results")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--n-ref", type=int, default=250)
    ap.add_argument("--n-control", type=int, default=100)
    a = ap.parse_args()
    out = Path(a.out)
    seed = a.seed
    (out / "seeds.json").write_text(json.dumps(
        {"base_seed": seed, "rng": "numpy.random.default_rng([base_seed, step, stratum, index])",
         "step_codes": STEP}, indent=2))

    targets, index = load_targets(a.targets)
    cls_all = targets["cls"].to_numpy()
    ident = pd.read_csv(a.identity, sep="\t")
    kin = pd.read_csv(a.kinship, sep="\t")
    kin_pairs = set(map(frozenset, zip(kin["#IID1"], kin["IID2"])))

    data, freq = {}, {}
    pooled_dose = []
    for c in COHORTS:
        samples, ref, alt = load_vcf(Path(a.vcf_dir) / f"{c}.vcf.gz", index, len(targets))
        data[c] = (samples, ref, alt)
        d = naive_dose(ref, alt)
        pooled_dose.append(d)
        freq[c] = alt_frequency(d)
        print(c, "libraries", len(samples), "sites with depth", int(((ref + alt).sum(axis=1) > 0).sum()), flush=True)
    p_all, _ = alt_frequency(np.concatenate(pooled_dose, axis=1))
    p_cohort = {c: np.where(freq[c][1] >= 20, freq[c][0], p_all) for c in COHORTS}
    pd.DataFrame({"chrom": targets["chrom"], "pos": targets["pos"], "ref": targets["ref"], "alt": targets["alt"],
                  "cls": cls_all, "p_all": p_all, **{f"p_{c}": p_cohort[c] for c in COHORTS},
                  **{f"n_{c}": freq[c][1] for c in COHORTS}}).to_csv(out / "site_table.tsv", sep="\t", index=False)

    # library table: depth, model-free summaries, STAR facts
    rows = []
    for c in COHORTS:
        samples, ref, alt = data[c]
        for s, run in enumerate(samples):
            r = {"run": run, "cohort": c, "panel_depth": int((ref[:, s] + alt[:, s]).sum()),
                 "usable_sites": int(((ref[:, s] + alt[:, s]) >= MIN_DEPTH).sum())}
            r.update(model_free(ref[:, s], alt[:, s], p_cohort[c], cls_all))
            rows.append(r)
    lib = pd.DataFrame(rows).merge(
        ident[["run", "expr_sex", "passes_min_sites", "individual_id", "other_allele_fraction_at_hom"]],
        on="run", how="left")
    facts = pd.DataFrame([star_facts(Path(a.results_root), r.cohort, r.run) for r in lib.itertuples()])
    lib = lib.merge(facts, on=["run", "cohort"], how="left")
    lib.to_csv(out / "library_table.tsv", sep="\t", index=False)
    lib.groupby("cohort").agg(
        libraries=("run", "size"), median_input_reads=("input_reads", "median"),
        median_read_length=("avg_input_read_length", "median"), median_unique_pct=("unique_pct", "median"),
        median_mismatch_rate_pct=("mismatch_rate_pct", "median"),
        median_strand_fwd_fraction=("strand_fwd_fraction", "median"),
        median_panel_depth=("panel_depth", "median"), median_usable_sites=("usable_sites", "median"),
        median_other_frac_hom=("other_frac_hom", "median"),
        median_other_frac_CT_direction=("other_frac_CT_direction", "median"),
        median_other_frac_TC_direction=("other_frac_TC_direction", "median"),
        median_other_frac_tv_direction=("other_frac_tv_direction", "median"),
        median_slope_on_other_allele_freq=("slope_on_other_allele_freq", "median"),
        median_intercept_at_freq0=("intercept_at_freq0", "median"),
    ).to_csv(out / "library_type_and_error_by_cohort.tsv", sep="\t")

    # site recurrence: share of other-allele reads at naive-hom sites from the top 5% of sites
    rec = []
    for c in COHORTS:
        _, ref, alt = data[c]
        n = ref + alt
        frac = np.where(n > 0, alt / np.maximum(n, 1), 0.0)
        other = np.where((n >= 20) & (frac < 0.15), alt, 0) + np.where((n >= 20) & (frac > 0.85), ref, 0)
        per_site = np.sort(other.sum(axis=1))[::-1]
        top = max(1, int(round(0.05 * np.sum(per_site >= 0))))
        rec.append({"cohort": c, "other_reads_total": int(per_site.sum()),
                    "share_from_top5pct_sites": per_site[:top].sum() / max(per_site.sum(), 1)})
    pd.DataFrame(rec).to_csv(out / "site_recurrence_by_cohort.tsv", sep="\t", index=False)

    # strata
    tg = lib[lib["cohort"] == TARGET].copy()
    tg["stratum"] = pd.qcut(tg["panel_depth"].rank(method="first"), 5, labels=False)
    strata_depth = tg.groupby("stratum")["panel_depth"].median().to_dict()
    ref_ok = lib[(lib["cohort"] != TARGET) & lib["passes_min_sites"].fillna(False).astype(bool)
                 & lib["expr_sex"].isin(["M", "F"])].copy()
    col = {c: {run: i for i, run in enumerate(data[c][0])} for c in COHORTS}

    def counts(c, run):
        s = col[c][run]
        return data[c][1][:, s], data[c][2][:, s]

    pool = mp.Pool(a.threads)

    # reference draw per stratum, evenly across cohorts
    ref_sets, control_sets = {}, {}
    for st, depth in strata_depth.items():
        elig = ref_ok[ref_ok["panel_depth"] >= depth]
        rng = np.random.default_rng([seed, STEP["pick_ref"], int(st), 0])
        per = int(np.ceil(a.n_ref / max(elig["cohort"].nunique(), 1)))
        chosen = []
        for c, grp in elig.groupby("cohort"):
            runs = grp["run"].to_numpy()
            chosen += list(rng.choice(runs, size=min(per, len(runs)), replace=False))
        chosen = list(rng.permutation(chosen))[: a.n_ref]
        ref_sets[st] = elig[elig["run"].isin(chosen)]
        control_sets[st] = elig[~elig["run"].isin(chosen)]

    thinned = {}
    for st, depth in strata_depth.items():
        for i, r in enumerate(ref_sets[st].itertuples()):
            rng = np.random.default_rng([seed, STEP["thin_ref"], int(st), i])
            rr, aa = counts(r.cohort, r.run)
            tr, ta = thin(rng, rr, aa, depth / r.panel_depth)
            d = usable(tr, ta, p_cohort[r.cohort], cls_all)
            if len(d[0]) >= MIN_SITES:
                thinned[(int(st), r.run, r.cohort)] = d
    print("reference tasks", len(thinned), flush=True)

    free = pd.DataFrame(pool.map(fit_free_rho, list(thinned.items()), chunksize=4))
    free[["stratum", "run", "cohort"]] = pd.DataFrame(free["key"].tolist(), index=free.index)
    rho = free.groupby("stratum").agg(rho_het=("free_rho_het", "median"), rho_hom=("free_rho_hom", "median"),
                                      n=("run", "size"))
    rho["stratum_median_panel_depth"] = [strata_depth[s] for s in rho.index]
    rho.to_csv(out / "reference_overdispersion_by_stratum.tsv", sep="\t")
    free.drop(columns="key").to_csv(out / "reference_free_rho_fits.tsv", sep="\t", index=False)
    print(rho, flush=True)

    def run_tasks(tasks):
        res = pd.DataFrame(pool.map(fit_models, tasks, chunksize=4))
        return res

    ref_res = run_tasks([(k, d, rho.loc[k[0], "rho_het"], rho.loc[k[0], "rho_hom"]) for k, d in thinned.items()])
    ref_res[["stratum", "run", "cohort"]] = pd.DataFrame(ref_res["key"].tolist(), index=ref_res.index)
    confirmed = set()
    same = kin[(kin["relation"] == "same_individual") & (kin["KINSHIP"] >= 0.45)]
    confirmed.update(same["#IID1"]); confirmed.update(same["IID2"])
    ref_res["confirmed_same_person_pair_member"] = ref_res["run"].isin(confirmed)
    ref_res.drop(columns="key").to_csv(out / "reference_null_fits.tsv", sep="\t", index=False)

    def empirical_p(stratum, llr, column):
        null = ref_res.loc[ref_res["stratum"] == stratum, column].to_numpy()
        return (1 + np.sum(null >= llr)) / (1 + len(null))

    def decide(df, stratum_col="stratum"):
        for tag in ("", "strict_"):
            df[f"{tag}p_emp"] = [empirical_p(s, l, f"{tag}LLR") for s, l in zip(df[stratum_col], df[f"{tag}LLR"])]
            df[f"{tag}called_mixed"] = (df[f"{tag}p_emp"] < 0.01) & (df[f"{tag}f_hat"] >= 0.02)
            df[f"{tag}source"] = np.select(
                [~df[f"{tag}called_mixed"], df[f"{tag}M2_minus_MP"] > 2, df[f"{tag}M2_minus_MP"] < -2],
                ["not_called", "two_person", "pooled_source"], "source_number_unresolved")
        return df

    # leave-one-out false-call rate of the rule within the reference itself
    ref_dec = decide(ref_res.copy())
    ref_dec.groupby(["stratum"]).agg(n=("run", "size"), called_mixed=("called_mixed", "sum"),
                                     strict_called_mixed=("strict_called_mixed", "sum")
                                     ).to_csv(out / "reference_rule_rate_by_stratum.tsv", sep="\t")

    # GSE193066
    tasks, untestable = [], []
    for r in tg.itertuples():
        rr, aa = counts(TARGET, r.run)
        d = usable(rr, aa, p_cohort[TARGET], cls_all)
        if len(d[0]) < MIN_SITES:
            untestable.append(r.run)
            continue
        st = int(r.stratum)
        tasks.append(((st, r.run, TARGET), d, rho.loc[st, "rho_het"], rho.loc[st, "rho_hom"]))
    tg_res = run_tasks(tasks)
    tg_res[["stratum", "run", "cohort"]] = pd.DataFrame(tg_res["key"].tolist(), index=tg_res.index)
    tg_res = decide(tg_res)
    for tag in ("", "strict_"):
        p = tg_res[f"{tag}p_emp"].to_numpy()
        order = np.argsort(p)
        q = np.empty(len(p))
        q[order] = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
        tg_res[f"{tag}bh_q"] = np.minimum(q, 1.0)
    full = tg.merge(tg_res.drop(columns=["key", "cohort", "stratum"]), on="run", how="left")
    full["existing_rule_flag"] = (full["expr_sex"] == "both") | (full["other_allele_fraction_at_hom"] > 0.02)
    full["n1_status"] = np.where(full["run"].isin(untestable), "untestable",
                                 np.where(full["called_mixed"].fillna(False).astype(bool), "remains_called_mixed",
                                          "not_supported_as_mixed"))
    full["n1_status_strict"] = np.where(full["run"].isin(untestable), "untestable",
                                        np.where(full["strict_called_mixed"].fillna(False).astype(bool),
                                                 "remains_called_mixed", "not_supported_as_mixed"))
    full.to_csv(out / "gse193066_library_calls.tsv", sep="\t", index=False)

    # controls at lowest and middle stratum
    st_low, st_mid = min(strata_depth, key=strata_depth.get), sorted(strata_depth, key=strata_depth.get)[2]
    pc_tasks = []
    for st in (st_low, st_mid):
        depth = strata_depth[st]
        cs = control_sets[st]
        cs = cs[cs["individual_id"].notna()]
        for fi, f in enumerate(FRACTIONS):
            rng = np.random.default_rng([seed, STEP["pc"], int(st), fi])
            made = 0
            while made < a.n_control:
                c = rng.choice(cs["cohort"].unique())
                grp = cs[cs["cohort"] == c]
                if len(grp) < 2:
                    continue
                i, j = rng.choice(len(grp), size=2, replace=False)
                A, B = grp.iloc[i], grp.iloc[j]
                if A["individual_id"] == B["individual_id"] or frozenset((A["run"], B["run"])) in kin_pairs:
                    continue
                ar, aa_ = counts(c, A["run"]); br, ba = counts(c, B["run"])
                r1, a1 = thin(rng, ar, aa_, (1 - f) * depth / A["panel_depth"])
                r2, a2 = thin(rng, br, ba, f * depth / B["panel_depth"])
                d = usable(r1 + r2, a1 + a2, p_cohort[c], cls_all)
                if len(d[0]) < MIN_SITES:
                    continue
                pc_tasks.append(((int(st), f"pc|{f}|{A['run']}|{B['run']}", c), d,
                                 rho.loc[st, "rho_het"], rho.loc[st, "rho_hom"]))
                made += 1
    pc = run_tasks(pc_tasks)
    pc[["stratum", "label", "cohort"]] = pd.DataFrame(pc["key"].tolist(), index=pc.index)
    pc["true_f"] = pc["label"].str.split("|").str[1].astype(float)
    pc = decide(pc)
    pc.drop(columns="key").to_csv(out / "positive_control_fits.tsv", sep="\t", index=False)
    pc_sum = pc.groupby(["stratum", "true_f"]).agg(
        n=("label", "size"), called_mixed=("called_mixed", "mean"), median_f_hat=("f_hat", "median"),
        two_person=("source", lambda s: np.mean(s == "two_person")),
        pooled_source=("source", lambda s: np.mean(s == "pooled_source")),
        strict_called_mixed=("strict_called_mixed", "mean"), strict_median_f_hat=("strict_f_hat", "median"))
    pc_sum["stratum_median_panel_depth"] = [strata_depth[s] for s in pc_sum.index.get_level_values(0)]
    pc_sum.to_csv(out / "positive_control_summary.tsv", sep="\t")

    nc_tasks = []
    cs = control_sets[st_mid]
    for kind_i, kind in enumerate(("deamination_3pct", "symmetric_flip_3pct")):
        rng = np.random.default_rng([seed, STEP["nc"], int(st_mid), kind_i])
        picks = cs.iloc[rng.choice(len(cs), size=min(a.n_control, len(cs)), replace=False)]
        for r in picks.itertuples():
            rr, aa_ = counts(r.cohort, r.run)
            tr, ta = thin(rng, rr, aa_, strata_depth[st_mid] / r.panel_depth)
            if kind == "deamination_3pct":
                up = np.where(cls_all == 0, rng.binomial(tr, 0.03), 0)    # C allele (REF) read as T (ALT)
                down = np.where(cls_all == 1, rng.binomial(ta, 0.03), 0)  # C allele (ALT) read as T (REF)
            else:
                up, down = rng.binomial(tr, 0.03), rng.binomial(ta, 0.03)
            tr2, ta2 = tr - up + down, ta + up - down
            d = usable(tr2, ta2, p_cohort[r.cohort], cls_all)
            if len(d[0]) >= MIN_SITES:
                nc_tasks.append(((int(st_mid), f"{kind}|{r.run}", r.cohort), d,
                                 rho.loc[st_mid, "rho_het"], rho.loc[st_mid, "rho_hom"]))
    nc = run_tasks(nc_tasks)
    nc[["stratum", "label", "cohort"]] = pd.DataFrame(nc["key"].tolist(), index=nc.index)
    nc["kind"] = nc["label"].str.split("|").str[0]
    nc = decide(nc)
    nc.drop(columns="key").to_csv(out / "negative_control_fits.tsv", sep="\t", index=False)
    nc_sum = nc.groupby("kind").agg(n=("label", "size"), called_mixed=("called_mixed", "mean"),
                                    median_f_hat=("f_hat", "median"), median_e_CT=("M0_e_CT", "median"),
                                    strict_called_mixed=("strict_called_mixed", "mean"))
    nc_sum.to_csv(out / "negative_control_summary.tsv", sep="\t")
    pool.close()

    # acceptance and headline
    sens10 = float(pc_sum.loc[(st_mid, 0.10), "called_mixed"])
    fp0 = float(pc_sum.loc[(st_mid, 0.0), "called_mixed"])
    accept = bool(sens10 >= 0.80 and fp0 <= 0.05 and (nc_sum["called_mixed"] <= 0.05).all())
    flagged = full[full["existing_rule_flag"]]
    summary = {
        "method_accepted": accept,
        "acceptance": {"sensitivity_f0.10_mid_depth": sens10, "false_call_f0_mid_depth": fp0,
                       "false_call_negative_controls": nc_sum["called_mixed"].to_dict()},
        "strict_method_accepted": bool(
            pc_sum.loc[(st_mid, 0.10), "strict_called_mixed"] >= 0.80
            and pc_sum.loc[(st_mid, 0.0), "strict_called_mixed"] <= 0.05
            and (nc_sum["strict_called_mixed"] <= 0.05).all()),
        "gse193066_libraries": int(len(full)),
        "existing_rule_flagged": int(len(flagged)),
        "flagged_by_status_primary": flagged["n1_status"].value_counts().to_dict(),
        "flagged_by_status_strict": flagged["n1_status_strict"].value_counts().to_dict(),
        "unflagged_by_status_primary": full[~full["existing_rule_flag"]]["n1_status"].value_counts().to_dict(),
        "flagged_remaining_by_source_primary": flagged[flagged["n1_status"] == "remains_called_mixed"]["source"].value_counts().to_dict(),
        "flagged_remaining_bh_q_lt_0.05_primary": int(((flagged["n1_status"] == "remains_called_mixed") & (flagged["bh_q"] < 0.05)).sum()),
        "testable_fraction_called_mixed_primary": float((full["n1_status"] == "remains_called_mixed").sum()
                                                        / max((full["n1_status"] != "untestable").sum(), 1)),
        "gse193066_medians": {k: float(full[k].median()) for k in
                              ("M0_e_CT", "M0_e_TC", "M0_e_tv", "M0_b", "M0_rho_hom", "f_hat", "LLR",
                               "other_frac_hom", "other_frac_CT_direction", "other_frac_TC_direction",
                               "other_frac_tv_direction", "slope_on_other_allele_freq", "intercept_at_freq0")},
        "reference_medians": {k: float(ref_res[k].median()) for k in
                              ("M0_e_CT", "M0_e_TC", "M0_e_tv", "M0_b", "M0_rho_hom", "f_hat", "LLR")},
        "expr_sex_both_by_status": full[full["expr_sex"] == "both"]["n1_status"].value_counts().to_dict(),
        "spearman_existing_fraction_vs_e_CT": float(full[["other_allele_fraction_at_hom", "M0_e_CT"]].corr(method="spearman").iloc[0, 1]),
        "spearman_existing_fraction_vs_f_hat": float(full[["other_allele_fraction_at_hom", "f_hat"]].corr(method="spearman").iloc[0, 1]),
        "strata_median_panel_depth": {str(k): float(v) for k, v in strata_depth.items()},
        "reference_libraries_per_stratum": ref_res.groupby("stratum").size().to_dict(),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
