#!/usr/bin/env python
"""H10: cirrhosis log-OR per SD ALT effect at colocalized loci.

Implements PRESPEC_H10.md in the output directory. Modality: GWAS summary
statistics. Unit: independent locus. All results are provisional on the current
COLOC release.

Usage: h10_cirrhosis_per_alt.py <repo_root> <out_dir>
"""
import json
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import norm

SEED = 20261007
N_BOOT = 10000
WINDOW = 500_000
P_GWS = 5e-8
Z95 = norm.ppf(0.975)
MHC = (6, 25_000_000, 34_000_000)
GPT_WIN = (8, 145_728_356 - 100_000, 145_732_557 + 100_000)
ENZYME = {"ALT", "AST", "GGT"}
DIRECT = {"NAFLD", "NASH", "PDFF"}

REPO, OUT = sys.argv[1], sys.argv[2]
FM = f"{REPO}/GWAS/finemapping"
REL = f"{FM}/results/susie_coloc_r2_20260924T053244Z"
SS = f"{FM}/data/sumstats"

FILES = {
    "MVP_ALT_EUR": f"{FM}/data/sumstats_af/MVP_ALT_EUR_af.tsv",  # same rows as sumstats plus af
    "UKBB_ALT": f"{SS}/UKBB_ALT_reformatted_hg19.tsv",
    "MVP_Cirrhosis_EUR": f"{SS}/MVP_Cirrhosis_EUR_reformatted_hg19.tsv",
    "Ghouse_Cirrhosis": f"{SS}/Ghouse_Cirrhosis_reformatted_hg19.tsv",
    "Ghouse_HCC": f"{SS}/Ghouse_HCC_reformatted_hg19.tsv",
}
N_ALT = {"MVP_ALT_EUR": 406_803, "UKBB_ALT": 355_729}

# name, gene, chr, pos (hg19, checked against Ensembl GRCh37)
CONTROLS = [
    ("rs738409", "PNPLA3", 22, 44324727),
    ("rs58542926", "TM6SF2", 19, 19379549),
    ("rs72613567", "HSD17B13", 4, 88231393),
    ("rs6834314", "HSD17B13", 4, 88213808),
    ("rs2642438", "MTARC1", 1, 220970028),
    ("rs2792751", "GPAM", 10, 113940329),
]

DESIGNS = [
    ("D1_cirrhosis", "MVP_ALT_EUR", "UKBB_ALT", "MVP_Cirrhosis_EUR"),
    ("D2_cirrhosis", "UKBB_ALT", "MVP_ALT_EUR", "Ghouse_Cirrhosis"),
    ("D2_hcc", "UKBB_ALT", "MVP_ALT_EUR", "Ghouse_HCC"),
]


def log(*a):
    print(*a, flush=True)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return q
    idx = np.where(ok)[0][np.argsort(p[ok])]
    ranked = p[idx] * m / np.arange(1, m + 1)
    q[idx] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return q


def poskey(chrom, pos):
    return np.asarray(chrom, np.int64) * 10_000_000_000 + np.asarray(pos, np.int64)


# ----------------------------------------------------------------------------
# 1. Loci from the COLOC release
# ----------------------------------------------------------------------------
genes = pd.read_csv(f"{REL}/gene_level_coloc_tier12.csv")
pos_genes = genes[genes.susie_state_t12 == "positive"].copy()
n_enz_driver = int(pos_genes.driving_trait.isin(ENZYME).sum())
n_dir_driver = int(pos_genes.driving_trait.isin(DIRECT).sum())
log(f"positive genes {len(pos_genes)}; driver enzyme {n_enz_driver}; driver direct {n_dir_driver}")
assert len(pos_genes) == 436 and n_enz_driver == 423 and n_dir_driver == 13, "COLOC release counts changed"

tier = pd.read_csv(f"{FM}/config/gwas_trait_tier.tsv", sep="\t")
tier12 = tier[tier.tier <= 2].set_index("study_name").trait.to_dict()

pairs = pd.read_csv(f"{REL}/susie_coloc_signal_pairs.tsv.gz", sep="\t")
pairs = pairs[
    pairs.eligible.astype(bool)
    & (pairs["PP.H4"] > 0.5)
    & pairs.gwas.isin(tier12)
    & pairs.ensembl.isin(pos_genes.ensembl)
].copy()
pairs["trait"] = pairs.gwas.map(tier12)
pairs["hchr"] = pairs.hit1.str.split(":").str[0].astype(int)
pairs["hpos"] = pairs.hit1.str.split(":").str[1].astype(int)

sym = pairs.groupby("ensembl").gene.first()
support = pairs.groupby("ensembl").trait.agg(lambda s: set(s))


def support_class(traits):
    e, d = bool(traits & ENZYME), bool(traits & DIRECT)
    return "both" if e and d else "enzyme_only" if e else "direct_only" if d else "none"


pos_genes["gene"] = pos_genes.ensembl.map(sym)
pos_genes["support_class"] = pos_genes.ensembl.map(support).apply(
    lambda t: support_class(t) if isinstance(t, set) else "no_eligible_pair"
)
pos_genes["driver_class"] = np.where(pos_genes.driving_trait.isin(ENZYME), "enzyme", "direct")
pos_genes.to_csv(f"{OUT}/genes_classes.tsv", sep="\t", index=False)
log("gene support class:", pos_genes.support_class.value_counts().to_dict())

# single-linkage clustering of candidate positions; merge clusters sharing a gene
cand = pairs[["hchr", "hpos"]].drop_duplicates().sort_values(["hchr", "hpos"]).reset_index(drop=True)
parent = list(range(len(cand)))


def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def union(i, j):
    parent[find(i)] = find(j)


for i in range(1, len(cand)):
    if cand.hchr[i] == cand.hchr[i - 1] and cand.hpos[i] - cand.hpos[i - 1] <= WINDOW:
        union(i, i - 1)
cidx = {(c, p): i for i, (c, p) in enumerate(zip(cand.hchr, cand.hpos))}
for _, g in pairs.groupby("ensembl"):
    ids = [cidx[(c, p)] for c, p in zip(g.hchr, g.hpos)]
    for j in ids[1:]:
        union(ids[0], j)
cand["root"] = [find(i) for i in range(len(cand))]
pairs["root"] = [find(cidx[(c, p)]) for c, p in zip(pairs.hchr, pairs.hpos)]

gene_info = pos_genes.set_index("ensembl")
loci = []
for k, (root, cg) in enumerate(cand.groupby("root", sort=False)):
    pg = pairs[pairs.root == root]
    ens = sorted(pg.ensembl.unique())
    classes = [gene_info.support_class[e] for e in ens]
    chrom, lo, hi = int(cg.hchr.iloc[0]), int(cg.hpos.min()), int(cg.hpos.max())
    loci.append(
        dict(
            locus_id=f"L{k + 1:03d}",
            chr=chrom,
            start=lo,
            end=hi,
            n_genes=len(ens),
            genes=",".join(sorted(str(sym[e]) for e in ens)),
            ensembl=",".join(ens),
            locus_class="enzyme_only" if all(c == "enzyme_only" for c in classes) else "direct_involved",
            driver_traits=",".join(sorted(set(gene_info.driving_trait[e] for e in ens))),
            driver_gwas=",".join(sorted(set(gene_info.susie_driving_gwas[e] for e in ens))),
            mhc=bool(chrom == MHC[0] and hi >= MHC[1] and lo <= MHC[2]),
            is_reference_locus=any(
                chrom == c and lo - WINDOW <= p <= hi + WINDOW for _, _, c, p in CONTROLS
            ),
            cand_positions=",".join(str(p) for p in cg.hpos),
        )
    )
loci = pd.DataFrame(loci)
loci.to_csv(f"{OUT}/loci.tsv", sep="\t", index=False)
log(f"loci {len(loci)}:", loci.locus_class.value_counts().to_dict(), "MHC", int(loci.mhc.sum()))
missing_genes = sorted(set(pos_genes.ensembl) - set(pairs.ensembl))
log("positive genes with no eligible tier-1/2 pair > 0.5:", len(missing_genes))

# ----------------------------------------------------------------------------
# 2. Pull needed rows from each summary-statistics file (chunked; low memory)
# ----------------------------------------------------------------------------
need = set(poskey(cand.hchr, cand.hpos).tolist()) | {int(poskey(c, p)) for _, _, c, p in CONTROLS}
need_arr = np.fromiter(need, np.int64)


def load(name):
    path = FILES[name]
    keep, samp = [], []
    cols = ["chromosome", "position", "allele1", "allele2", "beta", "se"] + (
        ["af"] if name == "MVP_ALT_EUR" else []
    )
    for ch in pd.read_csv(path, sep="\t", usecols=cols, chunksize=2_000_000, dtype={"allele1": str, "allele2": str}):
        ch = ch[pd.to_numeric(ch.chromosome, errors="coerce").between(1, 22)]
        ch["chromosome"] = ch.chromosome.astype(int)
        k = poskey(ch.chromosome, ch.position)
        gpt = (ch.chromosome == GPT_WIN[0]) & ch.position.between(GPT_WIN[1], GPT_WIN[2])
        keep.append(ch[np.isin(k, need_arr) | gpt.values])
        samp.append(ch[(ch.position.values % 23) == 7])  # fixed thinning for the calibration sample
    d, s = pd.concat(keep), pd.concat(samp)
    for x in (d, s):
        x["a1"] = x.allele1.str.upper()
        x["a2"] = x.allele2.str.upper()
        x["pair"] = np.where(x.a1 < x.a2, x.a1 + "/" + x.a2, x.a2 + "/" + x.a1)
        x["vid"] = x.chromosome.astype(str) + ":" + x.position.astype(str) + ":" + x["pair"]
        x["z"] = x.beta / x.se
    d = d.sort_values("z", key=np.abs, ascending=False).drop_duplicates("vid")
    s = s.drop_duplicates("vid")
    log(f"loaded {name}: locus rows {len(d)}, thinned rows {len(s)}")
    return d.set_index("vid"), s.set_index("vid")


DATA, SAMP = {}, {}
for nm in FILES:
    DATA[nm], SAMP[nm] = load(nm)
AF = pd.concat([DATA["MVP_ALT_EUR"].af, SAMP["MVP_ALT_EUR"].af])
AF = AF[~AF.index.duplicated()]


# ----------------------------------------------------------------------------
# 3. Per-variant ratio
# ----------------------------------------------------------------------------
def ratio_row(vid, sel, est, out):
    """Wald ratio and Fieller interval for one variant; None if absent from any file."""
    if not (vid in DATA[sel].index and vid in DATA[est].index and vid in DATA[out].index and vid in AF.index):
        return None
    rs, re, ro = DATA[sel].loc[vid], DATA[est].loc[vid], DATA[out].loc[vid]
    p = float(AF[vid])
    maf = min(p, 1 - p)
    if not (0 < p < 1):
        return None
    ze = re.z
    denom = np.sqrt(2 * p * (1 - p) * (N_ALT[est] + ze**2))
    b_a_signed, se_a = ze / denom, 1 / denom
    # orient the outcome to the estimation cohort's allele1, then to the ALT-raising allele
    flip_o = 1.0 if ro.a1 == re.a1 else -1.0
    flip_s = 1.0 if rs.a1 == re.a1 else -1.0
    s = 1.0 if b_a_signed >= 0 else -1.0
    b_a = abs(b_a_signed)
    b_c, se_c = s * flip_o * ro.beta, float(ro.se)
    theta = b_c / b_a
    g = (Z95 * se_a / b_a) ** 2
    if g < 1:
        half = (Z95 / b_a) * np.sqrt(se_c**2 * (1 - g) + theta**2 * se_a**2) / (1 - g)
        lo, hi = theta / (1 - g) - half, theta / (1 - g) + half
    else:
        lo, hi = -np.inf, np.inf
    zc = b_c / se_c
    return dict(
        lead=vid,
        alt_raising_allele=re.a1 if s > 0 else re.a2,
        other_allele=re.a2 if s > 0 else re.a1,
        maf=maf,
        palindromic=rs["pair"] in ("A/T", "C/G"),
        sel_z=s * flip_s * rs.z,
        sel_p=2 * norm.sf(abs(rs.z)),
        alt_native_beta=s * re.beta,
        alt_native_se=float(re.se),
        alt_z=abs(ze),
        alt_p=2 * norm.sf(abs(ze)),
        alt_sd=b_a,
        alt_sd_se=se_a,
        out_logor=b_c,
        out_se=se_c,
        out_z=zc,
        out_p=2 * norm.sf(abs(zc)),
        theta=theta,
        theta_lo=lo,
        theta_hi=hi,
        fieller_bounded=bool(g < 1),
    )


def run_design(tag, sel, est, out):
    log(f"\n=== {tag}: select {sel}; estimate {est}; outcome {out} ===")
    # --- positive controls
    ctl = []
    for rsid, gene, c, p in CONTROLS:
        hits = [v for v in DATA[sel].index[(DATA[sel].chromosome == c) & (DATA[sel].position == p)]]
        rows = [r for r in (ratio_row(v, sel, est, out) for v in hits) if r]
        if not rows:
            ctl.append(dict(control=rsid, gene=gene, present=False))
            continue
        r = min(rows, key=lambda x: x["sel_p"])
        ctl.append(dict(control=rsid, gene=gene, present=True, **r))
    ctl = pd.DataFrame(ctl)
    # HSD17B13: indel if present, else the SNV proxy
    if bool(ctl.loc[ctl.control == "rs72613567", "present"].iloc[0]):
        ctl = ctl[ctl.control != "rs6834314"]
    else:
        ctl = ctl[ctl.control != "rs72613567"]
    ctl = ctl.reset_index(drop=True)
    ctl["qualifies"] = ctl.present & (ctl.get("alt_p", np.nan) < P_GWS) & (ctl.get("theta", np.nan) > 0) & (
        ctl.get("out_p", np.nan) < 0.01
    )
    qual = ctl[ctl.qualifies]
    pc_pass = bool(
        len(qual) >= 3 and {"PNPLA3", "TM6SF2"} <= set(qual.gene)
    )
    L_ref = float(qual.theta.min()) if len(qual) else np.nan
    log(f"positive control pass: {pc_pass}; qualifying {list(qual.gene)}; L_ref {L_ref:.4g}")

    # --- loci
    rows = []
    for L in loci.itertuples():
        best = None
        for p in L.cand_positions.split(","):
            m = (DATA[sel].chromosome == L.chr) & (DATA[sel].position == int(p))
            for v in DATA[sel].index[m]:
                r = ratio_row(v, sel, est, out)
                if r and (best is None or (r["sel_p"], int(p)) < (best["sel_p"], int(best["lead"].split(":")[1]))):
                    best = r
        base = dict(
            locus_id=L.locus_id, chr=L.chr, start=L.start, end=L.end, genes=L.genes, n_genes=L.n_genes,
            locus_class=L.locus_class, mhc=L.mhc, is_reference_locus=L.is_reference_locus,
            driver_traits=L.driver_traits, driver_gwas=L.driver_gwas,
            coloc_driver_is_estimation_cohort=est in L.driver_gwas.split(","),
        )
        if best is None:
            rows.append(dict(**base, has_shared_variant=False, alt_evaluable=False))
        else:
            rows.append(dict(**base, has_shared_variant=True, alt_evaluable=bool(best["sel_p"] < P_GWS), **best))
    res = pd.DataFrame(rows)
    res["p_T2_theta_ge_Lref"] = norm.cdf(
        (res.out_logor - L_ref * res.alt_sd) / np.sqrt(res.out_se**2 + L_ref**2 * res.alt_sd_se**2)
    )
    res["q_T1"], res["q_T2"], res["class"] = np.nan, np.nan, "not_evaluable"
    res.loc[res.mhc, "class"] = "mhc_not_classified"
    for lc in ("enzyme_only", "direct_involved"):
        fam = (res.locus_class == lc) & res.alt_evaluable & ~res.mhc
        res.loc[fam, "q_T1"] = bh(res.loc[fam, "out_p"])
        res.loc[fam, "q_T2"] = bh(res.loc[fam, "p_T2_theta_ge_Lref"])
        t1 = fam & (res.q_T1 < 0.05)
        res.loc[t1 & (res.theta >= L_ref), "class"] = "PNPLA3_class"
        res.loc[t1 & (res.theta > 0) & (res.theta < L_ref), "class"] = "intermediate"
        res.loc[t1 & (res.theta < 0), "class"] = "opposite_direction"
        res.loc[fam & ~t1 & (res.q_T2 < 0.05), "class"] = "no_cirrhosis_effect"
        res.loc[fam & ~t1 & ~(res.q_T2 < 0.05), "class"] = "uninformative"
    if not pc_pass:
        res["class_interpretable"] = False
    else:
        res["class_interpretable"] = True

    # --- negative control: GPT locus
    d = DATA[sel]
    gw = d[(d.chromosome == GPT_WIN[0]) & d.position.between(GPT_WIN[1], GPT_WIN[2])]
    gpt = [r for r in (ratio_row(v, sel, est, out) for v in gw.index) if r and r["maf"] >= 0.01]
    neg = None
    if gpt:
        neg = min(gpt, key=lambda x: x["sel_p"])
        neg["n_shared_variants_in_window"] = len(gpt)
        neg["evaluable"] = bool(neg["sel_p"] < P_GWS)
        neg["passes"] = bool(not (neg["out_p"] < 0.05 and neg["theta"] >= L_ref)) if neg["evaluable"] else None
    log("negative control (GPT):", neg)

    # --- aggregate slope with LD-block bootstrap
    blocks = pd.read_csv(f"{FM}/data/ld_ref/1kg_eur/approx_LD_blocks.txt", sep=r"\s+")
    blocks.columns = ["chr", "start", "stop"]
    blocks["chr"] = blocks.chr.astype(str).str.replace("chr", "").astype(int)

    def block_of(vid):
        c, p = int(vid.split(":")[0]), int(vid.split(":")[1])
        b = blocks[(blocks.chr == c) & (blocks.start <= p) & (blocks.stop > p)]
        return f"{c}:{int(b.start.iloc[0])}" if len(b) else f"{c}:none:{p // 1_000_000}"

    agg = []
    rng = np.random.default_rng(SEED)

    def slope_of(x):
        w = 1 / x.out_se.values**2
        return float(np.sum(w * x.alt_sd.values * x.out_logor.values) / np.sum(w * x.alt_sd.values**2))

    for label, sub in (
        ("enzyme_only_evaluable", res[(res.locus_class == "enzyme_only") & res.alt_evaluable & ~res.mhc]),
        ("enzyme_only_evaluable_excl_reference_loci",
         res[(res.locus_class == "enzyme_only") & res.alt_evaluable & ~res.mhc & ~res.is_reference_locus]),
        ("direct_involved_evaluable", res[(res.locus_class == "direct_involved") & res.alt_evaluable & ~res.mhc]),
        ("qualifying_positive_controls", qual),
    ):
        sub = sub.copy()
        if len(sub) < 2:
            agg.append(dict(set=label, n_loci=len(sub)))
            continue
        sub["block"] = [block_of(v) for v in sub.lead]
        est_slope = slope_of(sub)
        q_stat = float(np.sum((sub.out_logor - est_slope * sub.alt_sd) ** 2 / sub.out_se**2))
        ub = sub.block.unique()
        groups = {b: sub[sub.block == b] for b in ub}
        boots = np.empty(N_BOOT)
        for i in range(N_BOOT):
            pick = rng.choice(ub, size=len(ub), replace=True)
            boots[i] = slope_of(pd.concat([groups[b] for b in pick]))
        agg.append(dict(set=label, n_loci=len(sub), n_blocks=len(ub), ivw_slope=est_slope,
                        boot_lo=float(np.quantile(boots, 0.025)), boot_hi=float(np.quantile(boots, 0.975)),
                        cochran_Q=q_stat, Q_df=len(sub) - 1, n_boot=N_BOOT, seed=SEED))
    agg = pd.DataFrame(agg)

    res.to_csv(f"{OUT}/{tag}_locus_results.tsv", sep="\t", index=False)
    ctl.to_csv(f"{OUT}/{tag}_positive_controls.tsv", sep="\t", index=False)
    agg.to_csv(f"{OUT}/{tag}_aggregate_slope.tsv", sep="\t", index=False)
    pd.DataFrame([neg] if neg else []).to_csv(f"{OUT}/{tag}_negative_control_GPT.tsv", sep="\t", index=False)
    return dict(tag=tag, selection=sel, estimation=est, outcome=out, positive_control_pass=pc_pass,
                L_ref=L_ref, qualifying_controls=list(qual.gene),
                negative_control=None if neg is None else {k: (v if not isinstance(v, (np.floating, np.bool_)) else v.item()) for k, v in neg.items()},
                class_counts={lc: res[res.locus_class == lc]["class"].value_counts().to_dict()
                              for lc in ("enzyme_only", "direct_involved")}), res


summary = dict(
    label="provisional on the current COLOC release",
    coloc_release=REL,
    n_positive_genes=len(pos_genes),
    driver_class={"enzyme": n_enz_driver, "direct": n_dir_driver},
    support_class=pos_genes.support_class.value_counts().to_dict(),
    n_loci=len(loci),
    locus_class=loci.locus_class.value_counts().to_dict(),
    n_mhc_loci=int(loci.mhc.sum()),
    designs={},
)
RES = {}
for tag, sel, est, out in DESIGNS:
    summary["designs"][tag], RES[tag] = run_design(tag, sel, est, out)

# D1 vs D2 agreement on loci evaluable in both
m = RES["D1_cirrhosis"][["locus_id", "genes", "locus_class", "class", "theta", "alt_evaluable"]].merge(
    RES["D2_cirrhosis"][["locus_id", "class", "theta", "alt_evaluable"]], on="locus_id", suffixes=("_D1", "_D2")
)
both = m[m.alt_evaluable_D1 & m.alt_evaluable_D2 & (m.locus_class == "enzyme_only")]
pd.crosstab(both.class_D1, both.class_D2).to_csv(f"{OUT}/D1_vs_D2_class_crosstab_enzyme_only.tsv", sep="\t")
m.to_csv(f"{OUT}/D1_vs_D2_per_locus.tsv", sep="\t", index=False)
summary["D1_D2_enzyme_only_both_evaluable"] = int(len(both))
summary["D1_D2_same_class"] = int((both.class_D1 == both.class_D2).sum())
if len(both) > 2:
    summary["D1_D2_theta_spearman"] = float(both[["theta_D1", "theta_D2"]].corr(method="spearman").iloc[0, 1])

# ----------------------------------------------------------------------------
# 4. Calibration: null z correlation (MVP ALT vs MVP cirrhosis) and null rate
# ----------------------------------------------------------------------------
a, u, c = SAMP["MVP_ALT_EUR"], SAMP["UKBB_ALT"], SAMP["MVP_Cirrhosis_EUR"]
ids = a.index.intersection(u.index).intersection(c.index)
a, u, c = a.loc[ids], u.loc[ids], c.loc[ids]
common = (np.minimum(a.af, 1 - a.af) >= 0.05).values
ids = ids[common]
rng = np.random.default_rng(SEED)
if len(ids) > 200_000:
    ids = pd.Index(rng.choice(ids.values, 200_000, replace=False))
a, u, c = a.loc[ids], u.loc[ids], c.loc[ids]
zc = np.where(a.a1.values == c.a1.values, c.z.values, -c.z.values)
ok = (np.abs(a.z.values) < 4) & (np.abs(zc) < 4)
null = (2 * norm.sf(np.abs(a.z.values)) > 0.5) & (2 * norm.sf(np.abs(u.z.values)) > 0.5)
cal = dict(
    n_sampled_common_variants=int(len(ids)),
    n_used_for_correlation=int(ok.sum()),
    null_z_correlation_MVP_ALT_vs_MVP_cirrhosis=float(np.corrcoef(a.z.values[ok], zc[ok])[0, 1]),
    n_ALT_null_variants=int(null.sum()),
    fraction_cirrhosis_p_lt_0p05_among_ALT_null=float(np.mean(2 * norm.sf(np.abs(zc[null])) < 0.05)),
    cirrhosis_lambda_median_among_ALT_null=float(np.median(zc[null] ** 2) / 0.4549364),
    note="variants thinned by position % 23 == 7, not LD-pruned",
)
summary["calibration_D1"] = cal
pd.DataFrame([cal]).to_csv(f"{OUT}/calibration_D1.tsv", sep="\t", index=False)

with open(f"{OUT}/summary.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)
log(json.dumps(summary, indent=2, default=str))
log("DONE")
