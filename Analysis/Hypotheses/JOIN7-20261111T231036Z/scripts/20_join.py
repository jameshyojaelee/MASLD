#!/usr/bin/env python3
"""JOIN7 step 2: genetic classes per gene and trait, mechanism x endpoint table, exact 2x2 tests.

Rules: PRESPEC_JOIN7.md. Usage: 20_join.py <project_root> <out_dir>
"""
import hashlib
import re
import sys
from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

T1 = "Analysis/Hypotheses/TRIAL1-20261110T205624Z"
HASH = {"trial_outcomes_v2.tsv": "0ff87dc41c074291f56b487ef226fa6d985f9a2b7237ef6885fb45085d20e2c5",
        "excluded_v2.tsv": "e40fa076297d87203124c988c60fbd8039c37e7a70ea90bf644a794cb1108c08",
        "CHANGES_v2.tsv": "d0cf36753977d5d353f2e13c9b9c4e2e059f1e36e737df05e4b20b577eea86eb"}
DIR6 = "Analysis/Hypotheses/DIR6-20261111T024808Z/out_23210381"
CATALOG = ("RNA-seq/results/evidence_passports/candidates/gene-catalog-r8-anchor-disease-sidecar-"
           "20261025T131731Z/catalog_r8_anchor_disease_liability.tsv")
FIG6 = "figures/main/fig6_gene_catalog/panels/fig6_drug_target_calibration_values.tsv"
TRAITS = ["pdff", "cirr", "alt_ukb", "alt_mvp"]
SUPPORT = {"directional_agree", "gene_level_support", "directional_disagree", "mixed_direction"}
RANK = {"directional_agree": 4, "gene_level_support": 3, "directional_disagree": 2, "no_support": 1,
        "untestable": 0}
IMPROVED = {"improved", "improved_unadjusted"}
NOT_IMPROVED = {"not_improved", "worsened"}
PRIMARY = [("pdff", "pdff"), ("cirr", "fibrosis_improvement"), ("alt_ukb", "alt")]
DESCRIPTIVE = [("alt_ukb", "fibrosis_improvement"), ("alt_ukb", "mash_resolution"),
               ("cirr", "mash_resolution")]
GENE_RE = re.compile(r"^[A-Z0-9]+$")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def role(note):
    n = str(note).lower()
    if "secondary" in n or "exploratory" in n:
        return "secondary_or_exploratory"
    if "primary" in n or "main efficacy" in n:
        return "primary"
    return "not stated"


def gene_classes(root, out):
    g = pd.read_csv(f"{out}/gene_trait_raw.tsv", sep="\t")
    win = pd.read_csv(f"{out}/gene_windows_grch37.tsv", sep="\t")
    test = set(TEST)
    genes = sorted(test | set(CONTEXT))
    full = pd.DataFrame(list(product(genes, TRAITS)), columns=["gene", "trait"])
    g = full.merge(g, on=["gene", "trait"], how="left")
    g["in_gene_file"] = g["gene"].isin(win["gene"])
    m = {t: int(g[(g.trait == t) & g.gene.isin(test)]["magma_p"].notna().sum()) for t in ("pdff", "cirr")}
    g["magma_m"] = g["trait"].map(m)
    g["magma_thr"] = 0.05 / g["magma_m"]
    g["magma_sig"] = g["magma_p"] < g["magma_thr"]
    g["var_sig"] = g["n_gws"].fillna(0) > 0
    g["scored"] = g["magma_p"].notna() | (g["n_var"].fillna(0) > 0)
    g["support"] = g["magma_sig"] | g["var_sig"]

    d6 = pd.read_csv(f"{root}/{DIR6}/dir6_gene_classes.tsv", sep="\t", dtype=str)
    d6 = d6[["symbol", "call", "defining_trait", "call_pdff", "call_cirr", "combined",
             "primary_evidence_class"]].drop_duplicates()
    cal = pd.read_csv(f"{root}/{DIR6}/calibration.tsv", sep="\t", dtype=str)

    def mapped(gene, trait):
        r = d6[d6.symbol == gene]
        for x in r.itertuples():
            if trait == "pdff" and x.call_pdff == "concordant":
                return x.call, "DIR6 H26 call, concordant on liver fat"
            if trait == "cirr" and x.call_cirr == "concordant":
                return x.call, "DIR6 H26 call, concordant on cirrhosis"
            if trait.startswith("alt") and x.defining_trait == "ALT":
                return x.call, "DIR6 H26 call, ALT-defined"
        if r.empty:
            c = cal[cal.symbol == gene]
            txt = "" if c.empty else str(c["published_coding_mechanism"].iloc[0]).lower()
            if ("loss of function" in txt or "loss-of-function" in txt) and "protect" in txt:
                return "inhibit", "DIR6 calibration text"
            if "gain of function" in txt or "gain-of-function" in txt:
                return "inhibit", "DIR6 calibration text"
        return None, "not mapped"

    mp = [mapped(a, b) for a, b in zip(g.gene, g.trait)]
    g["mapped_direction"] = [x[0] for x in mp]
    g["direction_source"] = [x[1] for x in mp]
    g["h26_call"] = g["gene"].map(d6.drop_duplicates("symbol").set_index("symbol")["call"])
    g["dir6_combined"] = g["gene"].map(d6.drop_duplicates("symbol").set_index("symbol")["combined"])
    g["gene_class"] = np.where(~g["scored"], "untestable", np.where(g["support"], "support", "no_support"))
    cat = pd.read_csv(f"{root}/{CATALOG}", sep="\t", dtype=str)
    cat = cat.drop_duplicates("symbol").set_index("symbol")
    g["catalog_r8_member"] = g["gene"].isin(cat.index)
    g["catalog_primary_evidence_class"] = g["gene"].map(cat["primary_evidence_class"])
    g["catalog_sidecar_class"] = g["gene"].map(cat["sidecar_class"])
    g["role"] = np.where(g["gene"].isin(test), "test", "context")
    g.to_csv(f"{out}/gene_trait_classes.tsv", sep="\t", index=False)
    return g, m


def mech_class(genes, action, gt, trait):
    if not genes:
        return "untestable", "no target"
    cls, parts = [], []
    act = str(action).split()[0].lower()
    for gene in genes:
        r = gt[(gt.gene == gene) & (gt.trait == trait)].iloc[0]
        if r.gene_class != "support":
            c = r.gene_class
        elif r.mapped_direction is None or pd.isna(r.mapped_direction):
            c = "gene_level_support"
        else:
            c = "directional_agree" if r.mapped_direction == act else "directional_disagree"
        cls.append(c)
        parts.append(f"{gene}:{c}")
    if "directional_agree" in cls and "directional_disagree" in cls:
        best = "mixed_direction"
    else:
        best = max(cls, key=lambda c: RANK[c])
    return best, ";".join(parts)


def min_attainable_p(r1, r2, c1):
    n = r1 + r2
    ps = []
    for a in range(max(0, c1 - r2), min(r1, c1) + 1):
        t = [[a, r1 - a], [c1 - a, r2 - c1 + a]]
        ps.append(fisher_exact(t).pvalue)
    return min(ps) if ps else np.nan


def two_by_two(df, gcol, endpoint, label, sub=None):
    d = df[df.endpoint == endpoint]
    if sub is not None:
        d = d[sub(d)]
    n_untest = int((d[gcol] == "untestable").sum())
    n_norep = int((~d.result_vs_control.isin(IMPROVED | NOT_IMPROVED)).sum())
    d = d[(d[gcol] != "untestable") & d.result_vs_control.isin(IMPROVED | NOT_IMPROVED)]
    s = d[gcol].isin(SUPPORT)
    imp = d.result_vs_control.isin(IMPROVED)
    a, b = int((s & imp).sum()), int((s & ~imp).sum())
    c, e = int((~s & imp).sum()), int((~s & ~imp).sum())
    t = [[a, b], [c, e]]
    row = dict(analysis=label, genetics=gcol.replace("class_", ""), endpoint=endpoint, n=a + b + c + e,
               support_improved=a, support_not=b, nosupport_improved=c, nosupport_not=e,
               excluded_untestable=n_untest, excluded_not_reported=n_norep)
    if min(a + b, c + e, a + c, b + e) == 0:
        row.update(or_cond=np.nan, ci_lo=np.nan, ci_hi=np.nan, p=1.0, min_attainable_p=1.0,
                   note="a margin is zero; no test possible")
        return row
    o = odds_ratio(t, kind="conditional")
    ci = o.confidence_interval(0.95)
    row.update(or_cond=o.statistic, ci_lo=ci.low, ci_hi=ci.high, p=fisher_exact(t).pvalue,
               min_attainable_p=min_attainable_p(a + b, c + e, a + c),
               note="CI spans 1: not informative at this number of mechanisms" if ci.low <= 1 <= ci.high
               else "CI excludes 1")
    row["mechanisms_support"] = ";".join(sorted(d[s].mechanism))
    return row


def main(root, out):
    for f, h in HASH.items():
        got = sha(f"{root}/{T1}/{f}")
        assert got == h, f"{f} changed: {got}"
    gt, m = gene_classes(root, out)

    tr = pd.read_csv(f"{root}/{T1}/trial_outcomes_v2.tsv", sep="\t", dtype=str)
    tr["row"] = np.arange(1, len(tr) + 1)
    tr["not_verified"] = tr.apply(lambda r: any("not verified" in str(v) for v in r.values), axis=1)
    rep = tr[tr.representative == "yes"].copy()

    targets = {}
    for mech, d in tr.groupby("mechanism", sort=False):
        toks = lambda s: [x.strip() for x in str(s).split(";") if GENE_RE.match(x.strip())]
        gs = sorted({x for s in rep[rep.mechanism == mech].target_genes for x in toks(s)})
        if not gs:
            gs = sorted({x for s in d.target_genes for x in toks(s)})
        targets[mech] = gs
    assert set().union(*targets.values()) == set(TEST), sorted(set().union(*targets.values()) ^ set(TEST))
    rep["targets_used"] = rep.mechanism.map(lambda k: ";".join(targets[k]))
    rep["n_targets"] = rep.mechanism.map(lambda k: len(targets[k]))
    rep["endpoint_role"] = rep.notes.map(role)
    alldrugs = tr.groupby("mechanism").apply(
        lambda d: "; ".join(sorted({f"{a} ({b})" for a, b in zip(d.drug, d.nct_id)})), include_groups=False)
    rep["all_drugs_nct"] = rep.mechanism.map(alldrugs)
    for t in TRAITS:
        res = [mech_class(targets[k], a, gt, t) for k, a in zip(rep.mechanism, rep.action_direction)]
        rep[f"class_{t}"] = [x[0] for x in res]
        rep[f"genes_{t}"] = [x[1] for x in res]
    rep["outcome_class"] = rep.result_vs_control.map(
        lambda r: "improved" if r in IMPROVED else "not improved" if r in NOT_IMPROVED else "not reported")
    keep = ["row", "mechanism", "drug", "nct_id", "phase", "endpoint", "endpoint_role", "result_vs_control",
            "outcome_class", "not_verified", "action_direction", "target_genes", "targets_used", "n_targets",
            "class_pdff", "class_cirr", "class_alt_ukb", "class_alt_mvp", "genes_pdff", "genes_cirr",
            "genes_alt_ukb", "genes_alt_mvp", "all_drugs_nct", "effect_size", "ci_or_p", "notes"]
    rep[keep].to_csv(f"{out}/mechanism_endpoint_table.tsv", sep="\t", index=False)

    rows = []
    for t, e in PRIMARY:
        gc = f"class_{t}"
        rows.append(two_by_two(rep, gc, e, "primary"))
        rows.append(two_by_two(rep, gc, e, "S1_drop_unadjusted",
                               lambda d: d.result_vs_control != "improved_unadjusted"))
        rows.append(two_by_two(rep, gc, e, "S2_drop_not_verified", lambda d: ~d.not_verified))
        rows.append(two_by_two(rep, gc, e, "S3_single_target", lambda d: d.n_targets == 1))
        r4 = rep.copy()
        r4[gc] = r4[gc].replace({"directional_disagree": "no_support"})
        rows.append(two_by_two(r4, gc, e, "S4_disagree_as_no"))
    rows.append(two_by_two(rep, "class_alt_mvp", "alt", "S5_mvp_alt"))
    for t, e in DESCRIPTIVE:
        rows.append(two_by_two(rep, f"class_{t}", e, "descriptive"))
    st = pd.DataFrame(rows)
    prim = st.analysis == "primary"
    p = st.loc[prim, "p"].values
    order = np.argsort(p)
    holm = np.empty(len(p))
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - k) * p[i]))
        holm[i] = run
    st.loc[prim, "p_holm"] = holm
    st.to_csv(f"{out}/two_by_two.tsv", sep="\t", index=False)
    print(st.drop(columns=["mechanisms_support"], errors="ignore").to_string(index=False))

    # controls
    ctl = gt[(gt.trait == "pdff") & gt.gene.isin(["PNPLA3", "TM6SF2", "GCKR"])]
    ctl[["gene", "magma_p", "magma_thr", "magma_sig", "n_gws", "min_p", "support"]].to_csv(
        f"{out}/control_genes_pdff.tsv", sep="\t", index=False)
    print(ctl[["gene", "magma_p", "magma_thr", "n_gws", "min_p", "support"]].to_string(index=False))
    print("MAGMA m per trait:", m)

    # old Figure 6 drug-target predictor, descriptive
    f6 = pd.read_csv(f"{root}/{FIG6}", sep="\t", dtype=str)
    f6 = f6[["gene", "drug_dev_status", "outcome", "max_phase_masld", "coloc_susie_pp4", "log2FC", "padj",
             "is_deg"]]
    wide = gt.pivot(index="gene", columns="trait", values="gene_class").add_prefix("class_").reset_index()
    tg = pd.DataFrame([(g, k) for k, gs in targets.items() for g in gs], columns=["gene", "mechanism"])
    oc = rep.pivot_table(index="mechanism", columns="endpoint", values="outcome_class",
                         aggfunc="first").add_prefix("trial_").reset_index()
    cmp_ = f6.merge(wide, on="gene", how="left").merge(tg, on="gene", how="left").merge(oc, on="mechanism",
                                                                                        how="left")
    cmp_.to_csv(f"{out}/fig6_predictor_comparison.tsv", sep="\t", index=False)
    if not ctl["support"].all():
        sys.exit("positive control failed")


CONTEXT = "TM6SF2 GCKR HSD17B13 CIDEB MTARC1 INHBE GPAM MOGAT2 NR3C1 RORA CNR1".split()
TEST = ("THRB GLP1R GIPR GCGR FGFR1 FGFR2 FGFR3 PPARA PPARD PPARG NR1H4 ACACA ACACB DGAT2 SCD FASN "
        "MAP3K5 CASP1 CASP2 CASP3 CASP4 CASP5 CASP6 CASP7 CASP8 CASP9 CASP10 CASP12 CASP14 CCR2 CCR5 "
        "LGALS3 LOXL2 SLC5A2 MPC1 MPC2 KHK FGF1 FGFR4 KLB ADORA3 AOC3 PNPLA3").split()

if __name__ == "__main__":
    np.random.seed(20261111)
    main(sys.argv[1], sys.argv[2])
