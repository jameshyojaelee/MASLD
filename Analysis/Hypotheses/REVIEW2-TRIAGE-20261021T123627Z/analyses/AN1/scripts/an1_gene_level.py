#!/usr/bin/env python3
"""AN1 gene-level sensitivity for the symmetric shared-posterior check (REVIEW2 B-G1).

Inputs: the r2 master (susie_coloc_all_gwas.csv) and the AN1 share-only row table(s).
A primary not_rerun_negative row whose every eQTL credible set keeps < 0.5 of its
posterior on the GWAS-fit variants becomes untestable (PP.H4/PP.H3 set to NA,
method = susie_untestable_insufficient_shared_posterior), exactly as the r2 rerun
marks such rows. Rows where only some sets fail keep their first-run PP.H4 (the
kept-pair maximum cannot exceed it and it is <= 0.5, so the call is unchanged);
their PP.H3 is unknown and is bounded.

Gene level, same code as the release:
  rule (a) adopted: rebuild_tier12.py build() + susie_state()
           positive > untestable > tested_le_0.5 > no_susie_result
  rule (b) PLAN D1: tested if any tier-1/2 row is eligible (finite SuSiE PP.H4)
Recomputed per rule: H57 primary statistic (functions copied verbatim from
Analysis/Hypotheses/H57-20261015T112421Z/scripts/01_h57_overlap.py), Fig 3G class
counts, the 3,918/3,335/3,319 PP.H3 counts, and Catalog r9d genes moving from tested
to untestable. Each baseline is reproduced first; the script stops if one does not.
Usage: an1_gene_level.py --out <new dir> --shares <row_shares.tsv> [<row_shares.tsv> ...]
"""
import argparse, collections, csv, hashlib, importlib.util, io, json, os, subprocess, sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
R2 = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z"
MASTER = R2 + "/susie_coloc_all_gwas.csv"
T12 = R2 + "/gene_level_coloc_tier12.csv"
BULK = ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
FIG3G = ROOT + "/figures/candidates/fig3g-working-layout-774-20261012T011006Z/source_tables/fig3g_raw_treat_interface.tsv"
CAT = ROOT + "/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet"
REB = ROOT + "/GWAS/finemapping/src/perf/rebuild_tier12.py"
TOL = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--shares", nargs="+", required=True)
args = ap.parse_args()
if os.path.lexists(args.out):
    sys.exit(f"refusing to reuse {args.out}")
os.makedirs(args.out)
OUT = args.out
np.random.seed(20261021)
LOG = open(f"{OUT}/log.txt", "w")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n"); LOG.flush()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stop(msg):
    log("STOP:", msg)
    sys.exit(1)


with open(f"{OUT}/input_sha256.tsv", "w") as fh:
    for p in [MASTER, T12, BULK, FIG3G, CAT, REB, __file__] + args.shares:
        fh.write(f"{p}\t{sha(p)}\n")

spec = importlib.util.spec_from_file_location("rebuild_tier12", REB)
rt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rt)
tier, trait, placement, label = rt.load_tiers()
UNTESTABLE = rt.UNTESTABLE

# ------------------------------------------------------------------ master and shares
COLS = ["gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H3.susie", "PP.H4.susie", "method",
        "eligibility_source"]
m = pd.read_csv(MASTER, usecols=COLS, dtype=str, keep_default_na=False)
m["tier"] = m.gwas_name.map(tier)
log("master rows", len(m))
sh = pd.concat([pd.read_csv(p, sep="\t", dtype={"ensembl": str, "gwas_name": str}) for p in args.shares],
               ignore_index=True)
if sh.duplicated(["gwas_name", "ensembl"]).any():
    stop("duplicate (study, gene) rows in the share tables")
log("share rows", len(sh), "| status", sh.status.value_counts(dropna=False).to_dict())
key = ["gwas_name", "ensembl"]
mk = m[key].copy(); mk["_i"] = np.arange(len(m))
j = sh.merge(mk, on=key, how="left")
if j._i.isna().any():
    stop(f"{int(j._i.isna().sum())} share rows not in the master")
if m.duplicated(key).any():
    dup = m[m.duplicated(key, keep=False)].merge(sh[key], on=key)
    if len(dup):
        stop(f"share rows map to duplicated master keys: {len(dup)}")
j["_i"] = j._i.astype(int)
assert (m.loc[j._i, "eligibility_source"].values == "not_rerun_negative").all()

j["tier_cfg"] = j.gwas_name.map(tier)
is_prim = j.tier_cfg.isin([1, 2])
if "primary" in j.columns and not (j.primary.astype(str).isin(["True", "TRUE", "1"]) == is_prim).all():
    stop("primary flag differs from tier 1/2")
prim = j[is_prim]
log("primary (tier-1/2) share rows", len(prim), "| status", prim.status.value_counts(dropna=False).to_dict())
n_full_fail = int((prim.status == "untestable_all_pairs_fail").sum())
n_partial = int((prim.status == "eligible_some_pairs_fail").sum())
n_ok = int((prim.status == "eligible_all_pairs").sum())
n_other = len(prim) - n_full_fail - n_partial - n_ok
log(f"primary rows: all sets fail {n_full_fail} | some fail {n_partial} | none fail {n_ok} | not reached {n_other}")

# row-level consistency of the rebuild with the first run (all share rows)
reached = j[j.status.isin(["untestable_all_pairs_fail", "eligible_some_pairs_fail", "eligible_all_pairs"])]
DIAG = {
    "share_rows": int(len(j)),
    "status_all_rows": j.status.value_counts(dropna=False).to_dict(),
    "status_non_primary_rows": j[~is_prim].status.value_counts(dropna=False).to_dict(),
    "n_snps_rebuilt_ne_master": int((j.n_snps_rebuilt != j.n_snps).sum()),
    "n_snps_rebuilt_missing": int(j.n_snps_rebuilt.isna().sum()),
    "n_cs_pairs_not_multiple_of_n_eqtl_cs": int((reached.n_cs_pairs % reached.n_eqtl_cs != 0).sum()),
    "gwas_ids_not_subset_of_eqtl": int((reached.gwas_ids_subset_of_eqtl.astype(str).str.upper() != "TRUE").sum()),
    "ld_route": j.ld_route.value_counts(dropna=False).to_dict(),
    "primary_status_by_ancestry": {f"{a}|{b}": int(n) for (a, b), n in
                                   prim.groupby(["ancestry.x" if "ancestry.x" in prim.columns else "ancestry",
                                                 "status"]).size().items()},
    "primary_genes_with_all_sets_fail_row": int(prim[prim.status == "untestable_all_pairs_fail"].ensembl.nunique()),
    "primary_genes": int(prim.ensembl.nunique()),
}
fin_all = pd.to_numeric(m["PP.H4.susie"], errors="coerce").notna()
DIAG["finite_susie_rows_master_all"] = int(fin_all.sum())
DIAG["finite_susie_rows_master_tier12"] = int((fin_all & m.tier.isin([1, 2])).sum())
DIAG["finite_susie_rows_after_all"] = DIAG["finite_susie_rows_master_all"] - \
    int((j.status == "untestable_all_pairs_fail").sum())
DIAG["finite_susie_rows_after_tier12"] = DIAG["finite_susie_rows_master_tier12"] - n_full_fail
for p_ in args.shares:
    idf = os.path.join(os.path.dirname(p_), "identity_check.tsv")
    if os.path.exists(idf):
        t = pd.read_csv(idf, sep="\t")
        DIAG.setdefault("identity_checks", []).append(
            dict(dir=os.path.basename(os.path.dirname(p_)), n=int(len(t)),
                 identical=int(t.identical.astype(str).str.upper().eq("TRUE").sum())))
log("row diagnostics", json.dumps(DIAG, default=str))

# sensitivity master: fully failing primary rows -> untestable
ms = m.copy()
fail_i = prim.loc[prim.status == "untestable_all_pairs_fail", "_i"].to_numpy()
ms.loc[fail_i, "PP.H4.susie"] = ""
ms.loc[fail_i, "PP.H3.susie"] = ""
ms.loc[fail_i, "method"] = UNTESTABLE
partial_i = set(prim.loc[prim.status == "eligible_some_pairs_fail", "_i"])
unreached_i = set(prim.loc[~prim.status.isin(["untestable_all_pairs_fail", "eligible_some_pairs_fail",
                                                "eligible_all_pairs"]), "_i"])

red_cols = ["gwas_name", "gene", "ensembl", "chr", "PP.H4.abf", "PP.H4.susie", "method"]
m[red_cols].to_csv(f"{OUT}/master_reduced_baseline.csv", index=False)
ms[red_cols].to_csv(f"{OUT}/master_reduced_sensitivity.csv", index=False)

# ------------------------------------------------------------------ gene level
def state_b(e):
    """PLAN D1: tested if any tier-1/2 row has a finite (eligible) SuSiE PP.H4."""
    if e["s"] is not None and e["s"] > 0.5:
        return "positive"
    if e["s"] is not None:
        return "tested_le_0.5"
    return "untestable" if e["u"] else "no_susie_result"


def gene_table(path):
    best, _ = rt.build(path, tier, trait, placement)
    rows = [dict(ensembl=g, coloc_best_susie_pp4=e["s"], coloc_best_abf_pp4=e["a"], susie_driving_gwas=e["s_g"],
                 susie_untestable_t12_studies=e["u"], state_a=rt.susie_state(e), state_b=state_b(e))
            for g, e in best.items()]
    df = pd.DataFrame(rows)
    for c in ["coloc_best_susie_pp4", "coloc_best_abf_pp4"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


gb = gene_table(f"{OUT}/master_reduced_baseline.csv")
gb_direct = gene_table(MASTER)
t12 = pd.read_csv(T12, float_precision="round_trip")   # exact Python float parse
cmp = t12.merge(gb_direct, on="ensembl", how="outer", suffixes=("_r2", ""), indicator=True)
ok_genes = (cmp._merge == "both").all()
s1 = pd.to_numeric(cmp.coloc_best_susie_pp4_r2, errors="coerce"); s2 = pd.to_numeric(cmp.coloc_best_susie_pp4)
ok_s = ((s1 == s2) | (s1.isna() & s2.isna())).all()
ok_state = (cmp.susie_state_t12 == cmp.state_a).all()
ok_u = (cmp.susie_untestable_t12_studies_r2 == cmp.susie_untestable_t12_studies).all()
same_reduced = gb.equals(gb_direct)
log(f"baseline gene table vs r2 gene_level_coloc_tier12.csv: genes match {ok_genes}, best SuSiE PP.H4 identical "
    f"{ok_s}, susie_state_t12 identical {ok_state}, untestable-study counts identical {ok_u}; reduced master gives "
    f"the same table {same_reduced}")
if not (ok_genes and ok_s and ok_state and ok_u and same_reduced):
    stop("baseline gene-level table does not reproduce r2")
gs = gene_table(f"{OUT}/master_reduced_sensitivity.csv")
G = gb.merge(gs, on="ensembl", suffixes=("_base", "_sens"))
assert len(G) == len(gb) == len(gs)
RES = {"diagnostics": DIAG, "rows": dict(primary=len(prim), all_sets_fail=n_full_fail, some_sets_fail=n_partial, none_fail=n_ok,
                    not_reached=n_other, status=prim.status.value_counts(dropna=False).to_dict())}
for r in ["a", "b"]:
    RES[f"states_rule_{r}"] = {"baseline": G[f"state_{r}_base"].value_counts().to_dict(),
                               "sensitivity": G[f"state_{r}_sens"].value_counts().to_dict(),
                               "transitions": {f"{a} -> {b}": int(n) for (a, b), n in
                                               G.groupby([f"state_{r}_base", f"state_{r}_sens"]).size().items()}}
    log(f"rule ({r}) states", json.dumps(RES[f"states_rule_{r}"]))
G.to_csv(f"{OUT}/gene_states.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ H57 (functions verbatim from 01_h57_overlap.py)
def make_bins(df, col, k):
    """df must already be sorted by Ensembl ID; ties in col broken by that order."""
    return pd.qcut(df[col].rank(method="first"), k, labels=False).to_numpy()


def bin_margins(bins, gen, deg):
    ks = np.unique(bins)
    n = np.array([(bins == k).sum() for k in ks])
    Gk = np.array([gen[bins == k].sum() for k in ks])
    Dk = np.array([deg[bins == k].sum() for k in ks])
    Ok = np.array([(gen & deg)[bins == k].sum() for k in ks])
    return n, Gk, Dk, Ok


def exact_null(n, Gk, Dk):
    """pmf of the total overlap: convolution of per-bin hypergeometrics (pop n_k, D_k successes, G_k draws)."""
    pmf = np.array([1.0])
    for nk, gk, dk in zip(n, Gk, Dk):
        x = np.arange(0, min(gk, dk) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, nk, dk, gk))
    return pmf


def exact_p(pmf, E, O):
    x = np.arange(len(pmf))
    return float(min(1.0, pmf[np.abs(x - E) >= abs(O - E) - TOL].sum()))


def mh(n, Gk, Dk, Ok):
    tabs = [np.array([[o, g - o], [d - o, nn - g - d + o]], dtype=float) for nn, g, d, o in zip(n, Gk, Dk, Ok)]
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint(alpha=0.05)
    cmh = st.test_null_odds(correction=False)
    bd = st.test_equal_odds(adjust=False)
    return dict(MH_OR=float(st.oddsratio_pooled), MH_lo=float(lo), MH_hi=float(hi),
                CMH_chi2=float(cmh.statistic), CMH_p=float(cmh.pvalue), BD_chi2=float(bd.statistic), BD_p=float(bd.pvalue))


def strat_test(bins, gen, deg, with_mh=True):
    n, Gk, Dk, Ok = bin_margins(bins, gen, deg)
    E = float((Gk * Dk / n).sum())
    O = int(Ok.sum())
    pmf = exact_null(n, Gk, Dk)
    r = dict(n_strata=len(n), O=O, E_strat=E, O_over_E=O / E if E > 0 else np.nan, p_exact=exact_p(pmf, E, O),
             null_mean=float((np.arange(len(pmf)) * pmf).sum()))
    if with_mh:
        r.update(mh(n, Gk, Dk, Ok))
    return r


def unstrat(gen, deg):
    N = len(gen); G_ = int(gen.sum()); D = int(deg.sum()); O = int((gen & deg).sum())
    orr, p = stats.fisher_exact([[O, G_ - O], [D - O, N - G_ - D + O]])
    return dict(N=N, G=G_, D=D, O=O, E=G_ * D / N, O_over_E=O / (G_ * D / N), fisher_OR=float(orr), fisher_p=float(p))


def is_deg(df):
    return (df["padj"] < 0.05) & (df["logFC"].abs() > 0.5)


bulk = pd.read_csv(BULK)
bulk["g"] = bulk["gene"].astype(str).str.split(".").str[0]


def h57(gt, state_col):
    # H57 read gene_level_coloc_tier12.csv with pandas' default float parser; pass the
    # table through the same CSV text (Python repr floats, as rebuild_tier12 writes it)
    # and the same parser so every PP.H4 carries the same last-digit rounding.
    coloc = pd.read_csv(io.StringIO(gt.rename(columns={state_col: "susie_state_t12"}).to_csv(index=False)))
    coloc["g"] = coloc["ensembl"].astype(str).str.split(".").str[0]
    coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
    jt = bulk.drop_duplicates("g").merge(
        coloc[["g", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_state_t12", "susie_driving_gwas"]],
        on="g", how="inner")
    jt["susie"] = pd.to_numeric(jt.coloc_best_susie_pp4, errors="coerce")
    jt["gen"] = (jt.susie >= 0.5).fillna(False).astype(bool)
    jt["deg"] = is_deg(jt).astype(bool)
    full = jt.sort_values("g").reset_index(drop=True)
    ue = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
    gen, deg = ue.gen.to_numpy(), ue.deg.to_numpy()
    b10 = make_bins(ue, "AveExpr", 10)
    _, Gk, Dk, _ = bin_margins(b10, gen, deg)
    K = 10 if int(((Gk >= 5) & (Dk >= 5)).sum()) >= 7 else 5
    bins = b10 if K == 10 else make_bins(ue, "AveExpr", K)
    r = strat_test(bins, gen, deg)
    r.update({f"unstrat_{k}": v for k, v in unstrat(gen, deg).items()})
    r.update(N=len(ue), G=int(gen.sum()), D=int(deg.sum()), K=K)
    # descriptive, not tested: Spearman rho of best PP.H4 with |log2FC| (fig2.md:23, fig3.md:45),
    # on U_eval and on every bulk-joined gene with a finite best PP.H4 (the Fig 3G set)
    rho, prho = stats.spearmanr(ue.susie, ue.logFC.abs())
    f7 = full[full.susie.notna()]
    rho7, prho7 = stats.spearmanr(f7.susie, f7.logFC.abs())
    r.update(rho_ueval=float(rho), rho_ueval_p=float(prho), n_finite=len(f7), rho_finite=float(rho7),
             rho_finite_p=float(prho7), n_joint=len(full),
             n_joint_untestable=int((full.susie_state_t12 == "untestable").sum()),
             n_joint_no_susie_result=int((full.susie_state_t12 == "no_susie_result").sum()),
             n_joint_untestable_deg=int(((full.susie_state_t12 == "untestable") & full.deg).sum()))
    return r, ue


H = {}
for r in ["a", "b"]:
    for ph in ["base", "sens"]:
        res, ue = h57(G.rename(columns={"coloc_best_susie_pp4_" + ph: "coloc_best_susie_pp4",
                                        "coloc_best_abf_pp4_" + ph: "coloc_best_abf_pp4",
                                        "susie_driving_gwas_" + ph: "susie_driving_gwas"}),
                      f"state_{r}_{ph}")
        H[f"{r}_{ph}"] = res
        ue[["g", "gene", "susie", "gen", "deg", "AveExpr"]].to_csv(f"{OUT}/h57_universe_rule_{r}_{ph}.tsv", sep="\t",
                                                                      index=False)
        log(f"H57 rule ({r}) {ph}: N {res['N']} G {res['G']} D {res['D']} O {res['O']} E {res['E_strat']:.4f} "
            f"MH OR {res['MH_OR']:.4f} ({res['MH_lo']:.4f}-{res['MH_hi']:.4f}) exact p {res['p_exact']:.4f} "
            f"K {res['K']} | unstrat E {res['unstrat_E']:.4f} Fisher OR {res['unstrat_fisher_OR']:.4f} "
            f"p {res['unstrat_fisher_p']:.4f} | rho U_eval {res['rho_ueval']:.4f} (p {res['rho_ueval_p']:.3f}); "
            f"finite {res['n_finite']} rho {res['rho_finite']:.4f} (p {res['rho_finite_p']:.3f}) | joint {res['n_joint']} "
            f"untestable {res['n_joint_untestable']} (DEG {res['n_joint_untestable_deg']}) "
            f"no_susie_result {res['n_joint_no_susie_result']}")
a0 = H["a_base"]
ok_a = (a0["N"], a0["G"], a0["D"], a0["O"]) == (3769, 426, 220, 26) and round(a0["E_strat"], 3) == 21.213 and \
    round(a0["MH_OR"], 4) == 1.2947 and round(a0["p_exact"], 4) == 0.2755
b0 = H["b_base"]
ok_b = (b0["N"], b0["O"]) == (3776, 26) and round(b0["E_strat"], 3) == 21.190 and round(b0["MH_OR"], 3) == 1.297
log(f"H57 baseline reproduces RESULT_H57 T1 (rule a): {ok_a}; D4 U_eval7 (rule b): {ok_b}")
if not (ok_a and ok_b):
    stop("H57 baseline does not reproduce")
RES["H57"] = H

# ------------------------------------------------------------------ Fig 3G
f3 = pd.read_csv(FIG3G, sep="\t")
f3["g"] = f3.gene_id.astype(str).str.split(".").str[0]
gg = G.assign(g=G.ensembl.str.split(".").str[0]).sort_values(
    ["coloc_best_susie_pp4_base", "coloc_best_abf_pp4_base"], ascending=False).drop_duplicates("g")
f3 = f3.merge(gg[["g", "coloc_best_susie_pp4_base", "coloc_best_susie_pp4_sens", "state_a_base", "state_a_sens",
                  "state_b_base", "state_b_sens"]], on="g", how="left")
ok_f3 = bool(np.allclose(f3.susie_pp4, f3.coloc_best_susie_pp4_base, rtol=0, atol=1e-12))


def f3_state(pp4, deg):
    gen = np.isfinite(pp4) & (pp4 >= 0.5)
    return np.where(gen & deg, "Both", np.where(gen, "Genetics sig.", np.where(deg, "Transcriptomics sig.", "Neither")))


deg3 = f3.canonical_deg.astype(str).str.upper().eq("TRUE").to_numpy()
f3["state_base_re"] = f3_state(f3.coloc_best_susie_pp4_base.to_numpy(float), deg3)
ok_f3 = ok_f3 and (f3.state_base_re == f3.state).all()
log(f"Fig 3G baseline: susie_pp4 equals the rebuilt best PP.H4 and states reproduce: {ok_f3}")
if not ok_f3:
    stop("Fig 3G baseline does not reproduce")
f3["finite_sens"] = np.isfinite(f3.coloc_best_susie_pp4_sens.to_numpy(float))
f3["state_sens"] = f3_state(f3.coloc_best_susie_pp4_sens.to_numpy(float), deg3)
F3 = {}
for ph, st, fin in [("base", "state", np.ones(len(f3), bool)), ("sens", "state_sens", f3.finite_sens.to_numpy())]:
    sub = f3[fin]
    F3[ph] = {"plotted": int(len(sub)),
              "by_state": sub[st].value_counts().to_dict(),
              "by_state_rule_a_tested": sub[sub[f"state_a_{ph}"] != "untestable"][st].value_counts().to_dict(),
              "by_state_rule_a_untestable": sub[sub[f"state_a_{ph}"] == "untestable"][st].value_counts().to_dict()}
    log(f"Fig 3G {ph}", json.dumps(F3[ph]))
F3["dropped_no_finite_pp4"] = f3[~f3.finite_sens].state.value_counts().to_dict()
RES["fig3g"] = F3
f3.to_csv(f"{OUT}/fig3g_states.tsv", sep="\t", index=False)

# ------------------------------------------------------------------ 3,918 / 3,335 / 3,319
cat = pd.read_parquet(CAT, columns=["ensembl_id", "symbol", "r8_genetics_call_state", "r8_genetics_testability_state",
                                    "gen_susie_state_t12", "primary_evidence_class", "class_display_r9"])
cat["g"] = cat.ensembl_id.astype(str).str.split(".").str[0]
if cat.g.duplicated().any():
    stop("duplicate catalog Ensembl ids")
callmap = dict(zip(cat.g, cat.r8_genetics_call_state))


def pp3_counts(mm, known_mask=None, assume_partial=None):
    t = mm[mm.tier.isin([1, 2])].copy()
    t["s4"] = pd.to_numeric(t["PP.H4.susie"], errors="coerce")
    t["s3"] = pd.to_numeric(t["PP.H3.susie"], errors="coerce")
    fin = t[t.s4.notna()]
    genes = set(fin.ensembl)
    best = fin.groupby("ensembl").s4.max()
    q = fin.s3 > 0.8
    if assume_partial is not None:
        part = fin.index.isin(list(partial_i))
        q = (q & ~part) if assume_partial == "no" else (q | part)
    hi3 = set(fin[q].ensembl)
    sel = {g for g in hi3 if best[g] <= 0.5}
    return genes, sel


g0, s0 = pp3_counts(m)
ind0 = {g for g in s0 if callmap.get(g.split(".")[0]) == "indeterminate"}
log(f"PP.H3 baseline: genes with finite tier-1/2 SuSiE {len(g0)}; PP.H3 > 0.8 and best PP.H4 <= 0.5 {len(s0)}; "
    f"indeterminate {len(ind0)}")
if (len(g0), len(s0), len(ind0)) != (3918, 3335, 3319):
    stop("PP.H3 baseline does not reproduce 3,918 / 3,335 / 3,319")
P = {"baseline": dict(finite=len(g0), pp3=len(s0), indeterminate=len(ind0),
                      not_indeterminate=collections.Counter(callmap.get(g.split('.')[0], 'not_in_catalog')
                                                            for g in s0 - ind0))}
gstate = {r: dict(zip(G.ensembl, G[f"state_{r}_sens"])) for r in "ab"}
gstate0 = {r: dict(zip(G.ensembl, G[f"state_{r}_base"])) for r in "ab"}
def call_after(g, st):
    """Catalog genetics call implied by a tier-1/2 SuSiE state (r9d: positive -> supported,
    tested_le_0.5 -> indeterminate, untestable -> untestable); genes absent from r9d stay absent."""
    c = callmap.get(g.split(".")[0])
    if c is None:
        return "not_in_catalog"
    return {"positive": "supported", "tested_le_0.5": "indeterminate", "untestable": "untestable"}.get(st, c)


chk = collections.Counter(call_after(g, gstate0["a"].get(g)) for g in s0)
log(f"PP.H3 baseline calls implied by rule (a) states: {dict(chk)}")
if chk["indeterminate"] != len(ind0):
    stop("implied catalog call does not reproduce 3,319")
P["baseline_rule_b"] = dict(collections.Counter(call_after(g, gstate0["b"].get(g)) for g in s0))
for variant, ap_ in [("carried_PP.H3", None), ("lower_partial_not_counted", "no"), ("upper_partial_counted", "yes")]:
    g1, s1_ = pp3_counts(ms, assume_partial=ap_)
    out = dict(finite=len(g1), pp3=len(s1_))
    for r in "ab":
        cc = collections.Counter(call_after(g, gstate[r].get(g)) for g in s1_)
        out[f"calls_rule_{r}"] = dict(cc)
        out[f"indeterminate_rule_{r}"] = cc["indeterminate"]
    P[variant] = out
    log(f"PP.H3 sensitivity ({variant})", json.dumps(out))
RES["pp3"] = P

# ------------------------------------------------------------------ Catalog r9d moves
C = G.assign(g=G.ensembl.str.split(".").str[0]).merge(cat, on="g", how="left")
MV = {}
for r in "ab":
    tested = ["positive", "tested_le_0.5"]
    mv = C[C[f"state_{r}_base"].isin(tested) & (C[f"state_{r}_sens"] == "untestable")]
    MV[r] = {"n_genes": int(len(mv)),
             "not_in_catalog": int(mv.ensembl_id.isna().sum()),
             "by_gen_susie_state_t12": mv.gen_susie_state_t12.fillna("NA").value_counts().to_dict(),
             "by_r8_genetics_call_state": mv.r8_genetics_call_state.fillna("NA").value_counts().to_dict(),
             "by_primary_evidence_class": mv.primary_evidence_class.fillna("NA").value_counts().to_dict(),
             "by_class_display_r9": mv.class_display_r9.fillna("NA").value_counts().to_dict()}
    mv[["ensembl", "symbol", f"state_{r}_base", f"state_{r}_sens", "coloc_best_susie_pp4_base",
        "coloc_best_susie_pp4_sens", "gen_susie_state_t12", "r8_genetics_call_state",
        "primary_evidence_class"]].to_csv(f"{OUT}/catalog_moves_rule_{r}.tsv", sep="\t", index=False)
    log(f"Catalog r9d moves rule ({r})", json.dumps(MV[r]))
# r9d rule check: every genetics-untestable gene is primary_evidence_class 'untested'
RES["catalog_r9d_untestable_to_class"] = cat[cat.r8_genetics_call_state == "untestable"] \
    .primary_evidence_class.value_counts().to_dict()
RES["catalog_moves"] = MV

with open(f"{OUT}/results.json", "w") as fh:
    json.dump(RES, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
with open(f"{OUT}/pip_freeze.txt", "w") as fh:
    fh.write(subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True).stdout)
log("DONE")
