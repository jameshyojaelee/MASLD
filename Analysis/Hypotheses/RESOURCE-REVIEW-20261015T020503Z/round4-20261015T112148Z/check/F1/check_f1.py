"""Blind re-derivation of F1 (Fig 1B) numbers from the release files.

Written independently of figures/candidates/review-r3-fixes-20261015T112148Z/F1/reproduce_f1.py and of
the H57 scripts: MH OR, RBG CI, exact convolution test, power and the planted-signal and
matched-list controls are coded here from their definitions. Inputs:
  COLOC r2 gene-level tier-1/2 table and all-rows table,
  adopted 774-unit refit DEG table (and the reproduction copy, for agreement only),
  Catalog r8 reviewed gene index (for the 1,129 split).
Unit: gene as a member of two lists. Seed 20260923.
"""
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = sys.argv[1]
C = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
CALL = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
B = ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
BR = ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z-reproduction/primary_source_controls/deg_results.csv"
CAT = ROOT + "/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet"
EXPECT = {C: "6360984750864e7ab56593bce90e8ce665d781ec72a1ba52790d6ba53e83b8d6",
          B: "0491ed10024388c3497dccb9f99bf669683edd5903c9d21c5308850d7c10bde6",
          BR: "e30509ac9c46449cfbc916191332890af7fef06e3208b34a91d7cdf90f9658e4",
          CAT: "ffd3c4650c77a899e9ebed7327697398c4b6462cfffcd7f0e635cd0d023db06b"}
res = {}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


for p, e in EXPECT.items():
    s = sha(p)
    print("sha256", p.split("/")[-3:], s == e)
    res["sha_ok_" + p.split("/")[-1]] = bool(s == e)

# ---- inputs
coloc = pd.read_csv(C)
coloc["g"] = coloc.ensembl.astype(str).str.replace(r"\.\d+$", "", regex=True)
print("coloc rows", len(coloc), "unique g", coloc.g.nunique())
bulk = pd.read_csv(B)
bulk["g"] = bulk.gene.astype(str).str.replace(r"\.\d+$", "", regex=True)
print("bulk rows", len(bulk), "unique g", bulk.g.nunique())
assert coloc.g.is_unique and bulk.g.is_unique

# BH across the full table
bh = stats.false_discovery_control(bulk["P.Value"].to_numpy(), method="bh")
print("padj equals BH(P.Value) over", len(bulk), "genes: max abs diff", float(np.max(np.abs(bh - bulk.padj))))
bulk["deg"] = (bulk.padj < 0.05) & (bulk.logFC.abs() > 0.5)
print("canonical DEG calls in full table", int(bulk.deg.sum()))
res["n_bulk"] = len(bulk); res["n_deg_all"] = int(bulk.deg.sum())

# agreement with reproduction copy
br = pd.read_csv(BR)
br["g"] = br.gene.astype(str).str.replace(r"\.\d+$", "", regex=True)
mm = bulk[["g", "deg", "AveExpr"]].merge(br.assign(deg2=(br.padj < 0.05) & (br.logFC.abs() > 0.5))[["g", "deg2", "AveExpr"]], on="g")
print("repro copy: genes", len(mm), "DEG disagreements", int((mm.deg != mm.deg2).sum()),
      "max |dAveExpr|", float((mm.AveExpr_x - mm.AveExpr_y).abs().max()))

# ---- universes
full = bulk.merge(coloc[["g", "coloc_best_susie_pp4", "susie_state_t12"]], on="g", how="inner")
full["pp4"] = pd.to_numeric(full.coloc_best_susie_pp4, errors="coerce")
full["gen"] = full.pp4.fillna(-1) > 0.5
assert ((full.pp4.fillna(-1) >= 0.5) == full.gen).all()
assert (full.gen == (full.susie_state_t12 == "positive")).all()
print("U_full", len(full))
ct = pd.crosstab(full.susie_state_t12, full.deg, margins=True)
print(ct)
print("finite PP.H4 by state:\n", full.groupby("susie_state_t12").pp4.agg(["count", "min", "max"]))
res["U_full"] = dict(N=len(full), G=int(full.gen.sum()), D=int(full.deg.sum()), O=int((full.gen & full.deg).sum()))
res["state_by_deg"] = {s: {"all": int((full.susie_state_t12 == s).sum()),
                           "deg": int(((full.susie_state_t12 == s) & full.deg).sum())}
                       for s in full.susie_state_t12.unique()}

ue = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].copy()
ue7 = full[full.pp4.notna()].copy()


def unstrat(df, name):
    N = len(df); G = int(df.gen.sum()); D = int(df.deg.sum()); O = int((df.gen & df.deg).sum())
    a, b, c, d = O, G - O, D - O, N - G - D + O
    E = G * D / N
    orr = a * d / (b * c)
    p = stats.fisher_exact([[a, b], [c, d]]).pvalue
    ci = stats.contingency.odds_ratio([[a, b], [c, d]]).confidence_interval()
    r = dict(N=N, G=G, D=D, O=O, E=E, OR=orr, p=p, ci_low=ci.low, ci_high=ci.high,
             mark_gen=E / G, mark_deg=E / D, obs_gen=O / G, obs_deg=O / D,
             nonoverlap_obs=(G - O) / G, nonoverlap_exp=1 - D / N,
             gen_only=b, deg_only=c, neither=d)
    print(name, json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()}))
    res[name] = r
    return r


unstrat(full, "U_full_unstrat")
unstrat(ue, "U_eval_unstrat")
unstrat(ue7, "U_eval7_unstrat")


def make_bins(df, k, how):
    order = np.lexsort((df.g.to_numpy(), df.AveExpr.to_numpy()))  # AveExpr, ties by Ensembl ID
    rank = np.empty(len(df), int); rank[order] = np.arange(1, len(df) + 1)
    if how == "qcut":
        return pd.qcut(rank, k, labels=False)
    b = np.empty(len(df), int)  # equal split of the sorted order
    for i, idx in enumerate(np.array_split(order, k)):
        b[idx] = i
    return b


def strat(gen, deg, bins):
    t = []
    for k in np.unique(bins):
        m = bins == k
        n = int(m.sum()); G = int(gen[m].sum()); D = int(deg[m].sum()); a = int((gen[m] & deg[m]).sum())
        t.append((n, G, D, a))
    return np.array(t)


def mh_rbg(t):
    n, G, D, a = t.T.astype(float)
    b = G - a; c = D - a; d = n - G - D + a
    R = a * d / n; S = b * c / n; P = (a + d) / n; Q = (b + c) / n
    orr = R.sum() / S.sum()
    v = ((P * R).sum() / (2 * R.sum() ** 2) + (P * S + Q * R).sum() / (2 * R.sum() * S.sum())
         + (Q * S).sum() / (2 * S.sum() ** 2))
    se = np.sqrt(v)
    return orr, np.exp(np.log(orr) - 1.959963985 * se), np.exp(np.log(orr) + 1.959963985 * se)


def null_pmf(t):
    pmf = np.array([1.0])
    for n, G, D, _ in t:
        x = np.arange(0, min(G, D) + 1)
        pk = stats.hypergeom.pmf(x, n, D, G)
        pmf = np.convolve(pmf, pk)
    return pmf


def exact_p(t, O):
    pmf = null_pmf(t)
    E = (t[:, 1] * t[:, 2] / t[:, 0]).sum()
    x = np.arange(len(pmf))
    return float(pmf[np.abs(x - E) >= abs(O - E) - 1e-9].sum()), E, pmf


def full_strat(df, k, how, name):
    bins = make_bins(df, k, how)
    t = strat(df.gen.to_numpy(), df.deg.to_numpy(), bins)
    O = int(t[:, 3].sum())
    p, E, pmf = exact_p(t, O)
    orr, lo, hi = mh_rbg(t)
    r = dict(O=O, E=E, OE=O / E, MH_OR=orr, lo=lo, hi=hi, p=p, mark_gen=E / t[:, 1].sum(),
             mark_deg=E / t[:, 2].sum(), bin_sizes=t[:, 0].tolist(), G_k=t[:, 1].tolist(), D_k=t[:, 2].tolist(),
             a_k=t[:, 3].tolist(), null_mean=float((np.arange(len(pmf)) * pmf).sum()))
    print(name, json.dumps({kk: (round(v, 6) if isinstance(v, float) else v) for kk, v in r.items()}))
    res[name] = r
    return r, bins, t


r10, bins10, t10 = full_strat(ue, 10, "qcut", "U_eval_10bin_qcut")
full_strat(ue, 10, "split", "U_eval_10bin_arraysplit")
full_strat(ue, 5, "qcut", "U_eval_5bin")
full_strat(ue, 20, "qcut", "U_eval_20bin")
full_strat(ue7, 10, "qcut", "U_eval7_10bin")
full_strat(full, 10, "qcut", "U_full_10bin")

# statsmodels cross-check of MH OR and RBG CI
try:
    from statsmodels.stats.contingency_tables import StratifiedTable
    tabs = [np.array([[a, G - a], [D - a, n - G - D + a]]) for n, G, D, a in t10]
    st = StratifiedTable(tabs)
    print("statsmodels MH OR", st.oddsratio_pooled, "CI", st.oddsratio_pooled_confint())
    res["statsmodels_MH"] = [float(st.oddsratio_pooled)] + [float(x) for x in st.oddsratio_pooled_confint()]
except Exception as ex:  # noqa
    print("statsmodels not available:", ex)

# compare against the candidate gene table bins
cand = pd.read_csv(ROOT + "/figures/candidates/review-r3-fixes-20261015T112148Z/F1/fig1b_U_eval_gene_table.tsv", sep="\t")
mine = ue.assign(mybin=bins10)[["g", "mybin", "gen", "deg", "AveExpr"]].merge(cand, on="g", suffixes=("", "_c"))
print("candidate gene table rows", len(cand), "matched", len(mine),
      "bin disagreements", int((mine.mybin != mine.aveexpr_bin10).sum()),
      "gen disagreements", int((mine.gen != mine.gen_c).sum()), "deg disagreements", int((mine.deg != mine.deg_c).sum()))

# H57 bins file
h57b = pd.read_csv(ROOT + "/Analysis/Hypotheses/H57-20261015T112421Z/run-20261015T113107Z/bins_aveexpr10.tsv", sep="\t")
print("H57 bins_aveexpr10.tsv:\n", h57b.to_string())

# ---- power (unstratified, Fisher noncentral hypergeometric), equal-tail two-sided 0.05
N, G, D = res["U_eval_unstrat"]["N"], res["U_eval_unstrat"]["G"], res["U_eval_unstrat"]["D"]
x = np.arange(0, min(G, D) + 1)
p0 = stats.hypergeom.pmf(x, N, D, G)
cdf = np.cumsum(p0); sf = np.cumsum(p0[::-1])[::-1]
cL = int(x[cdf <= 0.025].max()); cU = int(x[sf <= 0.025].min())
from scipy.optimize import brentq
pw_up = lambda o: stats.nchypergeom_fisher.sf(cU - 1, N, D, G, o) - 0.8
pw_dn = lambda o: stats.nchypergeom_fisher.cdf(cL, N, D, G, o) - 0.8
or_up = brentq(pw_up, 1.0, 5.0); or_dn = brentq(pw_dn, 0.05, 1.0)
print(f"power: cL {cL} cU {cU} size {cdf[cL] + sf[cU]:.4f} OR_up80 {or_up:.3f} OR_dn80 {or_dn:.3f}")
res["power"] = dict(cL=cL, cU=cU, OR_up=or_up, OR_dn=or_dn)
# power of the stratified exact test at OR 0.46 and 1.77 is checked through PC-B below

# ---- PC-B planted signal in U_eval, 10 bins (my implementation)
rng = np.random.default_rng(20260923)
p_null, E10, pmf10 = exact_p(t10, 26)
xs = np.arange(len(pmf10))
rej = lambda O: pmf10[np.abs(xs - E10) >= abs(O - E10) - 1e-9].sum() < 0.05
genv = ue.gen.to_numpy(); degv = ue.deg.to_numpy()


def pcb(s, reps=1000):
    det = 0
    for _ in range(reps):
        g = np.zeros(len(genv), bool)
        for k in np.unique(bins10):  # null labelling within bins
            idx = np.flatnonzero(bins10 == k)
            g[rng.choice(idx, size=int(genv[idx].sum()), replace=False)] = True
        ok = True
        for _ in range(abs(s)):
            if s > 0:
                A = g & ~degv; Bm = ~g & degv
            else:
                A = g & degv; Bm = ~g & ~degv
            w = np.array([A[bins10 == k].sum() * Bm[bins10 == k].sum() for k in range(10)], float)
            if w.sum() == 0:
                ok = False; break
            k = rng.choice(10, p=w / w.sum())
            i = rng.choice(np.flatnonzero(A & (bins10 == k))); j = rng.choice(np.flatnonzero(Bm & (bins10 == k)))
            g[i] = False; g[j] = True
        if ok and rej(int((g & degv).sum())):
            det += 1
    return det / reps


res["PCB_enrich_s14"] = pcb(14); res["PCB_deplete_s12"] = pcb(-12)
print("PC-B detection: enrichment s=+14", res["PCB_enrich_s14"], "depletion s=-12", res["PCB_deplete_s12"])

# ---- NC-1 style matched list: nearest unused non-genetic gene by AveExpr (ties to lower AveExpr)
ae = ue.AveExpr.to_numpy()
pool0 = np.flatnonzero(~genv)
cases = np.flatnonzero(genv)
ratios, Os = [], []
for rep in range(50):
    used = np.zeros(len(ae), bool)
    lst = []
    for i in rng.permutation(cases):
        cand_ = pool0[~used[pool0]]
        dd = np.abs(ae[cand_] - ae[i])
        m = dd.min()
        tie = cand_[dd == m]
        j = tie[np.argmin(ae[tie])]
        used[j] = True; lst.append(j)
    L = np.zeros(len(ae), bool); L[lst] = True
    tl = strat(L, degv, bins10)
    E_l = (tl[:, 1] * tl[:, 2] / tl[:, 0]).sum(); O_l = int((L & degv).sum())
    ratios.append(O_l / E_l); Os.append(O_l)
print("NC-1 style lists (50 random orders): O values", sorted(set(Os)), "mean strat O/E", np.mean(ratios),
      "E_strat range", None)
res["NC1_mean_OE"] = float(np.mean(ratios)); res["NC1_O"] = sorted(set(Os))
pl = exact_p(strat(L, degv, bins10), Os[-1])[0]
print("NC-1 last list exact p", pl)

# ---- 1,129 split by state (14,704) and Catalog r8 class check
only = full[full.deg & ~full.gen]
print("U_full DEG-only", len(only), only.susie_state_t12.value_counts().to_dict())
res["U_full_deg_only_by_state"] = only.susie_state_t12.value_counts().to_dict()
try:
    cat = pd.read_parquet(CAT)
    print("catalog columns", list(cat.columns)[:40])
    idcol = [c for c in cat.columns if "ensembl" in c.lower() or c in ("gene_id", "g")]
    print("catalog id columns", idcol)
    cc = "primary_evidence_class"
    cat["g"] = cat[idcol[0]].astype(str).str.replace(r"\.\d+$", "", regex=True)
    x_ = cat[[cc, "g"]].merge(full[["g", "susie_state_t12"]], on="g", how="left")
    tab = pd.crosstab(x_[cc], x_.susie_state_t12.fillna("not_in_U_full"))
    print(tab.loc[[i for i in tab.index if i in ("established_state_associated", "concordant", "genetically_anchored")]])
except Exception as ex:  # noqa
    print("catalog check failed:", repr(ex))

# ---- why the 935 DEGs have no SuSiE result (method values over tier-1/2 rows)
nev = set(full.loc[(full.susie_state_t12 == "no_susie_result") & full.deg, "g"])
tiers = coloc[["g", "driving_tier"]]
allr = pd.read_csv(CALL, usecols=["gwas_name", "ensembl", "method", "PP.H4.susie"])
allr["g"] = allr.ensembl.astype(str).str.replace(r"\.\d+$", "", regex=True)
tierf = ROOT + "/GWAS/finemapping/config/gwas_trait_tier.tsv"
tt = pd.read_csv(tierf, sep="\t")
print("tier table columns", list(tt.columns))
tcol = "tier"
ncol = "study_name"
t12 = set(tt.loc[tt[tcol].astype(str).isin(["1", "2"]), ncol])
sub = allr[allr.g.isin(nev) & allr.gwas_name.isin(t12)]
combo = sub.groupby("g").method.apply(lambda s: "+".join(sorted(set(s)))).value_counts()
print("935 never-evaluated DEGs: genes with tier-1/2 rows", sub.g.nunique(), "method combos", combo.to_dict())
res["nev_method_combos"] = combo.to_dict()
res["nev_any_finite_susie"] = int(sub["PP.H4.susie"].notna().sum())

json.dump(res, open(OUT + "/check_f1_results.json", "w"), indent=1, default=float)
print("done")
