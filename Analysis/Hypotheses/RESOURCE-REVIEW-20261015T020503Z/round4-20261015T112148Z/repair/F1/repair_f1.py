"""F1 repair: re-derive the disputed Figure 1B values from the release files.

Inputs (adopted releases only): COLOC r2 gene-level tier-1/2 table, COLOC r2 all-GWAS rows,
GWAS tier config, adopted 774-unit refit DEG table, evidence-classes r8 table (14,720 rows).
Implementations here are written independently of the H57 script (MH/RBG by formula,
exact convolution, own power solver, own greedy matcher).
Usage: python repair_f1.py <run_dir> <fig_dir>
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from scipy.stats import hypergeom, nchypergeom_fisher

SEED = 20260923
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
C = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/gene_level_coloc_tier12.csv"
A = ROOT + "/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv"
TIER = ROOT + "/GWAS/finemapping/config/gwas_trait_tier.tsv"
B = (ROOT + "/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/"
     "primary_source_controls/deg_results.csv")
EC = (ROOT + "/RNA-seq/results/manuscript_release/candidates/evidence-classes-corrected-r8-20261009T142027Z/"
      "evidence_class_table.tsv")
OLD_TABLE = ROOT + "/figures/candidates/review-r3-fixes-20261015T112148Z/F1/fig1b_U_eval_gene_table.tsv"
RUN, FIG = sys.argv[1], sys.argv[2]
R = {}


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- universes
coloc = pd.read_csv(C)
coloc["g"] = coloc.ensembl.astype(str).str.split(".").str[0]
R["coloc_rows"] = len(coloc)
R["coloc_dup_ids"] = int(coloc.g.duplicated().sum())
coloc = coloc.sort_values(["coloc_best_susie_pp4", "coloc_best_abf_pp4"], ascending=False).drop_duplicates("g")
bulk = pd.read_csv(B)
bulk["g"] = bulk.gene.astype(str).str.split(".").str[0]
R["bulk_rows"] = len(bulk)
R["bulk_dup_ids"] = int(bulk.g.duplicated().sum())
keep = ["g", "coloc_best_susie_pp4", "coloc_best_abf_pp4", "susie_state_t12"]
full = bulk.drop_duplicates("g").merge(coloc[keep], on="g", how="inner")
full["susie"] = pd.to_numeric(full.coloc_best_susie_pp4, errors="coerce")
full["abf"] = pd.to_numeric(full.coloc_best_abf_pp4, errors="coerce")
full["gen"] = (full.susie > 0.5).fillna(False).astype(bool)
assert int(full.gen.sum()) == int((full.susie >= 0.5).sum())
full["deg"] = ((full.padj < 0.05) & (full.logFC.abs() > 0.5)).astype(bool)
full = full.sort_values("g").reset_index(drop=True)
R["U_full_state_by_deg"] = pd.crosstab(full.susie_state_t12.fillna("NA"), full.deg).to_dict()
R["gen_equals_state_positive"] = bool((full.gen == (full.susie_state_t12 == "positive")).all())
ut = full[full.susie_state_t12 == "untestable"]
R["untestable_n"] = len(ut)
R["untestable_susie_range"] = [float(ut.susie.min()), float(ut.susie.max())]
ueval = full[full.susie_state_t12.isin(["positive", "tested_le_0.5"])].reset_index(drop=True)
ueval7 = full[full.susie.notna()].reset_index(drop=True)


def sizes(df, gen="gen"):
    g = df[gen].to_numpy(); d = df.deg.to_numpy()
    return len(df), int(g.sum()), int(d.sum()), int((g & d).sum())


for nm, df in [("U_full", full), ("U_eval", ueval), ("U_eval7", ueval7)]:
    N, G, D, O = sizes(df)
    t = [[O, G - O], [D - O, N - G - D + O]]
    orr, p = stats.fisher_exact(t)
    R[nm] = dict(N=N, G=G, D=D, O=O, E=G * D / N, fisher_OR=float(orr), fisher_p=float(p),
                 neither=N - G - D + O, nonoverlap_obs=(G - O) / G, nonoverlap_exp=1 - D / N,
                 marker_gen=(G * D / N) / G, marker_deg=(G * D / N) / D)
    log(nm, R[nm])


# ---------------------------------------------------------------- stratified test (own code)
def bins_of(df, col, k):
    # df sorted by Ensembl ID, so rank(method="first") breaks AveExpr ties by Ensembl ID
    return pd.qcut(df[col].rank(method="first"), k, labels=False).to_numpy()


def margins(b, gen, deg):
    ks = np.unique(b)
    n = np.array([(b == k).sum() for k in ks])
    G = np.array([gen[b == k].sum() for k in ks])
    D = np.array([deg[b == k].sum() for k in ks])
    O = np.array([(gen & deg)[b == k].sum() for k in ks])
    return n, G, D, O


def mh_rbg(n, G, D, O):
    a = O.astype(float); bb = (G - O).astype(float); c = (D - O).astype(float); d = (n - G - D + O).astype(float)
    Rk = a * d / n; Sk = bb * c / n; Pk = (a + d) / n; Qk = (bb + c) / n
    sR, sS = Rk.sum(), Sk.sum()
    orr = sR / sS
    v = (Pk * Rk).sum() / (2 * sR ** 2) + (Pk * Sk + Qk * Rk).sum() / (2 * sR * sS) + (Qk * Sk).sum() / (2 * sS ** 2)
    se = np.sqrt(v)
    return float(orr), float(np.exp(np.log(orr) - 1.959963984540054 * se)), float(np.exp(np.log(orr) + 1.959963984540054 * se))


def null_pmf(n, G, D):
    pmf = np.array([1.0])
    for nk, gk, dk in zip(n, G, D):
        x = np.arange(0, min(gk, dk) + 1)
        pmf = np.convolve(pmf, hypergeom.pmf(x, nk, dk, gk))
    return pmf


def strat(df, gen_col, col="AveExpr", k=10):
    b = bins_of(df, col, k)
    gen = df[gen_col].to_numpy(); deg = df.deg.to_numpy()
    n, G, D, O = margins(b, gen, deg)
    E = float((G * D / n).sum()); Ot = int(O.sum())
    pmf = null_pmf(n, G, D)
    x = np.arange(len(pmf))
    p = float(pmf[np.abs(x - E) >= abs(Ot - E) - 1e-9].sum())
    orr, lo, hi = mh_rbg(n, G, D, O)
    return b, dict(k=k, O=Ot, E_strat=E, O_over_E=Ot / E, MH_OR=orr, MH_lo=lo, MH_hi=hi, p_exact=min(1.0, p),
                   null_mean=float((x * pmf).sum()), marker_gen=E / G.sum(), marker_deg=E / D.sum())


bins_e, R["T1_U_eval_10bins"] = strat(ueval, "gen")
log("T1", R["T1_U_eval_10bins"])
for k in (5, 20):
    R[f"U_eval_{k}bins"] = strat(ueval, "gen", k=k)[1]
R["U_eval7_10bins"] = strat(ueval7, "gen")[1]
bins_f, R["U_full_10bins"] = strat(full, "gen")
log("U_full 10 bins", R["U_full_10bins"])

# ---------------------------------------------------------------- D8: ABF-only map in U_full
full["gen_abf"] = (full.abf >= 0.5).fillna(False).astype(bool)
N, G, D, O = sizes(full, "gen_abf")
R["D8_abf_U_full_unstrat"] = dict(N=N, G=G, D=D, O=O, E=G * D / N)
R["D8_abf_U_full_10bins"] = strat(full, "gen_abf")[1]
log("D8", R["D8_abf_U_full_unstrat"], R["D8_abf_U_full_10bins"])

# ---------------------------------------------------------------- power (unstratified exact, equal tails)
N, G, D, O = sizes(ueval)
X = hypergeom(N, D, G)
xs = np.arange(0, min(G, D) + 1)
c_lo = int(xs[X.cdf(xs) <= 0.025].max())
c_hi = int(xs[X.sf(xs - 1) <= 0.025].min())
size = float(X.cdf(c_lo) + X.sf(c_hi - 1))


def power(psi):
    Y = nchypergeom_fisher(N, D, G, psi)
    return float(Y.cdf(c_lo) + Y.sf(c_hi - 1))


or_up = brentq(lambda s: power(s) - 0.8, 1.0001, 10)
or_dn = brentq(lambda s: power(s) - 0.8, 0.05, 0.9999)
R["power"] = dict(c_lo=c_lo, c_hi=c_hi, size=size, OR_80_up=or_up, OR_80_down=or_dn,
                  power_at_MH_lo=power(R["T1_U_eval_10bins"]["MH_lo"]))
log("power", R["power"])

# ---------------------------------------------------------------- NC-1 style matched list (own greedy matcher)
ae = ueval.AveExpr.to_numpy(); gen = ueval.gen.to_numpy(); deg = ueval.deg.to_numpy()
gi = np.where(gen)[0]
pool = np.where(~gen)[0]
pool = pool[np.argsort(ae[pool], kind="stable")]
pv = ae[pool]
n_b, _, D_b, _ = margins(bins_e, gen, deg)
rng = np.random.default_rng(SEED)
nc = []
orders = [("ascending_AveExpr", gi[np.argsort(ae[gi], kind="stable")])] + \
         [(f"random_{i}", rng.permutation(gi)) for i in range(20)]
for name, order in orders:
    used = np.zeros(len(pool), bool); chosen = []
    for j in order:
        dist = np.abs(pv - ae[j]); dist[used] = np.inf
        m = dist.min()
        cand = np.where(dist == m)[0]
        c = cand[np.argmin(pv[cand])]  # tie: lower AveExpr
        used[c] = True; chosen.append(pool[c])
    lst = np.zeros(len(ueval), bool); lst[np.array(chosen)] = True
    _, Gk, _, Ok = margins(bins_e, lst, deg)
    E = float((Gk * D_b / n_b).sum())
    nc.append(dict(order=name, O=int(Ok.sum()), E_strat=E, strat_OE=int(Ok.sum()) / E))
nc = pd.DataFrame(nc)
nc.to_csv(os.path.join(RUN, "nc_matched_lists.tsv"), sep="\t", index=False)
R["NC_matched"] = dict(O_values=sorted(nc.O.unique().tolist()), E_strat_range=[float(nc.E_strat.min()), float(nc.E_strat.max())],
                       strat_OE_mean=float(nc.strat_OE.mean()), n_orders=len(nc))
log("NC", R["NC_matched"])

# ---------------------------------------------------------------- D7: what the 935 DEGs without a multi-signal result had
tier = pd.read_csv(TIER, sep="\t")
t12 = set(tier.study_name[tier.tier.astype(str).isin(["1", "2"])])
nores = full[(full.susie_state_t12 == "no_susie_result") & full.deg]
allg = pd.read_csv(A, usecols=["gwas_name", "ensembl", "method", "PP.H4.abf"])
allg["g"] = allg.ensembl.astype(str).str.split(".").str[0]
a12 = allg[allg.gwas_name.isin(t12) & allg.g.isin(set(nores.g))]
meth = a12.groupby("g").method.apply(lambda s: "+".join(sorted(set(s.astype(str)))))
R["no_result_DEGs"] = dict(n=len(nores), n_with_t12_rows=int(meth.shape[0]), method_combos=meth.value_counts().to_dict(),
                           n_with_finite_abf=int(a12.groupby("g")["PP.H4.abf"].apply(lambda s: s.notna().any()).sum()),
                           mean_AveExpr=float(nores.AveExpr.mean()),
                           mean_AveExpr_tested_nongenetic_DEG=float(ueval[ueval.deg & ~ueval.gen].AveExpr.mean()))
log("no-result DEGs", R["no_result_DEGs"])
del allg, a12

# ---------------------------------------------------------------- 14,720 symbol-grouped universe
ec = pd.read_csv(EC, sep="\t", low_memory=False)
jt = ec[ec.joint_testable.astype(str).str.lower() == "true"].copy()
jt["s"] = pd.to_numeric(jt.max_susie_pp4_all, errors="coerce")
jt["deg"] = (pd.to_numeric(jt.bulk_padj, errors="coerce") < 0.05) & (pd.to_numeric(jt.bulk_logFC, errors="coerce").abs() > 0.5)
R["U_14720"] = dict(N=len(jt), G=int((jt.s > 0.5).sum()), D=int(jt.deg.sum()), O=int(((jt.s > 0.5) & jt.deg).sum()),
                    no_susie=int(jt.s.isna().sum()), no_susie_deg=int((jt.s.isna() & jt.deg).sum()))
R["U_14720"]["E"] = R["U_14720"]["G"] * R["U_14720"]["D"] / R["U_14720"]["N"]
log("14,720", R["U_14720"])

# ---------------------------------------------------------------- panel source table
out = ueval[["g", "AveExpr", "logFC", "padj", "susie", "susie_state_t12", "gen", "deg"]].copy()
out["aveexpr_bin10"] = bins_e
out.to_csv(os.path.join(FIG, "fig1b_U_eval_gene_table.tsv"), sep="\t", index=False)
old = pd.read_csv(OLD_TABLE, sep="\t")
new = pd.read_csv(os.path.join(FIG, "fig1b_U_eval_gene_table.tsv"), sep="\t")
R["gene_table_equals_round1"] = bool(old.equals(new))
log("gene table equals round 1:", R["gene_table_equals_round1"])

with open(os.path.join(RUN, "rederive.json"), "w") as fh:
    json.dump(R, fh, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
log("done")
