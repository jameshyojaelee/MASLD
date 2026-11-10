#!/usr/bin/env python3
"""Rerun the adopted genetics-by-RNA overlap test (3,589-gene AN1 rule-(a) universe)
with the 90-gene eQTL parity refit applied. Test functions copied verbatim from
REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/scripts/an1_gene_level.py (from 01_h57_overlap.py).
Baseline must reproduce the adopted values before the refit is applied."""
import json, os, sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from statsmodels.stats.contingency_tables import StratifiedTable

TOL = 1e-9
AN1 = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN1/gene_level_20261021T165252Z/h57_universe_rule_a_sens.tsv"
RUN = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z"
OUT = sys.argv[1]


def make_bins(df, col, k):
    return pd.qcut(df[col].rank(method="first"), k, labels=False).to_numpy()


def bin_margins(bins, gen, deg):
    ks = np.unique(bins)
    n = np.array([(bins == k).sum() for k in ks])
    Gk = np.array([gen[bins == k].sum() for k in ks])
    Dk = np.array([deg[bins == k].sum() for k in ks])
    Ok = np.array([(gen & deg)[bins == k].sum() for k in ks])
    return n, Gk, Dk, Ok


def exact_null(n, Gk, Dk):
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
    return dict(MH_OR=float(st.oddsratio_pooled), MH_lo=float(lo), MH_hi=float(hi))


def run(ue):
    ue = ue.sort_values("g").reset_index(drop=True)
    gen, deg = ue.gen.to_numpy(bool), ue.deg.to_numpy(bool)
    b10 = make_bins(ue, "AveExpr", 10)
    _, Gk, Dk, _ = bin_margins(b10, gen, deg)
    K = 10 if int(((Gk >= 5) & (Dk >= 5)).sum()) >= 7 else 5
    bins = b10 if K == 10 else make_bins(ue, "AveExpr", K)
    n, Gk, Dk, Ok = bin_margins(bins, gen, deg)
    E = float((Gk * Dk / n).sum()); O = int(Ok.sum())
    pmf = exact_null(n, Gk, Dk)
    cdf = np.cumsum(pmf)
    r = dict(N=len(ue), G=int(gen.sum()), D=int(deg.sum()), O=O, E=E, K=K, p_exact=exact_p(pmf, E, O),
             null_q025=int(np.searchsorted(cdf, 0.025)), null_q975=int(np.searchsorted(cdf, 0.975)))
    r.update(mh(n, Gk, Dk, Ok))
    N = len(ue); G_ = int(gen.sum()); D = int(deg.sum())
    orr, fp = stats.fisher_exact([[O, G_ - O], [D - O, N - G_ - D + O]])
    r.update(unstrat_E=G_ * D / N, fisher_OR=float(orr), fisher_p=float(fp),
             neither=int((~gen & ~deg).sum()), rna_only=int((~gen & deg).sum()), gen_only=int((gen & ~deg).sum()))
    if "logFC" in ue:
        rho, prho = stats.spearmanr(ue.susie, ue.logFC.abs())
        r.update(rho=float(rho), rho_p=float(prho), gen_up=int((gen & (ue.logFC > 0)).sum()),
                 gen_down=int((gen & (ue.logFC < 0)).sum()))
    return r


BULK = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv"
bulk = pd.read_csv(BULK); bulk["g"] = bulk["gene"].astype(str).str.split(".").str[0]
bulk = bulk.drop_duplicates("g")[["g", "logFC"]]
u = pd.read_csv(AN1, sep="\t").merge(bulk, on="g", how="left")
assert u.logFC.notna().all()
u["gen"] = u.gen.astype(str).eq("True"); u["deg"] = u.deg.astype(str).eq("True")
base = run(u)
assert base["N"] == 3589 and base["G"] == 426 and base["D"] == 209 and base["O"] == 26, base
assert abs(base["E"] - 21.396) < 5e-4 and abs(base["MH_OR"] - 1.281506) < 1e-6 and abs(base["p_exact"] - 0.276987) < 1e-6, base

cmp = pd.read_csv(f"{RUN}/results/eqtl_dedup_gene_comparison.tsv", sep="\t")
unt = set(cmp.loc[cmp.new_state.str.startswith("untestable"), "ensembl"])
lost = set(cmp.loc[cmp.new_state.eq("lost"), "ensembl"])
assert len(unt) == 6 and len(lost) == 1 and unt | lost <= set(u.g)
new = u[~u.g.isin(unt)].copy()
new.loc[new.g.isin(lost), "gen"] = False
newpp = cmp.set_index("ensembl").new_best_susie.dropna()
hit = new.g.isin(newpp.index); new.loc[hit, "susie"] = new.loc[hit, "g"].map(newpp)
print("refit PP.H4 applied to", int(hit.sum()), "universe genes", file=sys.stderr); assert ((new.susie >= 0.5) == new.gen).all()
refit = run(new)
ub = pd.read_csv(AN1.replace("rule_a_sens", "rule_b_sens"), sep="\t")
ub["gen"] = ub.gen.astype(str).eq("True"); ub["deg"] = ub.deg.astype(str).eq("True")
base_b = run(ub); assert base_b["N"] == 3736 and base_b["O"] == 26, base_b
nb = ub[~ub.g.isin(unt)].copy(); nb.loc[nb.g.isin(lost), "gen"] = False
refit_b = run(nb)
os.makedirs(OUT, exist_ok=False)
pd.DataFrame([dict(set="adopted", **base), dict(set="eqtl_refit", **refit),
              dict(set="rule_b_adopted", **base_b), dict(set="rule_b_eqtl_refit", **refit_b)]).to_csv(f"{OUT}/overlap_test.tsv", sep="\t", index=False)
print(json.dumps(dict(adopted=base, refit=refit, rule_b=base_b, rule_b_refit=refit_b), indent=1))
