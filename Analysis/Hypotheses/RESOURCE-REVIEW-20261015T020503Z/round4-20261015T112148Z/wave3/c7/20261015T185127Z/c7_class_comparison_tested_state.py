#!/usr/bin/env python3
"""Conflict 7: hepatocyte-specificity class comparison restricted to state-only
genes that multi-signal COLOC evaluated.

Post hoc sensitivity analysis, outside any prespecified family.

Inputs (read only):
  - gene classes, AveExpr and hep_specificity from the published 774-unit /
    COLOC r2 class comparison (job 22236311):
    figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results/gene_classes.tsv.gz
  - published pairs and sign-flip diagnostic from the same directory
  - COLOC r2 evaluation state (positive / tested_le_0.5) from the F1 repair table:
    figures/candidates/review-r3-fixes-20261015T112148Z/F1/repair/fig1b_U_eval_gene_table.tsv

Matching function and sign-flip null are copied from 03_run_class_test.py and
04_signflip_diagnostic.py in that directory, unchanged.

Arms:
  A  reproduction of the published 281-pair comparison (sign-flip seed 20260907)
  B  same classes (COLOC r2, 400 / 1,129), sign-flip seed 20260923
  C  genetic-only vs state-only genes tested by multi-signal COLOC (PP.H4 <= 0.5),
     caliper recomputed by the published rule (0.05 x pooled SD), seed 20260923
  C2 as C, caliper fixed at the published arm-A caliper value, seed 20260923
  D  published 281 pairs split by whether the state partner was COLOC-tested
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

PROJ = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
CD = PROJ/"figures/candidates/bulk-context-774-20261014T140904Z_22236311/class_difference/results"
F1 = PROJ/"figures/candidates/review-r3-fixes-20261015T112148Z/F1/repair/fig1b_U_eval_gene_table.tsv"
OUT = Path(sys.argv[1])
CALIPER_SD = 0.05; N_DRAWS = 100_000
SEED_PUB = 20260907; SEED_NEW = 20260923
OUTCOME = "hep_specificity"

# ---- copied unchanged from 03_run_class_test.py
def match_hardest_first(a,b,cal,rng=None):
    a=a.sort_values("g").reset_index(drop=True); b=b.sort_values("g").reset_index(drop=True)
    av=a.AveExpr.to_numpy(); bv=b.AveExpr.to_numpy()
    npart=np.array([np.sum(np.abs(bv-x)<=cal) for x in av])
    order=np.argsort(npart, kind="stable")
    used=np.zeros(len(b),bool); pa=[]; pb=[]; ea=[]; eb=[]
    for i in order:
        d=np.abs(bv-av[i]); d[used]=np.inf
        j=int(np.argmin(d))
        if d[j]<=cal:
            used[j]=True; pa.append(a.g[i]); pb.append(b.g[j]); ea.append(av[i]); eb.append(bv[j])
    return pd.DataFrame({"g_genetic":pa,"g_state":pb,"expr_genetic":ea,"expr_state":eb})

# ---- copied from 04_signflip_diagnostic.py (one arm)
def signflip(d, seed):
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(N_DRAWS, len(d)))
    null = (signs * d).mean(axis=1)
    obs = d.mean()
    n_pos, n_neg = int((d > 0).sum()), int((d < 0).sum())
    return {"observed_mean_difference": float(obs),
            "signflip_null_sd": float(null.std(ddof=1)),
            "signflip_p_two_sided": float((1 + np.sum(np.abs(null) >= abs(obs))) / (1 + N_DRAWS)),
            "wilcoxon_p": float(stats.wilcoxon(d).pvalue),
            "sign_test_positive": n_pos, "sign_test_negative": n_neg,
            "sign_test_p": float(stats.binomtest(n_pos, n_pos + n_neg, 0.5).pvalue),
            "seed": seed, "n_draws": N_DRAWS}

cls = pd.read_csv(CD/"gene_classes.tsv.gz", sep="\t")
f1 = pd.read_csv(F1, sep="\t")[["g", "susie_state_t12", "AveExpr"]].rename(columns={"AveExpr": "AveExpr_f1"})
cls = cls.merge(f1, on="g", how="left")
# the two tables are the same 774-unit fit; refuse to continue if they disagree
both = cls.dropna(subset=["AveExpr_f1"])
assert float((both.AveExpr - both.AveExpr_f1).abs().max()) == 0.0
cls["coloc_state"] = cls["susie_state_t12"].fillna("no_multisignal_result")

counts = pd.crosstab(cls.klass, cls.coloc_state)
counts_sc = pd.crosstab(cls.klass[cls[OUTCOME].notna()], cls.coloc_state[cls[OUTCOME].notna()])

def run(name, a_mask, b_mask, seed, caliper_abs=None):
    s = cls[cls.AveExpr.notna() & cls[OUTCOME].notna()]
    a = s.loc[a_mask.reindex(s.index), ["g", "AveExpr"]]
    b = s.loc[b_mask.reindex(s.index), ["g", "AveExpr"]]
    sd = np.sqrt((a.AveExpr.var(ddof=1) + b.AveExpr.var(ddof=1)) / 2)
    cal = CALIPER_SD * sd if caliper_abs is None else caliper_abs
    m = match_hardest_first(a, b, cal)
    val = s.set_index("g")[OUTCOME]
    m["y_genetic"] = m.g_genetic.map(val); m["y_state"] = m.g_state.map(val)
    m = m.dropna(subset=["y_genetic", "y_state"])
    d = (m.y_genetic - m.y_state).to_numpy()
    res = {"arm": name, "n_genetic_pool": len(a), "n_state_pool": len(b),
           "mean_AveExpr_genetic": float(a.AveExpr.mean()),
           "mean_AveExpr_state": float(b.AveExpr.mean()),
           "SMD_before_matching": float((a.AveExpr.mean() - b.AveExpr.mean()) / sd),
           "pooled_SD": float(sd), "caliper_abs": float(cal),
           "n_pairs": len(d),
           "SMD_after_matching": float((m.expr_genetic - m.expr_state).mean() / sd),
           "mean_hep_spec_genetic_in_pairs": float(m.y_genetic.mean()),
           "mean_hep_spec_state_in_pairs": float(m.y_state.mean()),
           "median_paired_difference": float(np.median(d)),
           "cohen_d_paired": float(d.mean() / d.std(ddof=1)),
           "MDE_d80": float((1.959964 + 0.8416212) / np.sqrt(len(d)))}
    res.update(signflip(d, seed))
    m.to_csv(OUT/f"pairs_{name}.tsv", sep="\t", index=False)
    return res, m, cal

gen = cls.klass.eq("genetic_only")
st_all = cls.klass.eq("state_only")
st_tested = st_all & cls.coloc_state.eq("tested_le_0.5")

rows = []
rA, mA, calA = run("A_reproduce_published_seed20260907", gen, st_all, SEED_PUB)
pub = pd.read_csv(CD/"pairs_singlecell_hep_specificity.tsv", sep="\t")
pairs_identical = bool(mA[["g_genetic", "g_state"]].reset_index(drop=True)
                       .equals(pub[["g_genetic", "g_state"]].reset_index(drop=True)))
pub_diag = pd.read_csv(CD/"corrected_null_diagnostic.tsv", sep="\t").set_index("arm").loc["singlecell_hep_specificity"]
rA["pairs_identical_to_published"] = pairs_identical
rA["published_signflip_p"] = float(pub_diag.signflip_p_two_sided)
rA["published_mean"] = float(pub_diag.observed_mean_difference)
rows.append(rA)
rB, _, _ = run("B_coloc_r2_classes_seed20260923", gen, st_all, SEED_NEW); rows.append(rB)
rC, mC, _ = run("C_tested_state_only_seed20260923", gen, st_tested, SEED_NEW); rows.append(rC)
rC2, _, _ = run("C2_tested_state_only_fixed_published_caliper_seed20260923", gen, st_tested, SEED_NEW, caliper_abs=calA)
rows.append(rC2)

# D: split the published pairs by the state partner's COLOC evaluation
state_of = cls.set_index("g").coloc_state
mA["state_partner_coloc"] = mA.g_state.map(state_of)
for lab, sub in mA.groupby("state_partner_coloc"):
    d = (sub.y_genetic - sub.y_state).to_numpy()
    r = {"arm": f"D_published_pairs_partner_{lab}", "n_pairs": len(d),
         "mean_AveExpr_state": float(sub.expr_state.mean()),
         "mean_hep_spec_genetic_in_pairs": float(sub.y_genetic.mean()),
         "mean_hep_spec_state_in_pairs": float(sub.y_state.mean()),
         "median_paired_difference": float(np.median(d))}
    if len(d) > 5:
        r.update(signflip(d, SEED_NEW))
    else:
        r["observed_mean_difference"] = float(d.mean())
    rows.append(r)
mA.to_csv(OUT/"pairs_A_with_partner_coloc_state.tsv", sep="\t", index=False)

R = pd.DataFrame(rows)
R["analysis_status"] = "post_hoc_sensitivity_outside_prespecified_family"
R["unit"] = "gene pairs (not donors; genes share loci and are not independent)"
R.to_csv(OUT/"c7_results.tsv", sep="\t", index=False)
counts.to_csv(OUT/"class_by_coloc_state_14704.tsv", sep="\t")
counts_sc.to_csv(OUT/"class_by_coloc_state_singlecell_observable.tsv", sep="\t")
json.dump({"inputs": [str(CD/"gene_classes.tsv.gz"), str(CD/"pairs_singlecell_hep_specificity.tsv"),
                      str(CD/"corrected_null_diagnostic.tsv"), str(F1)],
           "seeds": {"reproduction": SEED_PUB, "new": SEED_NEW}, "n_draws": N_DRAWS,
           "caliper_rule": "0.05 x sqrt((var_genetic + var_state)/2) of AveExpr within each arm"},
          open(OUT/"c7_params.json", "w"), indent=1)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(counts); print(counts_sc)
print(R.T.to_string())
