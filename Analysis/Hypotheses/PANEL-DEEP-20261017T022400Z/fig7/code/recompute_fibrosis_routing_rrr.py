"""Recompute the 7B graded-fibrosis, routing and RNA->chromatin rows, the 7D routing block,
and the 7E routing ratio, from per-sample predictions, saved bootstrap indices and per-fold skills.

Inputs (analysis outputs):
  pooled_graded_fibrosis_ladder-20260828T114911Z/predictions_per_sample.tsv, bootstrap_indices.npz
  routed_multihead_eval-20260829T230000Z/predictions_per_sample.tsv, results.json (floor only)
  gse267145-crossmodal-ceiling-20260830/results_rna2h3_v2_20260830T170000/rna_conditioned_results.json (per-fold skill)
  gse296875-crossmodal-recurrence-20260830T193000Z/results/recurrence_results.json (jackknife se, LOWO arm)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as tdist

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUT = Path(sys.argv[1])
COH = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE240729"]
res = {}


def sp(a, b):
    ra = rankdata(a); rb = rankdata(b); ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return np.nan if d == 0 else float((ra * rb).sum() / d)


# ---------------- graded fibrosis ladder
L = pd.read_csv(EX / "pooled_graded_fibrosis_ladder-20260828T114911Z/predictions_per_sample.tsv", sep="\t")
IDX = np.load(EX / "pooled_graded_fibrosis_ladder-20260828T114911Z/bootstrap_indices.npz")
rungs = ["a", "c", "d"]
per = {k: {c: sp(L.loc[L.cohort == c, f"pred_{k}"], L.loc[L.cohort == c, "fibrosis_grade"]) for c in COH} for k in rungs}
macro = {k: float(np.mean(list(per[k].values()))) for k in rungs}
boot = {k: np.zeros((10000, 5)) for k in rungs}
for j, c in enumerate(COH):
    sub = L[L.cohort == c].reset_index(drop=True)
    I = IDX[c]
    y = sub["fibrosis_grade"].to_numpy()[I]
    ry = rankdata(y, axis=1); ry -= ry.mean(1, keepdims=True)
    for k in rungs:
        p = sub[f"pred_{k}"].to_numpy()[I]
        rp = rankdata(p, axis=1); rp -= rp.mean(1, keepdims=True)
        den = np.sqrt((rp ** 2).sum(1) * (ry ** 2).sum(1))
        boot[k][:, j] = np.where(den == 0, 0.0, (rp * ry).sum(1) / np.where(den == 0, 1, den))
bm = {k: boot[k].mean(1) for k in rungs}
fib = {"n_samples": int(len(L)), "n_per_cohort": L.cohort.value_counts().reindex(COH).astype(int).to_dict(),
       "macro": macro, "per_cohort": per}
for cid, (x, y_, bh_rank) in {"C1_c_minus_a": ("c", "a", 2), "C2_c_minus_d": ("c", "d", 1)}.items():
    dv = bm[x] - bm[y_]
    lvl = 0.05 * bh_rank / 3
    fib[cid] = {"point": macro[x] - macro[y_], "ci95": [float(np.quantile(dv, .025)), float(np.quantile(dv, .975))],
                "ci_bh": [float(np.quantile(dv, lvl / 2)), float(np.quantile(dv, 1 - lvl / 2))],
                "per_cohort": {c: per[x][c] - per[y_][c] for c in COH},
                "cohorts_positive": int(sum(per[x][c] - per[y_][c] > 0 for c in COH))}
res["graded_fibrosis"] = fib

# ---------------- routing
R = pd.read_csv(EX / "routed_multihead_eval-20260829T230000Z/predictions_per_sample.tsv", sep="\t")
rr = json.load(open(EX / "routed_multihead_eval-20260829T230000Z/results.json"))
arms = ["ref", "routed", "d", "c", "cont", "a"]
perR = {k: {c: sp(R.loc[R.cohort == c, f"pred_{k}"], R.loc[R.cohort == c, "fibrosis_grade"]) for c in COH} for k in arms}
macR = {k: float(np.mean(list(perR[k].values()))) for k in arms}
oracle = float(np.mean([max(perR[k][c] for k in ["ref", "c", "d", "cont"]) for c in COH]))
res["routing"] = {"n_samples": int(len(R)), "macro": macR, "routed_equals_ref_bitwise": bool(np.array_equal(R.pred_routed.to_numpy(), R.pred_ref.to_numpy())),
                  "routed_minus_ref": macR["routed"] - macR["ref"], "oracle_macro_over_ref_c_d_cont": oracle,
                  "oracle_gain_over_ref": oracle - macR["ref"], "floor": rr["floor"],
                  "oracle_gain_over_floor_ratio": (oracle - macR["ref"]) / rr["floor"]}

# ---------------- RNA -> H3K27ac (GSE267145) and transfer to ATAC (GSE296875)
h3 = json.load(open(EX / "gse267145-crossmodal-ceiling-20260830/results_rna2h3_v2_20260830T170000/rna_conditioned_results.json"))
pf = pd.DataFrame(h3["per_fold"])
res["rna_to_h3k27ac"] = {"skill": h3["skill"], "n": h3["n_scored"], "fold_min": float(pf.skill.min()), "fold_max": float(pf.skill.max()),
                         "fold_mean": float(pf.skill.mean()), "fold_n_weighted_mean": float((pf.skill * pf.n).sum() / pf.n.sum()),
                         "fold_t_ci95": [float(pf.skill.mean() - tdist.ppf(.975, 4) * pf.skill.std(ddof=1) / np.sqrt(5)),
                                         float(pf.skill.mean() + tdist.ppf(.975, 4) * pf.skill.std(ddof=1) / np.sqrt(5))]}
at = json.load(open(EX / "gse296875-crossmodal-recurrence-20260830T193000Z/results/recurrence_results.json"))
jk = at["primary"]["donor_jackknife"]
res["rna_to_atac_transfer"] = {"skill": at["primary"]["skill"], "donors": at["donors"], "jk_se": jk["se"], "jk_pseudomean": jk["pseudovalue_mean"],
                               "ci_z": [jk["pseudovalue_mean"] - 1.959964 * jk["se"], jk["pseudovalue_mean"] + 1.959964 * jk["se"]],
                               "ci_t38": [jk["pseudovalue_mean"] - tdist.ppf(.975, 38) * jk["se"], jk["pseudovalue_mean"] + tdist.ppf(.975, 38) * jk["se"]],
                               "recorded_ci": [jk["ci95_low"], jk["ci95_high"]],
                               "lowo_skill": at["secondary_leave_one_well_out"]["skill"],
                               "lowo_ci": [at["secondary_leave_one_well_out"]["donor_jackknife"]["ci95_low"], at["secondary_leave_one_well_out"]["donor_jackknife"]["ci95_high"]],
                               "fold_min": at["primary"]["fold_skill_min"], "fold_max": at["primary"]["fold_skill_max"]}
json.dump(res, open(OUT / "fibrosis_routing_rrr_recomputed.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
