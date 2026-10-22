"""Independent recomputation of the endogenous eQTL rows (Figure 7B |z| and signed rows,
7D signed-eQTL AUC block and LD-block join, 7E signed-eQTL ratio).

Inputs (analysis outputs, not figure tables):
  dnalm_heads-20260831T183647Z/tables/analysis_rows.parquet, oof_{hyenadna,nucleotide_transformer}.npz,
  bootstrap_draws.npz (producer's paired draws, used only as a second check)
  eqtl_sign_reliability-20260831T195702Z/tables/oof_*.npz (stratum refits)
  variant_program_substrate_rebuild-20260830T160000Z/tables/training_matrix.parquet + splits (LD-block join)
Held-out set: ld_fold == 4 (55 LD blocks). Bootstrap: own seed 7001, B = 2000, LD block unit.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

EX = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
OUT = Path(sys.argv[1])
B = 2000
rng = np.random.default_rng(7001)


def spear(x, y):
    a, b = rankdata(x), rankdata(y)
    a -= a.mean(); b -= b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else np.nan


def auc(score, lab):
    r = rankdata(score)
    n1 = lab.sum(); n0 = len(lab) - n1
    return float((r[lab == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def q(v, lo=0.025, hi=0.975):
    v = v[np.isfinite(v)]
    return float(np.quantile(v, lo)), float(np.quantile(v, hi))


A = pd.read_parquet(EX / "dnalm_heads-20260831T183647Z/tables/analysis_rows.parquet")
te = np.nonzero((A["ld_fold"] == 4).to_numpy())[0]
blk = A["block"].to_numpy()[te]
ub = np.unique(blk)
by = {b: np.nonzero(blk == b)[0] for b in ub}
draws = [np.concatenate([by[b] for b in rng.choice(ub, size=len(ub), replace=True)]) for _ in range(B)]
yz = A["abs_z"].to_numpy()[te]
yb = A["beta_alt"].to_numpy()[te]
lab = (yb > 0).astype(int)
res = {"heldout_rows": int(len(te)), "heldout_blocks": int(len(ub)), "heldout_variants": int(A["variant_key"].iloc[te].nunique())}

oof = {M: np.load(EX / f"dnalm_heads-20260831T183647Z/tables/oof_{M}.npz") for M in ["hyenadna", "nucleotide_transformer"]}
rows_7b = []
for M, z in oof.items():
    b1 = z["b1_dist__abs_z"][te]; m1 = z["m1_dist_emb__abs_z"][te]
    pt = spear(m1, yz) - spear(b1, yz)
    bd = np.array([spear(m1[d], yz[d]) - spear(b1[d], yz[d]) for d in draws])
    lo, hi = q(bd)
    rows_7b.append({"row": f"{M} |z| m1_dist_emb - b1_dist", "point": pt, "lo": lo, "hi": hi,
                    "b1_spearman": spear(b1, yz), "m1_spearman": spear(m1, yz)})
    # signed direction: sign AUC m1 vs b1, and vs chance 0.5
    s1 = z["m1_dist_emb__beta_alt"][te]; sb = z["b1_dist__beta_alt"][te]
    a1 = auc(s1, lab); ab = auc(sb, lab)
    bd2 = np.array([auc(s1[d], lab[d]) - auc(sb[d], lab[d]) for d in draws])
    bd3 = np.array([auc(s1[d], lab[d]) for d in draws])
    rows_7b.append({"row": f"{M} signAUC m1_dist_emb - b1_dist", "point": a1 - ab, "lo": q(bd2)[0], "hi": q(bd2)[1],
                    "m1_auc": a1, "b1_auc": ab, "m1_auc_lo": q(bd3)[0], "m1_auc_hi": q(bd3)[1]})

# producer's saved draws as a second check of the interval step
sd = np.load(EX / "dnalm_heads-20260831T183647Z/tables/bootstrap_draws.npz")
for M in ["hyenadna", "nucleotide_transformer"]:
    dv = sd[f"abs_z__{M}::m1_dist_emb"] - sd[f"abs_z__{M}::b1_dist"]
    rows_7b.append({"row": f"{M} |z| from producer draws", "point": np.nan, "lo": q(dv)[0], "hi": q(dv)[1]})
    if f"sign_auc__{M}::m1_dist_emb" in sd.files:
        dv = sd[f"sign_auc__{M}::m1_dist_emb"] - sd[f"sign_auc__{M}::b1_dist"]
        rows_7b.append({"row": f"{M} signAUC m1-b1 from producer draws", "point": np.nan, "lo": q(dv)[0], "hi": q(dv)[1]})
pd.DataFrame(rows_7b).to_csv(OUT / "eqtl_7b_rows_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(pd.DataFrame(rows_7b).to_string())

# ---------------- strata (7D signed eQTL AUC block, 7E ratio)
absz_all = A["abs_z"].to_numpy(); beta_all = A["beta_alt"].to_numpy(); fold_all = A["ld_fold"].to_numpy()
blk_all = A["block"].to_numpy(); ensg_all = A["ENSG"].to_numpy()
SO = {M: np.load(EX / f"eqtl_sign_reliability-20260831T195702Z/tables/oof_{M}.npz") for M in ["hyenadna", "nucleotide_transformer"]}
r_flat_all = norm.cdf(absz_all)


def ceiling_flat(r, y):
    Aq = float(r[y == 1].sum()); Bq = float((1 - r)[y == 1].sum())
    Cq = float((1 - r)[y == 0].sum()); Dq = float(r[y == 0].sum())
    return float((Aq * Dq + 0.5 * (Aq * Cq + Bq * Dq)) / ((Aq + Bq) * (Cq + Dq)))


strata_rows = []
for sname, thr in [("all", 0), ("z_ge_1", 1), ("z_ge_2", 2), ("z_ge_3", 3), ("z_ge_5", 5)]:
    pos = SO["hyenadna"][f"__index__{sname}"]
    assert np.array_equal(pos, SO["nucleotide_transformer"][f"__index__{sname}"])
    tl = np.nonzero(fold_all[pos] == 4)[0]
    g = pos[tl]
    assert np.all(absz_all[g] >= thr)
    lb = (beta_all[g] > 0).astype(int)
    bk = blk_all[g]; u = np.unique(bk); byb = {b: np.nonzero(bk == b)[0] for b in u}
    dr = [np.concatenate([byb[b] for b in rng.choice(u, size=len(u), replace=True)]) for _ in range(B)]
    heads = {}
    for M, z in SO.items():
        for fs in ["e_emb", "m1_dist_emb"]:
            heads[f"{M}::{fs}"] = z[f"{sname}__{fs}__beta_alt"][tl]
    pts = {k: auc(v, lb) for k, v in heads.items()}
    best = max(pts, key=pts.get)
    bd = np.array([auc(heads[best][d], lb[d]) for d in dr])
    lo, hi = q(bd)
    # within-block label permutation chance with a fixed random score (200 draws)
    order = np.argsort(bk, kind="stable"); segs = np.split(order, np.nonzero(np.diff(bk[order]))[0] + 1)
    rk = rankdata(rng.standard_normal(len(g))); n1 = lb.sum(); n0 = len(lb) - n1
    pa = []
    for _ in range(200):
        yp = lb.copy()
        for sg in segs:
            if len(sg) > 1:
                yp[sg] = lb[rng.permutation(sg)]
        pa.append((rk[yp == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    chance = float(np.median(pa))
    cflat = ceiling_flat(r_flat_all[g], lb)
    grp = pd.DataFrame({"g": ensg_all[g], "b": bk}).drop_duplicates()
    f_causal = len(grp) / len(g)
    seq_bound = 0.5 + f_causal * (cflat - 0.5)
    res_half = (hi - lo) / 2
    strata_rows.append({"stratum": sname, "heldout_rows": len(g), "heldout_blocks": len(u), "best_head": best,
                        "best_auc": pts[best], "auc_lo": lo, "auc_hi": hi, "resolution_half_width": res_half,
                        "chance_median": chance, "C_flat": cflat, "headroom_flat": cflat - chance,
                        "seq_bound": seq_bound, "seq_bound_minus_chance": seq_bound - chance,
                        "ratio_best_minus_chance_over_halfwidth": (pts[best] - chance) / res_half,
                        **{f"auc::{k}": v for k, v in pts.items()}})
st = pd.DataFrame(strata_rows)
st.to_csv(OUT / "eqtl_strata_recomputed.tsv", sep="\t", index=False, float_format="%.6f")
print(st.T.to_string())

# ---------------- LD-block join (7D): gene-mean predictor with own block left in vs held out
FIX = EX / "variant_program_substrate_rebuild-20260830T160000Z"
full = pd.read_parquet(FIX / "tables/training_matrix.parquet", columns=["variant_key", "ENSG", "beta_alt"]).merge(
    pd.read_csv(FIX / "splits/variant_splits.tsv", sep="\t")[["variant_key", "block", "ld_fold"]], on="variant_key", how="inner")
ftr = full[(full["ld_fold"] != 4) & full["ENSG"].notna()].copy()
gs = ftr.groupby("ENSG")["beta_alt"].agg(["sum", "count"])
gb = ftr.groupby(["ENSG", "block"])["beta_alt"].agg(["sum", "count"]).rename(columns={"sum": "bsum", "count": "bcount"})
x = ftr.join(gs, on="ENSG").join(gb, on=["ENSG", "block"])
left_in = x["sum"] / x["count"]
held = (x["sum"] - x["bsum"]) / (x["count"] - x["bcount"])
ok = held.notna() & np.isfinite(held)
r_in = float(np.corrcoef(left_in, x["beta_alt"])[0, 1])
r_out = float(np.corrcoef(held[ok], x.loc[ok, "beta_alt"])[0, 1])
res.update({"ld_join_left_in_pearson": r_in, "ld_join_held_out_pearson": r_out,
            "ld_join_rows_left_in": int(len(x)), "ld_join_rows_held_out": int(ok.sum())})
json.dump(res, open(OUT / "eqtl_misc_recomputed.json", "w"), indent=1)
print(res)
