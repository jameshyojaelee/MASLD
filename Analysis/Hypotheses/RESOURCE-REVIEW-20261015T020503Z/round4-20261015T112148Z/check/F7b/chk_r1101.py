"""Blind check of R11-01. (1) label-only Spearman, participation ratio, C1/C2 on the 99 GSE267145 grades.
(2) independent simulation of the association-vector statistics (PR of the 3x3 Spearman matrix of
per-feature association vectors, and the C3 rule 'any pair |rho| < 0.80') under: no signal, one shared
axis, three programs each tracking one observed grade, three programs each tracking the part of its grade
not shared with the other two. 20,000 features, 4 reps per cell, seed 20260923."""
import json, sys
import numpy as np, pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
J = f"{ROOT}/Analysis/MASLD_Model_Benchmark/executions/gse267145-authoritative-join-21064930/participant_join.tsv"
OUT = f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b"
j = pd.read_csv(J, sep="\t")
L = j[["steatosis", "ballooning", "lobular_inflammation"]].apply(pd.to_numeric, errors="coerce")
assert L.notna().all().all() and j.participant_id.is_unique
n = len(L)
R = L.corr(method="spearman").to_numpy()
ev = np.linalg.eigvalsh(R)
pr_lab = ev.sum() ** 2 / (ev ** 2).sum()
off = R[np.triu_indices(3, 1)]
res = dict(n=n, rho_sb=off[0], rho_si=off[1], rho_bi=off[2], label_pr=pr_lab,
           closed_form=9 / (3 + 2 * (off ** 2).sum()), C1_pass=bool((np.abs(off) < 0.80).all()), C2_pass=bool(pr_lab >= 2.0))

def zs(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std()
Lz = np.column_stack([zs(stats.rankdata(L.iloc[:, k])) for k in range(3)])
orth = []
for k in range(3):
    o = [m for m in range(3) if m != k]
    B = np.column_stack([np.ones(n), Lz[:, o]])
    res_k = Lz[:, k] - B @ np.linalg.lstsq(B, Lz[:, k], rcond=None)[0]
    orth.append(zs(res_k))
Lo = np.column_stack(orth)
shared = zs(Lz.sum(1))

def stats_of(X):
    Xr = np.apply_along_axis(stats.rankdata, 0, X); Xr = (Xr - Xr.mean(0)) / Xr.std(0)
    A = Xr.T @ Lz / n
    C = stats.spearmanr(A).statistic
    e = np.linalg.eigvalsh(C)
    offc = np.abs(C[np.triu_indices(3, 1)])
    return e.sum() ** 2 / (e ** 2).sum(), offc.min(), bool((offc < 0.80).any())

rng = np.random.default_rng(20260923)
G = 20000
rows = []
for model in ["noise", "one_axis", "grade_programs", "orthogonal_part_programs"]:
    for frac in ([0.0] if model == "noise" else [0.02, 0.05, 0.1, 0.2, 0.4]):
        for beta in ([0.0] if model == "noise" else [0.5, 1.0]):
            for rep in range(4):
                X = rng.standard_normal((n, G)); k = int(frac * G)
                if model == "one_axis":
                    X[:, :k] += beta * shared[:, None]
                elif model in ("grade_programs", "orthogonal_part_programs"):
                    S = Lz if model == "grade_programs" else Lo
                    for a in range(3):
                        X[:, a * k // 3:(a + 1) * k // 3] += beta * S[:, [a]]
                p, mn, c3 = stats_of(X)
                rows.append(dict(model=model, frac=frac, beta=beta, rep=rep, pr=p, min_pair_rho=mn, c3_pass=c3))
t = pd.DataFrame(rows)
t.to_csv(f"{OUT}/chk_r1101_sim_reps.tsv", sep="\t", index=False)
s = t.groupby("model").agg(pr_min=("pr", "min"), pr_max=("pr", "max"), c3_pass_frac=("c3_pass", "mean"))
res["sim_summary"] = s.reset_index().to_dict(orient="records")
res["sim_by_cell"] = t.groupby(["model", "frac", "beta"]).pr.mean().round(4).reset_index().to_dict(orient="records")
json.dump(res, open(f"{OUT}/chk_r1101.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
