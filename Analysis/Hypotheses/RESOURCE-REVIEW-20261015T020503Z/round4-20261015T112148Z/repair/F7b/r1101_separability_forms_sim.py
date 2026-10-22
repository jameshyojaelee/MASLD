"""F7b repair, R11-01: which form of NASH-CRN component separability can the GSE267145
statistics detect?

Two statistics, both computed on per-feature association vectors with the 99 deposited grades:
  C3   (RNA prespecification, decisive criterion): Spearman between the per-feature
       association vectors for each pair of components; C3 passes if at least one pair has
       |rho| < 0.80.
  PR   (aspect-chromatin analysis): participation ratio of the 3x3 Spearman matrix of the
       association vectors.

Simulated forms of signal (expression = noise + signal on a fraction of features):
  noise           no signal
  one_axis        every signal feature follows the summed grade (one shared program)
  grade_following each third of the signal features follows one recorded grade
  grade_specific  each third follows the part of its grade not shared with the other two
                  (grade z-rank residualized on the other two)
  mixture         half the signal features one_axis, half grade_specific

Independent of review job 4041 (own code, seed 20260923). Feature counts 30,605 (RNA universe
of the prespecified evaluator) and 96,460 (H3K27ac regions). 4 replicates per cell.
"""
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions"
JOIN = f"{ROOT}/gse267145-authoritative-join-21064930/participant_join.tsv"
RNA_VEC = f"{ROOT}/model-check-1044-21181334-aspect-separability/separability/gene_association_vectors.tsv"
SEP = f"{ROOT}/model-check-1044-21181334-aspect-separability/separability/aspect_separability.json"
OUT = sys.argv[1]
ASP = ["steatosis", "ballooning", "lobular_inflammation"]
SEED, REPS = 20260923, 4
FRACS = [0.02, 0.05, 0.10, 0.20, 0.40]
BETAS = [0.5, 1.0]
MODELS = ["one_axis", "grade_following", "grade_specific", "mixture"]


def pr(m):
    ev = np.linalg.eigvalsh(m)
    return float(ev.sum() ** 2 / (ev ** 2).sum())


def zs(v):
    return (v - v.mean()) / v.std()


lab = pd.read_csv(JOIN, sep="\t")[ASP].astype(float).to_numpy()
n = lab.shape[0]
assert n == 99
Lr = np.column_stack([stats.rankdata(lab[:, k]) for k in range(3)])
Lz = (Lr - Lr.mean(0)) / Lr.std(0)
R = stats.spearmanr(lab).statistic
print("label rho", np.round(R[np.triu_indices(3, 1)], 4), "label PR", round(pr(R), 6))

# observed values, read back from the producing files
sep = json.load(open(SEP))
obs_c3 = sep["criteria"]["c3_aspects_rank_genes_differently"] if "c3_aspects_rank_genes_differently" in sep["criteria"] else None
rv = pd.read_csv(RNA_VEC, sep="\t")
rr = stats.spearmanr(rv[ASP].to_numpy()).statistic
print("RNA association vectors: features", len(rv), "pair rho", np.round(rr[np.triu_indices(3, 1)], 4),
      "PR", round(pr(rr), 4))

sev = zs(stats.rankdata(lab.sum(1)))
spec = np.zeros_like(Lz)
for k in range(3):
    X = np.column_stack([np.ones(n), np.delete(Lz, k, axis=1)])
    coef, *_ = np.linalg.lstsq(X, Lz[:, k], rcond=None)
    spec[:, k] = zs(Lz[:, k] - X @ coef)


def stats_of(expr):
    """expr: features x donors, continuous. Returns min pair |rho|, the three pairs, PR."""
    rk = np.argsort(np.argsort(expr, axis=1), axis=1).astype(float)
    rk = (rk - rk.mean(1, keepdims=True)) / rk.std(1, keepdims=True)
    assoc = rk @ Lz / n  # Spearman of each feature with each grade (grades rank-averaged)
    cm = stats.spearmanr(assoc).statistic
    pairs = cm[np.triu_indices(3, 1)]
    return float(np.abs(pairs).min()), pairs, pr(cm)


rng = np.random.default_rng(SEED)
rows = []
t0 = time.time()
for G in [30605, 96460]:
    for rep in range(REPS):
        e = rng.standard_normal((G, n))
        mn, pairs, p = stats_of(e)
        rows.append(dict(G=G, model="noise", frac=0.0, beta=0.0, rep=rep, c3_min_abs_rho=mn,
                         rho_sb=pairs[0], rho_si=pairs[1], rho_bi=pairs[2], pr=p))
    for model in MODELS:
        for f in FRACS:
            for b in BETAS:
                for rep in range(REPS):
                    e = rng.standard_normal((G, n))
                    k = int(round(f * G))
                    sign = rng.choice([-1.0, 1.0], size=k)
                    which = np.arange(k) % 3
                    if model == "one_axis":
                        drv = np.tile(sev, (k, 1))
                    elif model == "grade_following":
                        drv = Lz[:, which].T
                    elif model == "grade_specific":
                        drv = spec[:, which].T
                    else:  # mixture
                        h = k // 2
                        drv = np.vstack([np.tile(sev, (h, 1)), spec[:, which[h:]].T])
                    e[:k] += b * sign[:, None] * drv
                    mn, pairs, p = stats_of(e)
                    rows.append(dict(G=G, model=model, frac=f, beta=b, rep=rep, c3_min_abs_rho=mn,
                                     rho_sb=pairs[0], rho_si=pairs[1], rho_bi=pairs[2], pr=p))
        print(f"G={G} {model} done {time.time() - t0:.0f}s", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/r1101_forms_sim_draws.tsv", sep="\t", index=False, float_format="%.6f")
summ = (df.groupby(["G", "model", "frac", "beta"])
          .agg(c3_min=("c3_min_abs_rho", "min"), c3_max=("c3_min_abs_rho", "max"),
               pr_min=("pr", "min"), pr_max=("pr", "max"), pr_mean=("pr", "mean"))
          .reset_index())
summ.to_csv(f"{OUT}/r1101_forms_sim_summary.tsv", sep="\t", index=False, float_format="%.4f")
pd.set_option("display.width", 220)
print(summ.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print("\nRange over the whole grid, by feature count and model:")
rng_t = df.groupby(["G", "model"]).agg(c3_lo=("c3_min_abs_rho", "min"), c3_hi=("c3_min_abs_rho", "max"),
                                         pr_lo=("pr", "min"), pr_hi=("pr", "max"),
                                         c3_pass_frac=("c3_min_abs_rho", lambda x: float((x < 0.80).mean())))
print(rng_t.to_string(float_format=lambda x: f"{x:.4f}"))
rng_t.to_csv(f"{OUT}/r1101_forms_sim_ranges.tsv", sep="\t", float_format="%.4f")
print(f"total {time.time() - t0:.0f}s")
