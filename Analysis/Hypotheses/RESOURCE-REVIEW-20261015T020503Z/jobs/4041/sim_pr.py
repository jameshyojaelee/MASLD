"""R11-01 check: can the association-vector participation ratio (PR) tell
separable NAS components from one shared axis, on the deposited 99 labels?

Reads only. Writes to this job directory.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
G = int(sys.argv[2]) if len(sys.argv) > 2 else 30605
REPS = int(sys.argv[3]) if len(sys.argv) > 3 else 5
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions")
JOIN = ROOT / "gse267145-authoritative-join-21064930/participant_join.tsv"
RNA_VEC = ROOT / "model-check-1044-21181334-aspect-separability/separability/gene_association_vectors.tsv"
ASPECTS = ["steatosis", "ballooning", "lobular_inflammation"]
SEED = 20261015


def pr(m):
    ev = np.linalg.eigvalsh(m)
    return float(ev.sum() ** 2 / (ev ** 2).sum()), sorted(ev.tolist(), reverse=True)


def zcols(a):
    return (a - a.mean(axis=0)) / a.std(axis=0, ddof=0)


def assoc_pr(expr, label_z):
    """expr: genes x donors (continuous). label_z: donors x K z-scored label ranks."""
    ranks = np.argsort(np.argsort(expr, axis=1), axis=1).astype(float)  # no ties: continuous
    rz = (ranks - ranks.mean(axis=1, keepdims=True)) / ranks.std(axis=1, ddof=0, keepdims=True)
    assoc = rz @ label_z / expr.shape[1]
    ranked = np.column_stack([rankdata(assoc[:, k]) for k in range(assoc.shape[1])])
    cm = np.corrcoef(ranked, rowvar=False)
    ratio, ev = pr(cm)
    pairs = [cm[0, 1], cm[0, 2], cm[1, 2]]
    return ratio, ev, pairs


def main():
    t0 = time.time()
    join = pd.read_csv(JOIN, sep="\t", dtype=str)
    lab = join[ASPECTS].astype(float).to_numpy()
    n = lab.shape[0]
    print(f"participants={n} unique={join['participant_id'].nunique()}")
    m = np.eye(3)
    for i in range(3):
        for j in range(i + 1, 3):
            m[i, j] = m[j, i] = spearmanr(lab[:, i], lab[:, j])[0]
    lpr, lev = pr(m)
    print("label Spearman: s|b={:.4f} s|i={:.4f} b|i={:.4f}".format(m[0, 1], m[0, 2], m[1, 2]))
    print(f"label-only PR={lpr:.6f} eigen={[round(x, 4) for x in lev]}")
    print("C1 pass (all |rho|<0.80):", bool(np.all(np.abs([m[0, 1], m[0, 2], m[1, 2]]) < 0.80)))
    print("C2 pass (PR>=2.0):", lpr >= 2.0)

    rv = pd.read_csv(RNA_VEC, sep="\t")
    rr = np.column_stack([rankdata(rv[a].to_numpy()) for a in ASPECTS])
    rpr, rev = pr(np.corrcoef(rr, rowvar=False))
    print(f"RNA deposited gene association vectors: genes={len(rv)} PR={rpr:.4f} eigen={[round(x, 4) for x in rev]}")

    label_z = zcols(np.column_stack([rankdata(lab[:, k]) for k in range(3)]))
    severity = rankdata(lab.sum(axis=1))
    sev_z = (severity - severity.mean()) / severity.std()
    # orthogonal parts: each aspect's z-rank residualized on the other two
    orth = np.zeros_like(label_z)
    for k in range(3):
        others = np.column_stack([np.ones(n), np.delete(label_z, k, axis=1)])
        coef, *_ = np.linalg.lstsq(others, label_z[:, k], rcond=None)
        r = label_z[:, k] - others @ coef
        orth[:, k] = (r - r.mean()) / r.std()

    rng = np.random.default_rng(SEED)
    rows = []
    fracs = [0.02, 0.05, 0.10, 0.20, 0.40]
    betas = [0.25, 0.5, 1.0, 2.0]
    models = ["one_axis", "separable_labels", "separable_orthogonal_parts"]
    # pure-noise reference (this is the joint donor-label null geometry)
    for rep in range(REPS):
        expr = rng.standard_normal((G, n))
        ratio, ev, pairs = assoc_pr(expr, label_z)
        rows.append(dict(model="noise_only", frac=0.0, beta=0.0, rep=rep, pr=ratio,
                         rho_sb=pairs[0], rho_si=pairs[1], rho_bi=pairs[2]))
    for model in models:
        for f in fracs:
            for b in betas:
                for rep in range(REPS):
                    expr = rng.standard_normal((G, n))
                    n_sig = int(round(f * G))
                    sign = rng.choice([-1.0, 1.0], size=n_sig)
                    if model == "one_axis":
                        expr[:n_sig] += b * sign[:, None] * sev_z[None, :]
                    else:
                        drivers = label_z if model == "separable_labels" else orth
                        which = np.arange(n_sig) % 3  # disjoint, equal-size aspect-specific sets
                        expr[:n_sig] += b * sign[:, None] * drivers[:, which].T
                    ratio, ev, pairs = assoc_pr(expr, label_z)
                    rows.append(dict(model=model, frac=f, beta=b, rep=rep, pr=ratio,
                                     rho_sb=pairs[0], rho_si=pairs[1], rho_bi=pairs[2]))
        print(f"{model} done at {time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sim_pr_draws.tsv", sep="\t", index=False)
    df["min_pair_abs_rho"] = df[["rho_sb", "rho_si", "rho_bi"]].abs().min(axis=1)
    summ = df.groupby(["model", "frac", "beta"]).agg(
        pr_mean=("pr", "mean"), pr_min=("pr", "min"), pr_max=("pr", "max"),
        min_pair_abs_rho_mean=("min_pair_abs_rho", "mean")).reset_index()
    summ.to_csv(OUT / "sim_pr_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 200)
    print(summ.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nPR range by model over the whole grid:")
    print(df.groupby("model")["pr"].agg(["min", "max", "mean"]).to_string(float_format=lambda x: f"{x:.4f}"))
    print("\nC3 (any vector pair |rho|<0.80) reachable, fraction of draws by model:")
    print((df.assign(c3=df["min_pair_abs_rho"] < 0.80).groupby("model")["c3"].mean()).to_string())
    print("\nmin pair |rho| range by model:")
    print(df.groupby("model")["min_pair_abs_rho"].agg(["min", "max"]).to_string(float_format=lambda x: f"{x:.4f}"))
    json.dump(dict(label_pr=lpr, label_rho=[m[0, 1], m[0, 2], m[1, 2]], rna_pr=rpr, G=G, reps=REPS, seed=SEED),
              open(OUT / "sim_pr_meta.json", "w"), indent=1)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
