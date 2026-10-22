"""Validity checks for candidate L03-04 (login node, small files only)."""
import numpy as np
import pandas as pd
from scipy.linalg import helmert

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/"
P = ROOT + "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/"
NEW = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = ROOT + "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv"
MAN = ROOT + "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
rng = np.random.default_rng(20261015)

dz = pd.read_csv(P + "model_design.tsv", sep="\t").set_index("sample_id")
man = pd.read_csv(MAN, sep="\t").set_index("sample_id")
new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]
old = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS])
ids = dz.index
new = new.loc[ids, ct]
old = old.loc[ids, ct]
man = man.loc[ids]
print("design columns:", list(dz.columns))
print("manifest columns:", list(man.columns))

# 1. near-zero parts
print("\n1. units (of 774) below floor, rebuilt / earlier")
for c in ct:
    print(f"  {c:25s} <1e-6 {int((new[c] < 1e-6).sum()):4d} / {int((old[c] < 1e-6).sum()):4d}"
          f"   <1e-4 {int((new[c] < 1e-4).sum()):4d} / {int((old[c] < 1e-4).sum()):4d}")


def ilr(f, floor=1e-6):
    x = np.clip(f.values, floor, None)
    x = x / x.sum(1, keepdims=True)
    lx = np.log(x)
    clr = lx - lx.mean(1, keepdims=True)
    H = helmert(x.shape[1])  # (D-1) x D orthonormal rows
    return clr @ H.T


def auc(score, y):
    order = np.argsort(score)
    r = np.empty(len(score)); r[order] = np.arange(1, len(score) + 1)
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def logit_fit(X, y, l2=1e-2, iters=200):
    """Ridge-penalized logistic regression by Newton steps (tiny penalty for separation)."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b))
        W = p * (1 - p)
        g = X.T @ (y - p) - l2 * b
        Hm = X.T @ (X * W[:, None]) + l2 * np.eye(X.shape[1])
        step = np.linalg.solve(Hm, g)
        b += step
        if np.abs(step).max() < 1e-8:
            break
    return b


def overlap(label, keep, y, dset, frac):
    Z = ilr(frac.loc[keep])
    Z = (Z - Z.mean(0)) / Z.std(0)
    D = pd.get_dummies(dset[keep]).values.astype(float)
    X = np.column_stack([D, Z])
    b = logit_fit(X, y)
    ps = 1 / (1 + np.exp(-X @ b))
    a_in = auc(X @ b, y)
    # 5-fold cross-validated AUC, folds stratified by dataset x group
    folds = np.zeros(len(y), int)
    strata = pd.Series(list(zip(dset[keep], y)))
    for _, idx in strata.groupby(strata).groups.items():
        idx = np.array(list(idx)); rng.shuffle(idx)
        folds[idx] = np.arange(len(idx)) % 5
    cv = np.zeros(len(y))
    for k in range(5):
        tr = folds != k
        bk = logit_fit(X[tr], y[tr])
        cv[~tr] = X[~tr] @ bk
    a_cv = auc(cv, y)
    lo1, hi1 = np.quantile(ps[y == 1], [0.05, 0.95])
    lo0, hi0 = np.quantile(ps[y == 0], [0.05, 0.95])
    n1_out = int((ps[y == 1] > ps[y == 0].max()).sum())
    n0_out = int((ps[y == 0] < ps[y == 1].min()).sum())
    print(f"  {label}: n1={int(y.sum())} n0={int((1-y).sum())} AUC in-sample {a_in:.3f}, 5-fold CV {a_cv:.3f}; "
          f"propensity 5-95% group1 [{lo1:.3f},{hi1:.3f}] group0 [{lo0:.3f},{hi0:.3f}]; "
          f"group1 above all group0 {n1_out}; group0 below all group1 {n0_out}")


dset = dz["dataset"] if "dataset" in dz.columns else man["dataset"]
grp = (man["group_binary"] == "Disease").astype(int).values
print("\n2. composition overlap (dataset dummies + 15 standardized ILR, ridge 1e-2 logistic)")
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    overlap(f"disease vs control, {nm}", np.ones(len(ids), bool), grp, dset, fr)
st = man["fibrosis_stage"]
keep = st.isin([0, 4]).values & dset.isin(["GSE130970", "GSE135251", "GSE162694"]).values
y4 = (st[keep] == 4).astype(int).values
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    overlap(f"F4 vs F0 (3 cohorts), {nm}", keep, y4, dset, fr)

# 3. ILR floor sensitivity of the group-on-composition association (dataset-adjusted)
print("\n3. partial R2 of disease group on 15 ILR | dataset (unweighted), by floor")
D = pd.get_dummies(dset).values.astype(float)
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    for fl in [1e-6, 1e-4, 1e-3]:
        Z = ilr(fr, fl)
        g = grp.astype(float)
        r0 = g - D @ np.linalg.lstsq(D, g, rcond=None)[0]
        X = np.column_stack([D, Z])
        r1 = g - X @ np.linalg.lstsq(X, g, rcond=None)[0]
        print(f"  {nm} floor {fl:g}: partial R2 {1 - (r1 @ r1) / (r0 @ r0):.3f}")

# 4. what is the largest single-coordinate leverage driven by floor values
print("\n4. ILR coordinate range (floor 1e-6), rebuilt: max |z| per coordinate")
Z = ilr(new)
Zs = (Z - Z.mean(0)) / Z.std(0)
print("  ", np.round(np.abs(Zs).max(0), 1).tolist())
hat_X = np.column_stack([D, Z])
H = hat_X @ np.linalg.pinv(hat_X)
lev = np.diag(H)
print(f"  leverage: mean {lev.mean():.4f}, units with leverage > 3x mean: {(lev > 3 * lev.mean()).sum()}, max {lev.max():.3f}")
