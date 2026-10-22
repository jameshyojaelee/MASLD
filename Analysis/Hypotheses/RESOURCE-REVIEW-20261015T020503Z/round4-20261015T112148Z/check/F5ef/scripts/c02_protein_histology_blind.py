"""Blind, independent Python re-derivation of the Fig 5E protein-histology column.

Raw inputs only: PXD051911 liver_protein_quant.txt and meta_data.txt, the v2
program registry and the v2 membership. Nothing is read from the proposer's
workstream. Random draws use numpy, so permutation P and bootstrap CIs agree
with the proposer only to Monte Carlo error.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

P = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = P / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F5ef/out"
RAW = P / "data/PXD051911/liver_protein_quant.txt"
META = P / "data/PXD051911/meta_data.txt"
REG = P / "figures/main/fig4_singlecell_programs/source_tables/source_diagnosis_correction_20261014/program_registry_v2.tsv"
MEM = P / ("Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/"
           "hotspot/program_membership_v2.tsv")
N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 9999
N_BOOT = int(sys.argv[2]) if len(sys.argv) > 2 else 4999
PERM_BASE, BOOT_BASE = 20260923, 30260923


def limma_quantile(A):
    """limma::normalizeQuantiles(ties=TRUE) with NA handling."""
    n, k = A.shape
    i = np.arange(n) / (n - 1)
    S = np.empty((n, k))
    for j in range(k):
        x = np.sort(A[~np.isnan(A[:, j]), j])
        m = len(x)
        S[:, j] = x if m == n else np.interp(i, np.arange(m) / (m - 1), x)
    mean = S.mean(axis=1)
    out = np.full_like(A, np.nan)
    for j in range(k):
        ok = ~np.isnan(A[:, j])
        m = ok.sum()
        r = rankdata(A[ok, j], method="average")
        out[ok, j] = np.interp((r - 1) / (m - 1), i, mean)
    return out


# ---- protein matrix ------------------------------------------------------
raw = pd.read_csv(RAW, sep="\t")
meta = pd.read_csv(META, sep="\t")
samples_all = [c for c in raw.columns if c not in ("ProteinAccessions", "Genes", "ProteinDescriptions")]
meta = meta[meta.liver_proteomics_filename.isin(samples_all)].drop_duplicates("liver_proteomics_filename")
meta = meta.set_index("liver_proteomics_filename").loc[samples_all].reset_index()
assert len(meta) == 58
raw = raw[raw.Genes.notna() & (raw.Genes != "") & ~raw.Genes.astype(str).str.contains(";", regex=False)]
X = raw[samples_all].to_numpy(dtype=float)
X[~np.isfinite(X) | (X <= 0)] = np.nan
X = np.log2(X)
keep = np.isnan(X).mean(axis=1) <= 0.5
X = X[keep]
genes = raw.Genes.to_numpy()[keep]
X = limma_quantile(X)
df = pd.DataFrame(X, columns=samples_all)
df["gene"] = genes
G = df.groupby("gene", sort=True)[samples_all].median()  # NaN-skipping median
mu = G.mean(axis=1, skipna=True)
sd = G.std(axis=1, skipna=True, ddof=1)
Z = G.sub(mu, axis=0).div(sd, axis=0)
Z = Z.replace([np.inf, -np.inf], np.nan)
print(f"single-gene groups after filter {keep.sum()}; symbols {len(Z)}", flush=True)

# ---- participants ----------------------------------------------------------
meta["masld"] = meta.saf_diagnosis != "No_MASLD"
meta["batch2020"] = (~meta.liver_proteomics_filename.str.startswith("2019")).astype(float)
meta["male"] = (meta.gender == "Male").astype(float)
meta["age"] = meta.alder.astype(float)
meta["bmi_n"] = meta.bmi.astype(float)
meta["Fibrosis"] = meta.kleiner_fibrosis_grade.str.replace("F", "", regex=False).astype(float)
meta["Steatosis"] = meta.steatosis_score.astype(float)
meta["Ballooning"] = meta.hepatocellular_ballooning_score.astype(float)
meta["Inflammation"] = meta.lobular_inflammation_score.astype(float)
COMP = ["Steatosis", "Ballooning", "Inflammation", "Fibrosis"]
mas = meta[meta.masld & meta[COMP].notna().all(axis=1)].reset_index(drop=True)
print("MASLD complete histology n =", len(mas), "unique patients", mas.patient_name.nunique(),
      "batch counts", mas.batch2020.value_counts().to_dict(), flush=True)


def burden(C):
    z = (C - C.mean(0)) / C.std(0, ddof=1)
    if not np.all(np.isfinite(z)):
        return None, None, None
    U, s, Vt = np.linalg.svd(z, full_matrices=False)
    load = Vt[0]
    score = z @ load
    if load.sum() < 0:
        load, score = -load, -score
    return score, load, s[0] ** 2 / (s ** 2).sum()


B, LOAD, VE = burden(mas[COMP].to_numpy())
print("burden loadings", dict(zip(COMP, np.round(LOAD, 6))), "PC1 var", round(VE, 5), flush=True)


def design(d, with_batch=True):
    cols = [np.ones(len(d))]
    if with_batch and d.batch2020.nunique() > 1:
        cols.append(d.batch2020.to_numpy())
    cols += [d.age.to_numpy(), d.bmi_n.to_numpy()]
    if d.male.nunique() > 1:
        cols.append(d.male.to_numpy())
    return np.column_stack(cols)


def resid(Xd, y):
    beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    return y - Xd @ beta


def pspear(score, bur, Xd):
    a = resid(Xd, rankdata(score))
    b = resid(Xd, rankdata(bur))
    return np.corrcoef(a, b)[0, 1]


# ---- program scores --------------------------------------------------------
reg = pd.read_csv(REG, sep="\t", usecols=["cell_type", "module", "program_uid", "module_name"])
mem = pd.read_csv(MEM, sep="\t")
mem = mem[mem.mapped_symbol_status != "unmapped_or_ambiguous"]
assert len(mem) == 6633
l1tot = pd.read_csv(MEM, sep="\t").groupby("program_uid").original_l1_weight.sum()
print("L1 weight totals over all source genes: min %.6f max %.6f" % (l1tot.min(), l1tot.max()), flush=True)
mcols = mas.liver_proteomics_filename.tolist()
Zm = Z[mcols]
Xfull = design(mas)
batches = mas.batch2020.to_numpy()
rows = []
scores_out = []
testable_index = 0
for _, rg in reg.iterrows():
    uid = rg.program_uid
    m = mem[mem.program_uid == uid].groupby("mapped_symbol", sort=False).original_l1_weight.sum()
    meas = m[m.index.isin(Z.index)]
    nm, ret = len(meas), meas.sum()
    row = dict(program_uid=uid, cell_type=rg.cell_type, module=rg.module, module_name=rg.module_name,
               n_members=len(m), n_measured=nm, retained_l1=ret, testable=(nm >= 8 and ret >= 0.20))
    if not row["testable"]:
        rows.append(row)
        continue
    testable_index += 1
    w = (meas / meas.sum())
    zz = Zm.loc[w.index].to_numpy()
    obs = np.isfinite(zz)
    ww = w.to_numpy()[:, None]
    sc = np.where(obs, zz, 0).T @ w.to_numpy() / (obs * ww).sum(0)
    eq = np.nanmean(zz, axis=0)
    top = w.idxmax()
    w2 = meas.drop(top)
    w2 = w2 / w2.sum()
    zz2 = Zm.loc[w2.index].to_numpy()
    obs2 = np.isfinite(zz2)
    lt = np.where(obs2, zz2, 0).T @ w2.to_numpy() / (obs2 * w2.to_numpy()[:, None]).sum(0)
    scores_out.append(pd.DataFrame(dict(program_uid=uid, sample_id=mcols, score=sc, equal_score=eq, leave_top_score=lt)))
    assert np.all(np.isfinite(sc))
    rho = pspear(sc, B, Xfull)
    rho_eq = pspear(eq, B, Xfull)
    rho_lt = pspear(lt, B, Xfull)
    rb = []
    for bval in (0.0, 1.0):
        sel = batches == bval
        d = mas[sel]
        rb.append(pspear(sc[sel], B[sel], design(d, with_batch=False)))
    # Freedman-Lane, batch-stratified, two-sided on |rho|
    rng = np.random.default_rng(PERM_BASE + testable_index)
    xr = resid(Xfull, rankdata(sc))
    yr = resid(Xfull, rankdata(B))
    H = Xfull @ np.linalg.pinv(Xfull)
    R = np.eye(len(yr)) - H
    idx = np.tile(np.arange(len(yr)), (N_PERM, 1))
    for bval in (0.0, 1.0):
        ii = np.where(batches == bval)[0]
        idx[:, ii] = ii[rng.permuted(np.tile(np.arange(len(ii)), (N_PERM, 1)), axis=1)]
    ystar = yr[idx] @ R.T  # residual of fitted + permuted residual
    xc = xr - xr.mean()
    yc = ystar - ystar.mean(1, keepdims=True)
    pstat = (yc @ xc) / np.sqrt((yc ** 2).sum(1) * (xc ** 2).sum())
    obs_r = np.corrcoef(xr, yr)[0, 1]
    k = int((np.abs(pstat) >= abs(obs_r) - 1e-12).sum())
    p = (k + 1) / (N_PERM + 1)
    # bootstrap within batch, re-deriving the burden axis
    rngb = np.random.default_rng(BOOT_BASE + testable_index)
    C = mas[COMP].to_numpy()
    bvals = []
    i0 = np.where(batches == 0.0)[0]
    i1 = np.where(batches == 1.0)[0]
    for _b in range(N_BOOT):
        bi = np.concatenate([rngb.choice(i0, len(i0), replace=True), rngb.choice(i1, len(i1), replace=True)])
        bb, ld, _ = burden(C[bi])
        if bb is None or np.any(ld <= 0):
            bvals.append(np.nan)
            continue
        bvals.append(pspear(sc[bi], bb, design(mas.iloc[bi])))
    bvals = np.array(bvals)
    ci = np.nanquantile(bvals, [0.025, 0.975])
    row.update(n=len(sc), rho=rho, perm_k=k, pvalue=p, ci_low=ci[0], ci_high=ci[1],
               n_boot_finite=int(np.isfinite(bvals).sum()), equal_rho=rho_eq, leave_top_rho=rho_lt,
               rho_2019=rb[0], rho_2020=rb[1], top_gene=top)
    rows.append(row)
    print(f"{testable_index:3d} {rg.cell_type}:{rg.module} rho={rho:.4f} p={p:.4f}", flush=True)

res = pd.DataFrame(rows)
t = res.testable
pv = res.loc[t, "pvalue"].to_numpy()
order = np.argsort(pv)
m = len(pv)
q = np.empty(m)
q[order] = np.minimum.accumulate((pv[order] * m / np.arange(1, m + 1))[::-1])[::-1]
q = np.minimum(q, 1)
res.loc[t, "padj"] = q
res["sens_agree"] = t & (np.sign(res.rho) == np.sign(res.equal_rho)) & (np.sign(res.rho) == np.sign(res.leave_top_rho))
res["batch_agree"] = t & (np.sign(res.rho) == np.sign(res.rho_2019)) & (np.sign(res.rho) == np.sign(res.rho_2020))
res["robust"] = t & (res.padj < 0.05) & res.sens_agree & res.batch_agree
res.to_csv(OUT / f"c02_protein_histology_blind_perm{N_PERM}_boot{N_BOOT}.tsv", sep="\t", index=False)
pd.concat(scores_out).to_csv(OUT / "c02_protein_scores_blind.tsv", sep="\t", index=False)
print("testable", int(t.sum()), "q<0.05", int((res.padj < 0.05).sum()), "robust", int(res.robust.sum()),
      "robust negative", int((res.robust & (res.rho < 0)).sum()), "robust positive", int((res.robust & (res.rho > 0)).sum()))
