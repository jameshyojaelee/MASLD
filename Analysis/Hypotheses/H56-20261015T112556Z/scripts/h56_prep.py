"""H56 stage 1 (inputs and design only; no expression values, no adjusted fit).

S1 input hashes, S3 fraction coverage, composition codings and ILR coordinates for every cell,
S4 dynamic-range recheck, diagnostic 2 (design-level SE inflation, S10), diagnostic 3
(propensity overlap), diagnostic 7 (specificity classes of the calls).
Functions for ILR, amalgamation and SE inflation are copied from
RESOURCE-REVIEW-20261015T020503Z/scratch/plan-H56/n_and_design.py so that S10 compares like with like.
"""
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import helmert

RUN = Path(os.environ["H56_RUN"])
OUT = Path(os.environ.get("H56_OUT", RUN / "out"))
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
PRESPEC_SHA = "cd786c6d17f60f2f88a7e3e11389820264cfa09e51e478d610516f38a847cf6d"
SEED = 20260923
OUT.mkdir(exist_ok=True)
LOG = open(OUT / "prep_diagnostics.txt", "w")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


def stop(code, msg):
    say(f"STOP {code}: {msg}")
    (OUT / "STOP_PREP.json").write_text(json.dumps({"stop": code, "message": msg}))
    sys.exit(3)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- prespec and S1
for p in (RUN / "PRESPEC_H56.md",
          ROOT / "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/PRESPEC_H56.md"):
    if sha256(p) != PRESPEC_SHA:
        stop("PRESPEC", f"prespec hash mismatch at {p}")
say("prespec sha256 OK", PRESPEC_SHA)
tab = pd.read_csv(RUN / "input_sha256.tsv", sep="\t", header=None, names=["sha", "bytes", "path"])
prespec_txt = (RUN / "PRESPEC_H56.md").read_text()
bad = []
for r in tab.itertuples():
    if r.sha not in prespec_txt:
        bad.append(("not_in_prespec", r.path))
        continue
    got = sha256(ROOT / r.path)
    if got != r.sha:
        bad.append((got, r.path))
say(f"S1 input hashes: {len(tab) - len(bad)} of {len(tab)} OK")
if len(tab) != 31 or bad:
    stop("S1", f"hash mismatch {bad}")

# ---------------------------------------------------------------- inputs
P = ROOT / "RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls"
NEW = ROOT / "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz"
OLD = str(ROOT / "Analysis/Deconvolution/results/{g}/{g}_bayesprism_proportions.tsv")
MAN = ROOT / "Analysis/Deconvolution/corrected_774_20261009T141331Z/bulk/sample_manifest.tsv"
DS = ["GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621"]
dz = pd.read_csv(P / "model_design.tsv", sep="\t").set_index("sample_id")
w = pd.read_csv(P / "lvqw_sample_weights.tsv", sep="\t").set_index("sample_id").sample_weight.loc[dz.index]
man = pd.read_csv(MAN, sep="\t").set_index("sample_id").loc[dz.index]
new = pd.read_csv(NEW, sep="\t").set_index("sample_id")
ct = [c for c in new.columns if c != "dataset"]
old_all = pd.concat([pd.read_csv(OLD.format(g=g), sep="\t", index_col=0) for g in DS])
have = dz.index.isin(new.index) & dz.index.isin(old_all.index)
say(f"S3 units with both fraction inputs: {int(have.sum())} of {len(dz)}; earlier rows {len(old_all)}")
if have.sum() < 774 or len(dz) != 774:
    stop("S3", "fewer than 774 units carry both fraction inputs")
new = new.loc[dz.index, ct]
old = old_all.loc[dz.index, ct]
ds = man.dataset
grp = dz.group_binaryDisease
assert man.analysis_unit_id.nunique() == 774

LIN7 = {"hep_epi": ["Hepatocytes", "Cholangiocytes"], "endo": ["Endothelial cells"], "fib": ["Fibroblasts"],
        "mac_mono": ["Macrophages", "Mono+mono derived cells"], "oth_my": ["cDC1s", "cDC2s", "pDCs", "Neutrophils", "Basophils"],
        "T_NK": ["T cells", "Circulating NK/NKT", "Resident NK"], "B_pl": ["B cells", "Plasma cells"]}
DR5 = {"hep_epi": LIN7["hep_epi"], "stromal": LIN7["fib"] + LIN7["endo"], "mac_mono": LIN7["mac_mono"],
       "oth_my": LIN7["oth_my"], "lymph": LIN7["T_NK"] + LIN7["B_pl"]}
DR4 = {"hep_epi": LIN7["hep_epi"], "mac_mono": LIN7["mac_mono"], "oth_my": LIN7["oth_my"],
       "npc_nonmy": LIN7["fib"] + LIN7["endo"] + LIN7["T_NK"] + LIN7["B_pl"]}
assert sorted(sum(LIN7.values(), [])) == sorted(ct)


def amalg(f, parts):
    a = np.column_stack([f[v].sum(1) for v in parts.values()])
    return a / a.sum(1, keepdims=True)


def ilr(x, floor):
    x = np.maximum(np.asarray(x, float), floor)
    D = x.shape[1]
    return np.log(x) @ (-helmert(D).T)


X0 = dz.drop(columns=["group_binaryDisease"]).to_numpy(float)
X0cols = [c for c in dz.columns if c != "group_binaryDisease"]
y = grp.to_numpy(float)
wt = w.to_numpy(float)


def infl(X0, Z, y, wt):
    def rv(X):
        sw = np.sqrt(wt)
        b, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
        r = (y - X @ b) * sw
        return r @ r
    return np.sqrt(rv(X0) / rv(np.hstack([X0, Z])))


main = ds.isin(["GSE130970", "GSE162694", "GSE213621"]).to_numpy()
sep = ds.isin(["GSE126848", "GSE135251"]).to_numpy()
allu = np.ones(774, bool)

# gate table (PRESPEC, "Power and minimum detectable change", SE inflation column)
GATE = {"A16_rebuilt": 1.157, "A16_earlier": 1.177, "C7_rebuilt": 1.087, "C7_earlier": 1.088,
        "D1_rebuilt": 1.148, "D1_earlier": 1.140, "DRmain_rebuilt": 1.116, "DRmain_earlier": 1.135,
        "DRsep_rebuilt": 1.001, "DRsep_earlier": 1.000}
FLAG = ["A16_rebuilt", "A16_earlier", "C7_rebuilt", "C7_earlier", "D1_rebuilt", "D1_earlier",
        "DRmain_rebuilt", "DRmain_earlier"]


def design_cols(keep, arm):
    """Comparator covariate columns (without group) for a cell's units."""
    nz = np.abs(X0[keep]).sum(0) > 0
    cols = [c for c, j in zip(X0cols, nz) if j]
    if arm == "DRmain":
        cols = [c for c in cols if c != "datasetGSE130970"]
    return cols


# ---------------------------------------------------------------- S4: dynamic-range recheck
say("\nS4 dynamic-range recheck (5th percentile of each DR part within each included cohort, both inputs)")
s4 = []
for arm, parts, coh in [("DRmain", DR5, ["GSE130970", "GSE162694", "GSE213621"]),
                        ("DRsep", DR4, ["GSE126848", "GSE135251"])]:
    for nm, fr in [("rebuilt", new), ("earlier", old)]:
        a = pd.DataFrame(amalg(fr, parts), index=fr.index, columns=list(parts))
        for c in coh:
            q = a.loc[(ds == c).to_numpy()].quantile(0.05)
            for part, v in q.items():
                s4.append([arm, nm, c, part, v])
s4 = pd.DataFrame(s4, columns=["arm", "input", "cohort", "part", "p05"])
s4.to_csv(OUT / "diag_S4_dynamic_range.tsv", sep="\t", index=False)
s4min = s4.groupby(["arm", "part"]).p05.min()
say(s4min.to_string())
drmain_ok = bool((s4[s4.arm == "DRmain"].p05 >= 1e-4).all())
drsep_ok = bool((s4[s4.arm == "DRsep"].p05 >= 1e-4).all())
say(f"S4 DR-main passes: {drmain_ok}; DR-sep passes: {drsep_ok}")

# ---------------------------------------------------------------- coordinates
cells = {}
for nm, fr in [("rebuilt", new), ("earlier", old)]:
    cells[f"A16_{nm}"] = ("A16", ilr(fr.to_numpy(), 1e-6), allu, fr.to_numpy(), 1e-6)
    cells[f"C7_{nm}"] = ("C7", ilr(amalg(fr, LIN7), 1e-6), allu, amalg(fr, LIN7), 1e-6)
    cells[f"D1_{nm}"] = ("D1", ilr(fr.to_numpy(), 1e-4), allu, fr.to_numpy(), 1e-4)
    cells[f"DRmain_{nm}"] = ("DRmain", ilr(amalg(fr, DR5), 1e-6), main, amalg(fr, DR5), 1e-6)
    cells[f"DRsep_{nm}"] = ("DRsep", ilr(amalg(fr, DR4), 1e-6), sep, amalg(fr, DR4), 1e-6)

# ---------------------------------------------------------------- diagnostic 2 and S10
say("\nDiagnostic 2: per cell design summary (stored sample weights)")
rows = []
for cell, (arm, Z, keep, parts_mat, floor) in cells.items():
    cols = design_cols(keep, arm)
    Xc = dz.loc[keep, cols].to_numpy(float)
    Zk = Z[keep]
    full = np.column_stack([Xc, Zk, y[keep]])
    rank = np.linalg.matrix_rank(full)
    k = infl(Xc, Zk, y[keep], wt[keep])
    pr2 = 1 - 1 / k**2
    nfloor = int((parts_mat[keep] < floor).sum())
    share = nfloor / parts_mat[keep].size
    sw = np.sqrt(wt[keep])
    Xw = full * sw[:, None]
    lev = np.einsum("ij,jk,ik->i", Xw, np.linalg.pinv(Xw.T @ Xw), Xw)
    rows.append(dict(cell=cell, arm=arm, n_units=int(keep.sum()), n_parts=parts_mat.shape[1], n_coords=Zk.shape[1],
                     design_cols=len(cols) + Zk.shape[1] + 1, design_rank=int(rank),
                     weighted_partial_R2_group_on_coords=round(pr2, 4), SE_inflation=round(k, 6),
                     SE_inflation_3dp=round(k, 3), gate_value=GATE[cell],
                     floor_replaced_n=nfloor, floor_replaced_share=round(share, 4), max_leverage=round(float(lev.max()), 4)))
d2 = pd.DataFrame(rows)
d2["S10_match"] = (d2.SE_inflation_3dp - d2.gate_value).abs() < 5e-4
d2["rank_full"] = d2.design_rank == d2.design_cols
say(d2.to_string(index=False))
d2.to_csv(OUT / "diag2_design.tsv", sep="\t", index=False)
if not d2.S10_match.all():
    stop("S10", "SE inflation differs from the gate table: " + d2.loc[~d2.S10_match, ["cell", "SE_inflation", "gate_value"]].to_string())
say("S10 passes: SE inflation matches the gate table to 3 decimals in all 10 cells")

# write coordinates (all 774 units; R subsets to each cell's units)
for cell, (arm, Z, keep, _, _) in cells.items():
    pd.DataFrame(Z, index=dz.index, columns=[f"ilr{i + 1}" for i in range(Z.shape[1])]).to_csv(
        OUT / f"coords_{cell}.tsv.gz", sep="\t")
pd.DataFrame({"cell": list(cells), "arm": [v[0] for v in cells.values()],
              "n_units": [int(v[2].sum()) for v in cells.values()],
              "comparator": ["U774" if v[0] in ("A16", "C7", "D1") else ("U537" if v[0] == "DRmain" else "U237")
                             for v in cells.values()],
              "covariate_cols": [";".join(design_cols(v[2], v[0])) for v in cells.values()],
              "flag_cell": [c in FLAG for c in cells],
              "estimable": [(r.rank_full and (drmain_ok if r.arm == "DRmain" else True)) for r in d2.itertuples()]}
             ).to_csv(OUT / "cells.tsv", sep="\t", index=False)
pd.DataFrame({"sample_id": dz.index, "dataset": ds.values, "group": y.astype(int), "sample_weight": wt,
              "in_DRmain": main, "in_DRsep": sep}).to_csv(OUT / "units.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- diagnostic 3 (propensity overlap)
rng = np.random.default_rng(SEED)


def auc(score, yy):
    order = np.argsort(score, kind="mergesort")
    r = np.empty(len(score))
    r[order] = np.arange(1, len(score) + 1)
    n1 = yy.sum()
    n0 = len(yy) - n1
    return (r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def logit_fit(X, yy, l2=1e-2, iters=200):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-X @ b))
        Wd = p * (1 - p)
        g = X.T @ (yy - p) - l2 * b
        Hm = X.T @ (X * Wd[:, None]) + l2 * np.eye(X.shape[1])
        step = np.linalg.solve(Hm, g)
        b += step
        if np.abs(step).max() < 1e-8:
            break
    return b


say("\nDiagnostic 3: 5-fold CV AUC of group on dataset + coordinates (ridge 1e-2 logistic), and overlap")
prop = pd.DataFrame(index=dz.index)
rows = []
for cell, (arm, Z, keep, _, _) in cells.items():
    yy = y[keep].astype(int)
    Zk = Z[keep]
    Zs = (Zk - Zk.mean(0)) / Zk.std(0)
    Dm = pd.get_dummies(ds[keep]).to_numpy(float)
    X = np.column_stack([Dm, Zs])
    b = logit_fit(X, yy)
    ps = 1 / (1 + np.exp(-X @ b))
    folds = np.zeros(len(yy), int)
    strata = pd.Series(list(zip(ds[keep], yy)))
    for _, idx in sorted(strata.groupby(strata).groups.items()):
        idx = np.array(list(idx))
        rng.shuffle(idx)
        folds[idx] = np.arange(len(idx)) % 5
    cv = np.zeros(len(yy))
    for k in range(5):
        tr = folds != k
        cv[~tr] = X[~tr] @ logit_fit(X[tr], yy[tr])
    above = (ps > ps[yy == 0].max()) & (yy == 1)
    below = (ps < ps[yy == 1].min()) & (yy == 0)
    rows.append(dict(cell=cell, n_units=int(keep.sum()), auc_in_sample=round(auc(X @ b, yy), 4), auc_cv5=round(auc(cv, yy), 4),
                     disease_above_all_controls=int(above.sum()), controls_below_all_disease=int(below.sum())))
    col = pd.Series(np.nan, index=dz.index)
    col[keep] = ps
    prop[f"ps_{cell}"] = col
    t = pd.Series(np.nan, index=dz.index)
    t[keep] = (~(above | below)).astype(float)
    prop[f"keep_{cell}"] = t
d3 = pd.DataFrame(rows)
say(d3.to_string(index=False))
d3.to_csv(OUT / "diag3_overlap.tsv", sep="\t", index=False)
prop.to_csv(OUT / "propensity.tsv.gz", sep="\t")

# ---------------------------------------------------------------- diagnostic 7 (specificity classes)
SC = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
PAIR = {d: ROOT / f"data/{d}/metadata/donor_pairing.csv" for d in ("GSE244832", "GSE185477", "GSE202379", "GSE136103")}
LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]
s2d = {}
for d, fp in PAIR.items():
    dp = pd.read_csv(fp, dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                s2d[s.strip()] = f"{d}_{r['donor_id']}"
prof, nd = {}, {}
for ln in LIN:
    m = pd.read_csv(SC / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    m = m.groupby(level=0).sum()
    cols = pd.Series(m.columns).map(s2d).fillna(pd.Series(m.columns)).values
    m = m.T.groupby(cols).sum().T
    cpm = m.div(m.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
    prof[ln] = np.log2(cpm + 1).mean(axis=1)
    nd[ln] = m.shape[1]
say("\nDiagnostic 7: donors per lineage profile", nd)
Pm = pd.DataFrame(prof).dropna(how="all")
Zp = Pm.sub(Pm.mean(1), axis=0).div(Pm.std(1).replace(0, np.nan), axis=0)
srt = np.sort(Pm.values, axis=1)
top_gap = pd.Series(srt[:, -1] - srt[:, -2], index=Pm.index)
rngv = pd.Series(srt[:, -1] - srt[:, 0], index=Pm.index)
top = Pm.idxmax(1)
deg = pd.read_csv(P / "deg_results.csv")
calls = deg[(deg.padj < 0.05) & (deg.logFC.abs() > 0.5)].copy()
sym_n = deg.symbol.value_counts()
calls["in_sc"] = calls.symbol.isin(Pm.index) & calls.symbol.map(sym_n).eq(1)
c = calls[calls.in_sc].set_index("symbol")
zmax = Zp.loc[c.index].max(1)
spec = (zmax >= 1.5) & (top_gap.loc[c.index] >= 1)
cls = pd.Series("other", index=c.index)
cls[spec] = "specific_" + top.loc[c.index][spec]
uni = (rngv.loc[c.index] < 1.5) & (Pm.loc[c.index].min(1) >= 1) & ~spec
cls[uni] = "uniform"
CLS = {"specific_Fibroblasts": "stromal", "specific_Endothelial_cells": "stromal", "specific_Hepatocytes": "hepatocyte",
       "specific_Macrophages": "immune", "specific_T_cells": "immune", "specific_Cholangiocytes": "cholangiocyte",
       "uniform": "uniform", "other": "other"}
c["sc_class_detail"] = cls
c["sc_class"] = cls.map(CLS)
out = calls[["gene", "symbol", "logFC"]].merge(c.reset_index()[["gene", "sc_class", "sc_class_detail"]], on="gene", how="left")
out["sc_class"] = out.sc_class.fillna("not_joined")
MEMB = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/program_membership_v2.tsv"
mb = pd.read_csv(MEMB, sep="\t")
ecm = set(mb[mb.program_uid == "hotspot_hepatocytes_f05c535ae5bbc0b9"].source_gene.unique())
out["stromal_ecm"] = out.symbol.isin(ecm)
out.to_csv(OUT / "call_classes.tsv", sep="\t", index=False)
cnt = out.sc_class.value_counts()
say(f"calls {len(out)}; joined {int(calls.in_sc.sum())}")
say(cnt.to_string())
say("Stromal ECM calls", int(out.stromal_ecm.sum()), sorted(out.symbol[out.stromal_ecm]))
exp = {"stromal": 159, "hepatocyte": 183, "immune": 164, "cholangiocyte": 167, "uniform": 52}
ok7 = int(calls.in_sc.sum()) == 1597 and all(cnt.get(k, 0) == v for k, v in exp.items()) and int(out.stromal_ecm.sum()) == 12
say("Diagnostic 7 class sizes match the prespec:", ok7)
if not ok7:
    say("WARNING: diagnostic 7 class sizes differ from the prespec (reported; not a stop rule)")
say("PREP_DONE")
