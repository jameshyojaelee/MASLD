"""H90 job 6300 (b): discovery bridge, gate G4 and diagnostic D4.

Computes the frozen-weight projection of Ductular injury (hepatocyte module 20)
in the discovery hepatocyte atlas on the 62-donor Fig 4C universe, standardizing
genes within each discovery dataset, and fits score ~ stage + dataset.
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import h90_common as C

C.check_inputs()
out = C.OUT


def ols(y, X):
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    e = y - X @ b
    return b, e, XtX_inv


def design_stage_dataset(df):
    ds = pd.get_dummies(df["dataset"], drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(df)), df["stage"].to_numpy(float), ds.to_numpy(float)])
    return X


# ---- (a) reproduction of the Fig 4C 62-donor universe --------------------------
sc = pd.read_csv(C.FIG4C_SCORES, sep="\t", usecols=["cell_type", "module", "donor", "score", "dataset"])
sc = sc[(sc.cell_type == C.PRIMARY[0]) & (sc.module == C.PRIMARY[1])]
n_scores = sc.donor.nunique()
meta = pd.read_csv(C.DISC_META, sep="\t")
d = sc.merge(meta[["sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"]],
             left_on=["donor", "dataset"], right_on=["sample", "dataset"], how="inner")
d = d[(~d.exclude_stage_analysis.astype(str).str.upper().eq("TRUE"))
      & d.disease_stage_coarse.isin(["Healthy", "Steatosis", "Steatohepatitis"])].copy()
d["stage"] = d.disease_stage_coarse.map({"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2})
d = d.sort_values(["dataset", "donor"]).reset_index(drop=True)
stage_counts = d.stage.value_counts().sort_index().to_dict()
X = design_stage_dataset(d)
b, e, _ = ols(d.score.to_numpy(float), X)
beta_fig4c = float(b[1])
g4a = (n_scores == 64 and len(d) == 62 and stage_counts == {0: 15, 1: 11, 2: 36}
       and abs(beta_fig4c - 0.402086) < 1e-6)
print(f"G4a: scores donors {n_scores}, universe {len(d)}, stages {stage_counts}, Fig 4C refit beta {beta_fig4c:.9f}")

# ---- gene sets ---------------------------------------------------------------
feat = pd.read_csv(C.FEATURES_FL514, sep="\t", header=None, names=["id", "symbol", "type"])
feat = feat[feat.type == "Gene Expression"]
f = h5py.File(C.DISC_H5AD, "r")
var = pd.Index(f["var/_index"][:].astype(str))
lib_genes = C.library_gene_set(pd.Index(feat.symbol), var)
(out / "library_genes.txt").write_text("\n".join(lib_genes) + "\n")
print(f"library gene set: {len(lib_genes)} genes (GSE289173 2020-A stripped x discovery hepatocyte var)")

mem = C.program_table()
w = C.program_weights(mem, *C.PRIMARY)
w = w[w.mapped_symbol.isin(set(feat.symbol)) & w.mapped_symbol.isin(set(var))].reset_index(drop=True)
print(f"Ductular injury genes in 2020-A features and discovery var: {len(w)}, weight {w.original_l1_weight.sum():.4f}")

# ---- read discovery counts for universe cells, chunked -----------------------
def cat(name):
    g = f[f"obs/{name}"]
    return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))

obs = pd.DataFrame({"sample": cat("sample"), "dataset": cat("dataset")})
s2d = C.build_sample_to_donor()
obs["donor"] = obs["sample"].astype(str).map(s2d).fillna(obs["sample"].astype(str))
uni = set(zip(d.donor, d.dataset))
keep = np.array([(a, b2) in uni for a, b2 in zip(obs.donor, obs.dataset.astype(str))])
print(f"discovery hepatocyte cells {len(obs):,}; in 62-donor universe {keep.sum():,}")
missing = uni - set(zip(obs.donor[keep], obs.dataset.astype(str)[keep]))
if missing:
    sys.exit(f"STOP G4: universe donors without cells: {sorted(missing)}")

lib_mask = np.asarray(var.isin(lib_genes))
gcol = np.full(len(var), -1, dtype=np.int64)
for j, g in enumerate(w.mapped_symbol):
    gcol[var.get_loc(g)] = j
G = len(w)
indptr = f["X/indptr"][:]
rows_keep = np.flatnonzero(keep)
lib = np.zeros(len(obs))
cnt = np.zeros((len(obs), G), dtype=np.float64)
step = 20000
for a in range(0, len(obs), step):
    bnd = min(a + step, len(obs))
    if not keep[a:bnd].any():
        continue
    s, t = indptr[a], indptr[bnd]
    data = f["X/data"][s:t].astype(np.float64)
    ind = f["X/indices"][s:t]
    rid = np.repeat(np.arange(a, bnd), np.diff(indptr[a:bnd + 1]))
    m = lib_mask[ind]
    lib[a:bnd] = np.bincount(rid[m] - a, weights=data[m], minlength=bnd - a)
    gc = gcol[ind]
    m2 = gc >= 0
    np.add.at(cnt, (rid[m2], gc[m2]), data[m2])
f.close()
lib = lib[rows_keep]
cnt = cnt[rows_keep]
cobs = obs.iloc[rows_keep].reset_index(drop=True)
x = np.log1p(np.divide(cnt * 1e4, lib[:, None], out=np.zeros_like(cnt), where=lib[:, None] > 0))
det = (cnt > 0).mean(axis=0)
scored = det >= 0.01
ws = w[scored].copy()
ws["weight"] = ws.original_l1_weight / ws.original_l1_weight.sum()
print(f"bridge scored genes (detection >= 1% in universe cells): {scored.sum()} of {G}; retained frozen weight {ws.original_l1_weight.sum():.4f}")
z = np.zeros((len(cobs), scored.sum()))
for ds in cobs.dataset.astype(str).unique():
    i = np.flatnonzero(cobs.dataset.astype(str).to_numpy() == ds)
    xs = x[np.ix_(i, np.flatnonzero(scored))]
    mu = xs.mean(axis=0)
    sd = xs.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd <= 0)] = 1.0
    z[i] = (xs - mu) / sd
cobs["proj"] = z @ ws.weight.to_numpy()
dproj = cobs.groupby(["donor", "dataset"], observed=True).agg(proj=("proj", "mean"), n_cells=("proj", "size")).reset_index()
dproj["dataset"] = dproj.dataset.astype(str)
br = d.merge(dproj, on=["donor", "dataset"], how="left")
if br.proj.isna().any() or len(br) != 62:
    sys.exit("STOP G4: bridge donor set differs from the 62 Fig 4C donors")
b2, e2, _ = ols(br.proj.to_numpy(float), design_stage_dataset(br))
slope_proj = float(b2[1])
br["proj_c"] = br.proj - br.groupby("dataset").proj.transform("mean")
br["fig4c_c"] = br.score - br.groupby("dataset").score.transform("mean")
rho = float(stats.spearmanr(br.proj_c, br.fig4c_c).statistic)
br.to_csv(out / "bridge_donor_scores.tsv", sep="\t", index=False)
ws.assign(detection=det[scored]).to_csv(out / "bridge_scored_genes.tsv", sep="\t", index=False)

# ---- D4: standardized stage means in GSE202379 and GSE244832 (gate convention) --
def stdmeans(h, col):
    grp = h.stage.to_numpy()
    yv = h[col].to_numpy(float)
    mu = {k: yv[grp == k].mean() for k in (0, 1, 2)}
    resid = yv - np.array([mu[k] for k in grp])
    sdr = np.sqrt((resid ** 2).sum() / (len(yv) - 3))
    return {k: (mu[k] - mu[0]) / sdr for k in (1, 2)}, len(yv)

d4 = []
for col, lab in (("proj", "projection"), ("score", "fig4c_hotspot")):
    r = {ds: stdmeans(br[br.dataset == ds], col) for ds in ("GSE202379", "GSE244832")}
    wsum = sum(n for _, n in r.values())
    m1 = sum(v[1] * n for v, n in r.values()) / wsum
    m2 = sum(v[2] * n for v, n in r.values()) / wsum
    ns = np.array([23, 6, 19])
    xx = np.repeat([0, 1, 2], ns)
    mux = np.repeat([0, m1, m2], ns)
    bimp = np.cov(xx, mux)[0, 1] / np.var(xx, ddof=1)
    res = np.var(mux - bimp * xx, ddof=1) * (len(xx) - 1) / (len(xx) - 2)
    se = np.sqrt((1 + res) / ((xx - xx.mean()) ** 2).sum())
    df = len(xx) - 2
    tc = stats.t.ppf(0.975, df)
    ncp = bimp / se
    powr = 1 - stats.nct.cdf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)
    for ds, (v, n) in r.items():
        d4.append(dict(score=lab, dataset=ds, n=n, m1=v[1], m2=v[2]))
    d4.append(dict(score=lab, dataset="n-weighted", n=wsum, m1=m1, m2=m2,
                   implied_slope_planned_23_6_19=bimp, expected_power_alpha05_planned=powr))
d4 = pd.DataFrame(d4)
d4.to_csv(out / "D4_bridge.tsv", sep="\t", index=False)
print(d4.to_string(index=False))

g4b = slope_proj > 0
rows = [
    dict(check="G4a_scores_donors_hep20", value=n_scores, expected=64, pass_=n_scores == 64),
    dict(check="G4a_universe_donors", value=len(d), expected=62, pass_=len(d) == 62),
    dict(check="G4a_stage_counts", value=f"{stage_counts.get(0)}/{stage_counts.get(1)}/{stage_counts.get(2)}", expected="15/11/36",
         pass_=stage_counts == {0: 15, 1: 11, 2: 36}),
    dict(check="G4a_fig4c_refit_beta", value=f"{beta_fig4c:.12f}", expected="0.402086 +/- 1e-6", pass_=abs(beta_fig4c - 0.402086) < 1e-6),
    dict(check="G4a_bridge_donor_set_equal", value=int(br.proj.notna().sum()), expected=62, pass_=True),
    dict(check="G4b_projection_slope_positive", value=f"{slope_proj:.6f}", expected="> 0", pass_=g4b),
    dict(check="G4c_spearman_within_dataset (qualifier)", value=f"{rho:.4f}", expected=">= 0.50 else weak-agreement wording", pass_=True),
    dict(check="bridge_scored_genes", value=int(scored.sum()), expected="of 21", pass_=True),
    dict(check="bridge_retained_weight", value=f"{ws.original_l1_weight.sum():.4f}", expected="", pass_=True),
    dict(check="bridge_cells", value=len(cobs), expected="", pass_=True),
]
gt = pd.DataFrame(rows).rename(columns={"pass_": "pass"})
gt.to_csv(out / "gate_G4_bridge.tsv", sep="\t", index=False)
print(gt.to_string(index=False))
if not g4a:
    (out / "STOP_G4.txt").write_text("G4a failed: Fig 4C universe or refit beta not reproduced\n")
    sys.exit("STOP G4a")
if not g4b:
    (out / "STOP_G4.txt").write_text("untestable as a replication: the frozen-weight projection does not carry the Fig 4C direction in discovery\n")
    sys.exit("STOP G4b")
print("G4 pass")
