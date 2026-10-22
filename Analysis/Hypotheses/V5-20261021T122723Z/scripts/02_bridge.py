"""V5 job 2: discovery bridge, gate G4. Follows PRESPEC_V5.md (G4).

H90 02_bridge_discovery.py logic, with two changes named in the prespec: the
projection uses the Ductular injury genes that pass G3 in V5 job 1, and G4b
also requires a positive disease (stage 1 or 2 vs 0) beta. Paths and copied
helpers are loaded from the pinned H90 module, not edited.

Usage: python -I 02_bridge.py
"""
import importlib.util
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
import v5_common as C

H90_COMMON = C.ROOT / "Analysis/Hypotheses/H90-20261015T120337Z/scripts/h90_common.py"
spec = importlib.util.spec_from_file_location("h90_common", H90_COMMON)
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)

out = C.OUT
if (out / "gate_G4_bridge.tsv").exists():
    sys.exit("refusing to overwrite gate_G4_bridge.tsv")
print(f"prespec sha256 {C.check_prespec()}", flush=True)
gate = pd.read_csv(out / "gate_G0_G3.tsv", sep="\t", index_col=0)["value"]
if not all(str(gate[k]) == "True" for k in ("G0_pass", "G1_pass", "G2_pass", "G3_primary_pass")):
    sys.exit("STOP: job 1 gates did not pass")


def ols(y, X):
    XtX_inv = np.linalg.inv(X.T @ X)
    b = XtX_inv @ X.T @ y
    return b, y - X @ b


def design(df, term):
    ds = pd.get_dummies(df["dataset"], drop_first=True, dtype=float)
    return np.column_stack([np.ones(len(df)), df[term].to_numpy(float), ds.to_numpy(float)])


# ---- (a) Fig 4C universe and refit ------------------------------------------------
sc = pd.read_csv(H.FIG4C_SCORES, sep="\t", usecols=["cell_type", "module", "donor", "score", "dataset"])
sc = sc[(sc.cell_type == "hepatocytes") & (sc.module == 20)]
n_scores = sc.donor.nunique()
meta = pd.read_csv(H.DISC_META, sep="\t")
d = sc.merge(meta[["sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis"]],
             left_on=["donor", "dataset"], right_on=["sample", "dataset"], how="inner")
d = d[(~d.exclude_stage_analysis.astype(str).str.upper().eq("TRUE"))
      & d.disease_stage_coarse.isin(["Healthy", "Steatosis", "Steatohepatitis"])].copy()
d["stage"] = d.disease_stage_coarse.map({"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2})
d["disease"] = (d.stage > 0).astype(int)
d = d.sort_values(["dataset", "donor"]).reset_index(drop=True)
stage_counts = d.stage.value_counts().sort_index().to_dict()
b, _ = ols(d.score.to_numpy(float), design(d, "stage"))
beta_fig4c = float(b[1])
g4a = n_scores == 64 and len(d) == 62 and stage_counts == {0: 15, 1: 11, 2: 36} and abs(beta_fig4c - 0.402086) < 1e-6
print(f"G4a: {n_scores} score donors, universe {len(d)}, stages {stage_counts}, refit beta {beta_fig4c:.9f}", flush=True)

# ---- gene sets: V5 G3-scored Ductular genes; library genes from the Flex axis -----
det = pd.read_csv(out / "D3_confirmatory_gene_detection.tsv", sep="\t")
v5_genes = set(det[(det.program_uid == C.PRIMARY) & det.scored].gene)
feat = pd.read_csv(sorted(C.MAT.glob("*_features.tsv.gz"))[0], sep="\t", header=None, names=["id", "symbol", "type"])
f = h5py.File(H.DISC_H5AD, "r")
var = pd.Index(f["var/_index"][:].astype(str))
lib_genes = H.library_gene_set(pd.Index(feat.symbol), var)
(out / "bridge_library_genes.txt").write_text("\n".join(lib_genes) + "\n")
w = C.weights(C.membership(), C.PRIMARY)
w = w[w.mapped_symbol.isin(v5_genes) & w.mapped_symbol.isin(set(var))].reset_index(drop=True)
print(f"bridge genes (V5 scored and in discovery var): {len(w)} of {len(v5_genes)}; weight {w.original_l1_weight.sum():.4f}", flush=True)


def cat(name):
    g = f[f"obs/{name}"]
    return pd.Categorical.from_codes(g["codes"][:], g["categories"][:].astype(str))


obs = pd.DataFrame({"sample": cat("sample"), "dataset": cat("dataset")})
s2d = H.build_sample_to_donor()
obs["donor"] = obs["sample"].astype(str).map(s2d).fillna(obs["sample"].astype(str))
uni = set(zip(d.donor, d.dataset))
keep = np.array([(a, b2) in uni for a, b2 in zip(obs.donor, obs.dataset.astype(str))])
missing = uni - set(zip(obs.donor[keep], obs.dataset.astype(str)[keep]))
if missing:
    sys.exit(f"STOP G4: universe donors without cells: {sorted(missing)}")

lib_mask = np.asarray(var.isin(lib_genes))
gcol = np.full(len(var), -1, dtype=np.int64)
for j, g in enumerate(w.mapped_symbol):
    gcol[var.get_loc(g)] = j
indptr = f["X/indptr"][:]
lib = np.zeros(len(obs))
cnt = np.zeros((len(obs), len(w)))
step = 20000
for a in range(0, len(obs), step):
    e = min(a + step, len(obs))
    if not keep[a:e].any():
        continue
    s, t = indptr[a], indptr[e]
    data = f["X/data"][s:t].astype(np.float64)
    ind = f["X/indices"][s:t]
    rid = np.repeat(np.arange(a, e), np.diff(indptr[a:e + 1]))
    m = lib_mask[ind]
    lib[a:e] = np.bincount(rid[m] - a, weights=data[m], minlength=e - a)
    gc = gcol[ind]
    m2 = gc >= 0
    np.add.at(cnt, (rid[m2], gc[m2]), data[m2])
f.close()
rows = np.flatnonzero(keep)
lib, cnt = lib[rows], cnt[rows]
cobs = obs.iloc[rows].reset_index(drop=True)
x = np.log1p(np.divide(cnt * 1e4, lib[:, None], out=np.zeros_like(cnt), where=lib[:, None] > 0))
dd = (cnt > 0).mean(axis=0)
scored = dd >= C.DET_MIN
ws = w[scored].copy()
ws["weight"] = ws.original_l1_weight / ws.original_l1_weight.sum()
z = np.zeros((len(cobs), scored.sum()))
for ds in cobs.dataset.astype(str).unique():
    i = np.flatnonzero(cobs.dataset.astype(str).to_numpy() == ds)
    xs = x[np.ix_(i, np.flatnonzero(scored))]
    sd = xs.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd <= 0)] = 1.0
    z[i] = (xs - xs.mean(axis=0)) / sd
cobs["proj"] = z @ ws.weight.to_numpy()
dproj = cobs.groupby(["donor", "dataset"], observed=True).agg(proj=("proj", "mean"), n_cells=("proj", "size")).reset_index()
dproj["dataset"] = dproj.dataset.astype(str)
br = d.merge(dproj, on=["donor", "dataset"], how="left")
if br.proj.isna().any() or len(br) != 62:
    sys.exit("STOP G4: bridge donor set differs from the 62 Fig 4C donors")
slope_stage = float(ols(br.proj.to_numpy(float), design(br, "stage"))[0][1])
beta_disease = float(ols(br.proj.to_numpy(float), design(br, "disease"))[0][1])
br["proj_c"] = br.proj - br.groupby("dataset").proj.transform("mean")
br["fig4c_c"] = br.score - br.groupby("dataset").score.transform("mean")
rho = float(stats.spearmanr(br.proj_c, br.fig4c_c).statistic)
br.to_csv(out / "bridge_donor_scores.tsv", sep="\t", index=False)
ws.assign(detection=dd[scored]).to_csv(out / "bridge_scored_genes.tsv", sep="\t", index=False)

g4b = slope_stage > 0 and beta_disease > 0
gt = pd.DataFrame([
    dict(check="G4a_scores_donors_hep20", value=n_scores, expected=64, passed=n_scores == 64),
    dict(check="G4a_universe_donors", value=len(d), expected=62, passed=len(d) == 62),
    dict(check="G4a_stage_counts", value=f"{stage_counts.get(0)}/{stage_counts.get(1)}/{stage_counts.get(2)}",
         expected="15/11/36", passed=stage_counts == {0: 15, 1: 11, 2: 36}),
    dict(check="G4a_fig4c_refit_beta", value=f"{beta_fig4c:.12f}", expected="0.402086 +/- 1e-6",
         passed=abs(beta_fig4c - 0.402086) < 1e-6),
    dict(check="G4b_projection_stage_slope", value=f"{slope_stage:.6f}", expected="> 0", passed=slope_stage > 0),
    dict(check="G4b_projection_disease_beta", value=f"{beta_disease:.6f}", expected="> 0", passed=beta_disease > 0),
    dict(check="G4c_spearman_within_dataset (qualifier)", value=f"{rho:.4f}", expected=">= 0.50 else weak-agreement wording", passed=True),
    dict(check="bridge_genes_scored", value=int(scored.sum()), expected=f"of {len(w)}", passed=True),
    dict(check="bridge_retained_weight", value=f"{ws.original_l1_weight.sum():.4f}", expected="", passed=True),
    dict(check="bridge_cells", value=len(cobs), expected="", passed=True),
])
gt.to_csv(out / "gate_G4_bridge.tsv", sep="\t", index=False)
print(gt.to_string(index=False), flush=True)
if not g4a:
    (out / "STOP_G4.txt").write_text("G4a failed: Fig 4C universe or refit beta not reproduced\n")
    sys.exit("STOP G4a")
if not g4b:
    (out / "STOP_G4.txt").write_text("untestable as a replication: the projection on the V5 gene set does not carry the discovery direction\n")
    sys.exit("STOP G4b")
print("G4 pass")
