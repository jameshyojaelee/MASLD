"""Job 1041 (review R04).

Part A: within-dataset stage-label permutation null for the Fig 4C primary model
  score ~ stage_ordinal + factor(dataset), 117 programs, donors as units.
  One permutation of donor stage labels within dataset is applied to all programs
  jointly (preserves between-program correlation). Seed 20261015. 9,999 draws.
  Reports per-program permutation p (floor 1/10,000), BH over 117, and max-|t| FWER p.
  Statistic: OLS t (homoskedastic) and HC3 t, both regenerated end to end per draw.

Part B: donor bootstrap of the Fig 4D content-lineage argmax (80 match / 37 differ).
  Re-implements label_vs_content.py: donor-collapsed per-lineage mean log2(CPM+1),
  gene z across six lineages, program L1-weighted mean z, argmax lineage.
  Resamples donors with replacement within each lineage, 200 draws, seed 20261015.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
NPERM = int(sys.argv[2])
NBOOT = int(sys.argv[3])
SEED = 20261015

# ---------------- Part A
V3 = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot"
META = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/metadata/donor_metadata_stage_corrected_donor.tsv"
STAGES = {"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2}
scores = pd.read_csv(V3 / "donor_program_scores_primary.tsv", sep="\t", usecols=["cell_type", "module", "donor", "score"])
m = pd.read_csv(META, sep="\t", keep_default_na=False)
m = m[m.exclude_stage_analysis.astype(str).str.upper().isin(["FALSE"]) & m.disease_stage_coarse.isin(STAGES)]
m = m.assign(stage=m.disease_stage_coarse.map(STAGES).astype(float))
donors = m[["sample", "dataset", "stage"]].rename(columns={"sample": "donor"}).reset_index(drop=True)
didx = {d: i for i, d in enumerate(donors.donor)}
dd = scores.merge(donors, on="donor")
assert dd.donor.nunique() == 62


def bh(p, n=117):
    p = np.asarray(p, float); o = np.argsort(p)
    r = p[o] * n / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty_like(q); out[o] = np.minimum(q, 1); return out


progs = []
for (ct, mod), d in dd.groupby(["cell_type", "module"]):
    idx = np.array([didx[x] for x in d.donor])
    codes, ds = np.unique(d.dataset.to_numpy(), return_inverse=True)
    X_ds = np.eye(len(codes))[ds][:, 1:]
    progs.append(dict(ct=ct, mod=mod, idx=idx, ds=ds, k=len(codes), y=d.score.to_numpy(float), Xds=X_ds))


def tstats(stage):
    t_ols = np.empty(len(progs)); t_hc3 = np.empty(len(progs))
    for i, p in enumerate(progs):
        s = stage[p["idx"]]
        X = np.column_stack([np.ones(len(s)), s, p["Xds"]])
        inv = np.linalg.inv(X.T @ X)
        b = inv @ (X.T @ p["y"])
        r = p["y"] - X @ b
        df = len(s) - X.shape[1]
        t_ols[i] = b[1] / np.sqrt(r @ r / df * inv[1, 1])
        h = np.einsum("ij,jk,ik->i", X, inv, X)
        hc = inv @ (X.T @ ((r / (1 - h))[:, None] ** 2 * X)) @ inv
        t_hc3[i] = b[1] / np.sqrt(hc[1, 1])
    return t_ols, t_hc3


rng = np.random.default_rng(SEED)
base = donors.stage.to_numpy()
groups = [np.where(donors.dataset.to_numpy() == g)[0] for g in donors.dataset.unique()]
obs_o, obs_h = tstats(base)
ex_o = np.zeros(len(progs)); ex_h = np.zeros(len(progs)); mx_o = np.zeros(len(progs)); mx_h = np.zeros(len(progs))
for _ in range(NPERM):
    s = base.copy()
    for g in groups:
        s[g] = base[rng.permutation(g)]
    to, th = tstats(s)
    ex_o += np.abs(to) >= np.abs(obs_o) - 1e-12
    ex_h += np.abs(th) >= np.abs(obs_h) - 1e-12
    mx_o += np.abs(to).max() >= np.abs(obs_o) - 1e-12
    mx_h += np.abs(th).max() >= np.abs(obs_h) - 1e-12
A = pd.DataFrame(dict(cell_type=[p["ct"] for p in progs], module=[p["mod"] for p in progs],
                      t_ols=obs_o, t_hc3=obs_h,
                      perm_p_ols=(ex_o + 1) / (NPERM + 1), perm_p_hc3=(ex_h + 1) / (NPERM + 1),
                      maxT_fwer_ols=(mx_o + 1) / (NPERM + 1), maxT_fwer_hc3=(mx_h + 1) / (NPERM + 1)))
A["perm_q_ols_bh117"] = bh(A.perm_p_ols)
A["perm_q_hc3_bh117"] = bh(A.perm_p_hc3)
A.to_csv(OUT / "within_dataset_permutation_117.tsv", sep="\t", index=False)
print(f"[A] {NPERM} within-dataset permutations, seed {SEED}, p floor {1/(NPERM+1):.2e}")
print(A.sort_values("perm_p_ols").head(8).to_string(index=False))
print("[A] perm q<0.05 (OLS t):", int((A.perm_q_ols_bh117 < .05).sum()), " (HC3 t):", int((A.perm_q_hc3_bh117 < .05).sum()))

# ---------------- Part B
SCDIR = ROOT / "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_fstage_pseudobulk"
REG = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv"
MEM = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
PUB = ROOT / "Analysis/SingleCell/results_gpu_v2/program_label_vs_content/plc-20260907T105813Z/program_label_vs_content.tsv"
L2L = {"hepatocytes": "Hepatocytes", "fibroblasts": "Fibroblasts", "cholangiocytes": "Cholangiocytes",
       "macrophages": "Macrophages", "tcells": "T_cells"}
LIN = ["Hepatocytes", "Fibroblasts", "Cholangiocytes", "Macrophages", "T_cells", "Endothelial_cells"]
s2d = {}
for ds in ("GSE244832", "GSE185477", "GSE202379", "GSE136103"):
    dp = pd.read_csv(ROOT / f"data/{ds}/metadata/donor_pairing.csv", dtype=str)
    for _, r in dp.iterrows():
        for s in str(r["rna_srrs"]).split(";"):
            if s and s != "nan":
                s2d[s.strip()] = f"{ds}_{r['donor_id']}"
logcpm = {}
for ln in LIN:
    x = pd.read_csv(SCDIR / f"pseudobulk_{ln}.tsv.gz", sep="\t", index_col=0)
    x = x.groupby(level=0).sum()
    cols = pd.Series(x.columns).map(s2d).fillna(pd.Series(x.columns)).values
    x = x.T.groupby(cols).sum().T
    cpm = x.div(x.sum(axis=0).replace(0, np.nan), axis=1) * 1e6
    logcpm[ln] = np.log2(cpm + 1).astype(np.float32)
genes = sorted(set.intersection(*[set(v.index) for v in logcpm.values()]))
reg = pd.read_csv(REG, sep="\t"); mem = pd.read_csv(MEM, sep="\t")
W = {}
for uid, g in mem.groupby("program_uid"):
    w = g.dropna(subset=["mapped_symbol"]).copy()
    w["wt"] = pd.to_numeric(w["original_l1_weight"], errors="coerce")
    w = w.dropna(subset=["wt"])
    W[uid] = w.groupby("mapped_symbol")["wt"].sum()
declared = {r.program_uid: L2L[r.cell_type] for r in reg.itertuples()}


def content(prof):
    P = pd.DataFrame(prof).dropna(how="all")
    Z = P.sub(P.mean(axis=1), axis=0).div(P.std(axis=1).replace(0, np.nan), axis=0)
    out = {}
    for uid, wt in W.items():
        wt = wt[wt.index.isin(Z.index)]
        if len(wt) < 5:
            out[uid] = None; continue
        wt = wt / wt.sum()
        s = (Z.loc[wt.index].T * wt.values).T.sum(axis=0)
        out[uid] = s.idxmax()
    return pd.Series(out)


full = content({ln: v.mean(axis=1) for ln, v in logcpm.items()})
pub = pd.read_csv(PUB, sep="\t").set_index("program_uid").content_lineage
agree_pub = int((full.reindex(pub.index) == pub).sum())
match_full = int(sum(full[u] == declared[u] for u in full.index))
print(f"[B] re-derived content lineage agrees with published for {agree_pub}/117; matches declared {match_full}/117")
rng = np.random.default_rng(SEED)
calls = []
for b in range(NBOOT):
    prof = {}
    for ln, v in logcpm.items():
        j = rng.integers(0, v.shape[1], v.shape[1])
        prof[ln] = v.iloc[:, j].mean(axis=1)
    calls.append(content(prof))
C = pd.concat(calls, axis=1)
stab = (C.eq(full, axis=0)).mean(axis=1)
match_frac = pd.Series({u: float((C.loc[u] == declared[u]).mean()) for u in C.index})
B = pd.DataFrame(dict(declared=pd.Series(declared), content_full=full, boot_same_call=stab, boot_match_declared=match_frac))
B.to_csv(OUT / "content_lineage_bootstrap.tsv", sep="\t")
n_match_stable = int(((B.content_full == B.declared) & (B.boot_same_call >= .95)).sum())
n_diff_stable = int(((B.content_full != B.declared) & (B.boot_same_call >= .95)).sum())
print(f"[B] {NBOOT} donor bootstraps, seed {SEED}: stable (>=95% same call) matches {n_match_stable}, stable differs {n_diff_stable}, "
      f"unstable {117 - n_match_stable - n_diff_stable}")
dec = pd.Series(declared).reindex(C.index)
per_draw = [int((C[c] == dec).sum()) for c in C.columns]
print("[B] per-draw match count: median", int(np.median(per_draw)), "range", min(per_draw), "-", max(per_draw))
for uid in ["hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e"]:
    print("[B]", uid, B.loc[uid].to_dict())
print("[B] donors per lineage:", {ln: int(v.shape[1]) for ln, v in logcpm.items()})
json.dump(dict(seed=SEED, nperm=NPERM, nboot=NBOOT), open(OUT / "run_params.json", "w"))
