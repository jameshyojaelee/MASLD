"""Fig 4G (node-level) and 4H (system-level) independent recomputation on the 367 disease units.

Upstream analysis outputs only (not the panel source tables):
  system scores  : molecular_systems_inference/.../tables/participant_system_scores.tsv.gz
  feature scores : molecular_layers/.../pathway_tf/pathways/<collection>/donor_scores.tsv.gz
                   molecular_layers/.../programs/hotspot/hotspot_participant_scores.tsv.gz
  continuum axes : candidates/hac-continuum-corrected774.../projection/participant_scores.tsv (fixed_projection)
                   candidates/hac-continuum-corrected774.../unsupervised/participant_scores.tsv (signature_pc1)
  geometry       : molecular_systems_atlas_20260824/geometry/latent_coordinates.tsv.gz (frozen)
Models (as the 435-unit producer 70_infer_molecular_systems.R defines them; validated below by
reproducing its 435-unit cohort tables):
  stage     : y ~ stage_z + sex           (stage_z standardized within the fitted units)
  continuum : y ~ axis_z + factor(stage) + sex  (axis_z standardized within the fitted units)
  y = score standardized within cohort over all cohort units ("zall") or re-standardized within
      the fitted disease units ("zdis"); classical OLS SE; fixed-effect inverse-variance meta.
Disease units: GSE162694 stage 0-4 (109), GSE213621 stage codes 1-3 (258); GSE213621 code 0 = 68 controls.
"""
import numpy as np, pandas as pd, gzip, sys, os
from scipy import stats
root = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
out = f"{root}/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
HAC = f"{root}/RNA-seq/results/histology_anchored_continuum"
SYS = f"{HAC}/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables"
ML = f"{HAC}/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z"
CORE = f"{HAC}/candidates/hac-continuum-corrected774-20261009T171017Z"
GEO = f"{root}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/geometry"
EVAL = ["GSE162694", "GSE213621"]


def ols_coef(X, Y):
    """X: n x p design (with intercept), Y: n x m. returns beta (p x m), se (p x m)"""
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    R = Y - X @ B
    df = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (R ** 2).sum(0) / df
    se = np.sqrt(np.outer(np.diag(XtX_inv), s2))
    return B, se, df


def zs(x):
    return (x - x.mean()) / x.std(ddof=1)


def fixed_meta(b, s):
    w = 1 / s ** 2
    bb = (w * b).sum(0) / w.sum(0)
    return bb, np.sqrt(1 / w.sum(0))


# ---- covariates and axes
fp = pd.read_csv(f"{CORE}/projection/participant_scores.tsv", sep="\t")
fp = fp[fp.axis_id == "fixed_projection"][["sample_id", "dataset", "axis_raw", "inferred_sex", "fibrosis_stage"]]
sp = pd.read_csv(f"{CORE}/unsupervised/participant_scores.tsv", sep="\t")
sp = sp[sp.axis_id == "signature_pc1"][["sample_id", "dataset", "axis_raw"]].rename(columns={"axis_raw": "axis_sigpc1"})
cov = fp.rename(columns={"axis_raw": "axis_fp"}).merge(sp, on=["sample_id", "dataset"], how="left")
cov = cov[cov.dataset.isin(EVAL)].copy()
print("units per cohort (all):", cov.groupby("dataset").size().to_dict())
cov["stage_complete"] = cov.fibrosis_stage.notna() & cov.inferred_sex.notna()
cov["disease"] = cov.stage_complete & ~((cov.dataset == "GSE213621") & (cov.fibrosis_stage == 0))
print("stage-complete (435 basis):", cov[cov.stage_complete].groupby("dataset").size().to_dict(),
      "; disease:", cov[cov.disease].groupby("dataset").size().to_dict())


def fit_family(wide, units_mask_col, ymode, axis_col):
    """wide: DataFrame indexed by sample_id, columns = features, scores standardized per cohort over all units
    (zall) supplied; returns per-feature fixed-meta stage beta & continuum beta and cohort betas"""
    res = {}
    for ds in EVAL:
        c = cov[(cov.dataset == ds) & cov[units_mask_col]].set_index("sample_id")
        c = c[c.index.isin(wide.index)]
        Y = wide.loc[c.index].values.astype(float)
        if ymode == "zdis":
            Y = (Y - Y.mean(0)) / Y.std(0, ddof=1)
        sex = (c.inferred_sex == c.inferred_sex.unique()[0]).astype(float).values
        stz = zs(c.fibrosis_stage.astype(float).values)
        X1 = np.column_stack([np.ones(len(c)), stz, sex])
        B1, S1, _ = ols_coef(X1, Y)
        ax = zs(c[axis_col].values.astype(float))
        lev = sorted(c.fibrosis_stage.unique())
        D = np.column_stack([(c.fibrosis_stage.values == L).astype(float) for L in lev[1:]])
        X2 = np.column_stack([np.ones(len(c)), ax, D, sex])
        B2, S2, _ = ols_coef(X2, Y)
        res[ds] = dict(n=len(c), stage_b=B1[1], stage_se=S1[1], cont_b=B2[1], cont_se=S2[1])
    sb, ss = fixed_meta(np.vstack([res[d]["stage_b"] for d in EVAL]), np.vstack([res[d]["stage_se"] for d in EVAL]))
    cb, cs = fixed_meta(np.vstack([res[d]["cont_b"] for d in EVAL]), np.vstack([res[d]["cont_se"] for d in EVAL]))
    df = pd.DataFrame({"feature": wide.columns, "stage_beta": sb, "stage_se": ss, "cont_beta": cb, "cont_se": cs})
    for d in EVAL:
        df[f"stage_b_{d}"] = res[d]["stage_b"]; df[f"cont_b_{d}"] = res[d]["cont_b"]
        df[f"cont_se_{d}"] = res[d]["cont_se"]; df[f"n_{d}"] = res[d]["n"]
    return df


# ================= 4H: systems
ss = pd.read_csv(f"{SYS}/participant_system_scores.tsv.gz", sep="\t")
ss = ss[(ss.aggregation_method == "collection_balanced") & ss.dataset.isin(EVAL)]
W = ss.pivot_table(index="sample_id", columns="community_id", values="system_score_z")
# validation: reproduce the producer's 435-unit cohort tables
v = fit_family(W, "stage_complete", "zall", "axis_fp")
ref_c = pd.read_csv(f"{SYS}/system_continuum_cohort_models.tsv", sep="\t")
ref_s = pd.read_csv(f"{SYS}/system_fibrosis_stage_cohort_models.tsv", sep="\t")
rc = ref_c[(ref_c.aggregation_method == "collection_balanced") & (ref_c.axis_id == "fixed_projection")]
rs = ref_s[ref_s.aggregation_method == "collection_balanced"]
for d in EVAL:
    a = v.set_index("feature")[f"cont_b_{d}"]; b = rc[rc.dataset == d].set_index("community_id").beta
    a2 = v.set_index("feature")[f"stage_b_{d}"]; b2 = rs[rs.dataset == d].set_index("community_id").beta
    print(f"validation 435 {d}: max|d cont beta| {np.abs(a - b.reindex(a.index)).max():.2e}; max|d stage beta| {np.abs(a2 - b2.reindex(a2.index)).max():.2e}")
v_sig = fit_family(W, "stage_complete", "zall", "axis_sigpc1")
rcs = ref_c[(ref_c.aggregation_method == "collection_balanced") & (ref_c.axis_id == "signature_pc1")]
for d in EVAL:
    a = v_sig.set_index("feature")[f"cont_b_{d}"]; b = rcs[rcs.dataset == d].set_index("community_id").beta
    print(f"validation 435 signature_pc1 {d}: max|d| {np.abs(a - b.reindex(a.index)).max():.2e}")

grid = []
for ymode in ["zall", "zdis"]:
    for axis in ["axis_fp", "axis_sigpc1"]:
        f = fit_family(W, "disease", ymode, axis)
        f["ymode"] = ymode; f["axis"] = axis
        grid.append(f)
G = pd.concat(grid)
G.to_csv(f"{out}/r4h_systems_367_grid.tsv", sep="\t", index=False)
for (ym, ax), f in G.groupby(["ymode", "axis"]):
    rho = stats.spearmanr(f.stage_beta, f.cont_beta).correlation
    same = int((np.sign(f.stage_beta) == np.sign(f.cont_beta)).sum())
    bp = int(((f.stage_beta > 0) & (f.cont_beta > 0)).sum()); bn = int(((f.stage_beta < 0) & (f.cont_beta < 0)).sum())
    print(f"367 {ym} {ax}: n={len(f)} rho={rho:.4f} same={same} (both+ {bp}, both- {bn}); n units {f.iloc[0][[c for c in f.columns if c.startswith('n_')]].to_dict()}")
# 435 roster version for the display-choice comparison
v["ymode"] = "zall"; v["axis"] = "axis_fp"; v.to_csv(f"{out}/r4h_systems_435_fp.tsv", sep="\t", index=False)
print("435 fp: rho", round(stats.spearmanr(v.stage_beta, v.cont_beta).correlation, 4),
      "same", int((np.sign(v.stage_beta) == np.sign(v.cont_beta)).sum()))

# ================= 4G: nodes
geo = pd.read_csv(f"{GEO}/latent_coordinates.tsv.gz", sep="\t")
print("geometry nodes", len(geo), geo.component_role.value_counts().to_dict())
frames = []
for coll in ["hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc"]:
    ids = set(geo[geo.collection == coll].feature_id)
    it = pd.read_csv(f"{ML}/pathway_tf/pathways/{coll}/donor_scores.tsv.gz", sep="\t", chunksize=2_000_000)
    for ch in it:
        ch = ch[ch.dataset.isin(EVAL) & ch.set_id.isin(ids)]
        frames.append(ch.rename(columns={"set_id": "feature_id", "pathway_score": "score"})[["feature_id", "sample_id", "dataset", "score"]])
hs = pd.read_csv(f"{ML}/programs/hotspot/hotspot_participant_scores.tsv.gz", sep="\t",
                 usecols=["program_uid", "sample_id", "dataset", "score_raw"])
hs = hs[hs.dataset.isin(EVAL) & hs.program_uid.isin(set(geo[geo.collection == "hotspot"].feature_id))]
frames.append(hs.rename(columns={"program_uid": "feature_id", "score_raw": "score"}))
F = pd.concat(frames)
del frames
print("feature score rows", len(F), "features", F.feature_id.nunique())
# standardize within cohort over all units of that cohort (producer convention for feature scores)
F["z"] = F.groupby(["feature_id", "dataset"]).score.transform(lambda x: (x - x.mean()) / x.std(ddof=1))
WF = F.pivot_table(index="sample_id", columns="feature_id", values="z")
WFr = F.pivot_table(index="sample_id", columns="feature_id", values="score")
del F
node = []
for ymode, Wsrc in [("zall", WF), ("zdis", WFr)]:
    for axis in ["axis_fp", "axis_sigpc1"]:
        f = fit_family(Wsrc, "disease", ymode, axis)
        f["ymode"] = ymode; f["axis"] = axis
        node.append(f)
N = pd.concat(node).merge(geo[["feature_id", "collection", "community_id", "umap_1", "umap_2", "component_role"]],
                          left_on="feature", right_on="feature_id", how="left")
N.to_csv(f"{out}/r4g_nodes_367_grid.tsv.gz", sep="\t", index=False, compression="gzip")
for (ym, ax), f in N.groupby(["ymode", "axis"]):
    q = f.cont_beta.abs().quantile(0.98)
    print(f"4G {ym} {ax}: nodes {len(f)}; median beta {f.cont_beta.median():.4f}; 98th pct |beta| {q:.4f}; "
          f"n |beta|>q {int((f.cont_beta.abs() > q).sum())}; frac positive {(f.cont_beta > 0).mean():.3f}")
