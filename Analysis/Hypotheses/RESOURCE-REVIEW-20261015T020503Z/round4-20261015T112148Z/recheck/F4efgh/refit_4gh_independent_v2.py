#!/usr/bin/env python3
"""Independent re-derivation of Figure 4G/4H system coefficients and the S3
layer effect rho on two rosters (435 as published; 367 disease units).

Written from the model definitions in 62_render_stage_continuum_comparison.R,
20_pathway_continuum.R and 67_overlay_molecular_systems_atlas.R; does not read
any candidate output except for the final comparison step.
Stage model per cohort:     outcome_z ~ stage_z + sex       (fixed IV meta)
Continuum model per cohort: outcome_z ~ axis_z + factor(stage) + sex (fixed IV meta)
System beta: mean over collections of the per-collection median node beta.
"""
import os, sys, gzip
import numpy as np
import pandas as pd
from scipy import stats

np.random.seed(20260923)
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = sys.argv[1]
RESTD = (sys.argv[2] == "restd") if len(sys.argv) > 2 else True
os.makedirs(OUT, exist_ok=True)
A = f"{R}/RNA-seq/results/histology_anchored_continuum/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z"
MAN = f"{R}/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/five_cohort_sample_manifest.tsv"
PROJ = f"{R}/RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261009T171017Z/projection/participant_scores.tsv"
ATLAS = f"{R}/figures/main/fig4_singlecell_programs/source_tables/current_candidate/molecular_systems_atlas_20260824/overlays/all_node_stage_continuum_overlays.tsv.gz"
INST = f"{R}/figures/main/fig4_singlecell_programs/source_tables/corrected_774/stage_continuum_comparison"
COHORTS = ["GSE162694", "GSE213621"]
COLLS = ["hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc"]

man = pd.read_csv(MAN, sep="\t", keep_default_na=False, na_values=["", "NA"])
proj = pd.read_csv(PROJ, sep="\t")
proj = proj[proj.axis_id == "fixed_projection"][["sample_id", "axis_raw"]]
man = man[man.dataset.isin(COHORTS)].merge(proj, on="sample_id", how="left")
base = man[man.fibrosis_stage.notna() & man.inferred_sex.notna() & np.isfinite(man.axis_raw)]
rosters = {
    "as_published_435": base,
    "disease_only_367": base[base.group_binary == "Disease"],
}
census = []
for rn, d in rosters.items():
    for (ds, g, st), n in d.groupby(["dataset", "group_binary", "fibrosis_stage"]).size().items():
        census.append((rn, ds, g, st, n))
pd.DataFrame(census, columns=["roster", "dataset", "group", "stage", "n"]).to_csv(f"{OUT}/roster_census_recheck.tsv", sep="\t", index=False)


def zs(x):
    x = np.asarray(x, float)
    return (x - x.mean(axis=-1, keepdims=True)) / x.std(axis=-1, ddof=1, keepdims=True)


def ols_coef(X, Y, j):
    """Y: n x p outcomes. return beta, se for column j."""
    XtX_inv = np.linalg.inv(X.T @ X)
    B = XtX_inv @ X.T @ Y
    resid = Y - X @ B
    df = X.shape[0] - X.shape[1]
    s2 = (resid ** 2).sum(0) / df
    return B[j], np.sqrt(s2 * XtX_inv[j, j])


def fit_feature_matrix(wide, units, restd_stage=True, restd_cont=True):
    """wide: features x samples (cohort subset). Returns stage and continuum beta/se."""
    res = {}
    for ds in COHORTS:
        u = units[units.dataset == ds].set_index("sample_id")
        Y = wide.loc[:, wide.columns.intersection(u.index)]
        u = u.loc[Y.columns]
        Yz = zs(Y.values)  # features x n
        Ys = Yz if restd_stage else Y.values
        Yc = Yz if restd_cont else Y.values
        ok = np.isfinite(Yz).all(1)
        sex = (u.inferred_sex.values == u.inferred_sex.values[0]).astype(float)
        stage_z = zs(u.fibrosis_stage.values.astype(float))
        Xs = np.column_stack([np.ones(len(u)), stage_z, sex])
        bs, ss = ols_coef(Xs, Ys[ok].T, 1)
        lev = sorted(u.fibrosis_stage.unique())
        dummies = [(u.fibrosis_stage.values == l).astype(float) for l in lev[1:]]
        axis_z = zs(u.axis_raw.values)
        Xc = np.column_stack([np.ones(len(u))] + dummies + [sex, axis_z])
        bc, sc = ols_coef(Xc, Yc[ok].T, Xc.shape[1] - 1)
        idx = Y.index[ok]
        res[ds] = pd.DataFrame({"sb": bs, "ss": ss, "cb": bc, "cs": sc, "n": len(u)}, index=idx)
    a, b = res[COHORTS[0]], res[COHORTS[1]]
    j = a.join(b, lsuffix="_1", rsuffix="_2", how="inner")
    out = pd.DataFrame(index=j.index)
    for k in ["s", "c"]:
        w1, w2 = 1 / j[f"{k}s_1"] ** 2, 1 / j[f"{k}s_2"] ** 2
        out[f"{k}_beta"] = (w1 * j[f"{k}b_1"] + w2 * j[f"{k}b_2"]) / (w1 + w2)
    out["n1"] = j.n_1; out["n2"] = j.n_2
    return out


atlas = pd.read_csv(ATLAS, sep="\t", usecols=["feature_id", "collection", "community_id", "component_role"])
atlas = atlas[atlas.component_role == "main_system_component"]
inst_path = pd.read_csv(f"{INST}/pathway_stage_vs_continuum.tsv.gz", sep="\t")
inst_hot = pd.read_csv(f"{INST}/hotspot_stage_vs_continuum.tsv", sep="\t")

all_feat = []
for rn, units in rosters.items():
    for coll in COLLS:
        sc = pd.read_csv(f"{A}/pathway_tf/pathways/{coll}/donor_scores.tsv.gz", sep="\t")
        sc = sc[sc.dataset.isin(COHORTS)]
        wide = sc.pivot(index="set_id", columns="sample_id", values="pathway_score")
        f = fit_feature_matrix(wide, units, restd_stage=RESTD, restd_cont=True)
        f["collection"] = coll; f["roster"] = rn
        all_feat.append(f.rename_axis("feature_id").reset_index())
        print(rn, coll, len(f), flush=True)
    hs = pd.read_csv(f"{A}/programs/hotspot/hotspot_participant_scores.tsv.gz", sep="\t",
                     usecols=["program_uid", "sample_id", "dataset", "outcome_z"])
    hs = hs[hs.dataset.isin(COHORTS)]
    wide = hs.pivot(index="program_uid", columns="sample_id", values="outcome_z")
    f = fit_feature_matrix(wide, units, restd_stage=RESTD, restd_cont=RESTD)
    f["collection"] = "hotspot"; f["roster"] = rn
    all_feat.append(f.rename_axis("feature_id").reset_index())
feat = pd.concat(all_feat)
feat.to_csv(f"{OUT}/feature_betas_recheck.tsv.gz", sep="\t", index=False)

# Layer rho restricted to installed estimable feature sets
sys.exit(0)
rows = []
for (rn, coll), g in feat[feat.feature_id.isin(testable)].groupby(["roster", "collection"]):
    g = g.dropna(subset=["s_beta", "c_beta"])
    rho = stats.spearmanr(g.s_beta, g.c_beta).correlation
    rows.append((rn, coll, len(g), rho, int((np.sign(g.s_beta) == np.sign(g.c_beta)).sum())))
lay = pd.DataFrame(rows, columns=["roster", "collection", "n", "rho", "same_direction"])
lay.to_csv(f"{OUT}/layer_rho_recheck.tsv", sep="\t", index=False)
print(lay.to_string())

# System aggregation
sysrows = []
nodes_out = []
for rn in rosters:
    f = feat[feat.roster == rn].merge(atlas, on=["feature_id", "collection"], how="inner")
    nodes_out.append(f)
    fam = f.groupby(["community_id", "collection"]).agg(s=("s_beta", "median"), c=("c_beta", "median")).reset_index()
    sy = fam.groupby("community_id").agg(stage_beta=("s", "mean"), continuum_beta=("c", "mean")).reset_index()
    sy["roster"] = rn
    sysrows.append(sy)
    rho = stats.spearmanr(sy.stage_beta, sy.continuum_beta).correlation
    same = int((np.sign(sy.stage_beta) == np.sign(sy.continuum_beta)).sum())
    q98 = np.quantile(np.abs(f.c_beta.dropna()), 0.98, method="median_unbiased")
    print(f"SYSTEM {rn}: n_sys={len(sy)} rho={rho:.6f} same={same} n_nodes={len(f)} q98abs_cont={q98:.6f} "
          f"S01 stage={sy.loc[sy.community_id=='system_01','stage_beta'].item():.6f} cont={sy.loc[sy.community_id=='system_01','continuum_beta'].item():.6f}")
pd.concat(sysrows).to_csv(f"{OUT}/system_betas_recheck.tsv", sep="\t", index=False)
pd.concat(nodes_out).to_csv(f"{OUT}/node_betas_recheck.tsv.gz", sep="\t", index=False)
