#!/usr/bin/env python3
"""U10 independent re-derivation of the molecular-system numbers (release r9).

Built from METHODS.md and the figure legends, starting at the per-participant pathway and
Hotspot program scores (not re-scored from expression). Producer tables are read only to
compare against.

  A. S3L: Pearson r between the two co-primary axes per cohort.
  B. Testable nodes (5,435 of 5,439) and collection-balanced participant system scores
     (feature z within cohort over all units; median within system x collection; mean over
     collections; z within cohort).
  C. Figure 4H donor-first on 367 disease units (rho, same sign, both-interval split) with a
     stratified donor bootstrap (resample units within cohort, B = 1000) at seed 20261017
     and two alternate seeds. Figure 4G node continuum betas, system medians and intervals.
     Feature-first 367 comparison values.
  D. 435 stage-complete units: 33 primary continuum calls and 27 stage calls (system scores);
     S4S feature-first rho and same-sign count from node betas.
  E. Hotspot S4I (93 of 117) and S4J classes (54/39/10/13) on 435 units.

Usage: u10_systems.py <out_dir>
"""
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats

OUT = sys.argv[1].rstrip("/")
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
HAC = f"{R}/RNA-seq/results/histology_anchored_continuum"
CORE = f"{HAC}/candidates/hac-continuum-corrected774-20261009T171017Z"
ML = f"{HAC}/molecular_layers/hac-molecular-layers-corrected774-20261009T171017Z"
SYS = f"{HAC}/molecular_systems_inference/hac-system-inference-corrected774-20261009T171017Z/tables"
GEO = f"{R}/figures/candidates/histology-continuum-molecular-systems-atlas-20260824T195148Z/geometry"
OVL = f"{R}/figures/candidates/continuum-774-panels-20261011T152106Z/molecular_systems_atlas/overlays"
PANEL4 = f"{R}/figures/candidates/panel-deep-fixes-20261017T035633Z/fig4/work"
S4J_SRC = (f"{R}/figures/main/fig4_singlecell_programs/source_tables/corrected_774/"
           "stage_continuum_comparison/hotspot_stage_vs_continuum.tsv")
EVAL = ["GSE162694", "GSE213621"]
AXES = ["signature_pc1", "fixed_projection"]
COLLS = ["hallmark", "kegg", "reactome", "go_bp", "go_mf", "go_cc"]
BOOT_SEEDS = [20261017, 20261101, 20261202]
B = 1000

log_lines = []


def log(msg):
    print(msg, flush=True)
    log_lines.append(str(msg))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def zscore(x):
    x = np.asarray(x, dtype=float)
    return (x - x.mean()) / x.std(ddof=1)


def zcols(Y):
    return (Y - Y.mean(0)) / Y.std(0, ddof=1)


def ols_coef1(X, Y):
    """OLS coefficient on column 1 and its classical SE for each column of Y."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    coef = XtX_inv @ X.T @ Y
    resid = Y - X @ coef
    dof = X.shape[0] - np.linalg.matrix_rank(X)
    s2 = (resid ** 2).sum(0) / dof
    return coef[1], np.sqrt(XtX_inv[1, 1] * s2)


def fixed_meta(betas, ses):
    b = np.vstack(betas)
    w = 1.0 / np.vstack(ses) ** 2
    est = (w * b).sum(0) / w.sum(0)
    se = np.sqrt(1.0 / w.sum(0))
    p = 2 * stats.norm.sf(np.abs(est / se))
    return est, se, p


def bh(p, family=None):
    """Benjamini-Hochberg over a declared family size (NaN p values stay NaN)."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = family if family is not None else ok.sum()
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, len(pv) + 1)
    adj = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    out = np.empty_like(adj)
    out[order] = adj
    q[ok] = out
    return q


def sex_m(c):
    return (c.inferred_sex.values == "M").astype(float)


def design_stage(c):
    return np.column_stack([np.ones(len(c)), zscore(c.fibrosis_stage.values), sex_m(c)])


def design_cont(c, axis_col):
    st = c.fibrosis_stage.values
    levels = sorted(np.unique(st))
    dummies = [(st == L).astype(float) for L in levels[1:]]
    cols = [np.ones(len(c)), zscore(c[axis_col].values)] + dummies + [sex_m(c)]
    return np.column_stack(cols)


inputs = {}

# ---------------------------------------------------------------- axes and covariates
proj_path = f"{CORE}/projection/participant_scores.tsv"
unsup_path = f"{CORE}/unsupervised/participant_scores.tsv"
inputs["projection_scores"] = proj_path
inputs["unsupervised_scores"] = unsup_path
proj = pd.read_csv(proj_path, sep="\t")
unsup = pd.read_csv(unsup_path, sep="\t")
fp = proj[proj.axis_id == "fixed_projection"][
    ["sample_id", "dataset", "axis_raw", "inferred_sex", "fibrosis_stage"]].rename(
    columns={"axis_raw": "fixed_projection"})
pc1 = unsup[unsup.axis_id == "signature_pc1"][["sample_id", "dataset", "axis_raw"]].rename(
    columns={"axis_raw": "signature_pc1"})
cov = fp.merge(pc1, on=["sample_id", "dataset"], how="left", validate="one_to_one")
assert len(cov) == 774 and cov.signature_pc1.notna().all()
log(f"units per cohort: {cov.dataset.value_counts().to_dict()}")

# A. S3L
s3l = []
for d, g in cov.groupby("dataset"):
    r_all = stats.pearsonr(g.fixed_projection, g.signature_pc1)[0]
    s3l.append(dict(dataset=d, n_units=len(g), pearson_r=r_all))
s3l = pd.DataFrame(s3l)
s3l.to_csv(f"{OUT}/A_s3l_axis_pearson.tsv", sep="\t", index=False)
log("A. S3L axis Pearson r per cohort\n" + s3l.to_string())

# ---------------------------------------------------------------- features
geo_path = f"{GEO}/latent_coordinates.tsv.gz"
inputs["geometry_latent_coordinates"] = geo_path
geo = pd.read_csv(geo_path, sep="\t")
log(f"geometry nodes {len(geo)}; roles {geo.component_role.value_counts().to_dict()}")
geo = geo[geo.component_role == "main_system_component"][["feature_id", "collection", "community_id"]]
systems = sorted(geo.community_id.unique())
assert len(geo) == 5439 and len(systems) == 43

ev = cov[cov.dataset.isin(EVAL)].copy()
frames = []
for coll in COLLS:
    path = f"{ML}/pathway_tf/pathways/{coll}/donor_scores.tsv.gz"
    inputs[f"pathway_{coll}"] = path
    ids = set(geo.loc[geo.collection == coll, "feature_id"])
    for ch in pd.read_csv(path, sep="\t", chunksize=2_000_000):
        ch = ch[ch.dataset.isin(EVAL) & ch.set_id.isin(ids)]
        frames.append(ch.rename(columns={"set_id": "feature_id", "pathway_score": "score"})[
            ["feature_id", "sample_id", "dataset", "score"]])
hs_path = f"{ML}/programs/hotspot/hotspot_participant_scores.tsv.gz"
inputs["hotspot_scores"] = hs_path
hs_all = pd.read_csv(hs_path, sep="\t", usecols=["program_uid", "sample_id", "dataset", "score_raw", "outcome_z"])
hs_all = hs_all[hs_all.dataset.isin(EVAL)]
hs_ids = set(geo.loc[geo.collection == "hotspot", "feature_id"])
frames.append(hs_all[hs_all.program_uid.isin(hs_ids)].rename(
    columns={"program_uid": "feature_id", "score_raw": "score"})[["feature_id", "sample_id", "dataset", "score"]])
F = pd.concat(frames, ignore_index=True)
del frames

RAW = {}
for d in EVAL:
    unit_ids = ev.loc[ev.dataset == d, "sample_id"].values
    w = F[F.dataset == d].pivot(index="sample_id", columns="feature_id", values="score")
    assert set(w.index) == set(unit_ids), f"{d}: score units differ from axis units"
    RAW[d] = w.loc[unit_ids]
testable = [f for f in geo.feature_id
            if all(f in RAW[d].columns and np.isfinite(RAW[d][f].values).all() for d in EVAL)]
untestable = geo[~geo.feature_id.isin(testable)]
log(f"B. testable main-component nodes {len(testable)} of {len(geo)}; untestable:\n"
    + untestable.to_string())
gt = geo[geo.feature_id.isin(testable)].reset_index(drop=True)
feat = gt.feature_id.values
coll_of = gt.collection.values
sys_of = gt.community_id.values
is_hot = coll_of == "hotspot"
ZALL = {d: pd.DataFrame(zcols(RAW[d][feat].values), index=RAW[d].index, columns=feat) for d in EVAL}

# B. collection-balanced participant system score
groups = {s: [g.feature_id.values for _, g in gt[gt.community_id == s].groupby("collection")] for s in systems}
SYSZ, SYSRAW = {}, {}
for d in EVAL:
    raw = pd.DataFrame({s: np.mean([ZALL[d][cols].median(axis=1).values for cols in groups[s]], axis=0)
                        for s in systems}, index=ZALL[d].index)
    SYSRAW[d] = raw
    SYSZ[d] = (raw - raw.mean()) / raw.std(ddof=1)
pss_path = f"{SYS}/participant_system_scores.tsv.gz"
inputs["stored_participant_system_scores"] = pss_path
pss = pd.read_csv(pss_path, sep="\t")
pss = pss[(pss.aggregation_method == "collection_balanced") & pss.dataset.isin(EVAL)]
for d in EVAL:
    st_raw = pss[pss.dataset == d].pivot(index="sample_id", columns="community_id", values="system_score_raw")
    st_z = pss[pss.dataset == d].pivot(index="sample_id", columns="community_id", values="system_score_z")
    dr = np.abs(SYSRAW[d].values - st_raw.loc[SYSRAW[d].index, systems].values).max()
    dz = np.abs(SYSZ[d].values - st_z.loc[SYSZ[d].index, systems].values).max()
    log(f"B. {d}: system score vs stored, max|d| raw {dr:.2e}, z {dz:.2e}")
nf = pss.groupby("community_id").n_features.first()
nt = gt.community_id.value_counts()
log("B. stored n_features minus testable count, systems differing: "
    + str({s: int(nf[s] - nt[s]) for s in systems if nf[s] != nt[s]}))

# ---------------------------------------------------------------- C. 367 disease units
dis = ev[ev.fibrosis_stage.notna() & ev.inferred_sex.notna()
         & ~((ev.dataset == "GSE213621") & (ev.fibrosis_stage == 0))]
units = {d: dis[dis.dataset == d].reset_index(drop=True) for d in EVAL}
log(f"C. disease units {[len(units[d]) for d in EVAL]}")
YS = {d: SYSZ[d].loc[units[d].sample_id, systems].values for d in EVAL}
YZ = {d: ZALL[d].loc[units[d].sample_id].values for d in EVAL}
YR = {d: RAW[d].loc[units[d].sample_id, feat].values for d in EVAL}
sys_index = {s: k for k, s in enumerate(systems)}
sys_codes = np.array([sys_index[s] for s in sys_of])


def node_y_panel(d, idx):
    """Figure 4G as drawn: pathways z over the fitted units, Hotspot z over all cohort units."""
    y = zcols(YR[d][idx])
    y[:, is_hot] = YZ[d][idx][:, is_hot]
    return y


def fit367(idx, nodes=True):
    sb, ss, cb, cs, nb, ns, nsb, nss = [], [], [], [], [], [], [], []
    for d in EVAL:
        c = units[d].iloc[idx[d]]
        Xs, Xc = design_stage(c), design_cont(c, "fixed_projection")
        b, s = ols_coef1(Xs, YS[d][idx[d]]); sb.append(b); ss.append(s)
        b, s = ols_coef1(Xc, YS[d][idx[d]]); cb.append(b); cs.append(s)
        if nodes:
            b, s = ols_coef1(Xc, node_y_panel(d, idx[d])); nb.append(b); ns.append(s)
            b, s = ols_coef1(Xs, YZ[d][idx[d]]); nsb.append(b); nss.append(s)
    out = dict(sys_stage=fixed_meta(sb, ss)[0], sys_cont=fixed_meta(cb, cs)[0])
    if nodes:
        out["node_cont"] = fixed_meta(nb, ns)[0]
        out["node_stage"] = fixed_meta(nsb, nss)[0]
    return out


def sys_median(v):
    return np.array([np.median(v[sys_codes == k]) for k in range(len(systems))])


def feature_first(v):
    t = pd.DataFrame(dict(s=sys_of, c=coll_of, v=v))
    return t.groupby(["s", "c"]).v.median().groupby("s").mean().reindex(systems).values


obs_idx = {d: np.arange(len(units[d])) for d in EVAL}
obs = fit367(obs_idx)
rho_df = stats.spearmanr(obs["sys_stage"], obs["sys_cont"])[0]
same_df = int((np.sign(obs["sys_stage"]) == np.sign(obs["sys_cont"])).sum())
ffx, ffy = feature_first(obs["node_stage"]), feature_first(obs["node_cont"])
rho_ff = stats.spearmanr(ffx, ffy)[0]
same_ff = int((np.sign(ffx) == np.sign(ffy)).sum())
med_obs = sys_median(obs["node_cont"])
log(f"C. 4H donor-first rho {rho_df:.4f}, same sign {same_df}/43, both>0 "
    f"{int(((obs['sys_stage'] > 0) & (obs['sys_cont'] > 0)).sum())}, both<0 "
    f"{int(((obs['sys_stage'] < 0) & (obs['sys_cont'] < 0)).sum())}")
log(f"C. 4H feature-first rho {rho_ff:.4f}, same sign {same_ff}/43")
log(f"C. 4G system medians > 0: {int((med_obs > 0).sum())}/43; nodes {len(feat)}")
# all-unit z for every node (instead of the panel's mixed standardization): sign count only
alt = []
for d in EVAL:
    c = units[d]
    alt.append(ols_coef1(design_cont(c, "fixed_projection"), YZ[d]))
alt_cont = fixed_meta([a[0] for a in alt], [a[1] for a in alt])[0]
log(f"C. 4G medians > 0 with all-unit z for all nodes: {int((sys_median(alt_cont) > 0).sum())}/43")

stored_h = pd.read_csv(f"{PANEL4}/fig4h_systems_donor_first_bootstrap.tsv", sep="\t").set_index("community_id").loc[systems]
stored_g = pd.read_csv(f"{PANEL4}/fig4g_system_medians_bootstrap.tsv", sep="\t").set_index("community_id").loc[systems]
inputs["panel_4h_table"] = f"{PANEL4}/fig4h_systems_donor_first_bootstrap.tsv"
inputs["panel_4g_table"] = f"{PANEL4}/fig4g_system_medians_bootstrap.tsv"
log(f"C. vs panel tables max|d|: stage {np.abs(obs['sys_stage'] - stored_h.stage_beta_donor_first.values).max():.2e}, "
    f"cont {np.abs(obs['sys_cont'] - stored_h.continuum_beta_donor_first.values).max():.2e}, "
    f"ff stage {np.abs(ffx - stored_h.stage_beta_feature_first.values).max():.2e}, "
    f"ff cont {np.abs(ffy - stored_h.continuum_beta_feature_first.values).max():.2e}, "
    f"4G median {np.abs(med_obs - stored_g.median_beta.values).max():.2e}")


def classify(slo, shi, clo, chi):
    s_res = (slo > 0) | (shi < 0)
    c_res = (clo > 0) | (chi < 0)
    res = s_res & c_res
    both_pos = res & (slo > 0) & (clo > 0)
    both_neg = res & (shi < 0) & (chi < 0)
    opp = res & ~both_pos & ~both_neg
    cls = np.where(both_pos, "both_positive", np.where(both_neg, "both_negative",
                   np.where(opp, "opposite", "not_resolved")))
    return cls


boot_rows, boot_sys = [], []
for seed in BOOT_SEEDS:
    rng = np.random.default_rng(seed)
    bx = np.empty((B, 43)); by = np.empty((B, 43)); bm = np.empty((B, 43)); br = np.empty(B)
    for b in range(B):
        idx = {d: rng.integers(0, len(units[d]), len(units[d])) for d in EVAL}
        f = fit367(idx, nodes=True)
        bx[b], by[b] = f["sys_stage"], f["sys_cont"]
        bm[b] = sys_median(f["node_cont"])
        br[b] = stats.spearmanr(f["sys_stage"], f["sys_cont"])[0]
    slo, shi = np.percentile(bx, 2.5, axis=0), np.percentile(bx, 97.5, axis=0)
    clo, chi = np.percentile(by, 2.5, axis=0), np.percentile(by, 97.5, axis=0)
    mlo, mhi = np.percentile(bm, 2.5, axis=0), np.percentile(bm, 97.5, axis=0)
    cls = classify(slo, shi, clo, chi)
    counts = pd.Series(cls).value_counts().to_dict()
    boot_rows.append(dict(seed=seed, B=B, rho_lo=np.percentile(br, 2.5), rho_hi=np.percentile(br, 97.5),
                          both_positive=counts.get("both_positive", 0), both_negative=counts.get("both_negative", 0),
                          opposite=counts.get("opposite", 0), not_resolved=counts.get("not_resolved", 0),
                          g_medians_interval_excludes_0=int(((mlo > 0) | (mhi < 0)).sum())))
    boot_sys.append(pd.DataFrame(dict(seed=seed, community_id=systems, stage_beta=obs["sys_stage"], stage_lo=slo,
                                      stage_hi=shi, cont_beta=obs["sys_cont"], cont_lo=clo, cont_hi=chi,
                                      class_4h=cls, g_median=med_obs, g_lo=mlo, g_hi=mhi)))
    if seed == 20261017:
        log(f"C. seed {seed} vs panel: max|d| stage_lo {np.abs(slo - stored_h.stage_lo.values).max():.2e}, "
            f"stage_hi {np.abs(shi - stored_h.stage_hi.values).max():.2e}, cont_lo {np.abs(clo - stored_h.cont_lo.values).max():.2e}, "
            f"cont_hi {np.abs(chi - stored_h.cont_hi.values).max():.2e}, g_lo {np.abs(mlo - stored_g.boot_lo.values).max():.2e}, "
            f"g_hi {np.abs(mhi - stored_g.boot_hi.values).max():.2e}")
    log(f"C. bootstrap seed {seed}: {boot_rows[-1]}")
boot_rows = pd.DataFrame(boot_rows)
boot_sys = pd.concat(boot_sys)
boot_rows.to_csv(f"{OUT}/C_4h_bootstrap_by_seed.tsv", sep="\t", index=False)
boot_sys.to_csv(f"{OUT}/C_4gh_systems_by_seed.tsv", sep="\t", index=False)
flip = boot_sys.pivot(index="community_id", columns="seed", values="class_4h")
flip = flip[flip.nunique(axis=1) > 1]
log("C. 4H systems whose class changes across seeds:\n" + (flip.to_string() if len(flip) else "none"))
pd.DataFrame(dict(feature_id=feat, collection=coll_of, community_id=sys_of, node_cont_panel=obs["node_cont"],
                  node_stage_allz=obs["node_stage"])).to_csv(f"{OUT}/C_4g_nodes_367.tsv.gz", sep="\t", index=False)

# ---------------------------------------------------------------- D. 435 stage-complete units
sc = ev[ev.fibrosis_stage.notna() & ev.inferred_sex.notna()]
U = {d: sc[sc.dataset == d].reset_index(drop=True) for d in EVAL}
log(f"D. stage-complete units {[len(U[d]) for d in EVAL]}; GSE213621 stage 0 n "
    f"{int((U['GSE213621'].fibrosis_stage == 0).sum())}")
st_b, st_s, ct = [], [], {a: ([], []) for a in AXES}
for d in EVAL:
    c = U[d]
    Y = SYSZ[d].loc[c.sample_id, systems].values
    b, s = ols_coef1(design_stage(c), Y); st_b.append(b); st_s.append(s)
    for a in AXES:
        b, s = ols_coef1(design_cont(c, a), Y); ct[a][0].append(b); ct[a][1].append(s)
stage_beta, _, stage_p = fixed_meta(st_b, st_s)
stage_q = bh(stage_p, 43)
stage_conc = np.sign(st_b[0]) == np.sign(st_b[1])
stage_call = (stage_q < 0.05) & stage_conc
cont_q, signs = {}, []
for a in AXES:
    cont_q[a] = bh(fixed_meta(*ct[a])[2], 43)
    signs += [np.sign(x) for x in ct[a][0]]
four_dir = np.all(np.vstack(signs) == signs[0], axis=0)
primary = (cont_q["signature_pc1"] < 0.05) & (cont_q["fixed_projection"] < 0.05) & four_dir
mem_path = f"{SYS}/system_inference_membership.tsv"
inputs["stored_membership"] = mem_path
mem = pd.read_csv(mem_path, sep="\t").set_index("community_id").loc[systems]
log(f"D. primary continuum calls {int(primary.sum())} (stored {int(mem.continuum_associated_system.sum())}, "
    f"mismatched systems {list(np.array(systems)[primary != mem.continuum_associated_system.values])}); "
    f"stage calls {int(stage_call.sum())} (stored {int(mem.fibrosis_stage_associated_system.sum())}, mismatched "
    f"{list(np.array(systems)[stage_call != mem.fibrosis_stage_associated_system.values])}); both {int((primary & stage_call).sum())}")
cont_fp_beta = fixed_meta(*ct["fixed_projection"])[0]
log(f"D. donor-first on 435: rho(stage, fixed_projection continuum) {stats.spearmanr(stage_beta, cont_fp_beta)[0]:.4f}, "
    f"same sign {int((np.sign(stage_beta) == np.sign(cont_fp_beta)).sum())}/43")
pd.DataFrame(dict(community_id=systems, stage_beta=stage_beta, stage_q=stage_q, stage_call=stage_call,
                  q_pc1=cont_q["signature_pc1"], q_fp=cont_q["fixed_projection"], four_dir=four_dir,
                  primary=primary, cont_fp_beta=cont_fp_beta)).to_csv(f"{OUT}/D_system_calls_435.tsv", sep="\t", index=False)

# S4S: node betas on 435, feature-first aggregation
nsb, nss, ncb, ncs = [], [], [], []
for d in EVAL:
    c = U[d]
    Y = ZALL[d].loc[c.sample_id].values
    b, s = ols_coef1(design_stage(c), Y); nsb.append(b); nss.append(s)
    b, s = ols_coef1(design_cont(c, "fixed_projection"), Y); ncb.append(b); ncs.append(s)
n_stage = fixed_meta(nsb, nss)[0]
n_cont = fixed_meta(ncb, ncs)[0]
ov_path = f"{OVL}/all_node_stage_continuum_overlays.tsv.gz"
inputs["overlay_nodes_774"] = ov_path
ov = pd.read_csv(ov_path, sep="\t").set_index("feature_id").loc[feat]
log(f"D. node betas vs 774 overlay: max|d| stage {np.nanmax(np.abs(n_stage - ov.stage_beta.values)):.2e}, "
    f"continuum {np.nanmax(np.abs(n_cont - ov.continuum_beta.values)):.2e}")
sx, sy = feature_first(n_stage), feature_first(n_cont)
fbx, fby = sys_median(n_stage), sys_median(n_cont)
log(f"D. S4S feature-first collection-balanced on 435: rho {stats.spearmanr(sx, sy)[0]:.4f}, same sign "
    f"{int((np.sign(sx) == np.sign(sy)).sum())}/43; feature-balanced median rho {stats.spearmanr(fbx, fby)[0]:.4f}, "
    f"same {int((np.sign(fbx) == np.sign(fby)).sum())}/43")

# ---------------------------------------------------------------- E. Hotspot S4I / S4J on 435
progs = sorted(hs_all.program_uid.unique())
log(f"E. Hotspot programs {len(progs)}")
H = {}
for d in EVAL:
    w = hs_all[hs_all.dataset == d].pivot(index="sample_id", columns="program_uid", values="score_raw")[progs]
    oz = hs_all[hs_all.dataset == d].pivot(index="sample_id", columns="program_uid", values="outcome_z")[progs]
    allu = ev.loc[ev.dataset == d, "sample_id"].values
    w, oz = w.loc[allu], oz.loc[allu]
    sd = w.std(ddof=1)
    z = (w - w.mean()) / sd
    diff = np.nanmax(np.abs(z.values - oz.values)) if np.isfinite(z.values).any() else np.nan
    log(f"E. {d}: zero-variance programs {list(sd.index[~(sd > 0)])}; max|my z - stored outcome_z| {diff:.2e}")
    H[d] = z
hb, hsb, hc = [], [], {a: ([], []) for a in AXES}
for d in EVAL:
    c = U[d]
    Y = H[d].loc[c.sample_id, progs].values
    b, s = ols_coef1(design_stage(c), Y); hb.append(b); hsb.append(s)
    for a in AXES:
        b, s = ols_coef1(design_cont(c, a), Y); hc[a][0].append(b); hc[a][1].append(s)
h_stage_beta, _, h_stage_p = fixed_meta(hb, hsb)
estimable = np.isfinite(h_stage_p)
h_stage_q = bh(h_stage_p, 117)
h_stage_sup = estimable & (np.sign(hb[0]) == np.sign(hb[1])) & (h_stage_q < 0.05)
hq, hsg = {}, []
for a in AXES:
    hq[a] = bh(fixed_meta(*hc[a])[2], 117)
    hsg += [np.sign(x) for x in hc[a][0]]
h_four = np.all(np.vstack(hsg) == hsg[0], axis=0)
h_cont = (hq["signature_pc1"] < 0.05) & (hq["fixed_projection"] < 0.05) & h_four
h_cont = np.where(np.isfinite(hq["fixed_projection"]), h_cont, False)
cls = np.where(~estimable, "not_estimable", np.where(h_stage_sup & h_cont, "Both",
               np.where(h_cont, "Continuum only", np.where(h_stage_sup, "Stage only", "Neither"))))
log(f"E. S4I continuum-associated {int(h_cont.sum())}/117 (estimable {int(estimable.sum())}); "
    f"S4J classes {pd.Series(cls).value_counts().to_dict()}; not estimable {list(np.array(progs)[~estimable])}")
inputs["s4j_source"] = S4J_SRC
s4j = pd.read_csv(S4J_SRC, sep="\t").set_index("program_uid").loc[progs]
mism = np.array(progs)[cls != s4j.support_class.fillna("not_estimable").astype(str).values]
log(f"E. S4J source classes {s4j.support_class.value_counts(dropna=False).to_dict()}; programs whose class differs: {list(mism)}")
pd.DataFrame(dict(program_uid=progs, stage_beta=h_stage_beta, stage_q=h_stage_q, stage_supported=h_stage_sup,
                  q_pc1=hq["signature_pc1"], q_fp=hq["fixed_projection"], four_dir=h_four, continuum=h_cont,
                  support_class=cls)).to_csv(f"{OUT}/E_hotspot_435.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- provenance
with open(f"{OUT}/inputs_sha256.tsv", "w") as fh:
    fh.write("role\tpath\tsha256\n")
    for k, v in inputs.items():
        fh.write(f"{k}\t{v}\t{sha256(v)}\n")
with open(f"{OUT}/run_parameters.json", "w") as fh:
    json.dump(dict(bootstrap_seeds=BOOT_SEEDS, B=B, bh_family_systems=43, bh_family_hotspot=117,
                   resampling_unit="genetic-identity cluster (one row per unit), within cohort",
                   unit_order="file order of fixed_projection rows in projection/participant_scores.tsv"), fh, indent=1)
with open(f"{OUT}/u10_systems.log", "w") as fh:
    fh.write("\n".join(log_lines) + "\n")
