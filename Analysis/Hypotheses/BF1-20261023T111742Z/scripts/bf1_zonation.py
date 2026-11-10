#!/usr/bin/env python
"""BF1: lobular zonation along the continuum (PRESPEC_BF1.md)."""
import argparse
import hashlib
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--draws", type=int, default=2000)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z"
LCPM = f"{ROOT}/Analysis/Hypotheses/CGA2-20261022T180846Z/logcpm-20261022T181752Z"
R9D = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105059Z/build_out/catalog_r9d_gene_table.parquet")
META = (f"{ROOT}/RNA-seq/results/histology_anchored_continuum/molecular_systems_inference/"
        "hac-system-inference-corrected774-20261009T171017Z/tables/system_continuum_meta_analysis.tsv")
SP = f"{ROOT}/Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled"
DMAP = f"{ROOT}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv"
GPC = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "gate/GAP-tissue/out/positive_controls_per_donor.tsv")
SEED = 20261023
H88_MARKERS = {"CYP2E1", "GLUL", "OAT", "CYP1A2", "CYP3A4", "CPS1", "HAL", "SDS", "ASS1", "PCK1",
               "KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"}
LAYERS = {"protein": ("zon_protein_call_weiss2026", {"central", "portal", "not_zonated"}),
          "mrna": ("zon_mrna_call_yakubovsky2026", {"central", "portal", "not_zonated_detected"})}
RES = {}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


with open(f"{a.out}/input_sha256.txt", "w") as f:
    for p in [f"{CG0}/genes.tsv", f"{CG0}/units.tsv", f"{CG0}/gene_system_weights.tsv.gz", R9D, META, DMAP, GPC,
              f"{LCPM}/GSE162694_logcpm.tsv.gz", f"{LCPM}/GSE213621_logcpm.tsv.gz"]:
        f.write(f"{sha(p)}  {p}\n")

# ---------------------------------------------------------------- gene annotation
g = pd.read_csv(f"{CG0}/genes.tsv", sep="\t")
g["signature_gene"] = g.signature_gene.fillna(False).astype(bool)
r9 = pd.read_parquet(R9D, columns=["ensembl_id", "zon_protein_call_weiss2026", "zon_mrna_call_yakubovsky2026",
                                   "lineage_hep_specificity_log2cpm"]).drop_duplicates("ensembl_id")
g = g.merge(r9, left_on="g", right_on="ensembl_id", how="left")
g["expr_bin"] = pd.qcut(g.AveExpr, 10, labels=False)
hs = g.lineage_hep_specificity_log2cpm
g["hep_bin"] = pd.qcut(hs, 3, labels=False).astype("float").fillna(3).astype(int)
g["stratum"] = g.expr_bin.astype(str) + "_" + g.hep_bin.astype(str)
g = g.set_index("g")

# ---------------------------------------------------------------- units
u = pd.read_csv(f"{CG0}/units.tsv", sep="\t")
u = u[u.dataset.isin(["GSE162694", "GSE213621"]) & u.fibrosis_stage.notna()]
u = u[~((u.dataset == "GSE213621") & (u.fibrosis_stage == 0))].copy()
n_units = u.dataset.value_counts().to_dict()
assert n_units == {"GSE213621": 258, "GSE162694": 109}, n_units
u["stage"] = u.fibrosis_stage.astype(int)
hep = u.bp_Hepatocytes.clip(1e-4, 1 - 1e-4)
u["logit_hep"] = np.log(hep / (1 - hep))


def layer_sets(layer, hep_only=False):
    col = LAYERS[layer][0]
    ok = ~g.signature_gene
    if hep_only:
        ok &= g.lineage_hep_specificity_log2cpm > 0
    C = g.index[ok & (g[col] == "central")]
    P = g.index[ok & (g[col] == "portal")]
    pool = g.index[~g.signature_gene & ~g[col].isin(["central", "portal"])]
    return list(C), list(P), list(pool)


def design(d, comp=False, stage=True):
    f = "continuum_z + C(inferred_sex)" + (" + C(stage)" if stage else "") + (" + logit_hep" if comp else "")
    import patsy
    X = patsy.dmatrix(f, d, return_type="dataframe")
    j = list(X.columns).index("continuum_z")
    H = np.linalg.pinv(X.values)[j]
    return f, H


def draw_sets(rng, genes, pools_by_stratum):
    """One matched gene per input gene, same stratum, without replacement within the call."""
    st = g.loc[genes, "stratum"].values
    out = np.empty(len(genes), dtype=object)
    for s in np.unique(st):
        pos = np.where(st == s)[0]
        pool = pools_by_stratum.get(s, np.array([]))
        if len(pool) == 0:
            pool = pools_by_stratum["__all__"]
        out[pos] = rng.choice(pool, size=len(pos), replace=len(pool) < len(pos))
    return list(out)


def draw_pair(rng, C, P, pst):
    r = draw_sets(rng, C + P, pst)
    return r[:len(C)], r[len(C):]


rows, nullrows = [], []
for coh in ["GSE162694", "GSE213621"]:
    d = u[u.dataset == coh].copy()
    L = pd.read_csv(f"{LCPM}/{coh}_logcpm.tsv.gz", sep="\t", index_col=0)
    L = L[d.sample_id]
    Z = L.sub(L.mean(1), axis=0).div(L.std(1).replace(0, 1), axis=0)
    Zt = Z.T.values  # units x genes
    col_of = {x: i for i, x in enumerate(Z.index)}
    for layer in LAYERS:
        for variant in ["primary", "composition", "hep_specific", "no_stage"]:
            Cg, Pg, pool = layer_sets(layer, hep_only=(variant == "hep_specific"))
            Cg = [x for x in Cg if x in col_of]; Pg = [x for x in Pg if x in col_of]
            pool = [x for x in pool if x in col_of]
            f, H = design(d, comp=(variant == "composition"), stage=(variant != "no_stage"))
            sc = lambda gs: Zt[:, [col_of[x] for x in gs]].mean(1)
            zc, zp = sc(Cg), sc(Pg)
            obs = dict(C=H @ zc, P=H @ zp, D=H @ (zc - zp))
            d2 = d.assign(zc=zc, zp=zp, zd=zc - zp)
            fits = {k: smf.ols(f"{v} ~ {f}", d2).fit(cov_type="HC3") for k, v in [("C", "zc"), ("P", "zp"), ("D", "zd")]}
            if variant in ("primary", "composition", "hep_specific", "no_stage"):
                pst = pd.Series(pool).groupby(g.loc[pool, "stratum"].values).apply(np.array).to_dict()
                pst["__all__"] = np.array(pool)
                rng = np.random.default_rng(SEED)
                B = a.draws if variant == "primary" else 500
                nd = np.empty((B, 3))
                for b in range(B):
                    rc, rp = draw_pair(rng, Cg, Pg, pst)
                    nc, npp = sc(rc), sc(rp)
                    nd[b] = [H @ nc, H @ npp, H @ (nc - npp)]
                for j, k in enumerate(["C", "P", "D"]):
                    nv = nd[:, j]; mu, sd = nv.mean(), nv.std(ddof=1)
                    p = (1 + np.sum(np.abs(nv - mu) >= abs(obs[k] - mu))) / (B + 1)
                    rows.append(dict(cohort=coh, layer=layer, variant=variant, score=k, n_C=len(Cg), n_P=len(Pg),
                                     n_units=len(d), beta=obs[k], se_hc3=fits[k].bse["continuum_z"],
                                     p_param=fits[k].pvalues["continuum_z"], null_mean=mu, null_sd=sd,
                                     z_vs_null=(obs[k] - mu) / sd, p_null=p, draws=B))
                if variant == "primary":
                    nullrows.append(pd.DataFrame(nd, columns=["C", "P", "D"]).assign(cohort=coh, layer=layer))
            print(coh, layer, variant, {k: round(v, 4) for k, v in obs.items()}, flush=True)
p1 = pd.DataFrame(rows)
p1.to_csv(f"{a.out}/p1_zone_scores.tsv", sep="\t", index=False)
pd.concat(nullrows).to_csv(f"{a.out}/p1_null_draws.tsv.gz", sep="\t", index=False)
prim = p1[(p1.layer == "protein") & (p1.variant == "primary") & (p1.score == "D")].set_index("cohort")
RES["P1_decision_zonally_biased"] = bool((prim.p_null < 0.05).all() and np.sign(prim.beta).nunique() == 1)
RES["P1_protein_D"] = prim[["beta", "se_hc3", "z_vs_null", "p_null", "n_C", "n_P", "n_units"]].to_dict(orient="index")

# ---------------------------------------------------------------- S1 gene level
s1 = []
for layer in LAYERS:
    Cg, Pg, pool = layer_sets(layer)
    pst = pd.Series(pool).groupby(g.loc[pool, "stratum"].values).apply(np.array).to_dict(); pst["__all__"] = np.array(pool)
    rng = np.random.default_rng(SEED)
    cz = g.cont_z
    for nm, gs in [("central", Cg), ("portal", Pg)]:
        obs = cz.loc[gs].mean()
        nv = np.array([cz.loc[draw_sets(rng, gs, pst)].mean() for _ in range(a.draws)])
        s1.append(dict(layer=layer, zone=nm, n=len(gs), mean_cont_z=obs, null_mean=nv.mean(), null_sd=nv.std(ddof=1),
                       p_null=(1 + np.sum(np.abs(nv - nv.mean()) >= abs(obs - nv.mean()))) / (a.draws + 1),
                       frac_negative=(cz.loc[gs] < 0).mean()))
pd.DataFrame(s1).to_csv(f"{a.out}/s1_gene_level.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- S2 system zone map
w = pd.read_csv(f"{CG0}/gene_system_weights.tsv.gz", sep="\t")
meta = pd.read_csv(META, sep="\t")
meta = meta[(meta.aggregation_method == "collection_balanced") & (meta.axis_id == "signature_pc1")].set_index("community_id")
s2 = []
for layer, (col, measured) in LAYERS.items():
    gm = g[g[col].isin(measured)]
    univ = set(gm.index) & set(w.gene_id_base)
    for s, dd in w.groupby("community_id"):
        mem = set(dd.gene_id_base) & univ
        r = dict(layer=layer, system=s, n_measured=len(mem))
        for z in ["central", "portal"]:
            inz = set(gm.index[gm[col] == z]) & univ
            a11 = len(mem & inz); a12 = len(mem - inz); a21 = len(inz - mem); a22 = len(univ - mem - inz)
            orr, p = stats.fisher_exact([[a11, a12], [a21, a22]])
            r[f"n_{z}"] = a11
            r[f"log_or_{z}"] = np.log((a11 + .5) * (a22 + .5) / ((a12 + .5) * (a21 + .5)))
            r[f"p_{z}"] = p
        s2.append(r)
s2 = pd.DataFrame(s2)
for layer in LAYERS:
    for z in ["central", "portal"]:
        m = s2.layer == layer
        s2.loc[m, f"bh_{z}"] = multipletests(s2.loc[m, f"p_{z}"], method="fdr_bh")[1]
s2["zone_bias"] = s2.log_or_central - s2.log_or_portal
s2["continuum_random_beta"] = s2.system.map(meta.random_beta)
s2["continuum_random_q"] = s2.system.map(meta.random_q_value)
s2.to_csv(f"{a.out}/s2_system_zone_map.tsv", sep="\t", index=False)
RES["S2_spearman_zone_bias_vs_continuum_beta"] = {
    layer: dict(zip(["rho", "p"], stats.spearmanr(s2[s2.layer == layer].zone_bias, s2[s2.layer == layer].continuum_random_beta)))
    for layer in LAYERS}

# ---------------------------------------------------------------- S3 genetic risk genes
s3 = []
cg = g.coloc_gen.fillna(False).astype(bool)
for layer, (col, measured) in LAYERS.items():
    gm = g[g[col].isin(measured)]
    cc = cg.loc[gm.index]
    for z in ["central", "portal"]:
        t = pd.crosstab(cc, gm[col] == z).reindex(index=[True, False], columns=[True, False]).fillna(0)
        orr, p = stats.fisher_exact(t.values)
        s3.append(dict(layer=layer, zone=z, coloc_in_zone=int(t.loc[True, True]), coloc_measured=int(t.loc[True].sum()),
                       other_in_zone=int(t.loc[False, True]), other_measured=int(t.loc[False].sum()), odds_ratio=orr, p=p))
pd.DataFrame(s3).to_csv(f"{a.out}/s3_coloc_zone.tsv", sep="\t", index=False)
named = ["PNPLA3", "TM6SF2", "GCKR", "MTARC1", "HSD17B13", "APOE", "GPAM", "MBOAT7"]
g.reset_index()[g.reset_index().gene_name.isin(named)][
    ["gene_name", "g", "cont_z", "coloc_gen", "zon_protein_call_weiss2026", "zon_mrna_call_yakubovsky2026",
     "lineage_hep_specificity_log2cpm"]].to_csv(f"{a.out}/s3_named_risk_genes.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- S4/S5 Visium
sym = g.gene_name
called_any = set(sym[g.zon_protein_call_weiss2026.isin(["central", "portal"]) |
                     g.zon_mrna_call_yakubovsky2026.isin(["central", "portal"])].dropna())
Cp = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "central")].dropna())
Pp = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "portal")].dropna())
U_arm = set(sym[g.signature_gene & (g.cont_z > 0)].dropna()) - called_any - H88_MARKERS
F_arm = set(sym[g.signature_gene & (g.cont_z < 0)].dropna()) - called_any - H88_MARKERS
dmap = pd.read_csv(DMAP, sep="\t")
need = Cp | Pp | U_arm | F_arm


def resid(s, cov):
    b, *_ = np.linalg.lstsq(cov, s, rcond=None)
    return s - cov @ b


vrows, pb = [], {}
for arr, grp in zip(dmap.array_id, dmap.group):
    A = ad.read_h5ad(f"{SP}/{arr}.h5ad")
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    tot = np.asarray(X.sum(1)).ravel()
    s_ = A.var.gene_symbol.astype(str).values
    keep = np.isin(s_, list(need))
    uu, inv = np.unique(s_[keep], return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uu)))
    Cm = (X[:, keep] @ M).toarray().astype(np.float32)
    pb[arr] = pd.Series(Cm.sum(0), index=uu)
    pb[arr]["__total__"] = tot.sum()
    det = (Cm > 0).mean(0)
    Lg = np.log1p(Cm / tot[:, None] * 1e4)
    Zs = (Lg - Lg.mean(0)) / np.where(Lg.std(0) > 0, Lg.std(0), 1)
    idx = {x: i for i, x in enumerate(uu)}
    cov = np.column_stack([np.ones(len(tot)), np.log1p(A.obs.total_counts.values), np.log1p(A.obs.n_genes.values)])
    sc = lambda gs: [x for x in gs if x in idx and det[idx[x]] >= 0.01]
    gc, gp, gu, gf = sc(Cp), sc(Pp), sc(U_arm), sc(F_arm)
    c_s = resid(Zs[:, [idx[x] for x in gc]].mean(1), cov); p_s = resid(Zs[:, [idx[x] for x in gp]].mean(1), cov)
    u_s = resid(Zs[:, [idx[x] for x in gu]].mean(1), cov); f_s = resid(Zs[:, [idx[x] for x in gf]].mean(1), cov)
    zeta = c_s - p_s
    vrows.append(dict(array_id=arr, group=grp, n_spots=len(tot), median_total=float(np.median(tot)),
                      n_central=len(gc), n_portal=len(gp), n_U=len(gu), n_F=len(gf),
                      r_zeta_U=np.corrcoef(zeta, u_s)[0, 1], r_zeta_F=np.corrcoef(zeta, f_s)[0, 1],
                      zonation_strength=-np.corrcoef(c_s, p_s)[0, 1]))
    print(arr, vrows[-1], flush=True)
v = pd.DataFrame(vrows)
# S5 pseudobulk continuum position
pbm = pd.DataFrame(pb).fillna(0)
tot = pbm.loc["__total__"]; pbm = pbm.drop("__total__")
lc = np.log2(pbm / tot * 1e6 + 1)
prim_ids = list(v.array_id[v.group == "primary"])
mu_p, sd_p = lc[prim_ids].mean(1), lc[prim_ids].std(1).replace(0, 1)
zz = lc.sub(mu_p, axis=0).div(sd_p, axis=0)  # all arrays scaled to the 32 primary arrays
ug = [x for x in U_arm if x in zz.index and pbm.loc[x, prim_ids].sum() > 0]
fg = [x for x in F_arm if x in zz.index and pbm.loc[x, prim_ids].sum() > 0]
pos = zz.loc[ug].mean() - zz.loc[fg].mean()
v["pseudobulk_position"] = v.array_id.map(pos)
gpc = pd.read_csv(GPC, sep="\t").set_index("array_id")
v["five_marker_strength"] = -v.array_id.map(gpc.pc1_pericentral_vs_periportal_r)
v.to_csv(f"{a.out}/s4_s5_visium_arrays.tsv", sep="\t", index=False)

rng = np.random.default_rng(SEED)
vp = v[v.group == "primary"].copy()
for k in ["r_zeta_U", "r_zeta_F"]:
    zf = np.arctanh(vp[k].values)
    flips = rng.choice([-1, 1], size=(100000, len(zf))) @ zf / len(zf)
    boots = [zf[rng.integers(0, len(zf), len(zf))].mean() for _ in range(10000)]
    RES[f"S4_{k}"] = dict(n_arrays=len(zf), mean_fisher_z=zf.mean(), n_positive=int((zf > 0).sum()),
                         signflip_p=(1 + np.sum(np.abs(flips) >= abs(zf.mean()))) / 100001,
                         boot_ci=list(np.percentile(boots, [2.5, 97.5])))
for k in ["zonation_strength", "five_marker_strength"]:
    rho = stats.spearmanr(vp.pseudobulk_position, vp[k])[0]
    perm = np.array([stats.spearmanr(rng.permutation(vp.pseudobulk_position.values), vp[k])[0] for _ in range(10000)])
    rk = vp[["pseudobulk_position", k]].rank()
    rk["lm"] = np.log(vp.median_total).rank()
    res_x = resid(rk.pseudobulk_position.values, np.column_stack([np.ones(len(rk)), rk.lm]))
    res_y = resid(rk[k].values, np.column_stack([np.ones(len(rk)), rk.lm]))
    RES[f"S5_{k}"] = dict(n_arrays=len(vp), spearman=rho, perm_p=(1 + np.sum(np.abs(perm) >= abs(rho))) / 10001,
                         partial_spearman_given_log_median_total=np.corrcoef(res_x, res_y)[0, 1])
RES["S5_controls"] = v[v.group != "primary"][["array_id", "pseudobulk_position", "zonation_strength"]].to_dict(orient="records")
json.dump(RES, open(f"{a.out}/results.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
