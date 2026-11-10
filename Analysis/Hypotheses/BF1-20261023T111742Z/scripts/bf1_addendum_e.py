#!/usr/bin/env python
"""BF1 addendum E: matched-gene null for S4, S4 in Vu 2025 and GSE192741, matched null for S3 (ADDENDUM_BF1_E.md)."""
import argparse
import hashlib
import json
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--draws", type=int, default=500)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
CG0 = f"{ROOT}/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z"
R9D = (f"{ROOT}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/"
       "catalog_r9d-20261016T105059Z/build_out/catalog_r9d_gene_table.parquet")
SP = f"{ROOT}/Analysis/Spatial/results/preprocessed"
DMAP = f"{ROOT}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv"
BF1_OUT = f"{ROOT}/Analysis/Hypotheses/BF1-20261023T111742Z/out-20261023T112740Z/s4_s5_visium_arrays.tsv"
VU = f"{SP}/merged_spatial_vu_raw.h5ad"
JBO = f"{SP}/merged_spatial_raw.h5ad"
SEED = 20261023
H88_MARKERS = {"CYP2E1", "GLUL", "OAT", "CYP1A2", "CYP3A4", "CPS1", "HAL", "SDS", "ASS1", "PCK1",
               "KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"}
RES = {}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


with open(f"{a.out}/input_sha256.txt", "w") as f:
    for p in [f"{CG0}/genes.tsv", R9D, DMAP, BF1_OUT, VU, JBO]:
        f.write(f"{sha(p)}  {p}\n")

# ---------------------------------------------------------------- gene annotation (as in bf1_zonation.py)
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

sym = g.gene_name
called_any = set(sym[g.zon_protein_call_weiss2026.isin(["central", "portal"]) |
                     g.zon_mrna_call_yakubovsky2026.isin(["central", "portal"])].dropna())
Cp = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "central")].dropna())
Pp = set(sym[(~g.signature_gene) & (g.zon_protein_call_weiss2026 == "portal")].dropna())
U_arm = set(sym[g.signature_gene & (g.cont_z > 0)].dropna()) - called_any - H88_MARKERS
sig_sym = set(sym[g.signature_gene].dropna())
excluded = sig_sym | called_any | H88_MARKERS
nonhep = set(g.dropna(subset=["gene_name"]).groupby("gene_name").lineage_hep_specificity_log2cpm.max()
             .loc[lambda s: s <= 0].index)


def resid_cols(Y, cov):
    B, *_ = np.linalg.lstsq(cov, Y, rcond=None)
    return Y - cov @ B


def array_test(X, syms, tot_obs, ngen_obs, rng, draws):
    """X: spots x genes raw counts (csr). Returns observed r(zeta, U) and the two matched-null distributions."""
    tot = np.asarray(X.sum(1)).ravel()
    uu, inv = np.unique(syms, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(uu)))
    Xs = (X @ M).tocsc()
    det = np.asarray((Xs > 0).sum(0)).ravel() / Xs.shape[0]  # count / n, as numpy mean in bf1_zonation.py
    keep = det >= 0.01
    uu, det = uu[keep], det[keep]
    Cm = Xs[:, keep].toarray().astype(np.float32)
    Lg = np.log1p(Cm / tot[:, None] * 1e4)
    del Cm
    sd = Lg.std(0)
    Zs = (Lg - Lg.mean(0)) / np.where(sd > 0, sd, 1)
    del Lg
    cov = np.column_stack([np.ones(len(tot)), np.log1p(tot_obs), np.log1p(ngen_obs)])
    R = resid_cols(Zs.astype(np.float64), cov).astype(np.float32)  # resid of a mean = mean of resids
    del Zs
    idx = {x: i for i, x in enumerate(uu)}
    ic = [idx[x] for x in Cp if x in idx]; ip = [idx[x] for x in Pp if x in idx]
    iu = np.array([idx[x] for x in U_arm if x in idx])
    zeta = R[:, ic].mean(1) - R[:, ip].mean(1)
    zc = (zeta - zeta.mean()) / zeta.std()

    def corr(cols):
        s = R[:, cols].mean(1)
        return float(np.mean(zc * (s - s.mean()) / s.std()))

    r_obs = corr(iu)
    dec = pd.qcut(pd.Series(det).rank(method="first"), 10, labels=False).values
    out = {"n_spots": len(tot), "n_detected": len(uu), "n_central": len(ic), "n_portal": len(ip), "n_U": len(iu),
           "r_obs": r_obs}
    for name, pool_set in [("pool_i", None), ("pool_ii", nonhep)]:
        ok = np.array([(x not in excluded) and (pool_set is None or x in pool_set) for x in uu])
        need = pd.Series(dec[iu]).value_counts()
        pools = {d: np.where(ok & (dec == d))[0] for d in need.index}
        short = {int(d): int(len(pools[d])) for d in need.index if len(pools[d]) < need[d]}
        null = np.empty(draws)
        for b in range(draws):
            cols = np.concatenate([rng.choice(pools[d], need[d], replace=len(pools[d]) < need[d]) for d in need.index])
            null[b] = corr(cols)
        zn = np.arctanh(null)
        out[name] = dict(pool_size=int(ok.sum()), short_deciles=short, null_mean_r=float(null.mean()),
                         diff_fisher_z=float(np.arctanh(r_obs) - zn.mean()),
                         pct_in_null=float((null < r_obs).mean()), below_p05=bool(r_obs < np.percentile(null, 5)))
    return out


def summarize(rows, key, rng):
    d = np.array([r[key] for r in rows])
    flips = rng.choice([-1, 1], size=(100000, len(d))) @ d / len(d)
    return dict(n=len(d), mean=float(d.mean()), n_negative=int((d < 0).sum()),
                signflip_p=float((1 + np.sum(np.abs(flips) >= abs(d.mean()))) / 100001))


def flat(name, grp, o):
    row = {k: v for k, v in o.items() if not isinstance(v, dict)}
    row.update(array_id=name, group=grp, fisher_z_obs=float(np.arctanh(o["r_obs"])))
    for p in ["pool_i", "pool_ii"]:
        for k, v in o[p].items():
            row[f"{p}_{k}"] = json.dumps(v) if isinstance(v, dict) else v
    return row


rng = np.random.default_rng(SEED)
# ---------------------------------------------------------------- E1 HMSMA, 32 primary arrays
dmap = pd.read_csv(DMAP, sep="\t")
dmap = dmap[dmap.group == "primary"]
bf1 = pd.read_csv(BF1_OUT, sep="\t").set_index("array_id")
rows = []
for arr in dmap.array_id:
    A = ad.read_h5ad(f"{SP}/HRA007511_starsolo_hexscaled/{arr}.h5ad")
    X = A.X.tocsr() if sp.issparse(A.X) else sp.csr_matrix(A.X)
    o = array_test(X, A.var.gene_symbol.astype(str).values, A.obs.total_counts.values, A.obs.n_genes.values,
                   rng, a.draws)
    assert abs(o["r_obs"] - bf1.loc[arr, "r_zeta_U"]) < 1e-3, (arr, o["r_obs"], bf1.loc[arr, "r_zeta_U"])
    rows.append(flat(arr, "HMSMA_primary", o))
    print(arr, rows[-1], flush=True)
e1 = pd.DataFrame(rows)
for p in ["pool_i", "pool_ii"]:
    RES[f"E1_HMSMA_{p}"] = dict(summarize(rows, f"{p}_diff_fisher_z", rng),
                               n_below_null_p05=int(e1[f"{p}_below_p05"].sum()),
                               median_pct_in_null=float(e1[f"{p}_pct_in_null"].median()))

# ---------------------------------------------------------------- E2 Vu 2025 and GSE192741
e2 = []
for path, label in [(VU, "Vu2025"), (JBO, "GSE192741")]:
    A = ad.read_h5ad(path)
    Xall = A.layers["counts"]
    Xall = Xall.tocsr() if sp.issparse(Xall) else sp.csr_matrix(Xall)
    for s in A.obs.sample_id.cat.categories if hasattr(A.obs.sample_id, "cat") else A.obs.sample_id.unique():
        m = (A.obs.sample_id == s).values
        ob = A.obs[m]
        o = array_test(Xall[m], A.var_names.astype(str).values, ob.total_counts.values,
                       ob.n_genes_by_counts.values, rng, a.draws)
        grp = f"{label}_{ob.condition.iloc[0]}"
        e2.append(dict(flat(s, grp, o), individual=ob.individual.iloc[0]))
        print(s, e2[-1], flush=True)
    del A, Xall
e2 = pd.DataFrame(e2)
vu = e2[e2.group.str.startswith("Vu2025")].to_dict("records")
RES["E2_Vu_arrays_observed"] = summarize(vu, "fisher_z_obs", rng)
RES["E2_Vu_arrays_observed"]["replicates_rule_met"] = bool(RES["E2_Vu_arrays_observed"]["mean"] < 0 and
                                                           RES["E2_Vu_arrays_observed"]["signflip_p"] < 0.05)
RES["E2_Vu_arrays_pool_i"] = dict(summarize(vu, "pool_i_diff_fisher_z", rng),
                                  n_below_null_p05=int(sum(r["pool_i_below_p05"] for r in vu)))
vp = (pd.DataFrame(vu).assign(prefix=lambda d: d.array_id.str.split("_").str[0])
      .groupby("prefix")[["fisher_z_obs", "pool_i_diff_fisher_z"]].mean())
signs = np.array(np.meshgrid(*[[-1, 1]] * len(vp))).reshape(len(vp), -1).T
for k in ["fisher_z_obs", "pool_i_diff_fisher_z"]:
    x = vp[k].values
    RES[f"E2_Vu_prefix_{k}"] = dict(n=len(x), mean=float(x.mean()), n_negative=int((x < 0).sum()),
                                    exact_signflip_p=float(np.mean(np.abs(signs @ x / len(x)) >= abs(x.mean()) - 1e-12)))
RES["E2_GSE192741_sections"] = e2[e2.group.str.startswith("GSE192741")][
    ["array_id", "group", "individual", "n_U", "r_obs", "pool_i_diff_fisher_z", "pool_i_pct_in_null"]].to_dict("records")
pd.concat([e1, e2]).to_csv(f"{a.out}/e1_e2_arrays.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- E3 matched null for S3
cg = g.coloc_gen.fillna(False).astype(bool)
e3 = []
for layer, col, measured in [("mrna", "zon_mrna_call_yakubovsky2026", {"central", "portal", "not_zonated_detected"}),
                             ("protein", "zon_protein_call_weiss2026", {"central", "portal", "not_zonated"})]:
    gm = g[g[col].isin(measured)]
    co, nc = gm[cg.loc[gm.index]], gm[~cg.loc[gm.index]]
    need = co.stratum.value_counts()
    pools = {s: nc[nc.stratum == s][col].values for s in need.index}
    empty = [s for s in need.index if len(pools[s]) == 0]
    for s in empty:  # no non-colocalized gene in this stratum: fall back to the same expression decile
        pools[s] = nc[nc.expr_bin == co[co.stratum == s].expr_bin.iloc[0]][col].values
    for z in ["central", "portal"]:
        obs = int((co[col] == z).sum())
        null = np.array([sum(int((rng.choice(pools[s], need[s], replace=len(pools[s]) < need[s]) == z).sum())
                             for s in need.index) for _ in range(2000)])
        dev = np.abs(null - null.mean())
        e3.append(dict(layer=layer, zone=z, prespecified=(z == "central"), n_coloc_measured=len(co), observed=obs,
                       null_mean=float(null.mean()), null_lo=float(np.percentile(null, 2.5)),
                       null_hi=float(np.percentile(null, 97.5)),
                       ratio_obs_to_null=float(obs / null.mean()),
                       p_two_sided=float((1 + np.sum(dev >= abs(obs - null.mean()))) / 2001),
                       short_strata=int(sum(len(pools[s]) < need[s] for s in need.index)),
                       empty_strata_fallback=";".join(f"{s}:{need[s]}" for s in empty)))
e3 = pd.DataFrame(e3)
e3.to_csv(f"{a.out}/e3_coloc_zone_matched.tsv", sep="\t", index=False)
RES["E3"] = e3.to_dict("records")
json.dump(RES, open(f"{a.out}/results.json", "w"), indent=1, default=float)
print(json.dumps(RES, indent=1, default=float))
