# GAP-tissue gate counting: HMSMA Visium donor counts, program measurability, positive controls, null-pair scale.
# Does NOT compute any statistic on frozen coordinated edges (only checks their measurability).
import os, sys, json, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
rng = np.random.default_rng(20260923)
ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = sys.argv[1]
SP = f"{ROOT}/Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled"
T4 = f"{ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
MEM = f"{ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
dmap = pd.read_csv(f"{ROOT}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv", sep="\t")
comm = pd.read_csv(f"{T4}/program_communities.tsv", sep="\t")
feats = set(comm.feature_id)
mem = pd.read_csv(MEM, sep="\t"); mem = mem[mem.program_uid.isin(feats)]
W = {p: g.groupby("mapped_symbol").original_l1_weight.sum() for p, g in mem.groupby("program_uid")}
coord = pd.read_csv(f"{T4}/coordinated_edges.tsv", sep="\t")
alle = pd.read_csv(f"{T4}/program_interaction_edges.tsv.gz", sep="\t", usecols=["program_a","program_b","n_shared","same_cell_type"])
ck = set(map(frozenset, zip(coord.program_a, coord.program_b)))
pool = alle[(~alle.same_cell_type) & (alle.n_shared == 0)]
pool = pool[[frozenset(x) not in ck for x in zip(pool.program_a, pool.program_b)]]
pool = pool[pool.program_a.isin(feats) & pool.program_b.isin(feats)]
nullpairs = pool.sample(n=min(300, len(pool)), random_state=20260923)
DUCT = "hotspot_hepatocytes_48f39dd4d817a10e"; COLL = "hotspot_fibroblasts_98c7913eef5f9315"
pc2_in_coord = frozenset([DUCT, COLL]) in ck
PERIC = ["CYP2E1","GLUL","OAT","CYP1A2","CYP3A4"]; PERIP = ["CPS1","HAL","SDS","ASS1","PCK1"]
need = set(mem.mapped_symbol) | set(PERIC) | set(PERIP)
meas_rows, pc_rows, null_rows = [], [], []
def resid(s, cov):
    b, *_ = np.linalg.lstsq(cov, s, rcond=None); r = s - cov @ b; return r
for arr in sorted(dmap.array_id):
    a = ad.read_h5ad(f"{SP}/{arr}.h5ad")
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    tot = np.asarray(X.sum(1)).ravel()
    sym = a.var.gene_symbol.astype(str).values
    keep = np.isin(sym, list(need))
    Xs = X[:, keep]; s2 = sym[keep]
    df = pd.DataFrame.sparse.from_spmatrix(Xs.T.tocsr() if False else Xs, columns=s2) if False else None
    # collapse duplicate symbols
    u, inv = np.unique(s2, return_inverse=True)
    M = sp.csr_matrix((np.ones(len(inv)), (np.arange(len(inv)), inv)), shape=(len(inv), len(u)))
    C = (Xs @ M).toarray().astype(np.float32)
    det = (C > 0).mean(0)
    L = np.log1p(C / tot[:, None] * 1e4)
    Z = (L - L.mean(0)) / np.where(L.std(0) > 0, L.std(0), 1)
    idx = {g: i for i, g in enumerate(u)}
    cov = np.column_stack([np.ones(len(tot)), np.log1p(a.obs.total_counts.values), np.log1p(a.obs.n_genes.values)])
    scores = {}
    for p, w in W.items():
        g = [x for x in w.index if x in idx and det[idx[x]] >= 0.01]
        ret = w[g].sum() / w.sum() if len(g) else 0.0
        meas_rows.append(dict(array_id=arr, program=p, n_genes=len(w), n_measured=len(g), weight_retained=ret, measurable_080=ret >= 0.80))
        if len(g) >= 8:
            ww = w[g].values / w[g].sum()
            s = Z[:, [idx[x] for x in g]] @ ww
            scores[p] = resid(s, cov)
    def mk(gl):
        g = [x for x in gl if x in idx and det[idx[x]] >= 0.01]
        return resid(Z[:, [idx[x] for x in g]].mean(1), cov), len(g)
    pc_s, npc = mk(PERIC); pp_s, npp = mk(PERIP)
    r1 = np.corrcoef(pc_s, pp_s)[0, 1]
    r2 = np.corrcoef(scores[DUCT], scores[COLL])[0, 1] if DUCT in scores and COLL in scores else np.nan
    pc_rows.append(dict(array_id=arr, n_spots=len(tot), pc1_pericentral_vs_periportal_r=r1, pc1_n_genes=f"{npc}/{npp}", pc2_ductular_vs_tgfb_collagen_r=r2))
    zs = [np.arctanh(np.corrcoef(scores[x], scores[y])[0, 1]) for x, y in zip(nullpairs.program_a, nullpairs.program_b) if x in scores and y in scores]
    null_rows.append(dict(array_id=arr, n_null_pairs=len(zs), null_mean_z=np.mean(zs), null_sd_z=np.std(zs)))
    print(arr, len(tot), round(r1, 3), round(r2, 3), flush=True)
meas = pd.DataFrame(meas_rows); meas.to_csv(f"{OUT}/program_measurability.tsv", sep="\t", index=False)
pc = pd.DataFrame(pc_rows).merge(dmap, on="array_id"); pc.to_csv(f"{OUT}/positive_controls_per_donor.tsv", sep="\t", index=False)
nl = pd.DataFrame(null_rows).merge(dmap, on="array_id"); nl.to_csv(f"{OUT}/null_pair_scale_per_donor.tsv", sep="\t", index=False)
prim = dmap[dmap.group == "primary"].array_id
mp = meas[meas.array_id.isin(prim)].groupby("program").measurable_080.all()
ok = set(mp[mp].index)
x = coord[~coord.same_cell_type]
summ = dict(n_arrays=len(dmap), n_individuals=dmap.individual.nunique(), group_counts=dmap.group.value_counts().to_dict(),
  n_spots_primary=int(pc[pc.array_id.isin(prim)].n_spots.sum()), spots_per_donor_min=int(pc[pc.array_id.isin(prim)].n_spots.min()),
  n_T4_programs=len(feats), programs_measurable_all32=len(ok),
  cross_lineage_coord_edges=len(x), cross_lineage_coord_edges_both_measurable=int(sum(a in ok and b in ok for a, b in zip(x.program_a, x.program_b))),
  null_pool_cross_lineage_disjoint=len(pool), pc2_pair_is_coordinated_edge=bool(pc2_in_coord))
for k, col in [("pc1", "pc1_pericentral_vs_periportal_r"), ("pc2", "pc2_ductular_vs_tgfb_collagen_r")]:
    v = np.arctanh(pc[pc.array_id.isin(prim)][col].dropna().values)
    summ[k] = dict(n=len(v), mean_z=float(v.mean()), sd_z=float(v.std(ddof=1)), dz=float(v.mean() / v.std(ddof=1)), n_neg=int((v < 0).sum()), n_pos=int((v > 0).sum()))
nv = nl[nl.array_id.isin(prim)].null_mean_z
summ["null_pair_mean_z_across_donors"] = dict(mean=float(nv.mean()), sd=float(nv.std(ddof=1)))
summ["within_donor_null_pair_sd_z_median"] = float(nl.null_sd_z.median())
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1, default=str); print(json.dumps(summ, indent=1, default=str))
