# H88 audit: independent recount of every n stated in PRESPEC_H88.md from metadata only.
# Reads edge metadata, the gate's measurability/PC tables and program membership. No spatial statistic.
import pandas as pd, numpy as np, itertools, json
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
T4 = f"{R}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
G = f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/GAP-tissue/out"
MEM = f"{R}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
out = {}
dm = pd.read_csv(f"{R}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv", sep="\t")
out["n_arrays"] = len(dm); out["n_individuals_unique"] = int(dm.individual.nunique())
out["groups"] = dm.group.value_counts().to_dict()
out["sex_by_group"] = {g: d.sex.value_counts().to_dict() for g, d in dm.groupby("group")}
out["inferred_control_individuals"] = sorted(dm[dm.group == "inferred_control"].individual.astype(str).tolist())
prim = dm[dm.group == "primary"].array_id.tolist()
comm = pd.read_csv(f"{T4}/program_communities.tsv", sep="\t"); feats = set(comm.feature_id)
out["n_T4_features"] = len(feats)
mem = pd.read_csv(MEM, sep="\t")
out["n_programs_in_membership"] = int(mem.program_uid.nunique())
# measurability recomputed from the gate table, by explicit column
m = pd.read_csv(f"{G}/program_measurability.tsv", sep="\t")
mp = m[m.array_id.isin(prim)]
out["meas_arrays_primary"] = int(mp.array_id.nunique())
allok = mp.assign(ok=mp.weight_retained >= 0.80).groupby("program").ok.agg(["all", "size"])
ok = sorted(allok[(allok["all"]) & (allok["size"] == 32)].index)
out["measurable_all32"] = len(ok)
out["measurable_flag_agrees_with_weight_rule"] = bool(((mp.weight_retained >= 0.80) == mp.measurable_080).all())
out["min_n_measured_among_ok"] = int(mp[mp.program.isin(ok)].n_measured.min())
# lineage from the edge table cell_type columns, not from the program name
e = pd.read_csv(f"{T4}/program_interaction_edges.tsv.gz", sep="\t")
ct = pd.concat([e[["program_a", "cell_type_a"]].set_axis(["p", "ct"], axis=1), e[["program_b", "cell_type_b"]].set_axis(["p", "ct"], axis=1)]).drop_duplicates()
assert ct.p.is_unique
ctm = dict(zip(ct.p, ct.ct))
out["lineage_counts_ok"] = pd.Series([ctm[p] for p in ok]).value_counts().to_dict()
out["edge_rows"] = len(e); out["C113_2"] = 113 * 112 // 2
c = pd.read_csv(f"{T4}/coordinated_edges.tsv", sep="\t")
out["coord_rows"] = len(c); out["coord_support_values"] = c.support.value_counts().to_dict()
key = lambda a, b: tuple(sorted((a, b)))
ck = {key(a, b) for a, b in zip(c.program_a, c.program_b)}
e["k"] = [key(a, b) for a, b in zip(e.program_a, e.program_b)]
out["edge_keys_unique"] = bool(e.k.is_unique)
out["coord_in_edge_table"] = int(e.k.isin(ck).sum())
e["coord"] = e.k.isin(ck)
out["edge_support_by_coord"] = pd.crosstab(e.support, e.coord).to_dict()
out["cross_lineage_coord_all113"] = int((e.coord & ~e.same_cell_type).sum())
out["coord_rows_cross_lineage_from_coord_file"] = int((~c.same_cell_type).sum())
# gate pool over 113 programs: cross-lineage, n_shared==0, not coordinated
out["gate_pool_113"] = int((~e.same_cell_type & (e.n_shared == 0) & ~e.coord).sum())
oks = set(ok)
x = e[e.program_a.isin(oks) & e.program_b.isin(oks)]
out["pairs_among_44"] = len(x); out["C44_2"] = 44 * 43 // 2
x = x[~x.same_cell_type]
out["cross_lineage_among_44"] = len(x)
out["cross_lineage_coord_44"] = int(x.coord.sum())
out["n_shared_gt0_44"] = int((x.n_shared > 0).sum()); out["entangled_44"] = int(x.membership_entangled.sum())
out["entangled_subset_of_shared"] = bool((x[x.membership_entangled].n_shared > 0).all())
out["coord_sharing_genes_44"] = int((x.coord & (x.n_shared > 0)).sum())
out["gate_pool_44"] = int((~x.coord & (x.n_shared == 0)).sum())
P = x[(x.n_shared == 0) & ~x.membership_entangled].copy()
out["P"] = [len(P), int(P.coord.sum()), int((~P.coord).sum())]
out["P_coord_rdisjoint_neg"] = int((P.coord & (P.r_disjoint < 0)).sum())
out["P_coord_rdisjoint_zero"] = int((P.coord & (P.r_disjoint == 0)).sum())
out["P_coord_rpartial_le0"] = int((P.coord & (P.r_partial <= 0)).sum())
Pp = P[P.r_disjoint > 0]
out["Pplus"] = [len(Pp), int(Pp.coord.sum()), int((~Pp.coord).sum())]
F2c = x[~x.membership_entangled]
out["F2c"] = [len(F2c), int(F2c.coord.sum()), int((~F2c.coord).sum())]
lab = {"cholangiocytes": "cholangiocyte", "fibroblasts": "fibroblast", "hepatocytes": "hepatocyte", "macrophages": "macrophage", "tcells": "T cell"}
P["lp"] = ["-".join(sorted([lab[a], lab[b]], key=lambda s: s.lower())) for a, b in zip(P.cell_type_a, P.cell_type_b)]
Pp = P[P.r_disjoint > 0]
st = pd.DataFrame({"P_coord": P[P.coord].lp.value_counts(), "P_non": P[~P.coord].lp.value_counts(),
                   "Pp_coord": Pp[Pp.coord].lp.value_counts(), "Pp_non": Pp[~Pp.coord].lp.value_counts()}).fillna(0).astype(int)
print(st.to_string())
out["strata_P_both"] = int(((st.P_coord > 0) & (st.P_non > 0)).sum()); out["strata_Pp_both"] = int(((st.Pp_coord > 0) & (st.Pp_non > 0)).sum())
out["G2_every_stratum_has_enough_noncoord"] = bool((st.P_non >= st.P_coord).all())
out["PC2_pool_size"] = int(st.P_non.sum() - st.P_coord.sum())
# hub structure of coordinated labels in P (relevant to the label-permutation null)
deg = pd.concat([P[P.coord].program_a, P[P.coord].program_b]).value_counts().reindex(ok).fillna(0).astype(int)
degall = pd.concat([P.program_a, P.program_b]).value_counts().reindex(ok).fillna(0).astype(int)
frac = (deg / degall)
out["coord_degree_in_P"] = dict(min=int(deg.min()), median=float(deg.median()), max=int(deg.max()), n_zero=int((deg == 0).sum()),
    top5_share_of_endpoints=float(deg.sort_values(ascending=False).head(5).sum() / (2 * deg.sum() / 2) / 2))
out["coord_fraction_per_program_in_P"] = dict(min=float(frac.min()), q25=float(frac.quantile(.25)), median=float(frac.median()), q75=float(frac.quantile(.75)), max=float(frac.max()))
# expected spread of per-program coordinated fraction if labels were random within strata: binomial approx
pbar = P.coord.mean()
exp_sd = np.sqrt(pbar * (1 - pbar) / degall).median()
out["coord_fraction_sd_observed"] = float(frac.std()); out["coord_fraction_sd_expected_random_median"] = float(exp_sd)
# chi-square-like dispersion of per-program coordinated counts vs random-label expectation (approx; programs overlap)
exp = degall * pbar
out["degree_dispersion_ratio"] = float((((deg - exp) ** 2) / (exp * (1 - pbar))).mean())
# permutation check of hub structure: shuffle labels within strata 2000x, compare SD of per-program coord degree
rng = np.random.default_rng(20260923)
lp = P.lp.values; co = P.coord.values.copy(); pa = P.program_a.values; pb = P.program_b.values
idx = {p: i for i, p in enumerate(ok)}
ia = np.array([idx[p] for p in pa]); ib = np.array([idx[p] for p in pb])
def degsd(lbl):
    d = np.bincount(ia[lbl], minlength=44) + np.bincount(ib[lbl], minlength=44)
    return (d / degall.values).std()
obs = degsd(co)
nulls = []
groups = [np.where(lp == s)[0] for s in np.unique(lp)]
for _ in range(2000):
    l2 = co.copy()
    for g in groups: l2[g] = rng.permutation(co[g])
    nulls.append(degsd(l2))
nulls = np.array(nulls)
out["hub_check_sd_of_program_coord_fraction"] = dict(observed=float(obs), null_mean=float(nulls.mean()), null_q95=float(np.quantile(nulls, .95)), p=float((1 + (nulls >= obs).sum()) / 2001))
# PC3 and PC1 genes: are they in the gate engine's `need` set (program members or the PC1 list)?
PC3 = ["KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]
memT4 = mem[mem.program_uid.isin(feats)]
out["PC3_genes_in_T4_membership_symbols"] = {g: bool(g in set(memT4.mapped_symbol)) for g in PC3}
# PC1 gate values
pc = pd.read_csv(f"{G}/positive_controls_per_donor.tsv", sep="\t")
pcp = pc[pc.array_id.isin(prim)]
z = np.arctanh(pcp.pc1_pericentral_vs_periportal_r)
out["PC1"] = dict(n=int(len(z)), mean=float(z.mean()), sd=float(z.std(ddof=1)), n_neg=int((z < 0).sum()))
out["spots_primary"] = int(pcp.n_spots.sum()); out["spots_min"] = int(pcp.n_spots.min()); out["spots_min_array"] = pcp.loc[pcp.n_spots.idxmin(), "array_id"]
z2 = np.arctanh(pcp.pc2_ductular_vs_tgfb_collagen_r.dropna())
out["gate_anchor_pair"] = dict(n=int(len(z2)), mean=float(z2.mean()))
nl = pd.read_csv(f"{G}/null_pair_scale_per_donor.tsv", sep="\t")
out["gate_null_pair_sd_median_all35"] = float(nl.null_sd_z.median()); out["gate_null_pair_mean_primary"] = float(nl[nl.array_id.isin(prim)].null_mean_z.mean())
print(json.dumps(out, indent=1, default=str))
