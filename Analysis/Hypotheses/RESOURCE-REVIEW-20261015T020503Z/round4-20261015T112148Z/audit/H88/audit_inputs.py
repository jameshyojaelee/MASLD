# H88 audit: recompute sha256 and bytes of every input named in PRESPEC_H88.md, check array headers,
# and check presence and detection (>=1% of spots) of the PC1 and PC3 marker genes in all 35 arrays.
# Computes no correlation, no program-pair statistic and no PC3 value.
import re, os, sys, hashlib, json
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
R = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = sys.argv[1]
PRE = f"{R}/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/audit/H88/PRESPEC_H88.before_audit.md"
txt = open(PRE).read()
T4 = "RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/systems/T4_interaction"
GATE = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/gate/GAP-tissue"
W = "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/write/H88"
meta = [f"{T4}/coordinated_edges.tsv", f"{T4}/program_interaction_edges.tsv.gz", f"{T4}/program_communities.tsv",
        "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv",
        "Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv", f"{GATE}/count_hmsma.py",
        f"{GATE}/out/program_measurability.tsv", f"{GATE}/out/summary.json", f"{GATE}/out/positive_controls_per_donor.tsv",
        f"{GATE}/out/null_pair_scale_per_donor.tsv", f"{W}/count_h88.py", f"{W}/count_h88.out",
        "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/plans/JOINT_RULE_H52_H54.md"]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
rows = []
for p in meta:
    bn = os.path.basename(p)
    m = re.search(r"\|[^|\n]*\|\s*`[^`]*" + re.escape(bn) + r"`\s*\|\s*([\d,]*)\s*\|\s*`([0-9a-f]{64})`", txt)
    exp_b = m.group(1).replace(",", "") if m else None; exp_h = m.group(2) if m else None
    full = f"{R}/{p}"; b = os.path.getsize(full); h = sha(full)
    rows.append(dict(file=p, bytes=b, sha256=h, prespec_bytes=exp_b or "", prespec_sha256=exp_h or "NOT_FOUND",
                     bytes_match=(exp_b in ("", None) or int(exp_b) == b), sha_match=(h == exp_h)))
SP = f"{R}/Analysis/Spatial/results/preprocessed/HRA007511_starsolo_hexscaled"
exp_arr = {a: (int(b.replace(",", "")), h) for a, b, h in re.findall(r"(HRA_\d+)\.h5ad\s*\|?\s*([\d,]+)\s*\|?\s*`([0-9a-f]{64})`", txt)}
dm = pd.read_csv(f"{R}/Analysis/Hypotheses/V4-20261008T015029Z/hmsma_gsa_donor_map.tsv", sep="\t")
PC1 = ["CYP2E1", "GLUL", "OAT", "CYP1A2", "CYP3A4", "CPS1", "HAL", "SDS", "ASS1", "PCK1"]
PC3 = ["KRT7", "KRT19", "EPCAM", "SOX9", "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN"]
arows, grows = [], []
for arr in sorted(dm.array_id):
    full = f"{SP}/{arr}.h5ad"; b = os.path.getsize(full); h = sha(full)
    eb, eh = exp_arr.get(arr, (None, None))
    a = ad.read_h5ad(full)
    X = a.X.tocsc() if sp.issparse(a.X) else sp.csc_matrix(a.X)
    sym = a.var.gene_symbol.astype(str).values
    arows.append(dict(array_id=arr, group=dm.set_index("array_id").group[arr], bytes=b, sha256=h, prespec_bytes=eb, prespec_sha256=eh,
                      bytes_match=(eb == b), sha_match=(eh == h), n_obs=a.n_obs, n_var=a.n_vars,
                      obs_cols_ok=all(c in a.obs.columns for c in ["array_row", "array_col", "total_counts", "n_genes"]),
                      obsm_spatial="spatial" in a.obsm, min_total_counts=float(a.obs.total_counts.min())))
    for g in PC1 + PC3:
        cols = np.where(sym == g)[0]
        if len(cols) == 0:
            grows.append(dict(array_id=arr, gene=g, set="PC1" if g in PC1 else "PC3", n_var_rows=0, det_frac=np.nan)); continue
        v = np.asarray(X[:, cols].sum(1)).ravel()
        grows.append(dict(array_id=arr, gene=g, set="PC1" if g in PC1 else "PC3", n_var_rows=len(cols), det_frac=float((v > 0).mean())))
    print(arr, a.n_obs, a.n_vars, flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/sha256_meta_recomputed.tsv", sep="\t", index=False)
A = pd.DataFrame(arows); A.to_csv(f"{OUT}/sha256_h5ad_recomputed.tsv", sep="\t", index=False)
Gd = pd.DataFrame(grows); Gd.to_csv(f"{OUT}/marker_gene_detection.tsv", sep="\t", index=False)
prim = A[A.group == "primary"]
summ = dict(meta_all_match=bool(pd.DataFrame(rows).sha_match.all() and pd.DataFrame(rows).bytes_match.all()),
            meta_mismatch=[r["file"] for r in rows if not (r["sha_match"] and r["bytes_match"])],
            h5ad_all_match=bool(A.sha_match.all() and A.bytes_match.all()), h5ad_mismatch=A[~(A.sha_match & A.bytes_match)].array_id.tolist(),
            n_h5ad_hashes_in_prespec=len(exp_arr), spots_primary=int(prim.n_obs.sum()), spots_min=int(prim.n_obs.min()),
            spots_min_array=prim.loc[prim.n_obs.idxmin(), "array_id"], HRA_55_n_var=int(A.set_index("array_id").n_var["HRA_55"]),
            all_obs_cols_ok=bool(A.obs_cols_ok.all()), all_obsm_spatial=bool(A.obsm_spatial.all()), min_total_counts=float(A.min_total_counts.min()),
            pc3_detection_primary=Gd[(Gd.set == "PC3") & Gd.array_id.isin(prim.array_id)].groupby("gene").det_frac.agg(["min", "median"]).round(4).to_dict(orient="index"),
            pc3_genes_below_1pct_in_any_primary=Gd[(Gd.set == "PC3") & Gd.array_id.isin(prim.array_id) & ~(Gd.det_frac >= 0.01)][["array_id", "gene", "det_frac"]].values.tolist(),
            pc1_genes_below_1pct_in_any_primary=Gd[(Gd.set == "PC1") & Gd.array_id.isin(prim.array_id) & ~(Gd.det_frac >= 0.01)][["array_id", "gene", "det_frac"]].values.tolist())
json.dump(summ, open(f"{OUT}/audit_inputs_summary.json", "w"), indent=1, default=str)
print(json.dumps(summ, indent=1, default=str))
