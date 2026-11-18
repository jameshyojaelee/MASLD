#!/usr/bin/env python3
"""All-16 Visium full-composition sensitivity (20261118 candidate, N08).

Uses the 20261118 within-section nonlinear ridge engine (04_spatial_projection_20261118.py)
with all 16 cell2location q05 factors simultaneously per physical section:
  X_s = [1, log1p(q05_16_s), q05_16_s, log1p(total_counts_s), log1p(n_genes_s)]
and evaluates weighted single-gene mean Moran's I across 9,999 matched-gene draws
for the two focal display programs (Stromal ECM & Ductular injury).
"""

from __future__ import annotations

import hashlib
import importlib.util
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ENGINE_PATH = BASE / "Analysis/Multimodal_Program_Projection/scripts/04_spatial_projection_20261118.py"
OUT = BASE / "Analysis/Multimodal_Program_Projection/candidates/spatial_20261118/full_composition_sensitivity"

FOCAL_UIDS = [
    "hotspot_hepatocytes_f05c535ae5bbc0b9",  # Stromal ECM (IGFBP7), hepatocytes::8
    "hotspot_hepatocytes_48f39dd4d817a10e",  # Ductular injury (BICC1), hepatocytes::20
]


def load_engine():
    spec = importlib.util.spec_from_file_location("spatial_engine_20261118", ENGINE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    eng = load_engine()

    gene_meta = pd.read_csv(eng.GENE_META_FILE, sep="\t", usecols=["gene_name", "gene_biotype", "chromosome"], compression="gzip")
    primary_chr = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}
    gene_meta_prim = gene_meta[gene_meta["chromosome"].isin(primary_chr)].drop_duplicates("gene_name")
    biotype = dict(zip(gene_meta_prim["gene_name"].astype(str), gene_meta_prim["gene_biotype"].astype(str)))

    reg = pd.read_csv(eng.V2_REGISTRY_FILE, sep="\t")
    reg = reg[reg["program_uid"].isin(FOCAL_UIDS)].copy()
    reg["program_id"] = reg["program_uid"]
    reg["program_name"] = reg["module_name"]
    reg = reg.sort_values(["cell_type", "module"]).reset_index(drop=True)
    reg["display_order"] = np.arange(1, len(reg) + 1)

    mem = pd.read_csv(eng.V2_MEMBERSHIP_FILE, sep="\t")
    mem["program_id"] = mem["program_uid"].astype(str)
    mem = mem[(mem["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed") & mem["mapped_symbol"].notna()].copy()
    mem["gene_symbol"] = mem["mapped_symbol"].astype(str).str.strip()
    union = set(mem["gene_symbol"])

    all_rows, unit_rows, match_rows, graph_rows = [], [], [], []
    for ds_idx, (dataset, path) in enumerate(eng.DATASETS.items()):
        adata = ad.read_h5ad(path)
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        counts = eng.dense_counts(adata)
        library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        norm = eng.log_normalize(counts, library)
        mean = np.asarray(norm.mean(axis=0)).ravel()
        detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        pool = eng.matched_pool_table(adata.var_names, mean, detection, biotype)

        factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
        q05_16 = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
        hep_col = q05_16[:, factor_names.index("Hepatocytes")]
        pool["lineage_corr"] = eng.abundance_correlation(norm, hep_col)
        pool["lineage_corr_bin"] = pd.qcut(pool["lineage_corr"].rank(method="first"), 10, labels=False, duplicates="drop")

        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        graphs, gaudit = eng.section_graphs(obs, coords, k=6)
        gaudit.insert(0, "dataset", dataset)
        graph_rows.append(gaudit)

        res_16 = eng.compute_within_section_single_gene_moran(norm, obs, graphs, q05_16, extra_cov=None, compute_raw=True)
        bar_res = res_16["bar_I_resid"]
        finite_m = np.isfinite(bar_res)
        g2i = {str(g): i for i, g in enumerate(adata.var_names)}
        ind_by_s = obs[["sample_id", "individual"]].drop_duplicates().set_index("sample_id")["individual"].astype(str).to_dict()

        rng = np.random.default_rng(eng.SEED + 1000 * ds_idx)
        for row in reg.itertuples(index=False):
            pid = str(row.program_id)
            m = mem[mem["program_id"] == pid].groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
            m["present"] = m["gene_symbol"].isin(g2i)
            m["detected"] = m["gene_symbol"].map(lambda g: (detection[g2i[g]] >= 0.01 and bool(finite_m[g2i[g]])) if g in g2i else False)
            measured = m[m["present"] & m["detected"]].copy()
            retained = float(measured["original_l1_weight"].sum())
            testable = len(measured) >= 8 and retained >= 0.20
            measured["weight"] = measured["original_l1_weight"] / retained
            idx = np.array([g2i[g] for g in measured["gene_symbol"]], dtype=int)
            w = measured["weight"].to_numpy(np.float64)

            sets, er, lr = eng.build_matched_sets(measured, pool, union, eng.N_NULL, rng, finite_mask=finite_m)
            sha = hashlib.sha256(sets.tobytes()).hexdigest()
            match_rows.append({
                "dataset": dataset, "program_id": pid, "n_measured": len(measured),
                "n_genes_expression_relaxed": int(np.sum(np.asarray(er) > 0)),
                "n_genes_lineage_relaxed": int(np.sum(np.asarray(lr) > 0)),
                "matched_set_sha256": sha,
            })
            stats = eng._score_weighted_mean_moran(bar_res, res_16["bar_I_raw"], None, idx, w, sets, eng.N_SENSITIVITY_NULL)
            stats.pop("null_res")
            all_rows.append({
                "program_id": pid, "dataset": dataset, "display_order": int(row.display_order),
                "cell_type": str(row.cell_type), "module": int(row.module), "program_name": str(row.program_name),
                "n_measured": len(measured), "retained_l1_weight": retained, "testable": testable,
                "n_null": eng.N_NULL, "matched_set_sha256": sha, **stats,
            })
            for s in res_16["arrays"]:
                si = eng.nan_weighted_dot(res_16["sec_resid"][s][idx], w)
                snull = eng.nan_weighted_mat(res_16["sec_resid"][s][sets], w)
                snull_mean = float(np.nanmean(snull))
                unit_rows.append({
                    "dataset": dataset, "program_id": pid, "sample_id": s, "individual": ind_by_s[s],
                    "residual_moran_i": si, "null_mean": snull_mean,
                    "centered_excess": si - snull_mean,
                })

    res_df = pd.DataFrame(all_rows)
    res_df["residual_padj"] = np.nan
    for ds in res_df["dataset"].unique():
        m = res_df["dataset"] == ds
        res_df.loc[m, "residual_padj"] = eng.bh(res_df.loc[m, "residual_pvalue"]).values
    res_df["robust"] = res_df["testable"] & (res_df["residual_padj"] < 0.05) & res_df["sensitivity_sign_agree"]
    res_df.to_csv(OUT / "full_composition_program_results.tsv", sep="\t", index=False)
    pd.DataFrame(unit_rows).to_csv(OUT / "full_composition_unit_results.tsv", sep="\t", index=False)
    pd.DataFrame(match_rows).to_csv(OUT / "full_composition_matching_audit.tsv", sep="\t", index=False)
    pd.concat(graph_rows, ignore_index=True).to_csv(OUT / "full_composition_graph_audit.tsv", sep="\t", index=False)
    print(res_df[["dataset", "program_name", "residual_moran_i", "residual_null_mean", "residual_moran_z", "residual_pvalue", "residual_padj", "robust"]].to_string(index=False))


if __name__ == "__main__":
    main()
