#!/usr/bin/env python3
"""Rebuild the observed donor-first residual Moran's I for Stromal ECM
(hepatocytes::8) under two memberships (v1 frozen, v2 mapped) and two designs
(matching-lineage abundance + 2 QC, all 16 cell2location factors + 2 QC).
Uses the copied pinned engine functions; no null draws."""
import importlib.util
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "eng", os.path.join(HERE, "engine_04_spatial_projection.py")
)
eng = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eng)

REAL = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
PATHS = {
    "GSE192741": REAL + "/Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": REAL + "/Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}

v1 = pd.read_csv(REAL + "/Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv", sep="\t")
v1 = v1[v1["mapped_symbol"].fillna(False) & v1["gene_symbol"].notna()]
v1 = v1[v1["program_id"] == "hepatocytes::8"]
v1 = v1.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()

v2 = pd.read_csv(
    REAL + "/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv",
    sep="\t",
)
v2 = v2[(v2["cell_type"] == "hepatocytes") & (v2["module"].astype(int) == 8) & v2["mapped_symbol"].notna()]
v2 = v2.rename(columns={"mapped_symbol": "gene_symbol_m"})
v2 = v2.groupby("gene_symbol_m", as_index=False)["original_l1_weight"].sum().rename(columns={"gene_symbol_m": "gene_symbol"})
MEMB = {"v1": v1, "v2": v2}
print("membership sizes:", {k: len(v) for k, v in MEMB.items()}, flush=True)

rows = []
for dataset, path in PATHS.items():
    adata = ad.read_h5ad(path)
    adata.var_names = adata.var_names.astype(str)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["individual"] = obs["individual"].astype(str)
    counts = eng.dense_counts(adata)
    library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
    norm = eng.log_normalize(counts, library)
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    g2i = {g: i for i, g in enumerate(adata.var_names)}
    factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
    abund = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
    hep = abund[:, factor_names.index("Hepatocytes")]
    coords = np.asarray(adata.obsm["spatial"])
    graphs, _ = eng.section_graphs(obs, coords, k=6)

    design_hep = eng.design_matrix(obs, hep)
    total = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
    detected = pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").to_numpy(float)
    cols = [np.ones(len(obs))] + [abund[:, j] for j in range(abund.shape[1])] + [np.log1p(total), np.log1p(detected)]
    design_all = np.column_stack(cols).astype(float)
    for j in range(1, design_all.shape[1]):
        mu = np.nanmean(design_all[:, j]); sd = np.nanstd(design_all[:, j], ddof=1)
        design_all[:, j] = (design_all[:, j] - mu) / sd if np.isfinite(sd) and sd > 0 else 0.0
    design_all[~np.isfinite(design_all)] = 0.0
    print(dataset, "n_spots", len(obs), "n_factors", abund.shape[1], "design_all cols", design_all.shape[1],
          "rank", np.linalg.matrix_rank(design_all), flush=True)

    for mname, m in MEMB.items():
        m = m.copy()
        m["present"] = m["gene_symbol"].isin(g2i)
        m["detected"] = m["gene_symbol"].map(lambda g: detection[g2i[g]] >= 0.01 if g in g2i else False)
        meas = m[m["present"] & m["detected"]].copy()
        retained = float(meas["original_l1_weight"].sum())
        w = (meas["original_l1_weight"] / retained).to_numpy(float)
        idx = np.array([g2i[g] for g in meas["gene_symbol"]], dtype=int)
        z = eng.extract_z(norm, idx)
        score = z @ w
        for dname, design in (("matching_lineage_hepatocytes_plus_2qc", design_hep),
                              ("all_16_cell2location_plus_2qc", design_all)):
            resid = eng.residualize(score, design)
            by = eng.moran_columns(resid, graphs)
            val = float(eng.donor_collapse(by, obs, dataset)[0])
            rows.append({
                "dataset": dataset, "membership": mname, "design": dname,
                "n_measured": len(meas), "retained_l1_weight": retained,
                "residual_moran_i": val,
                "measured_genes": ",".join(sorted(meas["gene_symbol"])),
            })
            print(dataset, mname, dname, len(meas), round(retained, 6), repr(val), flush=True)
    del adata, counts, norm

out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "observed_residual_moran_rebuild.tsv"), sep="\t", index=False)
print(out.drop(columns=["measured_genes"]).to_string(), flush=True)
