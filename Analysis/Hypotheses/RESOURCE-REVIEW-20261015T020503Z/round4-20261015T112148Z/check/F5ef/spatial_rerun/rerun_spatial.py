#!/usr/bin/env python3
"""Review G3: regenerate the Fig 5E/5F spatial Moran statistics from spot-level Visium data.

Uses the byte-pinned engine copy (04_spatial_projection.py, sha256 d1985837...) for every
scientific function (normalization, matched pools, graphs, Moran, donor collapse, matched sets).
One pass per dataset/program computes the residual Moran I under two nuisance designs on the
same matched-gene draws:
  lineage_only : hepatocyte q05 abundance + log1p total counts + log1p detected genes
                 (the v2 primary design; RESULTS "fell from 0.1390 / 0.1471"; section values)
  all16        : all 16 cell2location q05 factors + the two QC covariates (Fig 5F design)
Programs: v2 registry robust_display (Stromal ECM hepatocytes::8, Ductular injury hepatocytes::20).
Seeds: 42 (GSE192741), 1042 (Vu). 9,999 draws. Per-section and per-donor nulls are kept.
Also computes observed values with the v1 membership (Fig 5E source) for comparison.
Raw, zonation and the 999-draw equal-weight/leave-top sensitivity nulls are not recomputed.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np

JOB = Path(__file__).resolve().parent
PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUTDIR = JOB / "out"
OUTDIR.mkdir(exist_ok=True)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# Import the pinned engine with its output root redirected into the job directory.
os.environ["MASLD_PROJECT_ROOT"] = str(JOB / "fakeroot")
os.environ["FIG4_SPATIAL_SMOKE"] = "FALSE"
spec = importlib.util.spec_from_file_location("eng", JOB / "04_spatial_projection.py")
eng = importlib.util.module_from_spec(spec)
spec.loader.exec_module(eng)
eng.SMOKE = False
eng.N_NULL = 9999
eng.SEED = 42
eng.GRAPH_DISTANCE_MULTIPLIER = 2.5
eng.MIN_GRAPH_SPOTS = 8

import anndata as ad  # noqa: E402
import pandas as pd  # noqa: E402

N_NULL = 9999
DATASETS = {
    "GSE192741": PROJECT / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": PROJECT / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}
HOT = PROJECT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
RECORDED_SHA = {
    ("GSE192741", "hotspot_hepatocytes_f05c535ae5bbc0b9"): "f2ecc252d7142f2addbb4b2e2262ac890d40d8dd57061a2ede75811fd7adae01",
    ("GSE192741", "hotspot_hepatocytes_48f39dd4d817a10e"): "04c56766f81c38acdbef145fbb6cfaadda841f940e38c2cfc610460f25549727",
    ("Vu_et_al_2025", "hotspot_hepatocytes_f05c535ae5bbc0b9"): "9c2fbcafb172f83b78084ab89cdc42fb372cbe9f2a763dfd15c604b8995461d8",
    ("Vu_et_al_2025", "hotspot_hepatocytes_48f39dd4d817a10e"): "f0d724d33da6d7357ac2a2deff4218a22a52a34d1c572cd192f9261804fd1e17",
}


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# v2 registry and membership exactly as visium_rerun_lib.prepare_engine_inputs("v2")
reg_rows = [r for r in read_tsv(HOT / "program_registry_v2.tsv") if r["robust_display"] == "TRUE"]
reg_rows.sort(key=lambda r: (r["cell_type"], int(r["module"]), r["program_uid"]))
assert len(reg_rows) == 2, len(reg_rows)
registry = pd.DataFrame([
    {"program_id": r["program_uid"], "display_order": i + 1, "cell_type": r["cell_type"],
     "module": int(r["module"]), "program_name": r["module_name"]}
    for i, r in enumerate(reg_rows)
])
sel = {r["program_uid"]: r for r in reg_rows}
mem_rows = []
for r in read_tsv(HOT / "program_membership_v2.tsv"):
    if r["program_uid"] in sel and r["mapped_symbol_status"] == "gencode_v49_unique_symbol_confirmed":
        mem_rows.append({"program_id": r["program_uid"], "cell_type": r["cell_type"], "module": int(r["module"]),
                         "source_gene": r["source_gene"], "gene_symbol": r["mapped_symbol"].strip(),
                         "mapped_symbol": True, "original_l1_weight": float(r["original_l1_weight"])})
membership = pd.DataFrame(mem_rows)
log(f"v2 registry {registry[['program_id','program_name']].values.tolist()}; membership rows {len(membership)}")

# v1 registry/membership (Fig 5E source) for observed-only comparison
v1_mem = pd.read_csv(PROJECT / "Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv", sep="\t")
v1_mem = v1_mem[v1_mem["mapped_symbol"].fillna(False) & v1_mem["gene_symbol"].notna()].copy()

gene_meta = pd.read_csv(PROJECT / "data/gencode_v49_gene_metadata.tsv.gz", sep="\t",
                        usecols=["gene_name", "gene_biotype"], compression="gzip").drop_duplicates("gene_name")
biotype = dict(zip(gene_meta["gene_name"].astype(str), gene_meta["gene_biotype"].astype(str)))


def standardize_design(columns):
    design = np.column_stack(columns).astype(float, copy=False)
    for c in range(1, design.shape[1]):
        mu = np.nanmean(design[:, c])
        sd = np.nanstd(design[:, c], ddof=1)
        design[:, c] = (design[:, c] - mu) / sd if np.isfinite(sd) and sd > 0 else 0.0
    design[~np.isfinite(design)] = 0.0
    return design


def measured_table(mem, pid, gene_to_index, detection):
    m = mem[mem["program_id"] == pid].copy()
    m = m.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
    m["present"] = m["gene_symbol"].isin(gene_to_index)
    m["detected"] = m["gene_symbol"].map(
        lambda g: detection[gene_to_index[g]] >= 0.01 if g in gene_to_index else False)
    measured = m[m["present"] & m["detected"]].copy()
    retained = float(measured["original_l1_weight"].sum())
    measured["weight"] = measured["original_l1_weight"] / retained
    return measured, retained


def pval(obs, null):
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    k = int(np.sum(null >= obs))
    return k, (1.0 + k) / (len(null) + 1.0), float(np.mean(null)), float(np.std(null, ddof=1)), len(null)


prog_rows, unit_rows, v1_rows = [], [], []
npz_payload = {}

for dataset, path in DATASETS.items():
    log(f"loading {dataset}")
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
    gene_to_index = {g: i for i, g in enumerate(adata.var_names)}
    program_gene_union = set(membership["gene_symbol"].dropna().astype(str))
    factor_names = [str(x).replace("means_per_cluster_mu_fg_", "") for x in adata.uns["mod"]["factor_names"]]
    abundance_mat = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
    hep = abundance_mat[:, factor_names.index("Hepatocytes")]
    lineage_pool = pool.copy()
    lineage_pool["lineage_corr"] = eng.abundance_correlation(norm, hep)
    lineage_pool["lineage_corr_bin"] = pd.qcut(lineage_pool["lineage_corr"].rank(method="first"), 10,
                                               labels=False, duplicates="drop")
    coords = np.asarray(adata.obsm["spatial"])
    # all-16 design as in 31_run_full_composition_sensitivity.load_full_designs
    canonical = sorted(factor_names)
    assert len(canonical) == 16 and len(set(canonical)) == 16
    reorder = [factor_names.index(n) for n in canonical]
    ab16 = abundance_mat[:, reorder]
    detected_genes = pd.to_numeric(obs["n_genes_by_counts"], errors="coerce").to_numpy(float)
    design_all16 = standardize_design([np.ones(len(obs))] + [ab16[:, c] for c in range(16)]
                                      + [np.log1p(library), np.log1p(detected_genes)])
    design_lin = eng.design_matrix(obs, hep)
    designs = {"lineage_only": design_lin, "all16": design_all16}
    log(f"{dataset}: spots {len(obs)} genes {adata.n_vars} design ranks "
        f"{np.linalg.matrix_rank(design_lin)}/{np.linalg.matrix_rank(design_all16)} "
        f"cond all16 {np.linalg.cond(design_all16):.4f}")
    del adata
    graphs, graph_audit = eng.section_graphs(obs, coords, k=6)
    assert set(graphs) == set(obs["sample_id"])
    samples = list(graphs.keys())
    sample_ind = obs[["sample_id", "individual", "condition"]].drop_duplicates().set_index("sample_id")
    norm_csc = norm.tocsc()
    rng = np.random.default_rng(42 + (0 if dataset == "GSE192741" else 1000))

    # v1 membership observed values (no null)
    for pid in ["hepatocytes::8", "hepatocytes::20"]:
        measured, retained = measured_table(v1_mem, pid, gene_to_index, detection)
        idx = np.array([gene_to_index[g] for g in measured["gene_symbol"]], dtype=int)
        score = eng.extract_z(norm, idx) @ measured["weight"].to_numpy(float)
        for dname, des in designs.items():
            by = eng.moran_columns(eng.residualize(score, des), graphs)
            coll = float(eng.donor_collapse(by, obs, dataset)[0])
            v1_rows.append({"dataset": dataset, "program_v1": pid, "design": dname, "unit": "dataset_donor_collapsed",
                            "n_measured": len(measured), "retained_l1": retained, "residual_moran_i": coll})
            for s in samples:
                v1_rows.append({"dataset": dataset, "program_v1": pid, "design": dname, "unit": s,
                                "n_measured": len(measured), "retained_l1": retained,
                                "residual_moran_i": float(by[s][0])})

    for row in registry.itertuples(index=False):
        pid = row.program_id
        measured, retained = measured_table(membership, pid, gene_to_index, detection)
        indices = np.array([gene_to_index[g] for g in measured["gene_symbol"]], dtype=int)
        weights = measured["weight"].to_numpy(float)
        z = eng.extract_z(norm, indices)
        score = z @ weights
        observed_by = {}
        observed = {}
        for dname, des in designs.items():
            by = eng.moran_columns(eng.residualize(score, des), graphs)
            observed_by[dname] = {s: float(by[s][0]) for s in samples}
            observed[dname] = float(eng.donor_collapse(by, obs, dataset)[0])
        raw_by = eng.moran_columns(score, graphs)
        observed_raw = float(eng.donor_collapse(raw_by, obs, dataset)[0])
        log(f"{dataset} {row.program_name}: n={len(measured)} L1={retained:.4f} observed "
            f"lin={observed['lineage_only']:.6f} all16={observed['all16']:.6f} raw={observed_raw:.6f}")
        matched_sets, _, _ = eng.build_matched_sets(measured, lineage_pool, program_gene_union, N_NULL, rng)
        sha = hashlib.sha256(matched_sets.tobytes()).hexdigest()
        log(f"  matched sets built sha {sha[:16]} recorded match {sha == RECORDED_SHA[(dataset, pid)]}")
        null_sample = {d: {s: [] for s in samples} for d in designs}
        null_coll = {d: [] for d in designs}
        for start in range(0, N_NULL, 64):
            block = matched_sets[start:start + 64]
            uniq, inverse = np.unique(block, return_inverse=True)
            zz = eng.extract_z(norm_csc, uniq)
            mapped = inverse.reshape(block.shape)
            scores = np.zeros((norm.shape[0], len(block)), dtype=np.float32)
            for j, w in enumerate(weights):
                scores += zz[:, mapped[:, j]] * np.float32(w)
            for dname, des in designs.items():
                by = eng.moran_columns(eng.residualize(scores, des), graphs)
                for s in samples:
                    null_sample[dname][s].append(by[s])
                null_coll[dname].append(eng.donor_collapse(by, obs, dataset))
        log(f"  null loop done")
        for dname in designs:
            ns = {s: np.concatenate(v) for s, v in null_sample[dname].items()}
            nc = np.concatenate(null_coll[dname])
            k, p, mu, sd, n = pval(observed[dname], nc)
            prog_rows.append({"dataset": dataset, "program_id": pid, "program_name": row.program_name,
                              "design": dname, "n_measured": len(measured), "retained_l1": retained,
                              "residual_moran_i": observed[dname], "null_mean": mu, "null_sd": sd,
                              "excess": observed[dname] - mu, "k_null_ge_obs": k, "n_null": n, "p": p,
                              "matched_set_sha256": sha,
                              "sha_matches_recorded": sha == RECORDED_SHA[(dataset, pid)]})
            # per section
            for s in samples:
                k, p, mu, sd, n = pval(observed_by[dname][s], ns[s])
                unit_rows.append({"dataset": dataset, "program_name": row.program_name, "design": dname,
                                  "unit_type": "section", "unit_id": s,
                                  "individual": sample_ind.loc[s, "individual"],
                                  "condition": str(sample_ind.loc[s, "condition"]),
                                  "residual_moran_i": observed_by[dname][s], "null_mean": mu,
                                  "excess": observed_by[dname][s] - mu, "k_null_ge_obs": k, "p": p})
            # per donor (GSE) / per array (Vu, individual == array)
            table = obs[["sample_id", "individual"]].drop_duplicates()
            for ind, part in table.groupby("individual", observed=False):
                arr = [str(x) for x in part["sample_id"]]
                o = float(np.mean([observed_by[dname][s] for s in arr]))
                nd = np.mean(np.vstack([ns[s] for s in arr]), axis=0)
                k, p, mu, sd, n = pval(o, nd)
                unit_rows.append({"dataset": dataset, "program_name": row.program_name, "design": dname,
                                  "unit_type": "donor_or_array", "unit_id": ind, "individual": ind,
                                  "condition": str(sample_ind.loc[arr[0], "condition"]),
                                  "residual_moran_i": o, "null_mean": mu, "excess": o - mu,
                                  "k_null_ge_obs": k, "p": p})
                npz_payload[f"{dataset}|{row.program_name}|{dname}|{ind}"] = nd.astype(np.float64)
            npz_payload[f"{dataset}|{row.program_name}|{dname}|collapsed"] = nc.astype(np.float64)
        pd.DataFrame(prog_rows).to_csv(OUTDIR / "program_level.tsv", sep="\t", index=False)
        pd.DataFrame(unit_rows).to_csv(OUTDIR / "unit_level.tsv", sep="\t", index=False)
        pd.DataFrame(v1_rows).to_csv(OUTDIR / "v1_membership_observed.tsv", sep="\t", index=False)
        np.savez_compressed(OUTDIR / "null_donor_draws.npz", **npz_payload)
    del norm, norm_csc, counts

# BH within dataset across the two programs, per design
prog = pd.DataFrame(prog_rows)
prog["q_bh_two_program_within_dataset"] = np.nan
for (ds, dn), g in prog.groupby(["dataset", "design"]):
    prog.loc[g.index, "q_bh_two_program_within_dataset"] = eng.bh(g["p"]).to_numpy()
prog.to_csv(OUTDIR / "program_level.tsv", sep="\t", index=False)
log("done")
