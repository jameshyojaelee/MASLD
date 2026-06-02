#!/usr/bin/env python
"""
figS_scdrs_prep.py

One-shot preprocessor for the figS_scdrs supplementary figures.
The R figure scripts (figS_scdrs_bulk.R, figS_scdrs_gwas.R) need per-cell
scDRS scores joined with disease-stage metadata. The cell-score parquets
are huge (10.9M rows bulk; 15.5M rows GWAS) and the R rnaseq env lacks the
arrow package, so we do the join + subsampling here in Python where the
parquet readers exist, and dump compact CSVs the R scripts can fread().

Outputs (under disease_signatures/figure_prep/):
  bulk_cell_scores_top4ct_signed.csv.gz   — per-cell bulk-DEG scores
      cols: cell_id, cell_type, disease_stage, transition, norm_score
      restricted to: signed mode × {Hepatocytes, Macrophages, Fibroblasts,
                                    Cholangiocytes} × {Healthy, Steatosis,
                                    Steatohepatitis} × 4 F-transitions
      subsampled to <= MAX_PER_BIN cells per (CT × transition × stage) bin
  gwas_cell_scores_top4ct_polyfun.csv.gz  — per-cell GWAS-COLOC scores
      cols: cell_id, cell_type, disease_stage, trait_class, panel, norm_score
      restricted to: 4 phenotype panels × PolyFun LD × top 4 CTs × staged
"""
import os
import sys
import pandas as pd
import numpy as np

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SIG  = os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures")
OUT  = os.path.join(SIG, "figure_prep")
os.makedirs(OUT, exist_ok=True)

TOP_CT = ["Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes"]
STAGES = ["Healthy", "Steatosis", "Steatohepatitis"]
PHENO_TRAITS = [
    "gwas_coloc_nafld_specific_polyfun",
    "gwas_coloc_liver_enzymes_polyfun",
    "gwas_coloc_pdff_polyfun",
    "gwas_coloc_cirrhosis_hcc_polyfun",
]
MAX_PER_BIN = 6000  # subsample cap per (CT, transition/trait, stage) for plotting

# ---- Load obs metadata ----
print("[prep] loading obs metadata", flush=True)
obs = pd.read_csv(os.path.join(SIG, "scdrs_obs_metadata.csv"))
obs = obs[
    obs["disease_stage_coarse"].isin(STAGES)
    & obs["cell_type"].isin(TOP_CT)
][["cell_id", "cell_type", "disease_stage_coarse"]].rename(
    columns={"disease_stage_coarse": "disease_stage"}
)
print(f"[prep] obs rows (top4 CT × staged): {len(obs):,}", flush=True)
keep_ids = set(obs["cell_id"])

# ---- Bulk per-cell scores ----
print("[prep] reading bulk parquet (signed mode only)", flush=True)
cs = pd.read_parquet(
    os.path.join(SIG, "scdrs_cell_scores.parquet"),
    columns=["cell_id", "norm_score", "transition", "mode"],
    filters=[("mode", "==", "signed")],
)
cs = cs[cs["cell_id"].isin(keep_ids)]
cs = cs.merge(obs, on="cell_id", how="inner")
print(f"[prep] bulk per-cell joined rows: {len(cs):,}", flush=True)

# Subsample per (cell_type, transition, disease_stage)
rng = np.random.default_rng(42)
def _sub(df):
    if len(df) <= MAX_PER_BIN:
        return df
    return df.iloc[rng.choice(len(df), MAX_PER_BIN, replace=False)]

bulk_out = (
    cs.groupby(["cell_type", "transition", "disease_stage"], group_keys=False)
      .apply(_sub)
      .reset_index(drop=True)
)
bulk_path = os.path.join(OUT, "bulk_cell_scores_top4ct_signed.csv.gz")
bulk_out.to_csv(bulk_path, index=False, compression="gzip")
print(f"[prep] wrote {bulk_path}  ({len(bulk_out):,} rows)", flush=True)

# Median trajectory table (no subsampling — use full joined CS)
median_bulk = (
    cs.groupby(["cell_type", "transition", "disease_stage"])["norm_score"]
      .agg(median="median",
           q25=lambda x: x.quantile(0.25),
           q75=lambda x: x.quantile(0.75),
           n="size")
      .reset_index()
)
median_path = os.path.join(OUT, "bulk_cell_scores_median_by_stage.csv")
median_bulk.to_csv(median_path, index=False)
print(f"[prep] wrote {median_path}", flush=True)

# ---- GWAS per-cell scores ----
print("[prep] reading GWAS parquet (polyfun panels)", flush=True)
gw = pd.read_parquet(
    os.path.join(SIG, "scdrs_gwas_cell_scores.parquet"),
    columns=["cell_id", "norm_score", "trait", "trait_class", "panel"],
    filters=[("panel", "==", "polyfun")],
)
gw = gw[gw["trait"].isin(PHENO_TRAITS) & gw["cell_id"].isin(keep_ids)]
gw = gw.merge(obs, on="cell_id", how="inner")
print(f"[prep] GWAS per-cell joined rows: {len(gw):,}", flush=True)

gwas_out = (
    gw.groupby(["cell_type", "trait_class", "disease_stage"], group_keys=False)
      .apply(_sub)
      .reset_index(drop=True)
)
gwas_path = os.path.join(OUT, "gwas_cell_scores_top4ct_polyfun.csv.gz")
gwas_out.to_csv(gwas_path, index=False, compression="gzip")
print(f"[prep] wrote {gwas_path}  ({len(gwas_out):,} rows)", flush=True)

median_gwas = (
    gw.groupby(["cell_type", "trait_class", "disease_stage"])["norm_score"]
      .agg(median="median",
           q25=lambda x: x.quantile(0.25),
           q75=lambda x: x.quantile(0.75),
           n="size")
      .reset_index()
)
median_gwas_path = os.path.join(OUT, "gwas_cell_scores_median_by_stage.csv")
median_gwas.to_csv(median_gwas_path, index=False)
print(f"[prep] wrote {median_gwas_path}", flush=True)

# ---- Macrophage UMAP slice (full per-cell signed scores for UMAP) ----
print("[prep] preparing macrophage UMAP slice", flush=True)
mac_meta = pd.read_csv(
    os.path.join(BASE, "Analysis/SingleCell/results_gpu_v2/pseudotime/Macrophages_metadata.csv"),
    index_col=0,
)
if "cell_id" not in mac_meta.columns:
    mac_meta = mac_meta.reset_index().rename(columns={"index": "cell_id"})
mac_ids = set(mac_meta["cell_id"])

mac_cs = pd.read_parquet(
    os.path.join(SIG, "scdrs_cell_scores.parquet"),
    columns=["cell_id", "norm_score", "transition", "mode"],
    filters=[("mode", "==", "signed")],
)
mac_cs = mac_cs[mac_cs["cell_id"].isin(mac_ids)]
mac_join = mac_meta[["cell_id", "UMAP_1", "UMAP_2"]].merge(
    mac_cs, on="cell_id", how="inner"
)
mac_path = os.path.join(OUT, "macrophage_umap_signed_scores.csv.gz")
mac_join.to_csv(mac_path, index=False, compression="gzip")
print(f"[prep] wrote {mac_path}  ({len(mac_join):,} rows)", flush=True)

print("[prep] done.", flush=True)
