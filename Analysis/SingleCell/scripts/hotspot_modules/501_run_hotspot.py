"""Run Hotspot on a single cell-type subset of the integrated MASLD atlas.

Inputs:
  --cell-type {global,hepatocytes,macrophages,fibroblasts,endothelial_cells,cholangiocytes,tcells}
  --smoke    (subsample to 5K cells for fast iteration)
  --loo-dataset NAME   (exclude one dataset; used by 503)
  --exclude-datasets NAME [NAME ...]   (exclude one or more datasets; protocol-contamination remediation)

Env vars (tuning):
  HOTSPOT_MAX_GENES  Cap on autocorr-significant genes carried into local
                     correlations (default 5000). Local correlations scale
                     O(genes^2 * cells), so reducing to 2000 gives ~6x
                     speedup at minimal sensitivity loss for the top modules.
  HOTSPOT_JOBS       n_jobs passed to compute_autocorrelations and
                     compute_local_correlations (default 8).

Outputs (written to results_gpu_v2/hotspot_modules/<cell_type>/[loo/<dataset>/]):
  autocorr.tsv          gene,Z,FDR,C
  module_genes.tsv      gene,module,weight    (weight = mean local-corr Z within module)
  cell_scores.parquet   cell_id,module,score  (long form)
  donor_scores.tsv      sample,module,score   (mean across cells per donor)
  hotspot_obj.pkl       pickled Hotspot internals for downstream re-use
  run_metadata.json     n_cells, n_genes_kept, n_modules, params
"""
from __future__ import annotations
import argparse
import json
import os
import pickle
import sys
from pathlib import Path

import anndata as ad
import hotspot
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import (
    load_atlas, load_gene_biotypes, hotspot_outdir, RUN_ORDER,
)

CONFOUNDER_PREFIXES = ("RPS", "RPL", "MT-", "MRPS", "MRPL")
CONFOUNDER_LISTS = {
    "HB": ["HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBM", "HBQ1", "HBZ"],
    "XIST_ESCAPEES": ["XIST", "TSIX", "KDM6A", "DDX3X", "DDX3Y", "EIF1AY", "USP9Y", "UTY", "RPS4Y1", "ZFY", "KDM5D", "NLGN4Y"],
}


def strip_confounders(adata: ad.AnnData) -> ad.AnnData:
    """Drop ribosomal, mitochondrial, hemoglobin, sex-escape, chrY, chrM genes."""
    bio = load_gene_biotypes()
    keep_biotypes = {"protein_coding", "lncRNA"}
    keep_gene = set(bio.loc[bio["gene_biotype"].isin(keep_biotypes), "gene_name"])
    keep_gene -= set(bio.loc[bio["chromosome"].isin(["chrY", "chrM"]), "gene_name"])

    var_names = adata.var_names.astype(str)
    mask = (
        np.isin(var_names, list(keep_gene))
        & ~np.array([n.startswith(CONFOUNDER_PREFIXES) for n in var_names])
        & ~np.isin(var_names, CONFOUNDER_LISTS["HB"])
        & ~np.isin(var_names, CONFOUNDER_LISTS["XIST_ESCAPEES"])
    )
    out = adata[:, mask].copy()
    print(f"   stripped confounders: {adata.n_vars} -> {out.n_vars} genes")
    return out


def filter_detected(adata: ad.AnnData, min_frac: float = 0.01) -> ad.AnnData:
    """Keep genes detected in >=min_frac of cells."""
    X = adata.layers.get("counts", adata.X)
    n_cells_per_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    mask = n_cells_per_gene >= int(adata.n_obs * min_frac)
    out = adata[:, mask].copy()
    print(f"   detected >= {min_frac:.0%}: {adata.n_vars} -> {out.n_vars} genes")
    return out


def ensure_raw_layer(adata: ad.AnnData) -> None:
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if "n_counts" not in adata.obs:
        adata.obs["n_counts"] = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()


def resolve_latent(adata: ad.AnnData) -> str:
    """Pick an existing low-dim latent; if none present, compute PCA(30) on the fly."""
    for key in ("X_scvi", "X_scVI", "X_pca_harmony", "X_pca"):
        if key in adata.obsm:
            return key
    print("   no latent in obsm — computing PCA(30) on log-normalized data")
    norm = adata.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)
    sc.pp.pca(norm, n_comps=30)
    adata.obsm["X_pca"] = norm.obsm["X_pca"]
    return "X_pca"


def run_hotspot(
    adata: ad.AnnData,
    n_neighbors: int = 30,
    autocorr_z: float = 7.0,
    autocorr_fdr: float = 0.05,
    max_genes: int = 5000,
    min_gene_threshold: int = 20,
    fdr_threshold: float = 0.05,
    jobs: int = 8,
) -> tuple[hotspot.Hotspot, pd.DataFrame]:
    ensure_raw_layer(adata)
    latent_key = resolve_latent(adata)
    print(f"   using latent: {latent_key}")
    hs = hotspot.Hotspot(
        adata,
        layer_key="counts",
        model="danb",
        latent_obsm_key=latent_key,
        umi_counts_obs_key="n_counts",
    )
    print("   create_knn_graph...")
    hs.create_knn_graph(weighted_graph=False, n_neighbors=n_neighbors)
    print("   compute_autocorrelations...")
    hs_results = hs.compute_autocorrelations(jobs=jobs)
    autocorr_df = hs_results.copy()
    autocorr_df.index.name = "gene"
    autocorr_df = autocorr_df.reset_index()

    selected = (
        autocorr_df.query("FDR < @autocorr_fdr & Z >= @autocorr_z")
        .sort_values("Z", ascending=False)
        .head(max_genes)
    )
    print(f"   {len(selected)} genes pass autocorr Z>={autocorr_z}, FDR<{autocorr_fdr}")
    if len(selected) < min_gene_threshold:
        raise RuntimeError(
            f"Only {len(selected)} autocorr-significant genes - below min_gene_threshold."
        )

    hs_genes = selected["gene"].tolist()
    print("   compute_local_correlations...")
    hs.compute_local_correlations(hs_genes, jobs=jobs)
    print("   create_modules...")
    modules = hs.create_modules(
        min_gene_threshold=min_gene_threshold,
        core_only=True,
        fdr_threshold=fdr_threshold,
    )
    n_modules = int(modules.max())
    print(f"   created {n_modules} modules")
    return hs, autocorr_df


def aggregate_donor(cell_scores: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    """Mean module score per donor."""
    merged = cell_scores.merge(
        obs[["sample"]].reset_index().rename(columns={"index": "cell_id"}),
        on="cell_id",
    )
    return (
        merged.groupby(["sample", "module"], as_index=False)["score"].mean()
    )


def main(
    cell_type: str,
    smoke: bool,
    loo_dataset: str | None,
    exclude_datasets: list[str] | None = None,
) -> None:
    assert cell_type in RUN_ORDER, f"Unknown cell_type: {cell_type}"
    print(f"[1/6] Loading {cell_type} (smoke={smoke}, loo={loo_dataset}, exclude={exclude_datasets})")
    adata = load_atlas(cell_type, smoke=smoke)
    # Apply LOO filter (single dataset) and/or exclude-datasets filter (multiple)
    if loo_dataset is not None:
        n_before = adata.n_obs
        adata = adata[adata.obs["dataset"] != loo_dataset].copy()
        print(f"[501] LOO filter: removed {n_before - adata.n_obs} cells from {loo_dataset}")
    if exclude_datasets:
        n_before = adata.n_obs
        excl = set(exclude_datasets)
        adata = adata[~adata.obs["dataset"].isin(excl)].copy()
        print(f"[501] exclude_datasets filter: removed {n_before - adata.n_obs} cells from {sorted(excl)}")

    print(f"[2/6] Gene filtering")
    adata = strip_confounders(adata)
    adata = filter_detected(adata, min_frac=0.01)

    print(f"[3/6] Running Hotspot ({adata.n_obs:,} cells, {adata.n_vars:,} genes)")
    max_genes_env = int(os.environ.get("HOTSPOT_MAX_GENES", "5000"))
    jobs_env = int(os.environ.get("HOTSPOT_JOBS", "8"))
    # Macrophages: bump kNN to 50 because excluding contamination datasets (e.g., Liver_Atlas)
    # halves cell count (~42k -> ~22k) and makes the default k=30 graph too sparse.
    n_neighbors = 50 if cell_type == "macrophages" else 30
    print(f"   HOTSPOT_MAX_GENES={max_genes_env}, HOTSPOT_JOBS={jobs_env}, n_neighbors={n_neighbors}")
    hs, autocorr_df = run_hotspot(
        adata, n_neighbors=n_neighbors, max_genes=max_genes_env, jobs=jobs_env
    )

    print(f"[4/6] Calculating module scores")
    module_scores = hs.calculate_module_scores()
    modules_series = hs.modules
    local_corr = hs.local_correlation_z

    print(f"[5/6] Building output tables")
    rows = []
    for gene, mod in modules_series.items():
        if mod < 0:
            continue
        members = modules_series[modules_series == mod].index
        weight = float(local_corr.loc[gene, members].mean())
        rows.append({"gene": gene, "module": int(mod), "weight": weight})
    module_genes = pd.DataFrame(rows)

    cell_scores_long = (
        module_scores.reset_index()
        .melt(id_vars=module_scores.index.name or "index",
              var_name="module", value_name="score")
        .rename(columns={module_scores.index.name or "index": "cell_id"})
    )
    cell_scores_long["module"] = cell_scores_long["module"].astype(int)

    donor_scores = aggregate_donor(cell_scores_long, adata.obs)

    outdir = hotspot_outdir(cell_type)
    if loo_dataset:
        outdir = outdir / "loo" / loo_dataset
        outdir.mkdir(parents=True, exist_ok=True)
    # exclude_datasets alone overwrites the canonical dir (this is the protocol-remediation
    # production run; backup-restore on disk before invoking). Only nest under exclude/
    # when combined with --loo-dataset (sensitivity/audit runs).

    autocorr_df.to_csv(outdir / "autocorr.tsv", sep="\t", index=False)
    module_genes.to_csv(outdir / "module_genes.tsv", sep="\t", index=False)
    cell_scores_long.to_parquet(outdir / "cell_scores.parquet", index=False)
    donor_scores.to_csv(outdir / "donor_scores.tsv", sep="\t", index=False)

    with open(outdir / "hotspot_obj.pkl", "wb") as fh:
        pickle.dump({"modules": modules_series, "local_corr": local_corr,
                     "module_scores": module_scores}, fh)

    meta = {
        "cell_type": cell_type,
        "loo_dataset": loo_dataset,
        "exclude_datasets": exclude_datasets,
        "smoke": smoke,
        "n_cells": int(adata.n_obs),
        "n_genes_kept": int(adata.n_vars),
        "n_autocorr_genes": int(len(autocorr_df.query("FDR<0.05 & Z>=7"))),
        "n_modules": int(modules_series[modules_series >= 0].nunique()),
        "params": {"n_neighbors": n_neighbors, "autocorr_z": 7.0, "autocorr_fdr": 0.05,
                   "max_genes": max_genes_env, "min_gene_threshold": 20,
                   "fdr_threshold": 0.05, "jobs": jobs_env,
                   "model": "danb"},
    }
    with open(outdir / "run_metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)

    print(f"[6/6] Wrote outputs to {outdir}")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell-type", required=True, choices=RUN_ORDER)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--loo-dataset", default=None)
    ap.add_argument(
        "--exclude-datasets",
        nargs="+",
        default=None,
        help="Exclude one or more datasets from the cell selection (e.g., --exclude-datasets GSE136103 Liver_Atlas). "
             "Mutually compatible with --loo-dataset (both apply). Used for protocol-contamination remediation.",
    )
    args = ap.parse_args()
    main(args.cell_type, args.smoke, args.loo_dataset, exclude_datasets=args.exclude_datasets)
