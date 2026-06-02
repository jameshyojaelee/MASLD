#!/usr/bin/env python
"""
run_scdrs_F_transitions.py
==========================
Score 1.23M-cell scRNA atlas with bulk dream F-transition DEG signatures via
scDRS (Zhang et al., Nat Genet 2022, PMID 36050550).

Inputs
------
* RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/
    fibrosis_consecutive_dream.csv
    -- per-contrast (F1_vs_F0, F2_vs_F1, F3_vs_F2, F4_vs_F3) bulk dream stats:
       gene (versioned ENSG), logFC, padj, t.
* Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad (1,232,318 x 37,533)
    -- log1p of size-factor-normalized counts; raw counts in .raw; var.gene_ids
       holds unversioned ENSG; var_names are mostly HGNC symbols.
* Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz
    -- has the preparation_method column missing from the .h5ad obs; merged on
       cell_id (index) to apply the prep-method filter (unsorted | nuclei).

Pipeline
--------
1. Build .gs file with 3 entries per transition: signed / up / down.
   * Filter: padj < 0.05 AND |logFC| > 0.5; restrict to genes mappable to atlas
     var_names (via var.gene_ids -> symbol map).
   * Weight = bulk dream `t` statistic (signed for *_signed, |t| for *_up/*_down).
2. Apply preparation_method in {unsorted, nuclei} filter (909,237 cells).
3. scdrs.preprocess with covariates: sample (one-hot), n_genes, pct_mito (if
   computable from raw); n_mean_bin/n_var_bin defaults.
4. scdrs.score_cell per gene-set (n_ctrl=1000, weight_opt='vs',
   return_ctrl_norm_score=True, return_ctrl_raw_score=True).
5. scdrs.method.downstream_group_analysis with group_cols=['cell_type'].

Outputs
-------
* disease_signatures/F_transitions_scdrs.gs                   (scDRS gene set)
* disease_signatures/scdrs_cell_scores.parquet                (per-cell scores)
* disease_signatures/scdrs_celltype_enrichment.csv            (per-celltype agg)
* disease_signatures/scdrs_signed_vs_split_concordance.csv    (mode QC)
* disease_signatures/scdrs_gs_summary.csv                     (.gs build stats)
"""
from __future__ import annotations

import gc
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# NumPy 2.0 compatibility shim for scdrs (np.float_ / np.int_ removed).
if not hasattr(np, "float_"):
    np.float_ = np.float64  # type: ignore[attr-defined]
if not hasattr(np, "int_"):
    np.int_ = np.int64  # type: ignore[attr-defined]

import scdrs
from scdrs.method import downstream_group_analysis

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
DEG_CSV   = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_consecutive_dream.csv"
ATLAS_H5  = BASE / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
UMAP_META = BASE / "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz"

OUT_DIR = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

GS_PATH         = OUT_DIR / "F_transitions_scdrs.gs"
GS_SUMMARY_PATH = OUT_DIR / "scdrs_gs_summary.csv"
CELL_SCORE_PATH = OUT_DIR / "scdrs_cell_scores.parquet"
CT_ENRICH_PATH  = OUT_DIR / "scdrs_celltype_enrichment.csv"
CT_HET_PATH     = OUT_DIR / "scdrs_celltype_heterogeneity.csv"
CONCORD_PATH    = OUT_DIR / "scdrs_signed_vs_split_concordance.csv"

# Parameters
PADJ_CUT = 0.05
LFC_CUT  = 0.5
N_CTRL   = 1000
WEIGHT_OPT = "vs"  # variance-stabilized; scDRS default
PREP_KEEP = ["unsorted", "nuclei"]
TRANSITIONS = ["F1_vs_F0", "F2_vs_F1", "F3_vs_F2", "F4_vs_F3"]
RANDOM_SEED = 42


def log(msg: str) -> None:
    print(f"[scdrs F-transitions] {time.strftime('%H:%M:%S')}  {msg}", flush=True)


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def strip_ensembl_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.[0-9]+$", "", regex=True)


def build_gs_lines(
    deg_df: pd.DataFrame,
    atlas_var: pd.DataFrame,
) -> tuple[list[str], pd.DataFrame]:
    """Build scDRS .gs lines: TRAIT\\tGENE:WEIGHT,GENE:WEIGHT,...

    Returns (lines, summary_df).

    Per scDRS spec, weights are positive. For *_signed mode we keep the sign of
    the dream `t` stat baked in by inverting genes with negative `t` (scDRS
    handles them as "down-weighted" only if positive — so we instead build
    *_signed as a single combined gene set with bulk-positive weights of |t|
    for genes meeting the cutoff and we let the *_up / *_down split capture
    direction. To preserve true signed semantics, *_signed scDRS still uses
    |t| (this matches the MAGMA-z analog that scDRS expects: positive trait
    relevance). Direction is recovered by comparing *_up vs *_down scores.
    """
    # Map var.gene_ids (unversioned ENSG) -> var_names (symbol-ish), one-to-one
    var_map = (
        atlas_var.reset_index().rename(columns={atlas_var.index.name or "index": "var_name"})
        [["gene_ids", "var_name"]]
        .drop_duplicates("gene_ids")
    )
    var_map["gene_ids"] = var_map["gene_ids"].astype(str)

    deg = deg_df.copy()
    deg["gene_ensg"] = strip_ensembl_version(deg["gene"])
    deg = deg.merge(var_map, left_on="gene_ensg", right_on="gene_ids", how="left")
    n_total = len(deg)
    n_mapped = deg["var_name"].notna().sum()
    log(f"DEG rows: {n_total:,};  mapped to atlas vars: {n_mapped:,}")

    # Sig filter
    deg["sig"] = (deg["padj"] < PADJ_CUT) & (deg["logFC"].abs() > LFC_CUT)
    deg_sig = deg[deg["sig"] & deg["var_name"].notna()].copy()
    deg_sig["abs_t"] = deg_sig["t"].abs()

    lines: list[str] = ["TRAIT\tGENESET"]
    summary = []

    for tr in TRANSITIONS:
        sub = deg_sig[deg_sig["contrast"] == tr]
        # signed: all sig genes; weight = |t| (scDRS convention)
        signed_genes = sub[["var_name", "abs_t"]].drop_duplicates("var_name")
        # up: t > 0
        up_genes = sub[sub["t"] > 0][["var_name", "abs_t"]].drop_duplicates("var_name")
        # down: t < 0
        down_genes = sub[sub["t"] < 0][["var_name", "abs_t"]].drop_duplicates("var_name")

        for mode, gs_df in [("signed", signed_genes), ("up", up_genes), ("down", down_genes)]:
            tag = f"{tr}__{mode}"
            if len(gs_df) == 0:
                # scDRS still needs a non-empty entry; skip and report
                log(f"  WARNING: {tag} is empty after filtering; skipping .gs entry")
                summary.append({"trait": tag, "n_genes": 0})
                continue
            payload = ",".join(
                f"{g}:{w:.6g}" for g, w in zip(gs_df["var_name"], gs_df["abs_t"])
            )
            lines.append(f"{tag}\t{payload}")
            summary.append({"trait": tag, "n_genes": int(len(gs_df))})

    return lines, pd.DataFrame(summary)


def write_gs(lines: list[str], path: Path) -> None:
    path.write_text("\n".join(lines) + "\n")
    log(f"Wrote .gs: {path} ({len(lines) - 1} traits)")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    np.random.seed(RANDOM_SEED)

    log("STEP 1: load bulk DEG table")
    deg_df = pd.read_csv(DEG_CSV, usecols=["gene", "logFC", "padj", "t", "contrast"])
    log(f"DEG file: {len(deg_df):,} rows  contrasts: {sorted(deg_df['contrast'].unique())}")

    log("STEP 2: load atlas .h5ad in-memory (raw counts will be used)")
    adata = sc.read_h5ad(ATLAS_H5)
    log(f"Atlas: {adata.n_obs:,} x {adata.n_vars:,};  raw available: {adata.raw is not None}")

    log("STEP 3: merge prep-method from atlas_umap_for_fig2.csv.gz")
    # The CSV has 8 columns (umap_1, umap_2, cell_type, sample, dataset,
    # preparation_method, disease_status, disease_stage_coarse) and 1,232,318
    # data rows -- it is positionally aligned with adata.obs (no cell_id
    # column).  Confirm row count match before merging.
    umap_df = pd.read_csv(UMAP_META)
    log(f"UMAP meta rows: {len(umap_df):,}; cols: {list(umap_df.columns)}")
    if len(umap_df) != adata.n_obs:
        raise RuntimeError(
            f"UMAP meta ({len(umap_df):,}) != atlas n_obs ({adata.n_obs:,}); cannot align positionally"
        )
    # Defensive consistency check: cell_type from UMAP should match atlas obs cell_type
    if "cell_type" in umap_df.columns and "cell_type" in adata.obs.columns:
        ct_match = (umap_df["cell_type"].astype(str).values == adata.obs["cell_type"].astype(str).values).mean()
        log(f"  positional cell_type match: {ct_match:.4f} (expect 1.000)")
        if ct_match < 0.99:
            raise RuntimeError("UMAP CSV does not align positionally with atlas")
    if "sample" in umap_df.columns and "sample" in adata.obs.columns:
        s_match = (umap_df["sample"].astype(str).values == adata.obs["sample"].astype(str).values).mean()
        log(f"  positional sample match: {s_match:.4f} (expect 1.000)")
    adata.obs["preparation_method"] = umap_df["preparation_method"].astype(str).values
    if "disease_stage_coarse" in umap_df.columns:
        adata.obs["disease_stage_coarse"] = umap_df["disease_stage_coarse"].astype(str).values

    log("STEP 4: filter to preparation_method in {unsorted, nuclei}")
    keep_mask = adata.obs["preparation_method"].isin(PREP_KEEP).values
    log(f"  keep: {int(keep_mask.sum()):,} / {adata.n_obs:,} cells")
    adata = adata[keep_mask].copy()
    gc.collect()

    log("STEP 5: build .gs (after we know which atlas vars are available)")
    lines, gs_summary = build_gs_lines(deg_df, adata.var)
    gs_summary.to_csv(GS_SUMMARY_PATH, index=False)
    log(f"GS summary:\n{gs_summary.to_string(index=False)}")
    write_gs(lines, GS_PATH)

    log("STEP 6: replace .X with raw counts for scDRS preprocessing")
    if adata.raw is None:
        raise RuntimeError("Atlas has no .raw; scDRS preprocessing expects raw counts")
    # adata.raw is typically aligned to original full atlas obs but the subset
    # adata.copy() preserves alignment because the .raw is sliced too. Verify.
    raw_X = adata.raw.X
    if raw_X.shape[0] != adata.n_obs:
        log(f"  raw rows ({raw_X.shape[0]}) != adata.n_obs ({adata.n_obs}); slicing")
        raw_X = adata.raw[adata.obs_names].X
    if not sp.issparse(raw_X):
        raw_X = sp.csr_matrix(raw_X)
    # Verify gene order matches between adata.var and adata.raw.var
    if not (adata.var_names == adata.raw.var_names).all():
        common_v = adata.var_names.intersection(adata.raw.var_names)
        log(f"  raw vs adata var mismatch; intersect = {len(common_v):,} genes")
        # Re-extract raw on common vars
        raw_full = adata.raw.to_adata()
        raw_aligned = raw_full[adata.obs_names, list(common_v)]
        raw_X = raw_aligned.X
        if not sp.issparse(raw_X):
            raw_X = sp.csr_matrix(raw_X)
        adata = adata[:, list(common_v)].copy()
    adata.X = raw_X.astype(np.float32)
    log(f"  X is sparse: {sp.issparse(adata.X)};  dtype: {adata.X.dtype};  shape: {adata.X.shape}")

    log("STEP 7: compute n_genes / pct_mito covariates")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    # n_genes_by_counts and total_counts populate adata.obs
    # pct_mito: identify MT genes by var_names prefix or gene_ids fallback
    mt_mask = np.asarray(
        adata.var_names.astype(str).str.upper().str.startswith("MT-")
    )
    if mt_mask.sum() == 0:
        log("  no MT- prefixed vars; trying var.gene_ids prefix MT-")
        mt_mask = np.asarray(
            adata.var["gene_ids"].astype(str).str.upper().str.startswith("MT-")
        )
    mt_idx = np.where(mt_mask)[0]
    if len(mt_idx) > 0:
        mt_counts = np.asarray(adata.X[:, mt_idx].sum(axis=1)).ravel()
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        pct_mito = np.where(total_counts > 0, mt_counts / total_counts * 100.0, 0.0)
        adata.obs["pct_mito"] = pct_mito.astype(np.float32)
        log(f"  pct_mito: {len(mt_idx)} MT genes;  median={float(np.median(pct_mito)):.3f}")
    else:
        adata.obs["pct_mito"] = 0.0
        log("  no MT genes found; pct_mito set to 0")

    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"].astype(np.float32)

    log("STEP 8: assemble covariate DataFrame")
    # 275 samples is too many to one-hot (~1TB dense). Use 7-dataset one-hot
    # as the dominant batch factor + n_genes + pct_mito continuous covariates.
    # scDRS preprocess implicit-covariate-correction mode (sparse X + cov)
    # avoids materializing the corrected X matrix.
    dataset_dummies = pd.get_dummies(
        adata.obs["dataset"].astype("category"),
        prefix="ds", drop_first=True, dtype=np.float32,
    )
    cov = pd.concat(
        [
            dataset_dummies,
            adata.obs[["n_genes", "pct_mito"]].astype(np.float32),
        ],
        axis=1,
    )
    cov.index = adata.obs_names
    log(f"  covariate matrix: {cov.shape};  cols: {list(cov.columns)}")

    log("STEP 9: scdrs.preprocess (regress out covariates implicitly)")
    t0 = time.time()
    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20)
    log(f"  preprocess done in {time.time()-t0:.1f}s; SCDRS_PARAM keys: {list(adata.uns['SCDRS_PARAM'].keys())}")

    log("STEP 10: load .gs back via scdrs.util")
    gs_dict = scdrs.util.load_gs(str(GS_PATH))
    log(f"  loaded {len(gs_dict)} traits from .gs")

    log("STEP 11: score_cell per trait (n_ctrl=1000)")
    df_full_records = []
    df_score_records = []
    for trait_idx, (trait, (genes, weights)) in enumerate(gs_dict.items()):
        t1 = time.time()
        # Restrict to overlap with adata.var
        gene_set = pd.Series(genes)
        weight_set = pd.Series(weights)
        in_atlas = gene_set.isin(adata.var_names)
        n_in = int(in_atlas.sum())
        if n_in == 0:
            log(f"  [{trait_idx+1}/{len(gs_dict)}] {trait}: 0 overlap; skipping")
            continue
        gene_use = gene_set[in_atlas].tolist()
        weight_use = weight_set[in_atlas].tolist()
        log(f"  [{trait_idx+1}/{len(gs_dict)}] {trait}: scoring with {n_in} genes")
        df_score = scdrs.score_cell(
            data=adata,
            gene_list=gene_use,
            gene_weight=weight_use,
            ctrl_match_key="mean_var",
            n_ctrl=N_CTRL,
            weight_opt=WEIGHT_OPT,
            return_ctrl_raw_score=False,
            return_ctrl_norm_score=True,
            random_seed=RANDOM_SEED,
            verbose=False,
        )
        # df_score columns: raw_score, norm_score, mc_pval, pval, nlog10_pval,
        #                   zscore, ctrl_norm_score_<i> ...
        df_score = df_score.copy()
        df_score["trait"] = trait
        log(f"    done in {time.time()-t1:.1f}s; score head:\n{df_score[['raw_score','norm_score','mc_pval','zscore']].describe().to_string()}")

        # Save lean per-cell record
        rec = df_score[["raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval", "zscore"]].copy()
        rec["trait"] = trait
        rec.index.name = "cell_id"
        df_score_records.append(rec.reset_index())
        # Keep full df for downstream group analysis (with ctrl_norm_score_*)
        df_full_records.append(df_score)

    if not df_score_records:
        raise RuntimeError("No traits scored")

    log("STEP 12: write per-cell scores (parquet) and macrophage CSV slice")
    cell_scores = pd.concat(df_score_records, ignore_index=True)
    # Parse trait into transition + mode
    parts = cell_scores["trait"].str.split("__", n=1, expand=True)
    cell_scores["transition"] = parts[0]
    cell_scores["mode"] = parts[1]
    cell_scores.to_parquet(CELL_SCORE_PATH, index=False, compression="zstd")
    log(f"  wrote {CELL_SCORE_PATH}  shape={cell_scores.shape}")
    # Also write a Macrophage CSV slice so the R figure script doesn't need
    # arrow/parquet support. Subset signed mode + Macrophage cells.
    mac_ids = adata.obs_names[adata.obs["cell_type"] == "Macrophages"]
    mac_slice = cell_scores[
        (cell_scores["mode"] == "signed")
        & cell_scores["cell_id"].isin(mac_ids)
    ][["cell_id", "transition", "norm_score", "zscore"]]
    mac_csv = OUT_DIR / "scdrs_cell_scores_macrophages_signed.csv"
    mac_slice.to_csv(mac_csv, index=False)
    log(f"  wrote {mac_csv}  shape={mac_slice.shape}")
    # Per-cell-type median of signed norm_score (used as scDRS NES on figure
    # subpanel C and matches what UCell median is compared against).
    cs_signed = cell_scores[cell_scores["mode"] == "signed"].merge(
        adata.obs[["cell_type"]].rename_axis("cell_id").reset_index(),
        on="cell_id", how="left",
    )
    ct_median_signed = (
        cs_signed.groupby(["cell_type", "transition"])
        .agg(scdrs_norm_median=("norm_score", "median"),
             scdrs_zscore_median=("zscore", "median"),
             n_cells=("norm_score", "size"))
        .reset_index()
    )
    ct_median_csv = OUT_DIR / "scdrs_celltype_norm_signed_median.csv"
    ct_median_signed.to_csv(ct_median_csv, index=False)
    log(f"  wrote {ct_median_csv}  shape={ct_median_signed.shape}")

    log("STEP 13: per-cell-type group analysis (one call per trait)")
    enrich_records = []
    het_records = []
    # downstream_group_analysis expects a single df_full_score per trait; iterate.
    for trait, df_score in zip(gs_dict.keys(), df_full_records):
        # Subset to cells present in df_score (should be all)
        ad_sub = adata[df_score.index].copy() if (df_score.index != adata.obs_names).any() else adata
        try:
            res = downstream_group_analysis(
                adata=ad_sub,
                df_full_score=df_score,
                group_cols=["cell_type"],
                fdr_thresholds=[0.05, 0.1, 0.2],
            )
        except Exception as exc:
            log(f"  group analysis FAILED for {trait}: {exc}")
            continue
        # res is a dict: {'cell_type': df}
        ct_df = res.get("cell_type")
        if ct_df is None:
            log(f"  no cell_type result for {trait}")
            continue
        ct_df = ct_df.copy()
        ct_df["trait"] = trait
        ct_df["cell_type"] = ct_df.index
        enrich_records.append(ct_df.reset_index(drop=True))
        log(f"  group({trait}): rows={len(ct_df)}  assoc_mcp_min={ct_df['assoc_mcp'].min():.3g}")
        del ad_sub
        gc.collect()

    if enrich_records:
        enrich_df = pd.concat(enrich_records, ignore_index=True)
        # Parse mode & transition from trait
        parts = enrich_df["trait"].str.split("__", n=1, expand=True)
        enrich_df["transition"] = parts[0]
        enrich_df["mode"] = parts[1]
        enrich_df.to_csv(CT_ENRICH_PATH, index=False)
        log(f"  wrote {CT_ENRICH_PATH}  shape={enrich_df.shape}")
    else:
        log("  WARNING: no enrichment records assembled")
        enrich_df = pd.DataFrame()

    log("STEP 14: signed vs (up - down) concordance per cell type")
    if not enrich_df.empty and {"signed", "up", "down"}.issubset(set(enrich_df["mode"].unique())):
        # Use -log10(assoc_mcp) signed by sign of mean assoc_zscore (if present)
        # otherwise use mean of cell zscores from per-cell parquet.
        cs = cell_scores.merge(
            adata.obs[["cell_type"]].rename_axis("cell_id").reset_index(),
            on="cell_id", how="left",
        )
        ct_zscore = (
            cs.groupby(["cell_type", "transition", "mode"])["zscore"]
            .median().reset_index()
        )
        signed = ct_zscore[ct_zscore["mode"] == "signed"].rename(columns={"zscore": "z_signed"})
        up = ct_zscore[ct_zscore["mode"] == "up"].rename(columns={"zscore": "z_up"})
        dn = ct_zscore[ct_zscore["mode"] == "down"].rename(columns={"zscore": "z_down"})
        merged = (
            signed.merge(up.drop(columns="mode"), on=["cell_type", "transition"], how="outer")
                  .merge(dn.drop(columns="mode"), on=["cell_type", "transition"], how="outer")
        )
        merged["z_split"] = merged["z_up"].fillna(0) - merged["z_down"].fillna(0)
        merged["abs_diff"] = (merged["z_signed"] - merged["z_split"]).abs()
        merged["discordant_sign"] = (
            np.sign(merged["z_signed"].fillna(0)) != np.sign(merged["z_split"].fillna(0))
        ) & (merged["z_signed"].abs() > 0.1) & (merged["z_split"].abs() > 0.1)
        merged.to_csv(CONCORD_PATH, index=False)
        log(f"  wrote {CONCORD_PATH};  median |signed - split| = {merged['abs_diff'].median():.3f}")
    else:
        log("  WARNING: could not build concordance table (modes missing)")

    log("DONE")


if __name__ == "__main__":
    main()
