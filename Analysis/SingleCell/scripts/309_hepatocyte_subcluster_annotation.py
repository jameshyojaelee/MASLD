#!/usr/bin/env python3
"""
309: Hepatocyte Subcluster Annotation

Select optimal Leiden resolution, score biological axes, compute disease
enrichment, assign provisional subtype labels.

Input:
    hepatocyte_subtypes/hepatocyte_atlas.h5ad  (from Script 308)
    hepatocyte_subtypes/leiden_resolution_scores.csv
    data/GSE244832/metadata/donor_pairing.csv  (for GSE244832 per-donor staging)

Outputs (to hepatocyte_subtypes/):
    hepatocyte_atlas_annotated.h5ad
    hepatocyte_subtype_metadata.csv
    subtype_axis_scores.csv
    subtype_disease_enrichment.csv
    subtype_markers.csv                              VIZ_ONLY (do not publish)
    pseudobulk_markers/sample_metadata.csv
    pseudobulk_markers/cells_per_donor_subtype.csv
    pseudobulk_markers/subtype_{N}_pseudobulk.csv    raw-count input for limma-voom
    pseudobulk_markers/subtype_{N}_de.csv            PUBLISHED markers (written by
                                                     309b_pseudobulk_hepatocyte_markers.R)

Pachter / Squair 2021 P0 fix: subtype markers come from donor × subtype
pseudobulk + limma-voom (n=donors). The rank_genes_groups call at L319 is
retained for VIZ_ONLY label assignment but its output (subtype_markers.csv)
is anti-conservative (5–10× FPR per Squair et al. 2021 Nat Commun) and
MUST NOT enter the manuscript. Downstream scripts (S3 cross-modal, S4
figures, atlas integration) read pseudobulk_markers/subtype_{N}_de.csv.

Environment: spatial (CPU) for this Python script; rnaseq for the R step.
"""

import logging
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import fisher_exact, kruskal, spearmanr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SC_DIR = os.path.join(BASE, "Analysis/SingleCell")
RESULTS = os.path.join(SC_DIR, "results_gpu_v2")
OUT_DIR = os.path.join(RESULTS, "hepatocyte_subtypes")
DONOR_PAIRING_PATH = os.path.join(
    BASE, "data/GSE244832/metadata/donor_pairing.csv"
)

# ---------------------------------------------------------------------------
# Coarse disease stage mapping (from scRNA atlas condition labels)
# ---------------------------------------------------------------------------
CONDITION_TO_STAGE = {
    "Healthy": ("Healthy", 0),
    "NAFLD": ("Steatosis", 1),
    "NASH": ("Steatohepatitis", 2),
    "Cirrhotic": ("Cirrhosis", 3),
    "MASLD": (np.nan, np.nan),  # Unspecified — try to recover from GSE244832
}

GSE244832_CONDITION_MAP = {
    "NORMAL": ("Healthy", 0),
    "MASL": ("Steatosis", 1),
    "MASH": ("Steatohepatitis", 2),
}

# ---------------------------------------------------------------------------
# Biological axis gene sets
# ---------------------------------------------------------------------------
AXIS_GENES = {
    "periportal": ["ASS1", "HAL", "CPS1", "ALDOB", "PCK1", "SLC1A2"],
    "pericentral": ["GLUL", "CYP2E1", "CYP1A2", "CYP3A4", "AXIN2"],
    "lipid_accumulation": ["FASN", "SCD", "ACACA", "DGAT1", "DGAT2", "PLIN2"],
    "fatty_acid_oxidation": ["ACOX1", "HADHA", "CPT1A", "ACADM", "EHHADH"],
    "er_stress": ["DDIT3", "ATF4", "XBP1", "HSPA5", "ATF6"],
    "inflammatory": ["CXCL8", "IL6", "TNF", "NFKB1", "CCL2"],
    "ferroptosis": ["GPX4", "SLC7A11", "ACSL4", "PTGS2", "HMOX1"],
}


def write_pseudobulk_inputs(adata, subtypes, out_dir):
    """Aggregate raw counts per (sample, hepatocyte_subtype) and write a
    {Subtype}_pseudobulk.csv (genes × samples) + sample_metadata.csv for the
    downstream R limma-voom step (309b_pseudobulk_hepatocyte_markers.R).

    Squair et al. 2021 (Nat Commun) → marker discovery on n=cells is
    anti-conservative (5–10× FPR). This builds the n=donors design.

    Raw counts are pulled from `adata.layers["counts"]` if present, else
    `adata.raw.X` (preserved through Script 308 → 309). adata.X may be
    log-normalized so it MUST NOT be used.
    """
    # Locate raw-count matrix.
    if "counts" in adata.layers:
        X = adata.layers["counts"]
        var_names = adata.var_names
        log.info("Pseudobulk source: adata.layers['counts'] (%s)", type(X).__name__)
    elif adata.raw is not None:
        X = adata.raw.X
        var_names = adata.raw.var_names
        log.info("Pseudobulk source: adata.raw.X (%s)", type(X).__name__)
    else:
        # Fallback: log-normalized X is NOT a valid pseudobulk input but
        # log a loud warning and proceed so the pipeline doesn't deadlock.
        log.warning("No .layers['counts'] or .raw — falling back to adata.X. "
                   "Pseudobulk counts WILL BE INCORRECT if X is normalized. "
                   "Re-run Script 308 with raw-count preservation.")
        X = adata.X
        var_names = adata.var_names

    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    # Cell-level keys
    samples = adata.obs["sample"].astype(str).values
    subtype_vec = adata.obs["hepatocyte_subtype"].astype(str).values
    datasets = adata.obs["dataset"].astype(str).values if "dataset" in adata.obs.columns else np.array(["unknown"] * len(samples))
    conditions = adata.obs["condition_binary"].astype(str).values if "condition_binary" in adata.obs.columns else adata.obs["condition"].astype(str).values

    # Build sample-level metadata table (one row per sample).
    sample_meta_rows = {}
    for s, ds, cond in zip(samples, datasets, conditions):
        if s not in sample_meta_rows:
            sample_meta_rows[s] = (ds, cond)
    sample_meta_df = pd.DataFrame(
        [(s, ds, cond) for s, (ds, cond) in sample_meta_rows.items()],
        columns=["sample", "dataset", "condition"],
    )
    sample_meta_df.to_csv(os.path.join(out_dir, "sample_metadata.csv"), index=False)
    log.info("  Wrote sample_metadata.csv (%d samples)", len(sample_meta_df))

    # Cell-count contribution table (donor × subtype) for QC.
    cell_counts = (
        pd.DataFrame({"sample": samples, "subtype": subtype_vec})
        .groupby(["sample", "subtype"]).size()
        .unstack(fill_value=0)
    )
    cell_counts.to_csv(os.path.join(out_dir, "cells_per_donor_subtype.csv"))
    log.info("  Wrote cells_per_donor_subtype.csv (%s)", cell_counts.shape)

    # Per subtype: aggregate raw counts across cells per donor → genes × donors.
    for st in subtypes:
        st_mask = subtype_vec == st
        if st_mask.sum() < 100:
            log.warning("  Subtype %s: only %d cells — skipping pseudobulk", st, st_mask.sum())
            continue

        st_samples = samples[st_mask]
        unique_donors = sorted(set(st_samples))

        # Build sample × cell membership matrix → multiply with X to sum counts.
        # n_donors × n_cells boolean → sparse for efficiency.
        donor_to_idx = {d: i for i, d in enumerate(unique_donors)}
        row_ind = np.array([donor_to_idx[s] for s in st_samples])
        col_ind = np.where(st_mask)[0]
        data = np.ones(st_mask.sum(), dtype=np.float32)
        M = sparse.csr_matrix(
            (data, (row_ind, np.arange(len(col_ind)))),
            shape=(len(unique_donors), len(col_ind)),
        )
        # donor × gene = (donor × cell) @ (cell × gene)
        donor_gene = M @ X[col_ind, :]
        donor_gene = np.asarray(donor_gene.todense())  # n_donors × n_genes
        # Cast to integer counts (raw counts should already be integer-valued).
        donor_gene = np.rint(donor_gene).astype(np.int64)

        out_df = pd.DataFrame(
            donor_gene.T,
            index=var_names,
            columns=unique_donors,
        )
        out_df.index.name = "gene"
        out_file = os.path.join(out_dir, f"subtype_{st}_pseudobulk.csv")
        out_df.to_csv(out_file)
        log.info("  Subtype %s: %d cells → %d donors × %d genes → %s",
                 st, st_mask.sum(), len(unique_donors), len(var_names), out_file)


def main():
    log.info("=" * 60)
    log.info("309: Hepatocyte Subcluster Annotation")
    log.info("=" * 60)

    # ── 1. Load data ─────────────────────────────────────────────────────
    h5ad_path = os.path.join(OUT_DIR, "hepatocyte_atlas.h5ad")
    log.info("Loading %s ...", h5ad_path)
    adata = sc.read_h5ad(h5ad_path)
    log.info("Shape: %s", adata.shape)

    sil_df = pd.read_csv(os.path.join(OUT_DIR, "leiden_resolution_scores.csv"))
    log.info("Resolution scores:\n%s", sil_df.to_string(index=False))

    # ── 2. Select optimal Leiden resolution ──────────────────────────────
    n_total = len(adata)
    min_cluster_frac = 0.01  # 1%
    max_cluster_frac = 0.60  # 60%

    best_res = None
    best_sil = -1

    for _, row in sil_df.iterrows():
        res = row["resolution"]
        key = f"leiden_{res}"
        cluster_sizes = adata.obs[key].value_counts()
        min_frac = cluster_sizes.min() / n_total
        max_frac = cluster_sizes.max() / n_total
        n_clusters = row["n_clusters"]

        passes = (
            min_frac >= min_cluster_frac
            and max_frac <= max_cluster_frac
            and 6 <= n_clusters <= 10
        )

        log.info(
            "  res=%.1f: %d clusters, min_frac=%.3f, max_frac=%.3f, sil=%.4f → %s",
            res, n_clusters, min_frac, max_frac, row["silhouette"],
            "PASS" if passes else "FAIL",
        )

        if passes and row["silhouette"] > best_sil:
            best_sil = row["silhouette"]
            best_res = res

    # Fallback: if none pass all criteria, pick best silhouette with relaxed constraints
    if best_res is None:
        log.warning("No resolution passed all constraints. Relaxing cluster count range.")
        for _, row in sil_df.iterrows():
            res = row["resolution"]
            key = f"leiden_{res}"
            cluster_sizes = adata.obs[key].value_counts()
            min_frac = cluster_sizes.min() / n_total
            max_frac = cluster_sizes.max() / n_total
            passes = min_frac >= min_cluster_frac and max_frac <= max_cluster_frac
            if passes and row["silhouette"] > best_sil:
                best_sil = row["silhouette"]
                best_res = res

    if best_res is None:
        best_res = sil_df.loc[sil_df["silhouette"].idxmax(), "resolution"]
        log.warning("Using highest-silhouette resolution as final fallback: %.1f", best_res)

    leiden_key = f"leiden_{best_res}"
    log.info("Selected resolution: %.1f (%s)", best_res, leiden_key)
    adata.obs["hepatocyte_subtype"] = adata.obs[leiden_key].astype(str)
    subtypes = sorted(adata.obs["hepatocyte_subtype"].unique(), key=lambda x: int(x))
    log.info("Subtypes: %s", subtypes)

    # ── 3. Biological axis scoring ───────────────────────────────────────
    log.info("Scoring biological axes...")

    # Ensure we're working with log-normalized data for scoring
    # adata.X should already be log-normalized from 308 (post-HVG normalization)
    # But sc.tl.score_genes needs all genes, not just HVGs
    axis_scores = {}
    for axis_name, genes in AXIS_GENES.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            log.warning("Axis '%s': only %d/%d genes present — skipping", axis_name, len(present), len(genes))
            continue
        sc.tl.score_genes(adata, gene_list=present, score_name=f"axis_{axis_name}")
        log.info("  %s: %d/%d genes present", axis_name, len(present), len(genes))

    # Compute cluster-level mean axis scores
    axis_cols = [c for c in adata.obs.columns if c.startswith("axis_")]
    cluster_axis = adata.obs.groupby("hepatocyte_subtype")[axis_cols].mean()
    cluster_axis.index.name = "subtype"
    axis_path = os.path.join(OUT_DIR, "subtype_axis_scores.csv")
    cluster_axis.to_csv(axis_path)
    log.info("Saved axis scores to %s", axis_path)

    # ── 4. Disease enrichment ────────────────────────────────────────────
    log.info("Computing disease enrichment...")

    # Map condition to binary
    cond = adata.obs["condition"].copy()
    # Treat NAFLD and NASH as MASLD
    cond = cond.replace({"NAFLD": "MASLD", "NASH": "MASLD", "Cirrhotic": "MASLD"})
    adata.obs["condition_binary"] = cond

    # Per-cluster Fisher's exact (MASLD vs Healthy)
    enrichment_rows = []
    n_masld_total = (cond == "MASLD").sum()
    n_healthy_total = (cond == "Healthy").sum()

    for st in subtypes:
        mask = adata.obs["hepatocyte_subtype"] == st
        n_masld_in = ((cond == "MASLD") & mask).sum()
        n_healthy_in = ((cond == "Healthy") & mask).sum()
        n_masld_out = n_masld_total - n_masld_in
        n_healthy_out = n_healthy_total - n_healthy_in

        table = [[n_masld_in, n_healthy_in], [n_masld_out, n_healthy_out]]
        odds_ratio, pval = fisher_exact(table)

        enrichment_rows.append({
            "subtype": st,
            "n_cells": int(mask.sum()),
            "n_masld": int(n_masld_in),
            "n_healthy": int(n_healthy_in),
            "odds_ratio": odds_ratio,
            "fisher_pval": pval,
        })

    enrich_df = pd.DataFrame(enrichment_rows)

    # FDR correction
    _, enrich_df["fisher_padj"], _, _ = multipletests(enrich_df["fisher_pval"], method="fdr_bh")

    # Classify enrichment
    enrich_df["enrichment_class"] = "neutral"
    enrich_df.loc[
        (enrich_df["odds_ratio"] > 1.5) & (enrich_df["fisher_padj"] < 0.05),
        "enrichment_class",
    ] = "MASLD_enriched"
    enrich_df.loc[
        (enrich_df["odds_ratio"] < 0.67) & (enrich_df["fisher_padj"] < 0.05),
        "enrichment_class",
    ] = "Healthy_enriched"

    log.info("Enrichment:\n%s", enrich_df[["subtype", "odds_ratio", "fisher_padj", "enrichment_class"]].to_string(index=False))

    # ── 5. Coarse disease stage association ─────────────────────────────────
    log.info("Computing coarse disease stage association...")

    # Build SRR -> (stage_label, stage_numeric) lookup for GSE244832
    gse244832_srr_map = {}
    if os.path.exists(DONOR_PAIRING_PATH):
        donor_df = pd.read_csv(DONOR_PAIRING_PATH)
        for _, row in donor_df.iterrows():
            cond = str(row["condition"]).strip()
            if cond in GSE244832_CONDITION_MAP:
                stage_label, stage_num = GSE244832_CONDITION_MAP[cond]
                for srr in str(row["rna_srrs"]).split(";"):
                    srr = srr.strip()
                    if srr:
                        gse244832_srr_map[srr] = (stage_label, stage_num)
        log.info("GSE244832 donor pairing: %d SRR -> stage mappings", len(gse244832_srr_map))
    else:
        log.warning("GSE244832 donor_pairing.csv not found at %s", DONOR_PAIRING_PATH)

    # Map each cell to coarse disease stage via vectorized lookup
    stage_label_list = []
    stage_numeric_list = []
    for sample_id, dataset, condition in zip(
        adata.obs["sample"], adata.obs["dataset"], adata.obs["condition"]
    ):
        # GSE244832: recover per-donor stage from donor_pairing
        if dataset == "GSE244832" and sample_id in gse244832_srr_map:
            lbl, num = gse244832_srr_map[sample_id]
        elif condition in CONDITION_TO_STAGE:
            lbl, num = CONDITION_TO_STAGE[condition]
        else:
            lbl, num = np.nan, np.nan
        stage_label_list.append(lbl)
        stage_numeric_list.append(num)

    adata.obs["disease_stage_coarse"] = pd.Categorical(
        stage_label_list,
        categories=["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"],
        ordered=True,
    )
    # Plain numpy float64 (NaN for missing); pd.Float64 nullable type not serialisable by anndata HDF5 writer
    adata.obs["disease_stage_numeric"] = np.array(stage_numeric_list, dtype=float)

    has_stage = adata.obs["disease_stage_numeric"].notna()
    log.info(
        "Cells with disease stage: %d / %d (%.1f%%)",
        has_stage.sum(), len(adata), 100 * has_stage.mean(),
    )
    log.info(
        "Stage distribution:\n%s",
        adata.obs.loc[has_stage, "disease_stage_coarse"].value_counts().to_string(),
    )

    if has_stage.sum() > 100:
        stage_obs = adata.obs[has_stage].copy()

        # Sample-level proportions
        sample_props = (
            stage_obs.groupby(["sample", "hepatocyte_subtype"]).size()
            .unstack(fill_value=0)
        )
        sample_props = sample_props.div(sample_props.sum(axis=1), axis=0)

        # Get sample-level disease stage
        sample_stage_series = stage_obs.groupby("sample")["disease_stage_numeric"].first()

        for st in subtypes:
            if st not in sample_props.columns:
                continue
            prop = sample_props[st]
            stage = sample_stage_series.loc[prop.index]

            # Drop NAs for correlation
            valid = stage.notna()
            rho, rho_p = spearmanr(stage[valid], prop[valid])
            enrich_df.loc[enrich_df["subtype"] == st, "stage_spearman_rho"] = rho
            enrich_df.loc[enrich_df["subtype"] == st, "stage_spearman_p"] = rho_p

            # Kruskal-Wallis across stages (drop NA from unique stages)
            valid_stages = sorted([s for s in stage.unique() if pd.notna(s)])
            stage_groups = [prop[stage == s].values for s in valid_stages if (stage == s).sum() >= 3]
            if len(stage_groups) >= 3:
                kw_stat, kw_p = kruskal(*stage_groups)
                enrich_df.loc[enrich_df["subtype"] == st, "stage_kw_stat"] = kw_stat
                enrich_df.loc[enrich_df["subtype"] == st, "stage_kw_pval"] = kw_p

    enrich_path = os.path.join(OUT_DIR, "subtype_disease_enrichment.csv")
    enrich_df.to_csv(enrich_path, index=False)
    log.info("Saved enrichment to %s", enrich_path)

    # ── 6. Marker genes ──────────────────────────────────────────────────
    # VIZ_ONLY: rank_genes_groups is used here for provisional label assignment
    # (top-axis lookup, UMAP coloring, dotplot) ONLY. It is the n=cells fallacy
    # for DE inference (Squair et al. 2021, Nat Commun: 5–10× FPR vs pseudobulk).
    # The PUBLISHED hepatocyte subtype markers come from the pseudobulk + limma-voom
    # retrofit at Section 6b below (pseudobulk_markers/{Subtype}.csv).
    # Downstream consumers (S3 cross-modal, S4 figure builders) MUST read
    # pseudobulk_markers/, NOT subtype_markers.csv. subtype_markers.csv is
    # retained as a legacy/UMAP-overlay artifact, not a publishable table.
    log.info("Computing marker genes (Wilcoxon) — VIZ_ONLY for label assignment...")
    sc.tl.rank_genes_groups(adata, groupby="hepatocyte_subtype", method="wilcoxon", pts=True)

    marker_rows = []
    for st in subtypes:
        names = sc.get.rank_genes_groups_df(adata, group=st)
        names = names[names["pvals_adj"] < 0.05].head(50)
        names["subtype"] = st
        marker_rows.append(names)

    markers_df = pd.concat(marker_rows, ignore_index=True)
    markers_path = os.path.join(OUT_DIR, "subtype_markers.csv")
    markers_df.to_csv(markers_path, index=False)
    log.info("Saved VIZ_ONLY markers to %s (%d rows) — DO NOT publish; see pseudobulk_markers/", markers_path, len(markers_df))

    # ── 6b. Pseudobulk markers (PUBLISHED; Pachter P0 Squair-2021 fix) ───
    # Donor × subtype raw-count pseudobulk + limma-voom. Subtype-vs-rest
    # contrast per subtype. Filter padj<0.05, |logFC|>0.5. n=donors, not n=cells.
    log.info("Writing pseudobulk inputs for limma-voom (Pachter P0 retrofit)...")
    pseudobulk_dir = os.path.join(OUT_DIR, "pseudobulk_markers")
    os.makedirs(pseudobulk_dir, exist_ok=True)
    write_pseudobulk_inputs(adata, subtypes, pseudobulk_dir)
    log.info("Pseudobulk inputs written. Run 309b_pseudobulk_hepatocyte_markers.R "
            "to produce per-subtype limma-voom DE tables.")

    # ── 7. Assign provisional labels ─────────────────────────────────────
    log.info("Assigning provisional labels...")

    for st in subtypes:
        # Determine top axis
        if st in cluster_axis.index:
            top_axis = cluster_axis.loc[st].idxmax().replace("axis_", "")
        else:
            top_axis = "unknown"

        # Get enrichment class
        ec = enrich_df.loc[enrich_df["subtype"] == st, "enrichment_class"].values[0]

        label = f"Hep_{st}_{ec}_{top_axis}"
        adata.obs.loc[adata.obs["hepatocyte_subtype"] == st, "hepatocyte_subtype_label"] = label
        log.info("  Cluster %s → %s", st, label)

    # ── 8. Save outputs ──────────────────────────────────────────────────
    # Per-cell metadata CSV
    meta_cols = (
        ["sample", "dataset", "condition", "condition_binary",
         "hepatocyte_subtype", "hepatocyte_subtype_label",
         "disease_stage_coarse", "disease_stage_numeric"]
        + axis_cols
    )
    meta_out = adata.obs[meta_cols].copy()
    meta_out["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    meta_out["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    meta_path = os.path.join(OUT_DIR, "hepatocyte_subtype_metadata.csv")
    meta_out.to_csv(meta_path)
    log.info("Saved metadata to %s", meta_path)

    # Annotated h5ad
    out_h5ad = os.path.join(OUT_DIR, "hepatocyte_atlas_annotated.h5ad")
    adata.write_h5ad(out_h5ad)
    log.info("Saved annotated atlas to %s", out_h5ad)

    log.info("Done.")


if __name__ == "__main__":
    main()
