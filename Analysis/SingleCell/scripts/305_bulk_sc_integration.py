#!/usr/bin/env python3
"""
305: Bulk-SC Integration — Project bulk signatures onto single cells.

Projects bulk RNA-seq derived gene signatures (transition programs, NMF subtypes,
divergence genes) onto single cells, then validates single-cell pseudotime
against bulk-level findings.

Inputs:
    - Cell-type subsets: pseudotime/{celltype}_subset.h5ad
    - Palantir pseudotime: pseudotime/palantir_pseudotime_{celltype}.csv
    - Monocle 3 pseudotime: pseudotime/monocle3_pseudotime_{celltype}.csv
    - Bulk transition gene sets: RNA-seq/results/progression/transition_genesets.csv
    - S1/S2 NMF signatures: RNA-seq/results/subtypes/
    - Divergence genes: RNA-seq/results/progression/divergence_genes.csv

Outputs (to results_gpu_v2/pseudotime/):
    - bulk_signature_scores.csv
    - bulk_sc_concordance.csv
    - s2_substate_mapping.csv

Usage:
    sbatch run_pseudotime_pipeline.sh
"""

import os
import sys
import warnings
import logging
import gc

import numpy as np
import pandas as pd
from scipy import stats, sparse

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
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
PT_DIR = os.path.join(RESULTS, "pseudotime")

BULK_RESULTS = os.path.join(
    BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results"
)
PROGRESSION_DIR = os.path.join(BULK_RESULTS, "progression")
SUBTYPES_DIR = os.path.join(BASE, "RNA-seq/results/subtypes")

CELL_TYPES = [
    "Hepatocytes", "Macrophages", "Fibroblasts",
    "Endothelial_cells", "Cholangiocytes",
]

import scanpy as sc


# ---------------------------------------------------------------------------
# Load bulk signatures
# ---------------------------------------------------------------------------
def load_bulk_signatures():
    """Load all bulk-derived gene signatures."""
    signatures = {}

    # 1. Transition gene sets (Script 147)
    tgs_path = os.path.join(PROGRESSION_DIR, "transition_genesets.csv")
    if os.path.exists(tgs_path):
        tgs = pd.read_csv(tgs_path)
        log.info("Loaded transition gene sets: %d rows", len(tgs))
        # Detect gene column name
        gene_col = "gene_symbol" if "gene_symbol" in tgs.columns else "gene"
        # Group by transition + direction
        for name, grp in tgs.groupby(["transition", "direction"]):
            key = f"transition_{name[0]}_{name[1]}"
            signatures[key] = grp[gene_col].tolist()
    else:
        log.warning("Transition gene sets not found: %s", tgs_path)

    # 2. Divergence genes (Script 117)
    div_path = os.path.join(PROGRESSION_DIR, "divergence_genes.csv")
    if os.path.exists(div_path):
        div_df = pd.read_csv(div_path)
        log.info("Loaded divergence genes: %d genes", len(div_df))
        signatures["divergence_all"] = div_df["gene"].tolist()
        # S1-specific and S2-specific if available
        if "subtype_specific" in div_df.columns:
            for st in div_df["subtype_specific"].unique():
                if pd.notna(st):
                    mask = div_df["subtype_specific"] == st
                    signatures[f"divergence_{st}"] = div_df.loc[mask, "gene"].tolist()
    else:
        log.warning("Divergence genes not found: %s", div_path)

    # 3. NMF subtype signatures
    nmf_path = os.path.join(SUBTYPES_DIR, "nmf_basis_genes.csv")
    if os.path.exists(nmf_path):
        nmf_df = pd.read_csv(nmf_path)
        log.info("Loaded NMF basis genes: %d", len(nmf_df))
        for subtype in nmf_df.columns:
            if subtype != "gene":
                top_genes = nmf_df.nlargest(100, subtype)["gene"].tolist()
                signatures[f"nmf_{subtype}"] = top_genes
    else:
        # Try alternative: load from nmf_results_cache.rds summary
        alt_path = os.path.join(SUBTYPES_DIR, "nmf_assignments.csv")
        if os.path.exists(alt_path):
            log.info("NMF basis not found; using assignments only")

    # 4. Consensus pseudotime gene correlates (Script 116)
    cpt_path = os.path.join(PROGRESSION_DIR, "transport_gene_costs.csv")
    if os.path.exists(cpt_path):
        cpt_df = pd.read_csv(cpt_path)
        log.info("Loaded transport gene costs: %d genes", len(cpt_df))
        # Top transport genes (highest cost = most changed along trajectory)
        if "total_cost" in cpt_df.columns:
            top = cpt_df.nlargest(200, "total_cost")
            signatures["bulk_transport_top200"] = top["gene"].tolist()
    else:
        log.warning("Transport gene costs not found: %s", cpt_path)

    log.info("Total signatures loaded: %d", len(signatures))
    for name, genes in signatures.items():
        log.info("  %s: %d genes", name, len(genes))

    return signatures


# ---------------------------------------------------------------------------
# Score cells with gene signatures
# ---------------------------------------------------------------------------
def score_cells(adata, signatures):
    """Score each cell for each signature using scanpy.tl.score_genes."""
    scores = {}

    for name, genes in signatures.items():
        # Filter to genes present in adata
        available = [g for g in genes if g in adata.var_names]
        if len(available) < 5:
            log.warning("Skipping %s: only %d/%d genes available",
                        name, len(available), len(genes))
            continue

        try:
            sc.tl.score_genes(adata, gene_list=available, score_name=f"sig_{name}")
            scores[name] = adata.obs[f"sig_{name}"].values.copy()
            log.info("Scored %s: %d genes, mean=%.3f", name, len(available),
                     np.mean(scores[name]))
        except Exception as e:
            log.warning("Failed to score %s: %s", name, e)

    return scores


# ---------------------------------------------------------------------------
# Validation: gene overlap
# ---------------------------------------------------------------------------
def compute_gene_overlap(sc_genes, bulk_genes, universe_size=20000):
    """Jaccard index + hypergeometric enrichment."""
    sc_set = set(sc_genes)
    bulk_set = set(bulk_genes)
    intersection = sc_set & bulk_set
    union = sc_set | bulk_set

    jaccard = len(intersection) / len(union) if len(union) > 0 else 0

    # Hypergeometric test
    from scipy.stats import hypergeom
    M = universe_size  # total genes
    n = len(bulk_set)  # drawn from bulk
    N = len(sc_set)    # drawn from sc
    k = len(intersection)  # overlap
    pval = hypergeom.sf(k - 1, M, n, N) if k > 0 else 1.0

    return {
        "jaccard": jaccard,
        "overlap_n": len(intersection),
        "sc_n": len(sc_set),
        "bulk_n": len(bulk_set),
        "hypergeom_pval": pval,
        "fold_enrichment": (k / max(N, 1)) / (n / M) if n > 0 and N > 0 else 0,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("305: Bulk-SC Integration")
    log.info("=" * 60)

    # Load bulk signatures
    signatures = load_bulk_signatures()
    if not signatures:
        log.error("No bulk signatures loaded — exiting")
        return

    all_scores = []
    concordance_results = []

    for ct in CELL_TYPES:
        log.info("\n" + "=" * 50)
        log.info("Processing: %s", ct)
        log.info("=" * 50)

        h5ad_path = os.path.join(PT_DIR, f"{ct}_subset.h5ad")
        if not os.path.exists(h5ad_path):
            log.warning("Subset not found: %s — skipping", h5ad_path)
            continue

        adata = sc.read_h5ad(h5ad_path)
        log.info("Loaded %s: %s", ct, adata.shape)

        # Score cells
        scores = score_cells(adata, signatures)

        # Build scores DataFrame
        score_df = pd.DataFrame(index=adata.obs_names)
        score_df["cell_type"] = ct
        score_df["condition"] = adata.obs["condition"].values
        score_df["leiden_substate"] = adata.obs["leiden_substate"].values
        for name, vals in scores.items():
            score_df[f"sig_{name}"] = vals
        all_scores.append(score_df)

        # Load pseudotime results
        palantir_path = os.path.join(PT_DIR, f"palantir_pseudotime_{ct}.csv")
        dpt_path = os.path.join(PT_DIR, f"dpt_pseudotime_{ct}.csv")
        monocle3_path = os.path.join(PT_DIR, f"monocle3_pseudotime_{ct}.csv")

        pseudotime_sources = {}
        for name, path in [("palantir", palantir_path), ("dpt", dpt_path),
                            ("monocle3", monocle3_path)]:
            if os.path.exists(path):
                pt_df = pd.read_csv(path, index_col=0)
                pt_col = [c for c in pt_df.columns if "pseudotime" in c.lower()]
                if pt_col:
                    pseudotime_sources[name] = pt_df[pt_col[0]]

        # Correlate signatures with pseudotime
        for pt_name, pt_vals in pseudotime_sources.items():
            common = score_df.index.intersection(pt_vals.index)
            if len(common) < 30:
                continue

            pt = pt_vals.loc[common].values
            valid = np.isfinite(pt)

            for sig_name in scores.keys():
                sig_vals = score_df.loc[common, f"sig_{sig_name}"].values
                both_valid = valid & np.isfinite(sig_vals)
                if both_valid.sum() < 30:
                    continue

                rho, pval = stats.spearmanr(
                    pt[both_valid], sig_vals[both_valid]
                )
                concordance_results.append({
                    "cell_type": ct,
                    "pseudotime_method": pt_name,
                    "signature": sig_name,
                    "spearman_rho": rho,
                    "pval": pval,
                    "n_cells": both_valid.sum(),
                })

        # Gene overlap: single-cell trajectory genes vs bulk transport genes
        corr_path = os.path.join(PT_DIR, f"pseudotime_corr_{ct}.csv")
        if os.path.exists(corr_path) and "bulk_transport_top200" in signatures:
            corr_df = pd.read_csv(corr_path)
            if "padj" in corr_df.columns:
                sc_sig_genes = corr_df.loc[
                    (corr_df["padj"] < 0.05) & (corr_df["spearman_rho"].abs() > 0.1),
                    "gene"
                ].tolist()[:500]
                overlap = compute_gene_overlap(
                    sc_sig_genes,
                    signatures["bulk_transport_top200"],
                )
                overlap["cell_type"] = ct
                overlap["comparison"] = "sc_trajectory_vs_bulk_transport"
                concordance_results.append(overlap)
                log.info("Gene overlap (%s): Jaccard=%.3f, p=%.2e, fold=%.1f",
                         ct, overlap["jaccard"], overlap["hypergeom_pval"],
                         overlap["fold_enrichment"])

        del adata
        gc.collect()

    # Save all scores
    if all_scores:
        full_scores = pd.concat(all_scores, axis=0)
        out_path = os.path.join(PT_DIR, "bulk_signature_scores.csv")
        full_scores.to_csv(out_path)
        log.info("Saved bulk signature scores: %s", full_scores.shape)

    # Save concordance
    if concordance_results:
        conc_df = pd.DataFrame(concordance_results)
        out_path = os.path.join(PT_DIR, "bulk_sc_concordance.csv")
        conc_df.to_csv(out_path, index=False)
        log.info("Concordance results:\n%s", conc_df.to_string())

    # S2 substate mapping (if NMF signatures were scored)
    if all_scores:
        full_scores = pd.concat(all_scores, axis=0)
        s2_cols = [c for c in full_scores.columns if "nmf" in c.lower() and "s2" in c.lower()]
        if s2_cols:
            s2_col = s2_cols[0]
            mapping = full_scores.groupby(["cell_type", "leiden_substate"]).agg(
                mean_s2_score=(s2_col, "mean"),
                n_cells=(s2_col, "count"),
            ).reset_index()
            out_path = os.path.join(PT_DIR, "s2_substate_mapping.csv")
            mapping.to_csv(out_path, index=False)

    log.info("\n=== 305: Bulk-SC Integration COMPLETE ===")


if __name__ == "__main__":
    main()
