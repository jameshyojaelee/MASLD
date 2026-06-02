#!/usr/bin/env python3
"""
16h_commot_validation.py — Cross-dataset COMMOT validation using Vu et al.

Runs COMMOT on Vu et al. (JHEP Reports 2025) arrays using the same filtered
LR database and parameters from GSE192741. Compares pathway-level signaling
scores and signaling-regulated genes for independent replication.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=12:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config,
    save_checkpoint, save_csv, print_header, print_step,
)


def load_vu_data():
    """Load Vu et al. preprocessed spatial data."""
    path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"
    if not path.exists():
        print(f"  ERROR: {path} not found. Run 10_vu_et_al_validation.py first.")
        sys.exit(1)
    adata = sc.read_h5ad(path)
    print(f"  Loaded Vu et al.: {adata.n_obs} spots × {adata.n_vars} genes, "
          f"{adata.obs['sample_id'].nunique()} arrays")
    return adata


def load_filtered_lr_database():
    """Load filtered LR database from 16b."""
    path = RESULTS_DIR / "commot" / "lr_database_cellchat_filtered.csv"
    if not path.exists():
        print(f"  ERROR: {path} not found. Run 16b_prepare_commot.py first.")
        sys.exit(1)
    df = pd.read_csv(path, index_col=0)
    print(f"  Loaded filtered LR database: {len(df)} pairs")
    return df


def load_gse192741_pathway_scores():
    """Load GSE192741 differential communication results from 16f."""
    path = RESULTS_DIR / "commot" / "differential_communication.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def load_gse192741_signaling_genes():
    """Load GSE192741 signaling-regulated genes from 16e."""
    path = RESULTS_DIR / "commot" / "signaling_regulated_genes_summary.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def run_commot_on_vu(adata, lr_db, commot_config):
    """Run COMMOT on Vu et al. data per array."""
    import commot as ct

    arrays = sorted(adata.obs["sample_id"].unique())
    results_per_array = {}

    for i, array_id in enumerate(arrays):
        print_step(f"COMMOT for {array_id}", i + 1, len(arrays))
        sub = adata[adata.obs["sample_id"] == array_id].copy()

        if sub.n_obs < 100:
            print(f"    WARNING: {sub.n_obs} spots, skipping")
            continue

        # Ensure spatial coordinates
        if "spatial" not in sub.obsm:
            print(f"    WARNING: No spatial coordinates for {array_id}, skipping")
            continue

        # Ensure normalized expression
        if sub.X.max() > 50:
            sc.pp.normalize_total(sub, target_sum=1e4)
            sc.pp.log1p(sub)

        dis_thr = commot_config.get("dis_thr", 500)

        try:
            ct.tl.spatial_communication(
                sub,
                database_name="CellChat",
                df_ligrec=lr_db,
                dis_thr=dis_thr,
                heteromeric=True,
                pathway_sum=True,
            )
            results_per_array[array_id] = sub
            commot_keys = [k for k in sub.obsp.keys() if "commot" in k.lower()]
            print(f"    Success: {len(commot_keys)} communication matrices")
        except Exception as e:
            print(f"    WARNING: COMMOT failed for {array_id}: {e}")

    return results_per_array


def extract_array_pathway_scores(results_per_array):
    """Extract pathway-level scores per array."""
    records = []
    for array_id, adata_sub in results_per_array.items():
        for key in adata_sub.obsm.keys():
            if "commot" not in key.lower() or "sum" not in key:
                continue
            parts = key.split("-")
            if len(parts) < 4:
                continue
            pathway = parts[2]
            direction = parts[-1]

            scores = np.asarray(adata_sub.obsm[key])
            if scores.ndim > 1:
                total = float(np.nansum(scores))
                mean = float(np.nanmean(scores))
            else:
                total = float(np.nansum(scores))
                mean = float(np.nanmean(scores))

            records.append({
                "array_id": array_id,
                "pathway": pathway,
                "direction": direction,
                "total_score": total,
                "mean_score": mean,
                "n_spots": adata_sub.n_obs,
            })

    return pd.DataFrame(records)


def cross_dataset_pathway_correlation(vu_scores, gse_scores):
    """Compare pathway-level signaling between Vu and GSE192741."""
    results = []

    if len(vu_scores) == 0 or len(gse_scores) == 0:
        return pd.DataFrame()

    # Aggregate Vu scores to pathway level (mean across arrays)
    vu_agg = vu_scores.groupby("pathway").agg(
        mean_total=("total_score", "mean"),
        mean_mean=("mean_score", "mean"),
        n_arrays=("array_id", "nunique"),
    ).reset_index()

    # Get pathway-level scores from GSE192741
    # gse_scores has one row per pathway with various score columns
    score_cols = [c for c in gse_scores.columns if "total" in c.lower()]

    shared_pathways = sorted(
        set(vu_agg["pathway"]) & set(gse_scores["pathway"]
                                      if "pathway" in gse_scores.columns
                                      else gse_scores.index)
    )

    if len(shared_pathways) < 3:
        print(f"    Only {len(shared_pathways)} shared pathways, correlation unreliable")
        results.append({
            "metric": "n_shared_pathways",
            "value": len(shared_pathways),
            "note": "Too few for correlation",
        })
        return pd.DataFrame(results)

    # Build matched vectors
    vu_vals = []
    gse_vals = []
    for pw in shared_pathways:
        vu_row = vu_agg[vu_agg["pathway"] == pw]
        if len(vu_row) == 0:
            continue

        if "pathway" in gse_scores.columns:
            gse_row = gse_scores[gse_scores["pathway"] == pw]
        else:
            gse_row = gse_scores.loc[[pw]] if pw in gse_scores.index else pd.DataFrame()

        if len(gse_row) == 0:
            continue

        vu_vals.append(float(vu_row["mean_total"].iloc[0]))
        # Use first available score column from GSE
        for col in score_cols:
            if col in gse_row.columns and not pd.isna(gse_row[col].iloc[0]):
                gse_vals.append(float(gse_row[col].iloc[0]))
                break
        else:
            # Fallback to any numeric column
            gse_vals.append(0.0)

    if len(vu_vals) >= 3:
        rho, pval = spearmanr(vu_vals, gse_vals)
        results.append({
            "metric": "pathway_score_spearman_rho",
            "value": rho,
            "note": f"n={len(vu_vals)} pathways",
        })
        results.append({
            "metric": "pathway_score_spearman_pval",
            "value": pval,
            "note": "",
        })
        print(f"    Pathway score correlation: rho={rho:.3f}, p={pval:.3e} "
              f"({len(vu_vals)} pathways)")

    results.append({
        "metric": "n_shared_pathways",
        "value": len(shared_pathways),
        "note": "",
    })
    results.append({
        "metric": "n_vu_pathways",
        "value": vu_agg["pathway"].nunique(),
        "note": "",
    })

    return pd.DataFrame(results)


def cross_dataset_gene_replication(results_per_array, gse_sig_genes):
    """Compare signaling-regulated genes between Vu and GSE192741."""
    if len(gse_sig_genes) == 0:
        return pd.DataFrame()

    gene_col = "gene" if "gene" in gse_sig_genes.columns else gse_sig_genes.index.name
    if gene_col and gene_col in gse_sig_genes.columns:
        gse_gene_set = set(gse_sig_genes[gene_col].dropna())
    else:
        gse_gene_set = set(gse_sig_genes.index)

    # Collect genes present in Vu arrays
    vu_all_genes = set()
    for adata_sub in results_per_array.values():
        vu_all_genes.update(adata_sub.var_names)

    overlap = gse_gene_set & vu_all_genes
    results = [{
        "metric": "gse192741_signaling_genes",
        "value": len(gse_gene_set),
    }, {
        "metric": "vu_total_genes",
        "value": len(vu_all_genes),
    }, {
        "metric": "shared_genes",
        "value": len(overlap),
    }, {
        "metric": "jaccard_gene_overlap",
        "value": len(overlap) / max(len(gse_gene_set | vu_all_genes), 1),
    }]

    return pd.DataFrame(results)


def main():
    print_header("16h: COMMOT Validation (Vu et al.)")

    config = load_config()
    commot_config = config["commot"]
    output_dir = RESULTS_DIR / "commot" / "validation_vu"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load data ──
    print_step("Loading data", 1, 4)
    adata = load_vu_data()
    lr_db = load_filtered_lr_database()

    # ── Run COMMOT on Vu arrays ──
    print_step("Running COMMOT on Vu arrays", 2, 4)
    results_per_array = run_commot_on_vu(adata, lr_db, commot_config)
    print(f"\n  Successfully processed: {len(results_per_array)} arrays")

    if not results_per_array:
        print("  ERROR: No arrays processed successfully")
        sys.exit(1)

    # ── Extract and save per-array results ──
    print_step("Extracting pathway scores", 3, 4)
    vu_scores = extract_array_pathway_scores(results_per_array)
    save_csv(vu_scores, "vu_commot_results.csv", subdir="commot/validation_vu")

    # Pathway summary across arrays
    if len(vu_scores) > 0:
        pw_summary = vu_scores.groupby("pathway").agg(
            mean_total=("total_score", "mean"),
            std_total=("total_score", "std"),
            n_arrays=("array_id", "nunique"),
        ).sort_values("mean_total", ascending=False)
        print(f"\n  Top pathways in Vu et al.:")
        for pw, row in pw_summary.head(10).iterrows():
            print(f"    {pw}: mean={row['mean_total']:.2f} "
                  f"(sd={row['std_total']:.2f}, n={int(row['n_arrays'])} arrays)")

    # ── Cross-dataset validation ──
    print_step("Cross-dataset validation", 4, 4)

    # Pathway correlation
    print("\n  Pathway-level correlation:")
    gse_scores = load_gse192741_pathway_scores()
    pw_corr = cross_dataset_pathway_correlation(vu_scores, gse_scores)
    if len(pw_corr) > 0:
        save_csv(pw_corr, "cross_dataset_pathway_correlation.csv",
                 subdir="commot/validation_vu")

    # Gene replication
    print("\n  Signaling-regulated gene replication:")
    gse_sig = load_gse192741_signaling_genes()
    gene_repl = cross_dataset_gene_replication(results_per_array, gse_sig)
    if len(gene_repl) > 0:
        save_csv(gene_repl, "replicated_signaling_genes.csv",
                 subdir="commot/validation_vu")
        for _, row in gene_repl.iterrows():
            val = row["value"]
            if isinstance(val, float) and val < 1:
                print(f"    {row['metric']}: {val:.3f}")
            else:
                print(f"    {row['metric']}: {val}")

    # Top pathway ranking concordance
    if len(vu_scores) > 0 and len(gse_scores) > 0:
        vu_ranking = vu_scores.groupby("pathway")["total_score"].mean() \
            .sort_values(ascending=False)
        pw_col = "pathway" if "pathway" in gse_scores.columns else gse_scores.index.name
        score_cols = [c for c in gse_scores.columns if "total" in c.lower()]
        if pw_col and score_cols:
            if pw_col in gse_scores.columns:
                gse_ranking = gse_scores.set_index(pw_col)[score_cols[0]] \
                    .sort_values(ascending=False)
            else:
                gse_ranking = gse_scores[score_cols[0]].sort_values(ascending=False)

            shared = sorted(set(vu_ranking.index) & set(gse_ranking.index))
            if len(shared) >= 3:
                vu_ranks = [list(vu_ranking.index).index(p) for p in shared]
                gse_ranks = [list(gse_ranking.index).index(p) for p in shared]
                rho, pval = spearmanr(vu_ranks, gse_ranks)
                print(f"\n  Pathway ranking concordance: rho={rho:.3f}, p={pval:.3e}")

    # ── Summary ──
    print(f"\n  Summary:")
    print(f"    Vu arrays processed: {len(results_per_array)}")
    if len(vu_scores) > 0:
        print(f"    Vu pathways detected: {vu_scores['pathway'].nunique()}")
        print(f"    Vu pathway-array records: {len(vu_scores)}")

    print_header("16h: Complete")


if __name__ == "__main__":
    main()
