#!/usr/bin/env python3
"""
12_zonation_lr_validation.py — Zonation and LR interaction validation.

For each spatial dataset (Guilliams et al. + Vu et al.):
  1. Zonation: Score spots on periportal-pericentral axis using marker genes,
     classify DEGs by zonation, compare disease vs control zonation shift
  2. LR: Identify spatially co-localized ligand-receptor pairs involving DEGs,
     find disease-altered LR interactions
  3. Cross-dataset: Compare zonation and LR findings between datasets

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import spearmanr, mannwhitneyu, kruskal
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_dream_degs,
    load_conserved, build_ensembl_to_symbol_map,
    save_csv, save_checkpoint, print_header,
)

OUTPUT_DIR = RESULTS_DIR / "validation_bulk"


# ── Zonation ────────────────────────────────────────────────────────────────

def score_zonation(adata, pp_markers, pc_markers, ctrl_mult=5):
    """Compute continuous periportal–pericentral score per spot."""
    pp_present = [g for g in pp_markers if g in adata.var_names]
    pc_present = [g for g in pc_markers if g in adata.var_names]
    print(f"    PP markers: {len(pp_present)}/{len(pp_markers)} — {pp_present}")
    print(f"    PC markers: {len(pc_present)}/{len(pc_markers)} — {pc_present}")

    if len(pp_present) < 2 or len(pc_present) < 2:
        print("    WARNING: <2 markers in a category, zonation may be unreliable")
        return adata

    sc.tl.score_genes(adata, gene_list=pp_present, score_name="periportal_score",
                      ctrl_size=len(pp_present) * ctrl_mult)
    sc.tl.score_genes(adata, gene_list=pc_present, score_name="pericentral_score",
                      ctrl_size=len(pc_present) * ctrl_mult)

    adata.obs["zonation_score"] = (
        adata.obs["pericentral_score"] - adata.obs["periportal_score"]
    )

    # Quantile bins
    try:
        adata.obs["zonation_bin"] = pd.qcut(
            adata.obs["zonation_score"], q=5,
            labels=["PP1", "PP2", "Mid", "PC2", "PC1"],
            duplicates="drop",
        )
    except ValueError:
        adata.obs["zonation_bin"] = pd.qcut(
            adata.obs["zonation_score"], q=3,
            labels=["PP", "Mid", "PC"],
            duplicates="drop",
        )

    return adata


def classify_deg_zonation(adata, bulk_df, dataset_label):
    """For each DEG, compute expression correlation with zonation axis."""
    if "zonation_score" not in adata.obs.columns:
        return pd.DataFrame()

    zon = adata.obs["zonation_score"].values
    records = []

    # Get genes present in both spatial and bulk
    degs = bulk_df[bulk_df["dream_padj"] < 0.1]
    gene_col = "symbol"
    spatial_genes = set(adata.var_names)

    for _, row in degs.iterrows():
        gene = row[gene_col]
        if gene not in spatial_genes:
            continue

        # Get expression vector
        gene_idx = list(adata.var_names).index(gene)
        expr = np.array(adata.X[:, gene_idx].todense()).flatten() \
            if hasattr(adata.X, "todense") else adata.X[:, gene_idx].flatten()

        if expr.std() < 1e-10:
            continue

        rho, pval = spearmanr(zon, expr)

        # Classify
        if pval < 0.05:
            if rho > 0.1:
                zon_class = "Pericentral-enriched"
            elif rho < -0.1:
                zon_class = "Periportal-enriched"
            else:
                zon_class = "Pan-lobular"
        else:
            zon_class = "Non-zoned"

        # Mean expression per zonation bin
        bin_means = {}
        for b in adata.obs["zonation_bin"].cat.categories:
            mask = adata.obs["zonation_bin"] == b
            bin_means[f"mean_{b}"] = float(expr[mask].mean()) if mask.sum() > 0 else np.nan

        records.append({
            "gene": gene,
            "dataset": dataset_label,
            "dream_logFC": row["logFC"],
            "spearman_rho": rho,
            "spearman_pval": pval,
            "zonation_class": zon_class,
            **bin_means,
        })

    result = pd.DataFrame(records)
    if len(result) > 0:
        _, padj, _, _ = multipletests(result["spearman_pval"].values, method="fdr_bh")
        result["padj_bh"] = padj

        # Summary
        for cls in ["Pericentral-enriched", "Periportal-enriched", "Pan-lobular", "Non-zoned"]:
            n = (result["zonation_class"] == cls).sum()
            print(f"    {cls}: {n} DEGs")

    return result


def zonation_disruption(adata, dataset_label):
    """Compare zonation scores between conditions (if applicable)."""
    if "condition" not in adata.obs.columns or "zonation_score" not in adata.obs.columns:
        return None

    conditions = adata.obs["condition"].unique()
    if len(conditions) < 2:
        print(f"    Only one condition ({conditions[0]}), skipping disruption analysis")
        return None

    results = []
    for c in conditions:
        scores = adata.obs.loc[adata.obs["condition"] == c, "zonation_score"]
        results.append({
            "dataset": dataset_label,
            "condition": c,
            "n_spots": len(scores),
            "mean_zonation": scores.mean(),
            "median_zonation": scores.median(),
            "std_zonation": scores.std(),
        })

    return pd.DataFrame(results)


# ── Ligand-Receptor ─────────────────────────────────────────────────────────

def run_lr_analysis(adata, dataset_label, bulk_df, config):
    """Run squidpy LR analysis and cross-reference with DEGs."""
    print(f"\n  LR analysis for {dataset_label}...")

    # Ensure spatial neighbors
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

    # Run ligrec using permutation-based spatial co-expression
    try:
        sq.gr.ligrec(
            adata,
            n_perms=config.get("n_perms", 1000),
            cluster_key="zonation_bin" if "zonation_bin" in adata.obs.columns else None,
            interactions_params={"resources": "CellPhoneDB_v4.0"},
            use_raw=False,
        )
    except Exception as e:
        print(f"    WARNING: ligrec failed ({e}), trying spatial_autocorr on LR pairs")
        return _lr_coexpression_fallback(adata, dataset_label, bulk_df, config)

    # Extract results
    if "ligrec_test" in adata.uns:
        pvals = adata.uns["ligrec_test"]["pvalues"]
        means = adata.uns["ligrec_test"]["means"]

        # Flatten to per-pair results
        records = []
        for (ligand, receptor) in pvals.index:
            for cluster_pair in pvals.columns:
                p = pvals.loc[(ligand, receptor), cluster_pair]
                m = means.loc[(ligand, receptor), cluster_pair]
                if pd.notna(p) and p < 0.05:
                    records.append({
                        "ligand": ligand,
                        "receptor": receptor,
                        "cluster_pair": str(cluster_pair),
                        "pval": p,
                        "mean_expr": m,
                        "dataset": dataset_label,
                    })

        lr_df = pd.DataFrame(records)
        if len(lr_df) > 0:
            # Cross-reference with DEGs
            deg_symbols = set(bulk_df[bulk_df["dream_padj"] < 0.1]["symbol"])
            lr_df["ligand_is_deg"] = lr_df["ligand"].isin(deg_symbols)
            lr_df["receptor_is_deg"] = lr_df["receptor"].isin(deg_symbols)
            lr_df["either_deg"] = lr_df["ligand_is_deg"] | lr_df["receptor_is_deg"]
            lr_df["both_deg"] = lr_df["ligand_is_deg"] & lr_df["receptor_is_deg"]

            n_total = len(lr_df)
            n_deg = lr_df["either_deg"].sum()
            n_both = lr_df["both_deg"].sum()
            print(f"    {n_total} significant LR pairs, {n_deg} involving DEGs, "
                  f"{n_both} both DEGs")

        return lr_df

    return pd.DataFrame()


def _lr_coexpression_fallback(adata, dataset_label, bulk_df, config):
    """Fallback: compute spatial co-expression of known LR pairs."""
    # Use MASLD-relevant LR pairs from config
    fibrosis_lr = config.get("fibrosis_lr", [])
    inflammation_lr = config.get("inflammation_lr", [])
    sinusoidal_lr = config.get("sinusoidal_lr", [])

    all_pairs = []
    for pair_list, category in [(fibrosis_lr, "fibrosis"),
                                 (inflammation_lr, "inflammation"),
                                 (sinusoidal_lr, "sinusoidal")]:
        for pair_str in pair_list:
            if "_" in pair_str:
                parts = pair_str.split("_")
                if len(parts) == 2:
                    all_pairs.append((parts[0], parts[1], category))

    records = []
    deg_symbols = set(bulk_df[bulk_df["dream_padj"] < 0.1]["symbol"])
    spatial_genes = set(adata.var_names)

    for ligand, receptor, category in all_pairs:
        if ligand not in spatial_genes or receptor not in spatial_genes:
            continue

        # Get expression
        lig_idx = list(adata.var_names).index(ligand)
        rec_idx = list(adata.var_names).index(receptor)
        lig_expr = np.array(adata.X[:, lig_idx].todense()).flatten() \
            if hasattr(adata.X, "todense") else adata.X[:, lig_idx].flatten()
        rec_expr = np.array(adata.X[:, rec_idx].todense()).flatten() \
            if hasattr(adata.X, "todense") else adata.X[:, rec_idx].flatten()

        # Spatial co-expression (Spearman of expression vectors)
        if lig_expr.std() > 0 and rec_expr.std() > 0:
            rho, pval = spearmanr(lig_expr, rec_expr)
        else:
            rho, pval = 0, 1

        records.append({
            "ligand": ligand,
            "receptor": receptor,
            "category": category,
            "spatial_coexpr_rho": rho,
            "spatial_coexpr_pval": pval,
            "mean_ligand_expr": float(lig_expr.mean()),
            "mean_receptor_expr": float(rec_expr.mean()),
            "ligand_is_deg": ligand in deg_symbols,
            "receptor_is_deg": receptor in deg_symbols,
            "dataset": dataset_label,
        })

    lr_df = pd.DataFrame(records)
    if len(lr_df) > 0:
        n_sig = (lr_df["spatial_coexpr_pval"] < 0.05).sum()
        print(f"    {n_sig}/{len(lr_df)} LR pairs with significant spatial co-expression")
    return lr_df


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print_header("12: Zonation & LR Validation")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    config = load_config()
    zon_config = config.get("zonation", {})
    comm_config = config.get("communication", {})

    pp_markers = zon_config.get("periportal_markers",
                                ["HAL", "SDS", "ASS1", "CPS1", "ALB", "ASL", "HAMP", "HSD17B13"])
    pc_markers = zon_config.get("pericentral_markers",
                                ["GLS2", "CYP2E1", "CYP1A2", "GLUL", "CYP3A4", "CYP2A6"])

    # Load bulk
    print("  Loading bulk results...")
    dream = load_dream_degs(padj_thresh=1.0, lfc_thresh=0.0)
    gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
    padj_col = "padj" if "padj" in dream.columns else "adj.P.Val"
    bulk = dream.rename(columns={gene_col: "symbol", padj_col: "dream_padj"})
    bulk = bulk.dropna(subset=["symbol"]).drop_duplicates(subset=["symbol"], keep="first")

    cc_genes = set(load_conserved())

    # Load spatial datasets.
    # F008 fix (2026-06-01): the "Guilliams et al." arm must be GSE192741-ONLY.
    # merged_spatial.h5ad is the Vu-merged 24k-spot object (dataset ∈
    # {Vu_et_al_2025, GSE192741}); the old `species=='human'` filter was a no-op
    # (every spot is human), so the Guilliams arm was a Guilliams+Vu mixture and
    # the cross-dataset comparison below was largely a self-comparison. Prefer a
    # dedicated GSE192741 object if it exists, else subset on dataset=='GSE192741'.
    datasets = {}
    gse_path = RESULTS_DIR / "preprocessed" / "merged_spatial_gse192741.h5ad"
    merged_path = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"
    if gse_path.exists():
        datasets["Guilliams et al."] = sc.read_h5ad(gse_path)
    elif merged_path.exists():
        adata_g = sc.read_h5ad(merged_path)
        if "dataset" in adata_g.obs.columns:
            adata_g = adata_g[adata_g.obs["dataset"] == "GSE192741"].copy()
        elif "sample_id" in adata_g.obs.columns:
            # Fallback: GSE192741 sample_ids all start with 'JBO'.
            adata_g = adata_g[adata_g.obs["sample_id"].astype(str)
                              .str.startswith("JBO")].copy()
        print(f"  Guilliams (GSE192741-only): {adata_g.n_obs} spots")
        datasets["Guilliams et al."] = adata_g

    vu_path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"
    if vu_path.exists():
        datasets["Vu et al."] = sc.read_h5ad(vu_path)

    all_zonation = []
    all_disruption = []
    all_lr = []

    for ds_label, adata in datasets.items():
        print_header(f"Dataset: {ds_label}")

        # ── Zonation ──
        print("  Scoring zonation...")
        adata = score_zonation(adata, pp_markers, pc_markers)

        if "zonation_score" in adata.obs.columns:
            # Validate: Moran's I on zonation score
            sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
            # Store zonation score in a gene-like slot for Moran's I
            adata.obs["_zon_score"] = adata.obs["zonation_score"]
            print(f"    Zonation score range: {adata.obs['zonation_score'].min():.3f} "
                  f"to {adata.obs['zonation_score'].max():.3f}")

            bins = adata.obs["zonation_bin"].value_counts()
            for b, n in bins.items():
                print(f"      {b}: {n} spots")

            # Classify DEGs by zonation
            print("\n  Classifying DEGs by zonation...")
            zon_df = classify_deg_zonation(adata, bulk, ds_label)
            if len(zon_df) > 0:
                all_zonation.append(zon_df)
                save_csv(zon_df,
                         f"deg_zonation_{ds_label.replace(' ', '_').replace('.', '')}.csv",
                         subdir="validation_bulk")

            # Zonation disruption
            print("\n  Checking zonation disruption...")
            disrupt = zonation_disruption(adata, ds_label)
            if disrupt is not None:
                all_disruption.append(disrupt)

        # ── LR ──
        print("\n  Running LR analysis...")
        lr_df = _lr_coexpression_fallback(adata, ds_label, bulk, comm_config)
        if len(lr_df) > 0:
            all_lr.append(lr_df)
            save_csv(lr_df,
                     f"lr_pairs_{ds_label.replace(' ', '_').replace('.', '')}.csv",
                     subdir="validation_bulk")

    # ── Cross-dataset comparison ──
    if len(all_zonation) == 2:
        print_header("Cross-Dataset Zonation Comparison")
        z1 = all_zonation[0].set_index("gene")
        z2 = all_zonation[1].set_index("gene")
        shared = z1.index.intersection(z2.index)
        print(f"  Shared DEGs: {len(shared)}")

        rho1 = z1.loc[shared, "spearman_rho"]
        rho2 = z2.loc[shared, "spearman_rho"]
        mask = rho1.notna() & rho2.notna()
        if mask.sum() > 10:
            r, p = spearmanr(rho1[mask], rho2[mask])
            print(f"  Zonation rho correlation: r={r:.3f}, p={p:.2e}")

            # Class concordance
            c1 = z1.loc[shared[mask], "zonation_class"]
            c2 = z2.loc[shared[mask], "zonation_class"]
            concordant = (c1 == c2).sum()
            pct = concordant / len(c1) * 100
            print(f"  Zonation class concordance: {concordant}/{len(c1)} ({pct:.1f}%)")

    if len(all_lr) == 2:
        print_header("Cross-Dataset LR Comparison")
        lr1 = all_lr[0]
        lr2 = all_lr[1]
        # Match on ligand+receptor
        lr1["pair"] = lr1["ligand"] + "_" + lr1["receptor"]
        lr2["pair"] = lr2["ligand"] + "_" + lr2["receptor"]
        shared_pairs = set(lr1["pair"]) & set(lr2["pair"])
        print(f"  Shared LR pairs: {len(shared_pairs)}")

        if len(shared_pairs) > 3:
            m1 = lr1.set_index("pair").loc[list(shared_pairs)]
            m2 = lr2.set_index("pair").loc[list(shared_pairs)]
            if "spatial_coexpr_rho" in m1.columns:
                r, p = spearmanr(m1["spatial_coexpr_rho"], m2["spatial_coexpr_rho"])
                print(f"  LR co-expression rho correlation: r={r:.3f}, p={p:.2e}")

    # ── Save combined ──
    if all_zonation:
        combined_zon = pd.concat(all_zonation, ignore_index=True)
        save_csv(combined_zon, "deg_zonation_combined.csv", subdir="validation_bulk")
    if all_disruption:
        combined_dis = pd.concat(all_disruption, ignore_index=True)
        save_csv(combined_dis, "zonation_disruption.csv", subdir="validation_bulk")
    if all_lr:
        combined_lr = pd.concat(all_lr, ignore_index=True)
        save_csv(combined_lr, "lr_pairs_combined.csv", subdir="validation_bulk")

    print_header("12: Complete")


if __name__ == "__main__":
    main()
