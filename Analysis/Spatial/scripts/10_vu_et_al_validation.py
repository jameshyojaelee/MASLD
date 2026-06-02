#!/usr/bin/env python3
"""
10_vu_et_al_validation.py — Spatial validation using Vu et al. JHEP Reports 2025.

Detects SVGs per-array via Moran's I, computes consensus across arrays,
and cross-references with GSE192741 SVGs, dream DEGs, and Conserved
for independent spatial replication.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_dream_degs,
    load_conserved, build_ensembl_to_symbol_map,
    save_csv, print_header, print_step,
)

OUTPUT_DIR = RESULTS_DIR / "validation_vu"


def load_vu_data():
    """Load Vu et al. preprocessed checkpoint."""
    path = RESULTS_DIR / "preprocessed" / "merged_spatial_vu.h5ad"
    print(f"  Loading: {path}")
    adata = sc.read_h5ad(path)
    print(f"  {adata.n_obs} spots × {adata.n_vars} genes, "
          f"{adata.obs['sample_id'].nunique()} arrays")
    return adata


def load_gse192741_svgs():
    """Load GSE192741 SVG results for comparison."""
    results = {}
    svg_dir = RESULTS_DIR / "svg"
    for csv_path in svg_dir.glob("svgs_*.csv"):
        condition = csv_path.stem.replace("svgs_", "")
        df = pd.read_csv(csv_path, index_col=0)
        results[condition] = df
        n_svg = df["svg"].sum() if "svg" in df.columns else 0
        print(f"  GSE192741 {condition}: {n_svg} SVGs / {len(df)} tested")
    return results


def detect_svgs_per_array(adata, config):
    """Run Moran's I per array, return per-array SVG DataFrames."""
    arrays = sorted(adata.obs["sample_id"].unique())
    per_array = {}

    for i, array_id in enumerate(arrays):
        print_step(f"SVGs for {array_id}", i + 1, len(arrays))
        sub = adata[adata.obs["sample_id"] == array_id].copy()

        if sub.n_obs < 100:
            print(f"    WARNING: {sub.n_obs} spots, skipping")
            continue

        # Build spatial graph
        sq.gr.spatial_neighbors(sub, coord_type="generic", n_neighs=6)

        # Test HVGs
        hvgs = sub.var_names[sub.var["highly_variable"]].tolist()
        if len(hvgs) == 0:
            hvgs = sub.var_names.tolist()[:2000]  # fallback

        sq.gr.spatial_autocorr(
            sub, mode="moran", genes=hvgs,
            n_perms=config.get("n_perms", 100),
            n_jobs=config.get("n_jobs", 8),
        )

        svg_df = sub.uns["moranI"].copy()
        # BH FDR correction
        pvals = svg_df["pval_norm"].fillna(1.0).values
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        svg_df["padj_bh"] = padj
        svg_df["svg"] = (
            (svg_df["padj_bh"] < config.get("fdr_threshold", 0.05))
            & (svg_df["I"] > config.get("min_morans_i", 0.1))
        )
        svg_df = svg_df.sort_values("I", ascending=False)

        n_svg = svg_df["svg"].sum()
        top3 = ", ".join(f"{g}(I={svg_df.loc[g, 'I']:.2f})"
                         for g in svg_df.head(3).index)
        print(f"    {n_svg} SVGs / {len(hvgs)} tested. Top: {top3}")

        per_array[array_id] = svg_df

    return per_array


def build_consensus_svgs(per_array, min_arrays=3):
    """Build consensus SVG table across arrays."""
    all_genes = set()
    for df in per_array.values():
        all_genes.update(df.index)

    records = []
    for gene in sorted(all_genes):
        morans = []
        n_svg = 0
        n_tested = 0
        for array_id, df in per_array.items():
            if gene in df.index:
                n_tested += 1
                morans.append(df.loc[gene, "I"])
                if df.loc[gene, "svg"]:
                    n_svg += 1

        records.append({
            "gene": gene,
            "n_arrays_tested": n_tested,
            "n_arrays_svg": n_svg,
            "mean_morans_i": np.mean(morans) if morans else np.nan,
            "median_morans_i": np.median(morans) if morans else np.nan,
            "max_morans_i": np.max(morans) if morans else np.nan,
            "consensus_svg": n_svg >= min_arrays,
        })

    consensus = pd.DataFrame(records).set_index("gene")
    consensus = consensus.sort_values("mean_morans_i", ascending=False)

    n_consensus = consensus["consensus_svg"].sum()
    print(f"\n  Consensus SVGs (≥{min_arrays} arrays): {n_consensus} / {len(consensus)}")
    return consensus


def cross_reference_gse192741(consensus, gse_svgs):
    """Compare Vu et al. consensus SVGs with GSE192741 SVGs."""
    results = []

    vu_svg_set = set(consensus[consensus["consensus_svg"]].index)

    for condition, gse_df in gse_svgs.items():
        gse_svg_set = set(gse_df[gse_df["svg"]].index)
        shared_genes = set(consensus.index) & set(gse_df.index)

        overlap = vu_svg_set & gse_svg_set
        vu_only = vu_svg_set - gse_svg_set
        gse_only = gse_svg_set - vu_svg_set

        # Fisher's exact on shared gene universe
        a = len(overlap)
        b = len(vu_svg_set & shared_genes) - a
        c = len(gse_svg_set & shared_genes) - a
        d = len(shared_genes) - a - b - c
        if min(a + b, c + d, a + c, b + d) > 0:
            odds, pval = fisher_exact([[a, b], [c, d]])
        else:
            odds, pval = np.nan, np.nan

        # Moran's I correlation for shared genes
        shared_list = sorted(shared_genes)
        vu_i = consensus.loc[shared_list, "mean_morans_i"]
        gse_i = gse_df.loc[shared_list, "I"]
        mask = vu_i.notna() & gse_i.notna()
        if mask.sum() > 10:
            rho, rho_p = spearmanr(vu_i[mask], gse_i[mask])
        else:
            rho, rho_p = np.nan, np.nan

        results.append({
            "comparison": f"Vu_consensus_vs_GSE192741_{condition}",
            "n_vu_svgs": len(vu_svg_set),
            "n_gse_svgs": len(gse_svg_set),
            "n_shared_genes": len(shared_genes),
            "n_overlap": a,
            "n_vu_only": len(vu_only),
            "n_gse_only": len(gse_only),
            "jaccard": a / max(len(vu_svg_set | gse_svg_set), 1),
            "fisher_OR": odds,
            "fisher_pval": pval,
            "morans_i_spearman_rho": rho,
            "morans_i_spearman_pval": rho_p,
        })

        print(f"\n  {condition}:")
        print(f"    Overlap: {a} genes (Jaccard={a / max(len(vu_svg_set | gse_svg_set), 1):.3f})")
        print(f"    Fisher OR={odds:.2f}, p={pval:.2e}")
        print(f"    Moran's I correlation: rho={rho:.3f}, p={rho_p:.2e}")
        if a > 0:
            top_overlap = sorted(overlap, key=lambda g: consensus.loc[g, "mean_morans_i"]
                                 if g in consensus.index else 0, reverse=True)[:10]
            print(f"    Top replicated: {', '.join(top_overlap)}")

    return pd.DataFrame(results)


def cross_reference_bulk(consensus, min_arrays=3):
    """Cross-reference with dream DEGs and Conserved."""
    results = []
    vu_svg_set = set(consensus[consensus["consensus_svg"]].index)
    all_genes = set(consensus.index)

    # 1. Dream DEGs
    try:
        dream = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
        gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
        deg_set = set(dream[gene_col].dropna()) & all_genes

        a = len(vu_svg_set & deg_set)
        b = len(vu_svg_set - deg_set)
        c = len(deg_set - vu_svg_set)
        d = len(all_genes - vu_svg_set - deg_set)
        odds, pval = fisher_exact([[a, b], [c, d]])
        results.append({
            "gene_set": "Dream_DEGs",
            "n_overlap": a, "n_svg": len(vu_svg_set), "n_set": len(deg_set),
            "odds_ratio": odds, "fisher_pval": pval,
        })
        print(f"\n  Dream DEGs: {a}/{len(vu_svg_set)} SVGs are DEGs "
              f"(OR={odds:.2f}, p={pval:.2e})")
    except Exception as e:
        print(f"  WARNING: Could not load dream DEGs: {e}")

    # 2. Conserved
    cc = load_conserved()
    if cc:
        cc_set = set(cc) & all_genes
        a = len(vu_svg_set & cc_set)
        b = len(vu_svg_set - cc_set)
        c = len(cc_set - vu_svg_set)
        d = len(all_genes - vu_svg_set - cc_set)
        odds, pval = fisher_exact([[a, b], [c, d]])
        results.append({
            "gene_set": "Conserved",
            "n_overlap": a, "n_svg": len(vu_svg_set), "n_set": len(cc_set),
            "odds_ratio": odds, "fisher_pval": pval,
        })
        print(f"  Conserved: {a}/{len(vu_svg_set)} SVGs in CC "
              f"(OR={odds:.2f}, p={pval:.2e})")

    # 3. Up-regulated DEGs specifically
    try:
        dream_up = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.5)
        gene_col = "symbol" if "symbol" in dream_up.columns else dream_up.columns[0]
        up_set = set(dream_up[dream_up["logFC"] > 0][gene_col].dropna()) & all_genes

        a = len(vu_svg_set & up_set)
        b = len(vu_svg_set - up_set)
        c = len(up_set - vu_svg_set)
        d = len(all_genes - vu_svg_set - up_set)
        odds, pval = fisher_exact([[a, b], [c, d]])
        results.append({
            "gene_set": "Dream_Up_DEGs",
            "n_overlap": a, "n_svg": len(vu_svg_set), "n_set": len(up_set),
            "odds_ratio": odds, "fisher_pval": pval,
        })
        print(f"  Up-regulated DEGs: {a}/{len(vu_svg_set)} "
              f"(OR={odds:.2f}, p={pval:.2e})")
    except Exception as e:
        print(f"  WARNING: Could not load up DEGs: {e}")

    tests_df = pd.DataFrame(results)
    if len(tests_df) > 0:
        _, padj, _, _ = multipletests(tests_df["fisher_pval"].values, method="fdr_bh")
        tests_df["fisher_padj_bh"] = padj

    return tests_df


def main():
    print_header("10: Vu et al. Spatial Validation")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    config = load_config()
    svg_config = config.get("svg", {
        "n_perms": 100, "n_jobs": 8,
        "fdr_threshold": 0.05, "min_morans_i": 0.1,
    })

    # ── Load data ──
    print("  Loading Vu et al. data...")
    adata = load_vu_data()

    print("\n  Loading GSE192741 SVGs...")
    gse_svgs = load_gse192741_svgs()

    # ── Per-array SVG detection ──
    print_header("Per-Array SVG Detection")
    per_array = detect_svgs_per_array(adata, svg_config)

    # Save per-array results
    for array_id, df in per_array.items():
        save_csv(df, f"svgs_vu_{array_id}.csv", subdir="validation_vu")

    # ── Consensus SVGs ──
    print_header("Consensus SVGs")
    consensus = build_consensus_svgs(per_array, min_arrays=3)
    save_csv(consensus, "consensus_svgs_vu.csv", subdir="validation_vu")

    # Top consensus SVGs
    top = consensus[consensus["consensus_svg"]].head(20)
    print("\n  Top 20 consensus SVGs:")
    for gene, row in top.iterrows():
        print(f"    {gene}: mean_I={row['mean_morans_i']:.3f}, "
              f"{int(row['n_arrays_svg'])}/{int(row['n_arrays_tested'])} arrays")

    # ── Cross-dataset concordance ──
    print_header("Cross-Dataset Concordance (vs GSE192741)")
    if gse_svgs:
        concordance = cross_reference_gse192741(consensus, gse_svgs)
        save_csv(concordance, "concordance_gse192741.csv", subdir="validation_vu")

    # ── Bulk RNA-seq cross-reference ──
    print_header("Bulk RNA-seq Cross-Reference")
    bulk_tests = cross_reference_bulk(consensus)
    if len(bulk_tests) > 0:
        save_csv(bulk_tests, "bulk_enrichment_tests.csv", subdir="validation_vu")

    # ── Summary ──
    print_header("Validation Summary")
    n_consensus = consensus["consensus_svg"].sum()
    print(f"  Vu et al. arrays: {len(per_array)}")
    print(f"  Total genes tested: {len(consensus)}")
    print(f"  Consensus SVGs (≥3 arrays): {n_consensus}")
    if gse_svgs:
        for cond, gse_df in gse_svgs.items():
            gse_svg_set = set(gse_df[gse_df["svg"]].index)
            vu_svg_set = set(consensus[consensus["consensus_svg"]].index)
            overlap = vu_svg_set & gse_svg_set
            print(f"  Replicated in GSE192741 {cond}: {len(overlap)}")

    print_header("10: Complete")


if __name__ == "__main__":
    main()
