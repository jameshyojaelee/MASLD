#!/usr/bin/env python3
"""
09c_niche_drug_targets.py — Niche-specific drug target prioritization.

Identifies disease niches (spatial domains enriched in MASLD) and
cross-references niche markers with DGIdb druggable genes, LINCS
reversal compounds, and network proximity scores.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_multi_evidence_atlas,
    save_csv, print_header, print_step,
)


def identify_disease_niches(adata, condition_healthy, condition_masld):
    """Identify spatial domains whose proportion increases in MASLD.

    Returns dict of {domain_id: fold_change_in_masld}.
    """
    if "spatial_domain" not in adata.obs.columns:
        print("  WARNING: No spatial_domain column")
        return {}

    domains = adata.obs["spatial_domain"].unique()
    niche_fc = {}

    n_healthy = (adata.obs["condition"] == condition_healthy).sum()
    n_masld = (adata.obs["condition"] == condition_masld).sum()

    for d in domains:
        d_healthy = ((adata.obs["condition"] == condition_healthy) &
                     (adata.obs["spatial_domain"] == d)).sum()
        d_masld = ((adata.obs["condition"] == condition_masld) &
                   (adata.obs["spatial_domain"] == d)).sum()

        prop_h = d_healthy / max(n_healthy, 1)
        prop_m = d_masld / max(n_masld, 1)
        fc = prop_m / max(prop_h, 1e-6)
        niche_fc[d] = {
            "fold_change": fc,
            "prop_healthy": prop_h,
            "prop_masld": prop_m,
            "n_spots_healthy": d_healthy,
            "n_spots_masld": d_masld,
        }

    return niche_fc


def get_domain_markers(adata, domain_id, n_top=50):
    """Get top marker genes for a spatial domain using rank_genes_groups."""
    if "spatial_domain" not in adata.obs.columns:
        return []

    # Use Wilcoxon rank-sum for domain vs rest
    adata_copy = adata.copy()
    adata_copy.obs["spatial_domain"] = adata_copy.obs["spatial_domain"].astype(str).astype("category")

    try:
        sc.tl.rank_genes_groups(
            adata_copy, groupby="spatial_domain", method="wilcoxon",
            use_raw=False, n_genes=n_top,
        )
        markers = sc.get.rank_genes_groups_df(adata_copy, group=str(domain_id))
        return markers.head(n_top)["names"].tolist()
    except Exception as e:
        print(f"  WARNING: rank_genes_groups failed for domain {domain_id}: {e}")
        return []


def load_lincs_compounds():
    """Load LINCS annotated compound data."""
    path = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/lincs_annotated_compounds.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def load_network_proximity():
    """Load network proximity results."""
    path = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/network_proximity/network_proximity_scores.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def load_drug_targets():
    """Load multi-layer drug targets."""
    path = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/multi_layer_drug_targets.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main():
    print_header("09c: Niche-Specific Drug Target Prioritization")

    config = load_config()
    output_dir = RESULTS_DIR / "drug_repurposing"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load spatial data with domains
    dom_h5ad = RESULTS_DIR / "domains" / "spatial_with_domains.h5ad"
    if not dom_h5ad.exists():
        # Try alternative path
        dom_h5ad = RESULTS_DIR / "preprocessed" / "merged_spatial.h5ad"
    if not dom_h5ad.exists():
        print("  ERROR: No spatial domain data found")
        sys.exit(1)

    adata = sc.read_h5ad(dom_h5ad)
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    # Load domain assignments if separate
    dom_csv = RESULTS_DIR / "domains" / "spatial_domains.csv"
    if dom_csv.exists() and "spatial_domain" not in adata.obs.columns:
        domains = pd.read_csv(dom_csv, index_col=0)
        shared = adata.obs.index.intersection(domains.index)
        adata.obs.loc[shared, "spatial_domain"] = domains.loc[shared, "spatial_domain"]

    # Load zonation if available
    zon_csv = RESULTS_DIR / "zonation" / "zonation_scores.csv"
    if zon_csv.exists() and "zonation_score" not in adata.obs.columns:
        zon = pd.read_csv(zon_csv, index_col=0)
        shared = adata.obs.index.intersection(zon.index)
        for col in ["zonation_score", "zonation_bin"]:
            if col in zon.columns:
                adata.obs.loc[shared, col] = zon.loc[shared, col]

    # Identify conditions
    conditions = adata.obs["condition"].unique().tolist()
    cond_healthy = [c for c in conditions if "healthy" in c.lower() or c == "Healthy"]
    cond_masld = [c for c in conditions if c not in cond_healthy]
    if not cond_healthy or not cond_masld:
        print("  ERROR: Need both healthy and MASLD conditions")
        sys.exit(1)

    # Identify disease niches
    print("\n  Identifying disease niches...")
    niche_fc = identify_disease_niches(adata, cond_healthy[0], cond_masld[0])

    niche_summary = pd.DataFrame([
        {"domain": d, **v} for d, v in niche_fc.items()
    ])
    if len(niche_summary) > 0:
        niche_summary = niche_summary.sort_values("fold_change", ascending=False)
        save_csv(niche_summary, "disease_niche_proportions.csv", subdir="drug_repurposing")
        print(f"  Domain proportions:")
        for _, row in niche_summary.iterrows():
            direction = "↑" if row["fold_change"] > 1.2 else "↓" if row["fold_change"] < 0.8 else "→"
            print(f"    D{row['domain']}: FC={row['fold_change']:.2f} {direction} "
                  f"(H={row['prop_healthy']:.3f}, M={row['prop_masld']:.3f})")

    # Get marker genes for disease-enriched domains
    disease_domains = [d for d, v in niche_fc.items() if v["fold_change"] > 1.2]
    print(f"\n  Disease-enriched domains (FC>1.2): {disease_domains}")

    all_niche_markers = []
    for domain in disease_domains:
        markers = get_domain_markers(adata, domain, n_top=50)
        print(f"    Domain {domain}: {len(markers)} marker genes")
        for gene in markers:
            all_niche_markers.append({"gene": gene, "domain": domain})

    niche_markers_df = pd.DataFrame(all_niche_markers)
    if len(niche_markers_df) > 0:
        niche_markers_df = niche_markers_df.drop_duplicates(subset=["gene"])
        niche_marker_genes = set(niche_markers_df["gene"])
    else:
        niche_marker_genes = set()

    # Load multi-evidence atlas for druggability
    print("\n  Loading multi-evidence atlas...")
    atlas = load_multi_evidence_atlas()
    sym_col = "symbol" if "symbol" in atlas.columns else "human_symbol"

    # Cross-reference with druggable genes
    druggable_col = [c for c in atlas.columns if "dgidb" in c.lower() or "druggable" in c.lower()]
    if druggable_col:
        druggable = set(atlas[atlas[druggable_col[0]] == True][sym_col])
        niche_druggable = niche_marker_genes & druggable
        print(f"  Niche markers that are druggable: {len(niche_druggable)}/{len(niche_marker_genes)}")
    else:
        niche_druggable = set()

    # Load LINCS compounds
    lincs = load_lincs_compounds()
    if len(lincs) > 0:
        lincs_target_col = [c for c in lincs.columns if "target" in c.lower() or "gene" in c.lower()]
        if lincs_target_col:
            lincs_targets = set(lincs[lincs_target_col[0]].dropna())
            niche_lincs = niche_marker_genes & lincs_targets
            print(f"  Niche markers with LINCS compounds: {len(niche_lincs)}")

    # Load network proximity
    net_prox = load_network_proximity()
    drug_targets = load_drug_targets()

    # Build spatial drug target cards
    print("\n  Building spatial drug target cards...")
    target_cards = []
    for _, row in niche_markers_df.iterrows():
        gene = row["gene"]
        domain = row["domain"]

        card = {
            "gene": gene,
            "disease_niche_domain": domain,
            "niche_fc": niche_fc.get(domain, {}).get("fold_change", np.nan),
        }

        # Add zonation info
        if "zonation_bin" in adata.obs.columns:
            zon_class_path = RESULTS_DIR / "zonation" / "deg_zonation_classification.csv"
            if zon_class_path.exists():
                zon_class = pd.read_csv(zon_class_path)
                gene_zon = zon_class[zon_class["gene"] == gene]
                if len(gene_zon) > 0:
                    card["zonation_class"] = gene_zon.iloc[0].get("zonation_class", "Unknown")

        # Add atlas evidence
        gene_atlas = atlas[atlas[sym_col] == gene]
        if len(gene_atlas) > 0:
            row_atlas = gene_atlas.iloc[0]
            card["dream_logFC"] = row_atlas.get("dream_logFC", np.nan)
            card["dream_padj"] = row_atlas.get("dream_padj", np.nan)
            card["is_conserved"] = row_atlas.get("is_conserved", False)
            card["is_druggable"] = gene in niche_druggable
            # mr_pval removed 2026-04-22 — MR ditched from paper.
            for col in ["twas_pval", "coloc_pp4", "sex_class"]:
                if col in row_atlas.index:
                    card[col] = row_atlas[col]

        # Add drug info
        if len(drug_targets) > 0:
            dt_col = [c for c in drug_targets.columns if "gene" in c.lower() or "target" in c.lower()]
            if dt_col:
                gene_dt = drug_targets[drug_targets[dt_col[0]] == gene]
                if len(gene_dt) > 0:
                    card["known_drugs"] = True
                    drug_col = [c for c in drug_targets.columns if "drug" in c.lower()]
                    if drug_col:
                        card["drug_names"] = "; ".join(gene_dt[drug_col[0]].dropna().tolist()[:3])

        target_cards.append(card)

    if target_cards:
        cards_df = pd.DataFrame(target_cards)
        # Sort by druggability + niche FC
        cards_df["priority_score"] = (
            cards_df["is_druggable"].fillna(False).astype(int) * 2 +
            cards_df["niche_fc"].fillna(1.0)
        )
        cards_df = cards_df.sort_values("priority_score", ascending=False)
        save_csv(cards_df, "spatial_niche_targets.csv", subdir="drug_repurposing")

        print(f"\n  Total niche drug target cards: {len(cards_df)}")
        print(f"  Druggable: {cards_df['is_druggable'].sum() if 'is_druggable' in cards_df.columns else 0}")
        print(f"\n  Top 10 spatial niche targets:")
        for _, row in cards_df.head(10).iterrows():
            drug_info = row.get("drug_names", "")
            print(f"    {row['gene']} (D{row['disease_niche_domain']}, "
                  f"FC={row['niche_fc']:.2f}"
                  f"{', drugs: ' + drug_info if drug_info else ''})")

    print_header("09c: Complete")


if __name__ == "__main__":
    main()
