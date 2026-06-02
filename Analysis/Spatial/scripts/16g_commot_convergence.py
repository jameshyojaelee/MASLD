#!/usr/bin/env python3
"""
16g_commot_convergence.py — COMMOT convergence with multi-evidence atlas.

Tests enrichment of signaling-regulated genes in COLOC hits, drug targets,
Conserved, and dream DEGs. Builds drug-target signaling evidence cards.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config,
    load_coloc_results, load_drug_targets, load_conserved,
    load_dream_degs, load_multi_evidence_atlas,
    fisher_test, save_csv, print_header, print_step,
)


def load_signaling_regulated_genes():
    """Load signaling-regulated gene summary from 16e."""
    path = RESULTS_DIR / "commot" / "signaling_regulated_genes_summary.csv"
    if not path.exists():
        # Try per-condition files
        all_dfs = []
        for csv_path in (RESULTS_DIR / "commot").glob("signaling_regulated_genes_*.csv"):
            if "summary" in csv_path.name:
                continue
            all_dfs.append(pd.read_csv(csv_path, index_col=0))
        if all_dfs:
            combined = pd.concat(all_dfs, ignore_index=True)
            return combined
        print("  WARNING: No signaling-regulated gene results found. "
              "Producing empty convergence output.")
        return pd.DataFrame(columns=["gene", "pathway", "condition", "impact_score"])
    return pd.read_csv(path, index_col=0)


def get_signaling_gene_set(sig_df):
    """Extract unique signaling-regulated gene names."""
    gene_col = "gene" if "gene" in sig_df.columns else sig_df.index.name
    if gene_col and gene_col in sig_df.columns:
        return set(sig_df[gene_col].dropna())
    return set(sig_df.index)


def run_enrichment_tests(sig_genes, universe):
    """Run Fisher's exact tests against key gene sets."""
    results = []

    # 1. COLOC overlap
    print_step("COLOC enrichment", 1, 4)
    coloc_df = load_coloc_results(pp4_threshold=0.5)
    if len(coloc_df) > 0:
        gene_cols = ["gene_symbol", "symbol", "gene_name", "gene"]
        coloc_gene_col = next((c for c in gene_cols if c in coloc_df.columns), None)
        if coloc_gene_col:
            coloc_genes = set(coloc_df[coloc_gene_col].dropna().unique())
            odds, pval, n_overlap = fisher_test(sig_genes, coloc_genes, universe)
            results.append({
                "gene_set": "COLOC_hits",
                "n_sig_genes": len(sig_genes),
                "n_target_set": len(coloc_genes),
                "n_overlap": n_overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })
            print(f"    COLOC: {n_overlap} overlap, OR={odds:.2f}, p={pval:.2e}")

    # 2. Drug targets
    print_step("Drug target enrichment", 2, 4)
    drug_df = load_drug_targets()
    if len(drug_df) > 0:
        dt_col = [c for c in drug_df.columns if "gene" in c.lower() or "target" in c.lower()]
        if dt_col:
            drug_genes = set(drug_df[dt_col[0]].dropna().unique())
            odds, pval, n_overlap = fisher_test(sig_genes, drug_genes, universe)
            results.append({
                "gene_set": "Drug_targets",
                "n_sig_genes": len(sig_genes),
                "n_target_set": len(drug_genes),
                "n_overlap": n_overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })
            print(f"    Drug targets: {n_overlap} overlap, OR={odds:.2f}, p={pval:.2e}")

    # 3. Conserved
    print_step("Conserved enrichment", 3, 4)
    cc_genes = set(load_conserved()) & universe
    if cc_genes:
        odds, pval, n_overlap = fisher_test(sig_genes, cc_genes, universe)
        results.append({
            "gene_set": "Conserved",
            "n_sig_genes": len(sig_genes),
            "n_target_set": len(cc_genes),
            "n_overlap": n_overlap,
            "odds_ratio": odds,
            "fisher_pval": pval,
        })
        print(f"    Conserved: {n_overlap} overlap, OR={odds:.2f}, p={pval:.2e}")

    # 4. Dream DEGs
    print_step("Dream DEG enrichment", 4, 4)
    dream = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
    gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
    deg_genes = set(dream[gene_col].dropna()) & universe
    if deg_genes:
        odds, pval, n_overlap = fisher_test(sig_genes, deg_genes, universe)
        results.append({
            "gene_set": "Dream_DEGs",
            "n_sig_genes": len(sig_genes),
            "n_target_set": len(deg_genes),
            "n_overlap": n_overlap,
            "odds_ratio": odds,
            "fisher_pval": pval,
        })
        print(f"    Dream DEGs: {n_overlap} overlap, OR={odds:.2f}, p={pval:.2e}")

    # DGIdb druggable from atlas
    atlas = load_multi_evidence_atlas()
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
    if "dgidb_druggable" in atlas.columns:
        druggable = set(atlas[atlas["dgidb_druggable"] == True][sym_col]) & universe
        if druggable:
            odds, pval, n_overlap = fisher_test(sig_genes, druggable, universe)
            results.append({
                "gene_set": "DGIdb_druggable",
                "n_sig_genes": len(sig_genes),
                "n_target_set": len(druggable),
                "n_overlap": n_overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })
            print(f"    DGIdb druggable: {n_overlap} overlap, OR={odds:.2f}, p={pval:.2e}")

    return pd.DataFrame(results)


def build_drug_target_cards(sig_df, sig_genes):
    """Build evidence cards for drug targets that are signaling-regulated."""
    drug_df = load_drug_targets()
    if len(drug_df) == 0:
        return pd.DataFrame()

    dt_col = [c for c in drug_df.columns if "gene" in c.lower() or "target" in c.lower()]
    if not dt_col:
        return pd.DataFrame()

    drug_col = [c for c in drug_df.columns if "drug" in c.lower()]
    drug_genes = set(drug_df[dt_col[0]].dropna())
    overlap = sig_genes & drug_genes

    if not overlap:
        print("  No drug targets are signaling-regulated")
        return pd.DataFrame()

    atlas = load_multi_evidence_atlas()
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"

    cards = []
    gene_col = "gene" if "gene" in sig_df.columns else sig_df.index.name
    for gene in sorted(overlap):
        card = {"gene": gene, "is_signaling_regulated": True}

        # Add signaling info
        if gene_col and gene_col in sig_df.columns:
            gene_rows = sig_df[sig_df[gene_col] == gene]
        else:
            gene_rows = sig_df.loc[[gene]] if gene in sig_df.index else pd.DataFrame()

        if len(gene_rows) > 0:
            row = gene_rows.iloc[0]
            for col in ["mean_impact", "max_impact", "n_pathways", "top_pathway",
                         "impact_score", "pathway"]:
                if col in row.index:
                    card[f"commot_{col}"] = row[col]

        # Add drug info
        gene_drugs = drug_df[drug_df[dt_col[0]] == gene]
        if len(gene_drugs) > 0 and drug_col:
            card["drugs"] = "; ".join(gene_drugs[drug_col[0]].dropna().tolist()[:5])

        # Add atlas evidence
        gene_atlas = atlas[atlas[sym_col] == gene]
        if len(gene_atlas) > 0:
            r = gene_atlas.iloc[0]
            # mr_pval removed 2026-04-22 — MR ditched from paper.
            for col in ["dream_logFC", "dream_padj", "is_conserved",
                         "twas_pval"]:
                if col in r.index:
                    card[col] = r[col]

        cards.append(card)

    return pd.DataFrame(cards)


def main():
    print_header("16g: COMMOT Convergence Analysis")

    config = load_config()
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load signaling-regulated genes ──
    print("  Loading signaling-regulated genes...")
    sig_df = load_signaling_regulated_genes()
    sig_genes = get_signaling_gene_set(sig_df)
    print(f"  Signaling-regulated genes: {len(sig_genes)}")

    # Build universe from atlas
    atlas = load_multi_evidence_atlas()
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
    universe = set(atlas[sym_col].dropna())
    sig_genes_in_universe = sig_genes & universe
    print(f"  In atlas universe: {len(sig_genes_in_universe)}/{len(sig_genes)}")

    # ── Enrichment tests ──
    print_header("Enrichment Tests")
    enrichment = run_enrichment_tests(sig_genes_in_universe, universe)

    if len(enrichment) > 0:
        # BH correction
        from statsmodels.stats.multitest import multipletests
        pvals = enrichment["fisher_pval"].fillna(1.0).values
        _, padj, _, _ = multipletests(pvals, method="fdr_bh")
        enrichment["fisher_padj_bh"] = padj
        save_csv(enrichment, "coloc_signaling_enrichment.csv", subdir="commot")

    # ── Drug-target signaling cards ──
    print_header("Drug-Target Signaling Cards")
    cards = build_drug_target_cards(sig_df, sig_genes)
    if len(cards) > 0:
        save_csv(cards, "drug_target_signaling_cards.csv", subdir="commot")
        print(f"  Drug targets with signaling regulation: {len(cards)}")
        for _, row in cards.head(10).iterrows():
            drugs = row.get("drugs", "")
            pathway = row.get("commot_top_pathway", row.get("commot_pathway", ""))
            print(f"    {row['gene']}: pathway={pathway}"
                  f"{', drugs=' + drugs if drugs else ''}")

    # ── Convergence summary ──
    print_header("Convergence Summary")
    summary_records = []

    summary_records.append({
        "metric": "n_signaling_regulated_genes",
        "value": len(sig_genes),
    })
    summary_records.append({
        "metric": "n_in_atlas_universe",
        "value": len(sig_genes_in_universe),
    })

    if len(enrichment) > 0:
        for _, row in enrichment.iterrows():
            summary_records.append({
                "metric": f"overlap_{row['gene_set']}",
                "value": row["n_overlap"],
            })
            summary_records.append({
                "metric": f"OR_{row['gene_set']}",
                "value": row["odds_ratio"],
            })
            summary_records.append({
                "metric": f"pval_{row['gene_set']}",
                "value": row["fisher_pval"],
            })

    summary_records.append({
        "metric": "n_drug_target_cards",
        "value": len(cards),
    })

    summary_df = pd.DataFrame(summary_records)
    save_csv(summary_df, "convergence_summary.csv", subdir="commot")

    for _, row in summary_df.iterrows():
        val = row["value"]
        if isinstance(val, float) and val < 0.01:
            print(f"    {row['metric']}: {val:.2e}")
        else:
            print(f"    {row['metric']}: {val}")

    print_header("16g: Complete")


if __name__ == "__main__":
    main()
