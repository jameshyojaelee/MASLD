#!/usr/bin/env python3
"""
63c_combine_otters_results.py
---------------------------------------------------------------------------
Combine per-chromosome OTTERS TWAS results with ACAT p-value combination.
Run this if the OTTERS TWAS job (63_otters_twas.sh) completes chromosomes
but doesn't reach the ACAT combination step.

Input:
  - RNA-seq/results/causal_inference/otters_broadaway/{gwas}/chr*/P0.05.txt
  - RNA-seq/results/causal_inference/otters_broadaway/{gwas}/chr*/P0.001.txt
  - data/broadaway_eqtl/otters_format/gene_anno.txt

Output:
  - RNA-seq/results/causal_inference/otters_broadaway/{gwas}/otters_twas_combined.csv

Usage:
  python RNA-seq/63c_combine_otters_results.py [--gwas ghodsian]
---------------------------------------------------------------------------
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

ANNO_FILE = os.path.join(BASE, "data/broadaway_eqtl/otters_format/gene_anno.txt")
OTTERS_DIR = os.path.join(BASE, "RNA-seq/results/causal_inference/otters_broadaway")

GWAS_LIST = ["ghodsian", "ukbb_alt", "ukbb_ast", "ukbb_ggt"]


def acat_pval(pvals):
    """ACAT (Cauchy combination) of p-values."""
    pvals = np.array([p for p in pvals if np.isfinite(p) and 0 < p < 1])
    if len(pvals) == 0:
        return np.nan
    if len(pvals) == 1:
        return pvals[0]
    # Cauchy combination
    t_stats = np.tan((0.5 - pvals) * np.pi)
    t_acat = np.mean(t_stats)
    # Convert back to p-value
    p_acat = 0.5 - np.arctan(t_acat) / np.pi
    return p_acat


def combine_gwas(gwas_name):
    """Combine all chromosome results for one GWAS."""
    gwas_dir = os.path.join(OTTERS_DIR, gwas_name)
    if not os.path.isdir(gwas_dir):
        print(f"  {gwas_name}: directory not found, skipping")
        return None

    # Load gene annotations
    anno = pd.read_csv(ANNO_FILE, sep="\t")
    gene_map = dict(zip(anno["TargetID"], anno["GeneName"]))

    # Collect per-chromosome, per-model results
    model_results = {}
    for model in ["P0.05", "P0.001"]:
        all_chr = []
        for chrom in range(1, 23):
            fpath = os.path.join(gwas_dir, f"chr{chrom}", f"{model}.txt")
            if not os.path.exists(fpath):
                continue
            try:
                df = pd.read_csv(fpath, sep="\t")
                if len(df) > 0:
                    df["chr"] = chrom
                    df["model"] = model
                    all_chr.append(df)
            except Exception:
                pass
        if all_chr:
            model_results[model] = pd.concat(all_chr, ignore_index=True)

    if not model_results:
        print(f"  {gwas_name}: no results found")
        return None

    # For each gene, collect p-values from all models and ACAT-combine
    all_genes = set()
    for model, df in model_results.items():
        all_genes.update(df["TargetID"].unique())

    results = []
    for gene_id in sorted(all_genes):
        gene_pvals = {}
        gene_z = {}
        gene_n_snps = {}
        gene_chr = None

        for model, df in model_results.items():
            gene_df = df[df["TargetID"] == gene_id]
            if len(gene_df) > 0:
                row = gene_df.iloc[0]
                try:
                    pval = float(row["FUSION_PVAL"])
                    z = float(row["FUSION_Z"])
                    nsnps = int(float(row["n_snps"]))
                except (ValueError, TypeError):
                    continue
                if not np.isfinite(pval) or pval <= 0 or pval >= 1:
                    continue
                gene_pvals[model] = pval
                gene_z[model] = z
                gene_n_snps[model] = nsnps
                gene_chr = row["chr"] if gene_chr is None else gene_chr

        # ACAT combination across models
        pvals_list = list(gene_pvals.values())
        if len(pvals_list) == 0:
            continue  # skip genes with no valid p-values
        acat_p = acat_pval(pvals_list)

        # Best z-score (from best model)
        best_model = min(gene_pvals, key=gene_pvals.get)
        best_z = gene_z[best_model]

        results.append({
            "TargetID": gene_id,
            "gene_symbol": gene_map.get(gene_id, ""),
            "chr": gene_chr,
            "n_models": len(gene_pvals),
            "best_z": best_z,
            "best_pval": min(pvals_list),
            "acat_pval": acat_p,
            "P0.05_pval": gene_pvals.get("P0.05", np.nan),
            "P0.05_z": gene_z.get("P0.05", np.nan),
            "P0.001_pval": gene_pvals.get("P0.001", np.nan),
            "P0.001_z": gene_z.get("P0.001", np.nan),
            "n_snps": max(gene_n_snps.values()),
        })

    combined = pd.DataFrame(results)

    # FDR correction on ACAT p-values
    # Replace exact 0 with smallest positive float for BH correction
    combined["acat_pval_adj"] = combined["acat_pval"].clip(lower=np.finfo(float).tiny)
    valid = combined["acat_pval_adj"].notna() & (combined["acat_pval_adj"] > 0)
    combined.loc[valid, "fdr"] = multipletests(
        combined.loc[valid, "acat_pval_adj"].values, method="fdr_bh"
    )[1]
    combined = combined.drop(columns=["acat_pval_adj"])

    combined = combined.sort_values("acat_pval")

    # Save
    out_file = os.path.join(gwas_dir, "otters_twas_combined.csv")
    combined.to_csv(out_file, index=False)

    n_sig = (combined["fdr"] < 0.05).sum()
    n_total = len(combined)
    print(f"  {gwas_name}: {n_total} genes, {n_sig} FDR<0.05")
    print(f"    Saved: {out_file}")

    return combined


def main():
    parser = argparse.ArgumentParser(description="Combine OTTERS TWAS results")
    parser.add_argument("--gwas", default=None, help="Specific GWAS (default: all)")
    args = parser.parse_args()

    gwas_list = [args.gwas] if args.gwas else GWAS_LIST

    print("=== Combining OTTERS TWAS Results ===")
    for gwas_name in gwas_list:
        print(f"\n--- {gwas_name} ---")
        combine_gwas(gwas_name)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
