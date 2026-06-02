#!/usr/bin/env python3
"""
61_otters_format_broadaway.py
Phase 2A-Step1: Reformat Broadaway liver eQTL summary statistics into OTTERS input format.

Input:
  - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv  (Broadaway marginal eQTL, hg19)
  - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv  (lead SNPs for gene annotation)
  - data/1kg_eur/chr{1-22}_eur.pvar                              (1KG EUR Phase 1, hg19)

Output (to data/broadaway_eqtl/otters_format/):
  - chr{1-22}_broadaway.txt   — OTTERS sseQTL format (CHROM, POS, A1, A2, Zscore, TargetID, N)
  - gene_anno.txt             — Gene annotation (CHROM, GeneStart, GeneEnd, TargetID, GeneName)
  - snp_match_chr{1-22}.txt   — Variant ID matching (broadaway_id, otters_id, flip)

Author: James Lee
Date: 2026-03-24
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from pathlib import Path


def get_project_root():
    """Get project root from env or default."""
    return os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
    )


def load_pvar(pvar_path, chrom):
    """Load 1KG EUR pvar file, extracting SNP-level info for matching.

    Returns DataFrame with columns: POS, REF, ALT, pvar_id
    Only keeps biallelic SNPs (single-character REF and ALT).
    """
    # pvar has ##INFO comment lines, then a #CHROM header line, then data.
    # comment="#" skips both ## lines and the #CHROM header, so we assign names manually.
    df = pd.read_csv(
        pvar_path,
        sep="\t",
        comment="#",
        header=None,
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"],
        usecols=[1, 2, 3, 4],  # POS, ID, REF, ALT by position
    )
    df.columns = ["POS", "pvar_id", "REF", "ALT"]
    df["POS"] = df["POS"].astype(np.int64)
    df["pvar_id"] = df["pvar_id"].astype(str)
    df["REF"] = df["REF"].astype(str)
    df["ALT"] = df["ALT"].astype(str)

    # Filter to biallelic SNPs only (single nucleotide REF and ALT)
    snp_mask = (df["REF"].str.len() == 1) & (df["ALT"].str.len() == 1)
    df = df[snp_mask].copy()

    # Drop duplicate positions (rare multi-allelic sites collapsed by plink2)
    df = df.drop_duplicates(subset=["POS", "REF", "ALT"], keep="first")

    return df


def load_broadaway_chr(marginal_path, chrom, min_n=100):
    """Load Broadaway marginal summary stats for one chromosome.

    Applies QC filters:
      - N >= min_n
      - SE > 0
      - Biallelic SNPs only (single nucleotide EA and NEA)
    """
    df = pd.read_csv(
        marginal_path,
        sep="\t",
        dtype=str,  # read everything as string first to avoid mixed-type issues
    )
    # Cast numeric columns
    for col in ["POS", "N"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in ["EAF", "Beta", "SE", "PVAL"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    # Drop rows where critical columns are NA after coercion
    df = df.dropna(subset=["POS", "N", "Beta", "SE"])

    n_raw = len(df)

    # Filter: N >= min_n
    df = df[df["N"] >= min_n]
    n_after_n = len(df)

    # Filter: SE > 0 (avoid division by zero in Z-score)
    df = df[df["SE"] > 0]
    n_after_se = len(df)

    # Filter: biallelic SNPs only (single nucleotide alleles)
    snp_mask = (df["EA"].str.len() == 1) & (df["NEA"].str.len() == 1)
    df = df[snp_mask]
    n_after_snp = len(df)

    # Compute Z-score
    df["Zscore"] = df["Beta"] / df["SE"]

    # Drop any rows with NaN Z-score (shouldn't happen after SE>0 filter, but be safe)
    df = df.dropna(subset=["Zscore"])

    print(f"  chr{chrom} Broadaway: {n_raw:,} raw -> {n_after_n:,} (N>={min_n}) "
          f"-> {n_after_se:,} (SE>0) -> {n_after_snp:,} (SNPs) -> {len(df):,} (final)")

    return df


def match_variants(broadaway_df, pvar_df, chrom):
    """Match Broadaway variants to 1KG EUR by position and alleles.

    Handles allele flipping: if Broadaway EA/NEA match 1KG ALT/REF, no flip needed.
    If Broadaway EA/NEA match 1KG REF/ALT (reversed), flip the Z-score sign.

    Returns:
      - matched_df: Broadaway rows with matched 1KG IDs and possibly flipped Z-scores
      - match_df: mapping table (broadaway_id, otters_id, flip)
    """
    # Create a lookup from pvar: keyed by POS
    # For efficiency, merge on POS first, then check alleles
    merged = broadaway_df.merge(pvar_df, on="POS", how="inner", suffixes=("", "_1kg"))

    n_pos_match = len(merged)

    # --- Step 1: Direct allele matching (same strand) ---
    # Forward: Broadaway EA == 1KG ALT and NEA == REF (no Z-score flip)
    direct_fwd = (merged["EA"] == merged["ALT"]) & (merged["NEA"] == merged["REF"])

    # Flipped: Broadaway EA == 1KG REF and NEA == ALT (flip Z-score)
    direct_flip = (merged["EA"] == merged["REF"]) & (merged["NEA"] == merged["ALT"])

    direct_any = direct_fwd | direct_flip

    # --- Step 2: Complement matching (opposite strand) ---
    # Only attempt for variants NOT already matched directly.
    # Strand-ambiguous SNPs (A/T or C/G) are indeterminate under complement —
    # a direct forward match looks identical to a complement flip — so we
    # exclude them from complement matching entirely.
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def complement(s):
        return s.map(comp)

    # Identify strand-ambiguous SNPs (A/T or C/G pairs in 1KG)
    ambiguous = ((merged["REF"] == "A") & (merged["ALT"] == "T")) | \
                ((merged["REF"] == "T") & (merged["ALT"] == "A")) | \
                ((merged["REF"] == "C") & (merged["ALT"] == "G")) | \
                ((merged["REF"] == "G") & (merged["ALT"] == "C"))

    # Complement forward: comp(EA)==ALT and comp(NEA)==REF
    comp_fwd = (~direct_any) & (~ambiguous) & \
               (complement(merged["EA"]) == merged["ALT"]) & \
               (complement(merged["NEA"]) == merged["REF"])

    # Complement flipped: comp(EA)==REF and comp(NEA)==ALT
    comp_flip = (~direct_any) & (~ambiguous) & \
                (complement(merged["EA"]) == merged["REF"]) & \
                (complement(merged["NEA"]) == merged["ALT"])

    # --- Step 3: Combine ---
    is_fwd = direct_fwd | comp_fwd
    is_flip = direct_flip | comp_flip
    any_match = is_fwd | is_flip

    matched = merged[any_match].copy()

    # Assign flip flag: 1 = Z-score sign is negated
    matched["flip"] = is_flip[any_match].astype(int).values
    matched.loc[matched["flip"] == 1, "Zscore"] = -matched.loc[matched["flip"] == 1, "Zscore"]

    # Build the OTTERS-style A1/A2 columns (A1=effect allele in 1KG frame = ALT, A2=REF)
    matched["A1"] = matched["ALT"]
    matched["A2"] = matched["REF"]

    # Drop duplicates: same POS + same gene might appear if pvar has dups
    matched = matched.drop_duplicates(subset=["POS", "ENSG"], keep="first")

    n_direct_fwd = direct_fwd.sum()
    n_direct_flip = direct_flip.sum()
    n_comp_fwd = comp_fwd.sum()
    n_comp_flip = comp_flip.sum()
    n_ambiguous_skipped = (ambiguous & ~direct_any).sum()

    print(f"  chr{chrom} matching: {n_pos_match:,} position matches -> "
          f"{len(matched):,} allele matches "
          f"(direct_fwd={n_direct_fwd:,}, direct_flip={n_direct_flip:,}, "
          f"comp_fwd={n_comp_fwd:,}, comp_flip={n_comp_flip:,}, "
          f"ambiguous_unmatched={n_ambiguous_skipped:,})")

    # Build match table
    # Broadaway ID format: CHR_POS_NEA_EA
    match_records = pd.DataFrame({
        "broadaway_id": matched["Variant"],
        "otters_id": matched["pvar_id"],
        "flip": matched["flip"],
    })

    return matched, match_records


def build_gene_anno(all_matched, leads_path):
    """Build gene annotation file from matched data.

    Uses the leads file for gene name mapping, and computes gene boundaries
    from min/max POS of tested variants in the marginal data (cis-window proxy).
    For a more precise annotation, we also attempt to extract any TSS/TES
    from the leads file, but since it lacks those columns, we use the
    marginal cis-window boundaries.
    """
    # Load leads file for gene name <-> ENSG mapping
    leads = pd.read_csv(leads_path, sep="\t", usecols=["Gene", "Ensembl"])
    leads = leads.drop_duplicates(subset=["Ensembl"])
    leads_map = dict(zip(leads["Ensembl"], leads["Gene"]))

    # Aggregate from matched marginal data
    gene_agg = all_matched.groupby(["ENSG"]).agg(
        CHROM=("CHR", "first"),
        GeneStart=("POS", "min"),
        GeneEnd=("POS", "max"),
        GeneName=("GeneSymbol", "first"),
    ).reset_index()

    # Use leads gene name if available (more reliable than marginal GeneSymbol)
    gene_agg["GeneName"] = gene_agg["ENSG"].map(leads_map).fillna(gene_agg["GeneName"])

    # Rename for OTTERS format
    gene_agg = gene_agg.rename(columns={"ENSG": "TargetID"})

    # Filter out genes with no name
    gene_agg = gene_agg.dropna(subset=["GeneName"])

    # Ensure integer coordinates
    gene_agg["GeneStart"] = gene_agg["GeneStart"].astype(int)
    gene_agg["GeneEnd"] = gene_agg["GeneEnd"].astype(int)

    # Sort by CHROM (numeric) and GeneStart
    gene_agg["CHROM"] = gene_agg["CHROM"].astype(str)
    gene_agg["_chrom_num"] = gene_agg["CHROM"].astype(int)
    gene_agg = gene_agg.sort_values(["_chrom_num", "GeneStart"]).drop(columns=["_chrom_num"])

    # Select output columns
    gene_agg = gene_agg[["CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName"]]

    return gene_agg


def process_chromosome(chrom, eqtl_dir, ld_dir, out_dir, min_n=100):
    """Process a single chromosome: load, filter, match, write.

    Returns:
      - matched DataFrame (for gene_anno aggregation)
      - summary dict with stats
    """
    marginal_path = os.path.join(eqtl_dir, f"chr{chrom}_marginal_summary_results.tsv")
    pvar_path = os.path.join(ld_dir, f"chr{chrom}_eur.pvar")

    # Check files exist
    if not os.path.exists(marginal_path):
        print(f"  WARNING: {marginal_path} not found, skipping chr{chrom}")
        return None, None
    if not os.path.exists(pvar_path):
        print(f"  WARNING: {pvar_path} not found, skipping chr{chrom}")
        return None, None

    # Load data
    t0 = time.time()
    broadaway = load_broadaway_chr(marginal_path, chrom, min_n=min_n)
    pvar = load_pvar(pvar_path, chrom)
    t_load = time.time() - t0

    n_broadaway = len(broadaway)
    n_pvar = len(pvar)
    n_genes_in = broadaway["ENSG"].nunique()

    # Match variants
    t0 = time.time()
    matched, match_records = match_variants(broadaway, pvar, chrom)
    t_match = time.time() - t0

    n_matched = len(matched)
    n_genes_out = matched["ENSG"].nunique()
    match_rate = n_matched / n_broadaway * 100 if n_broadaway > 0 else 0

    # Write OTTERS sseQTL file
    otters_df = matched[["CHR", "POS", "A1", "A2", "Zscore", "ENSG", "N"]].copy()
    otters_df = otters_df.rename(columns={"CHR": "CHROM", "ENSG": "TargetID"})
    # Sort by TargetID then POS for clean output
    otters_df = otters_df.sort_values(["TargetID", "POS"])

    otters_path = os.path.join(out_dir, f"chr{chrom}_broadaway.txt")
    otters_df.to_csv(otters_path, sep="\t", index=False)

    # Write SNP match file
    match_path = os.path.join(out_dir, f"snp_match_chr{chrom}.txt")
    match_records.to_csv(match_path, sep="\t", index=False)

    stats = {
        "chrom": chrom,
        "n_broadaway_snps": n_broadaway,
        "n_pvar_snps": n_pvar,
        "n_matched": n_matched,
        "match_rate_pct": round(match_rate, 1),
        "n_genes_in": n_genes_in,
        "n_genes_out": n_genes_out,
        "load_time_s": round(t_load, 1),
        "match_time_s": round(t_match, 1),
    }

    print(f"  chr{chrom} DONE: {n_matched:,} variants, {n_genes_out} genes, "
          f"{match_rate:.1f}% match rate ({t_load + t_match:.1f}s)")

    return matched, stats


def main():
    parser = argparse.ArgumentParser(
        description="Reformat Broadaway liver eQTL summary stats to OTTERS input format"
    )
    parser.add_argument(
        "--project-root", type=str, default=None,
        help="Project root directory (default: MASLD_PROJECT_ROOT env or standard path)"
    )
    parser.add_argument(
        "--min-n", type=int, default=100,
        help="Minimum sample size to keep a variant (default: 100)"
    )
    parser.add_argument(
        "--chromosomes", type=str, default="1-22",
        help="Chromosome range to process (default: 1-22)"
    )
    args = parser.parse_args()

    # Resolve paths
    project_root = args.project_root or get_project_root()
    eqtl_dir = os.path.join(project_root, "data", "broadaway_eqtl")
    ld_dir = os.path.join(project_root, "data", "1kg_eur")
    leads_path = os.path.join(eqtl_dir, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv")
    out_dir = os.path.join(eqtl_dir, "otters_format")

    # Parse chromosome range
    if "-" in args.chromosomes:
        start, end = args.chromosomes.split("-")
        chroms = list(range(int(start), int(end) + 1))
    else:
        chroms = [int(c) for c in args.chromosomes.split(",")]

    print("=" * 80)
    print("OTTERS FORMAT — Broadaway Liver eQTL")
    print("=" * 80)
    print(f"Project root:  {project_root}")
    print(f"eQTL dir:      {eqtl_dir}")
    print(f"LD ref dir:    {ld_dir}")
    print(f"Output dir:    {out_dir}")
    print(f"Chromosomes:   {chroms}")
    print(f"Min N:         {args.min_n}")
    print("=" * 80)

    # Create output directory
    os.makedirs(out_dir, exist_ok=True)

    # Process each chromosome
    all_matched = []
    all_stats = []
    t_total = time.time()

    for chrom in chroms:
        print(f"\n--- Chromosome {chrom} ---")
        matched, stats = process_chromosome(chrom, eqtl_dir, ld_dir, out_dir, min_n=args.min_n)
        if matched is not None and len(matched) > 0:
            # Keep only columns needed for gene_anno to limit memory
            all_matched.append(matched[["CHR", "POS", "ENSG", "GeneSymbol"]].copy())
            all_stats.append(stats)

    # Build gene annotation from all chromosomes
    if all_matched:
        print("\n--- Building gene annotation ---")
        combined = pd.concat(all_matched, ignore_index=True)
        gene_anno = build_gene_anno(combined, leads_path)
        gene_anno_path = os.path.join(out_dir, "gene_anno.txt")
        gene_anno.to_csv(gene_anno_path, sep="\t", index=False)
        print(f"  Gene annotation: {len(gene_anno)} genes written to {gene_anno_path}")

    # Print summary
    t_elapsed = time.time() - t_total
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        total_broadaway = stats_df["n_broadaway_snps"].sum()
        total_matched = stats_df["n_matched"].sum()
        total_genes = stats_df["n_genes_out"].sum()
        overall_match_rate = total_matched / total_broadaway * 100 if total_broadaway > 0 else 0

        print(f"\n{'Chr':<6} {'Broadaway':>12} {'1KG':>10} {'Matched':>10} {'Rate':>8} {'Genes':>8}")
        print("-" * 60)
        for s in all_stats:
            print(f"{s['chrom']:<6} {s['n_broadaway_snps']:>12,} {s['n_pvar_snps']:>10,} "
                  f"{s['n_matched']:>10,} {s['match_rate_pct']:>7.1f}% {s['n_genes_out']:>8}")
        print("-" * 60)
        print(f"{'TOTAL':<6} {total_broadaway:>12,} {'':>10} "
              f"{total_matched:>10,} {overall_match_rate:>7.1f}% {total_genes:>8}")

        # Unique genes across all chromosomes
        if all_matched:
            unique_genes = combined["ENSG"].nunique()
            print(f"\nUnique genes (all chromosomes): {unique_genes}")

    print(f"\nTotal elapsed time: {t_elapsed:.1f}s")
    print(f"Output directory: {out_dir}")
    print("\nOutput files:")
    for f in sorted(os.listdir(out_dir)):
        fpath = os.path.join(out_dir, f)
        fsize = os.path.getsize(fpath) / (1024 * 1024)
        print(f"  {f}  ({fsize:.1f} MB)")

    print("\nDone.")


if __name__ == "__main__":
    main()
