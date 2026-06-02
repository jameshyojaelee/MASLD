#!/usr/bin/env python3
"""
09_finngen_credset_comparison.py
Compare our SuSiE/CARMA fine-mapping results against FinnGen R12 credible sets
retrieved from the public PheWeb API.

Strategy:
  1. Pull FinnGen R12 credible sets via https://r12.finngen.fi/api/autoreport/{pheno}
  2. Build hg38→hg19 position map from our raw FinnGen sumstats (hg38 'bp' column)
     and reformatted sumstats (hg19 'position' column), matched by alleles
  3. Match variants by chr:pos_hg19:alleles against our combined_finemapping.csv
  4. Compute PIP correlation, credible set overlap, and concordance metrics

Usage:
  python3 09_finngen_credset_comparison.py

Output:
  results/finngen_credset_comparison.csv   — per-variant matched PIPs
  results/finngen_credset_summary.csv      — per-locus concordance summary
  figures/finngen_credset_comparison.pdf    — correlation + CS overlap plots (if matplotlib)
"""

import json
import os
import sys
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BASE_DIR = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
FM_DIR = os.path.join(BASE_DIR, "GWAS/finemapping")
RESULTS_DIR = os.path.join(FM_DIR, "results")
FIG_DIR = os.path.join(FM_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Our raw FinnGen sumstats (hg38) and reformatted (hg19)
RAW_DIR = os.path.join(FM_DIR, "GWAS/finemapping")  # not used — see below
SUMSTATS_DIR = os.path.join(FM_DIR, "data/sumstats")

# FinnGen phenocodes mapped to our study names
PHENO_MAP = {
    "NAFLD": "FinnGen_NAFLD",
    "CHIRHEP_NAS": "FinnGen_NASH",
    "C3_HEPATOCELLU_CARC_EXALLC": "FinnGen_HCC",
}

# FinnGen raw sumstats paths (hg38, FinnGen-format columns)
FINNGEN_RAW = {
    "NAFLD": os.path.join(
        BASE_DIR, "GWAS/MR_Data/FinnGen/finngen_R12_NAFLD.gz"
    ),
    "CHIRHEP_NAS": os.path.join(
        BASE_DIR, "GWAS/MR_Data/FinnGen/finngen_R12_CHIRHEP_NAS.gz"
    ),
    "C3_HEPATOCELLU_CARC_EXALLC": os.path.join(
        BASE_DIR, "GWAS/MR_Data/FinnGen/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz"
    ),
}

# Our reformatted FinnGen sumstats (hg19)
FINNGEN_HG19 = {
    "NAFLD": os.path.join(SUMSTATS_DIR, "FinnGen_NAFLD_reformatted_hg19.tsv"),
    "CHIRHEP_NAS": os.path.join(SUMSTATS_DIR, "FinnGen_NASH_reformatted_hg19.tsv"),
    "C3_HEPATOCELLU_CARC_EXALLC": os.path.join(
        SUMSTATS_DIR, "FinnGen_HCC_reformatted_hg19.tsv"
    ),
}


# --------------------------------------------------------------------------
# 1. Pull FinnGen credible sets from API
# --------------------------------------------------------------------------
def fetch_finngen_credsets(phenocode):
    """Fetch credible sets from FinnGen R12 PheWeb API."""
    url = f"https://r12.finngen.fi/api/autoreport/{phenocode}"
    print(f"  Fetching {url} ...")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            data = json.loads(r.read())
    except Exception as e:
        print(f"  ERROR fetching {phenocode}: {e}")
        return []

    rows = []
    for locus in data:
        locus_id = locus.get("locus_id", "")
        gene = locus.get("lead_most_severe_gene", "NA")
        good_cs = locus.get("good_cs", False)
        cs_size = locus.get("cs_size", 0)
        cs_bf = locus.get("cs_log_bayes_factor", np.nan)
        pval = locus.get("pval", np.nan)

        cs_str = locus.get("credible_set_variants", "")
        if not cs_str:
            continue

        for entry in cs_str.split(";"):
            parts = entry.split("|")
            if len(parts) < 2:
                continue
            variant = parts[0]  # chrN_pos_ref_alt
            pip = float(parts[1])

            # Parse variant: chr19_19282905_G_A -> chr=19, pos=19282905, ref=G, alt=A
            vparts = variant.split("_")
            if len(vparts) < 4:
                continue
            chrom = vparts[0].replace("chr", "")
            pos_hg38 = int(vparts[1])
            ref = vparts[2]
            alt = vparts[3]

            rows.append(
                {
                    "phenocode": phenocode,
                    "locus_id": locus_id,
                    "gene": gene,
                    "good_cs": good_cs,
                    "cs_size": cs_size,
                    "cs_log_bf": cs_bf,
                    "locus_pval": pval,
                    "chromosome": int(chrom),
                    "position_hg38": pos_hg38,
                    "ref": ref,
                    "alt": alt,
                    "finngen_pip": pip,
                    "variant_hg38": variant,
                }
            )

    print(f"  {phenocode}: {len(data)} loci, {len(rows)} CS variants")
    return rows


# --------------------------------------------------------------------------
# 2. Build hg38 -> hg19 position map from our sumstats
# --------------------------------------------------------------------------
def build_position_map(phenocode, target_positions):
    """
    Build hg38 -> hg19 position map for specific target variants only.

    Instead of loading the full ~20M-row FinnGen sumstats, we stream through the
    raw file and only keep rows matching our target hg38 positions. This reduces
    memory from ~9GB to <100MB.

    Args:
        phenocode: FinnGen phenotype code
        target_positions: set of (chromosome, position_hg38) tuples to map
    """
    raw_path = FINNGEN_RAW.get(phenocode)
    hg19_path = FINNGEN_HG19.get(phenocode)

    if not raw_path or not hg19_path:
        print(f"  No sumstats paths for {phenocode}")
        return {}

    if not os.path.exists(raw_path):
        print(f"  Raw sumstats not found: {raw_path}")
        return {}
    if not os.path.exists(hg19_path):
        print(f"  hg19 sumstats not found: {hg19_path}")
        return {}

    print(f"  Building hg38->hg19 map for {phenocode} ({len(target_positions)} targets)...")

    # Step 1: Stream raw file (hg38), keep only target positions + their betas
    import gzip

    raw_matches = {}  # (chr, pos_hg38) -> (ref, alt, beta)
    opener = gzip.open if raw_path.endswith(".gz") else open
    with opener(raw_path, "rt") as f:
        header = f.readline().strip().split("\t")
        idx_chrom = header.index("#chrom")
        idx_pos = header.index("pos")
        idx_ref = header.index("ref")
        idx_alt = header.index("alt")
        idx_beta = header.index("beta")

        for line in f:
            parts = line.strip().split("\t")
            try:
                chrom = int(parts[idx_chrom])
                pos = int(parts[idx_pos])
            except (ValueError, IndexError):
                continue

            if (chrom, pos) in target_positions:
                raw_matches[(chrom, pos)] = (
                    parts[idx_ref],
                    parts[idx_alt],
                    round(float(parts[idx_beta]), 6),
                )

    print(f"    Found {len(raw_matches)} target variants in raw hg38 file")

    if not raw_matches:
        return {}

    # Step 2: Stream hg19 file, match by chr + beta + alleles
    # Build lookup: (chr, beta_round, allele1, allele2) -> pos_hg19
    hg19_lookup = {}
    with open(hg19_path, "rt") as f:
        header = f.readline().strip().split("\t")
        idx_chrom = header.index("chromosome")
        idx_pos = header.index("position")
        idx_a1 = header.index("allele1")
        idx_a2 = header.index("allele2")
        idx_beta = header.index("beta")

        for line in f:
            parts = line.strip().split("\t")
            try:
                chrom = int(parts[idx_chrom])
                pos_hg19 = int(parts[idx_pos])
                a1 = parts[idx_a1]
                a2 = parts[idx_a2]
                beta = round(float(parts[idx_beta]), 6)
            except (ValueError, IndexError):
                continue
            key = (chrom, beta, a1, a2)
            hg19_lookup[key] = pos_hg19

    print(f"    Loaded {len(hg19_lookup)} hg19 variants for matching")

    # Step 3: Match raw -> hg19
    pos_map = {}
    for (chrom, pos_hg38), (ref, alt, beta) in raw_matches.items():
        # In our convention: allele1=alt (effect), allele2=ref
        key = (chrom, beta, alt, ref)
        if key in hg19_lookup:
            pos_map[(chrom, pos_hg38)] = hg19_lookup[key]

    print(f"    Mapped {len(pos_map)}/{len(raw_matches)} positions (hg38->hg19)")
    return pos_map


# --------------------------------------------------------------------------
# 3. Load our finemapping results for FinnGen studies
# --------------------------------------------------------------------------
def load_our_finemapping():
    """Load combined_finemapping.csv filtered to FinnGen studies."""
    path = os.path.join(RESULTS_DIR, "combined_finemapping.csv")
    print(f"Loading our finemapping from {path}...")

    # Only load FinnGen rows to save memory
    chunks = []
    for chunk in pd.read_csv(path, chunksize=500_000):
        fg = chunk[chunk["study"].str.startswith("FinnGen_")]
        if len(fg) > 0:
            chunks.append(fg)

    if not chunks:
        print("  No FinnGen entries found!")
        return pd.DataFrame()

    df = pd.concat(chunks, ignore_index=True)
    print(f"  Loaded {len(df)} FinnGen variants across {df['study'].nunique()} studies")
    return df


# --------------------------------------------------------------------------
# 4. Match and compare
# --------------------------------------------------------------------------
def compare_pips(finngen_cs, pos_map, our_fm, study_name):
    """Match FinnGen CS variants to our finemapping by hg19 position + alleles."""
    our_study = our_fm[our_fm["study"] == study_name].copy()
    if our_study.empty:
        return pd.DataFrame()

    # Index our data by (chr, pos_hg19)
    our_idx = {}
    for _, row in our_study.iterrows():
        key = (int(row["chromosome"]), int(row["position"]))
        # Store as list to handle potential duplicates
        if key not in our_idx:
            our_idx[key] = []
        our_idx[key].append(row)

    matched = []
    for _, fg in finngen_cs.iterrows():
        chrom = fg["chromosome"]
        pos_hg38 = fg["position_hg38"]

        # Liftover hg38 -> hg19
        pos_hg19 = pos_map.get((chrom, pos_hg38))
        if pos_hg19 is None:
            continue

        # Look up in our data
        key = (chrom, pos_hg19)
        if key not in our_idx:
            continue

        for our_row in our_idx[key]:
            # Check allele match (FinnGen alt = our allele1, FinnGen ref = our allele2)
            allele_ok = (fg["alt"] == our_row["allele1"] and fg["ref"] == our_row["allele2"]) or \
                        (fg["ref"] == our_row["allele1"] and fg["alt"] == our_row["allele2"])
            if not allele_ok:
                continue

            matched.append(
                {
                    "chromosome": chrom,
                    "position_hg19": pos_hg19,
                    "position_hg38": pos_hg38,
                    "allele1": our_row["allele1"],
                    "allele2": our_row["allele2"],
                    "finngen_phenocode": fg["phenocode"],
                    "our_study": study_name,
                    "locus_gene": fg["gene"],
                    "locus_pval": fg["locus_pval"],
                    "good_cs": fg["good_cs"],
                    "finngen_pip": fg["finngen_pip"],
                    "our_susie_pip": our_row["susie_pip"],
                    "our_carma_pip": our_row["carma_pip"],
                    "our_recommended_pip": our_row["recommended_pip"],
                    "our_susie_cs": our_row["susie_cs"],
                    "our_carma_cs": our_row["carma_cs"],
                    "our_concordant_pip": our_row["concordant_pip"],
                    "finngen_in_cs": True,  # all FinnGen variants here are in a CS
                    "our_in_susie_cs": our_row["susie_cs"] > 0,
                    "our_in_carma_cs": our_row["carma_cs"] > 0,
                    "our_in_either_cs": our_row["either_in_cs"]
                    if "either_in_cs" in our_row.index
                    else (our_row["susie_cs"] > 0 or our_row["carma_cs"] > 0),
                }
            )

    return pd.DataFrame(matched)


# --------------------------------------------------------------------------
# 5. Compute summary stats
# --------------------------------------------------------------------------
def compute_summary(matched_df):
    """Per-locus concordance metrics."""
    if matched_df.empty:
        return pd.DataFrame()

    summaries = []
    for (pheno, gene), grp in matched_df.groupby(
        ["finngen_phenocode", "locus_gene"]
    ):
        n = len(grp)
        good_cs = grp["good_cs"].iloc[0]

        # Drop rows with NaN PIPs for correlation
        grp_clean = grp.dropna(subset=["finngen_pip", "our_susie_pip", "our_carma_pip", "our_recommended_pip"])
        n_clean = len(grp_clean)

        # PIP correlation
        rho_susie = (
            grp_clean["finngen_pip"].corr(grp_clean["our_susie_pip"]) if n_clean >= 3 else np.nan
        )
        rho_carma = (
            grp_clean["finngen_pip"].corr(grp_clean["our_carma_pip"]) if n_clean >= 3 else np.nan
        )
        rho_rec = (
            grp_clean["finngen_pip"].corr(grp_clean["our_recommended_pip"]) if n_clean >= 3 else np.nan
        )

        # Top variant agreement (use clean subset)
        # FIX (review A06#2): initialize before the conditional branches so the
        # summaries.append below never raises UnboundLocalError when PIPs are all-NaN
        # (a normal SuSiE/CARMA non-convergence condition).
        fg_top = None
        our_top_susie = None
        our_top_carma = None
        if grp_clean.empty:
            top_agree_susie = False
            top_agree_carma = False
        else:
            fg_top = grp_clean.loc[grp_clean["finngen_pip"].idxmax()]
            susie_valid = grp_clean["our_susie_pip"].notna()
            carma_valid = grp_clean["our_carma_pip"].notna()
            if susie_valid.any():
                our_top_susie = grp_clean.loc[grp_clean.loc[susie_valid, "our_susie_pip"].idxmax()]
                top_agree_susie = fg_top["position_hg19"] == our_top_susie["position_hg19"]
            else:
                top_agree_susie = False
            if carma_valid.any():
                our_top_carma = grp_clean.loc[grp_clean.loc[carma_valid, "our_carma_pip"].idxmax()]
                top_agree_carma = fg_top["position_hg19"] == our_top_carma["position_hg19"]
            else:
                top_agree_carma = False

        # CS overlap (Jaccard): FinnGen CS variants vs our CS variants
        fg_cs_set = set(grp["position_hg19"])  # all are in FinnGen CS
        our_cs_set = set(grp[grp["our_in_either_cs"]]["position_hg19"])
        if fg_cs_set or our_cs_set:
            jaccard = len(fg_cs_set & our_cs_set) / len(fg_cs_set | our_cs_set)
        else:
            jaccard = np.nan

        summaries.append(
            {
                "phenocode": pheno,
                "gene": gene,
                "good_cs": good_cs,
                "n_matched_variants": n,
                "finngen_top_pip": fg_top["finngen_pip"] if fg_top is not None else np.nan,
                "our_top_susie_pip": our_top_susie["our_susie_pip"] if our_top_susie is not None else np.nan,
                "our_top_carma_pip": our_top_carma["our_carma_pip"] if our_top_carma is not None else np.nan,
                "top_variant_agree_susie": top_agree_susie,
                "top_variant_agree_carma": top_agree_carma,
                "pip_rho_susie": rho_susie,
                "pip_rho_carma": rho_carma,
                "pip_rho_recommended": rho_rec,
                "cs_jaccard": jaccard,
                "locus_pval": grp["locus_pval"].iloc[0],
            }
        )

    return pd.DataFrame(summaries)


# --------------------------------------------------------------------------
# 6. Plot
# --------------------------------------------------------------------------
def make_plots(matched_df, summary_df, outpath):
    """Generate comparison figure."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print("  matplotlib not available, skipping plots")
        return

    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35)

    good = matched_df[matched_df["good_cs"]].copy()
    if good.empty:
        print("  No good_cs variants to plot")
        return

    # Panel A: FinnGen PIP vs our SuSiE PIP
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.scatter(good["finngen_pip"], good["our_susie_pip"], alpha=0.5, s=20, c="#1976D2")
    ax1.plot([0, 1], [0, 1], "k--", alpha=0.3)
    rho = good["finngen_pip"].corr(good["our_susie_pip"])
    ax1.set_xlabel("FinnGen PIP (SuSiE)")
    ax1.set_ylabel("Our SuSiE PIP")
    ax1.set_title(f"SuSiE PIP (rho={rho:.3f})")
    ax1.set_xlim(-0.05, 1.05)
    ax1.set_ylim(-0.05, 1.05)

    # Panel B: FinnGen PIP vs our CARMA PIP
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(good["finngen_pip"], good["our_carma_pip"], alpha=0.5, s=20, c="#C62828")
    ax2.plot([0, 1], [0, 1], "k--", alpha=0.3)
    rho = good["finngen_pip"].corr(good["our_carma_pip"])
    ax2.set_xlabel("FinnGen PIP (SuSiE)")
    ax2.set_ylabel("Our CARMA PIP")
    ax2.set_title(f"CARMA PIP (rho={rho:.3f})")
    ax2.set_xlim(-0.05, 1.05)
    ax2.set_ylim(-0.05, 1.05)

    # Panel C: FinnGen PIP vs our recommended PIP
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.scatter(
        good["finngen_pip"], good["our_recommended_pip"], alpha=0.5, s=20, c="#2E7D32"
    )
    ax3.plot([0, 1], [0, 1], "k--", alpha=0.3)
    rho = good["finngen_pip"].corr(good["our_recommended_pip"])
    ax3.set_xlabel("FinnGen PIP (SuSiE)")
    ax3.set_ylabel("Our Recommended PIP")
    ax3.set_title(f"Recommended PIP (rho={rho:.3f})")
    ax3.set_xlim(-0.05, 1.05)
    ax3.set_ylim(-0.05, 1.05)

    # Panel D: Per-locus Jaccard CS overlap
    ax4 = fig.add_subplot(gs[1, 0])
    good_summary = summary_df[summary_df["good_cs"]].copy()
    if not good_summary.empty:
        bars = good_summary.sort_values("cs_jaccard", ascending=False)
        ax4.barh(
            range(len(bars)),
            bars["cs_jaccard"],
            color="#FF8F00",
            edgecolor="k",
            linewidth=0.5,
        )
        ax4.set_yticks(range(len(bars)))
        ax4.set_yticklabels(bars["gene"], fontsize=8)
        ax4.set_xlabel("CS Jaccard Index")
        ax4.set_title("Credible Set Overlap")
        ax4.set_xlim(0, 1.05)

    # Panel E: Top variant agreement
    ax5 = fig.add_subplot(gs[1, 1])
    if not good_summary.empty:
        agree_susie = good_summary["top_variant_agree_susie"].sum()
        agree_carma = good_summary["top_variant_agree_carma"].sum()
        total = len(good_summary)
        ax5.bar(
            ["SuSiE", "CARMA"],
            [agree_susie / total, agree_carma / total],
            color=["#1976D2", "#C62828"],
            edgecolor="k",
        )
        ax5.set_ylabel("Fraction")
        ax5.set_title(f"Top Variant Agreement\n({total} loci)")
        ax5.set_ylim(0, 1.05)
        for i, v in enumerate([agree_susie, agree_carma]):
            ax5.text(i, v / total + 0.03, f"{v}/{total}", ha="center", fontsize=10)

    # Panel F: PIP correlation per locus
    ax6 = fig.add_subplot(gs[1, 2])
    if not good_summary.empty:
        x = np.arange(len(good_summary))
        w = 0.35
        bars1 = good_summary.sort_values("gene")
        ax6.bar(x - w / 2, bars1["pip_rho_susie"], w, label="SuSiE", color="#1976D2")
        ax6.bar(x + w / 2, bars1["pip_rho_carma"], w, label="CARMA", color="#C62828")
        ax6.set_xticks(x)
        ax6.set_xticklabels(bars1["gene"], rotation=45, ha="right", fontsize=8)
        ax6.set_ylabel("PIP Correlation (rho)")
        ax6.set_title("Per-Locus PIP Correlation")
        ax6.legend(fontsize=8)
        ax6.axhline(0, color="k", linewidth=0.5)

    fig.suptitle(
        "FinnGen R12 vs Our Fine-Mapping: Credible Set Sanity Check",
        fontsize=13,
        fontweight="bold",
    )

    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"  Saved figure: {outpath}")
    plt.close()


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("FinnGen R12 Credible Set Comparison")
    print("=" * 70)

    # Step 1: Fetch FinnGen credible sets
    print("\n--- Step 1: Fetching FinnGen credible sets ---")
    all_fg = []
    for phenocode in PHENO_MAP:
        rows = fetch_finngen_credsets(phenocode)
        all_fg.extend(rows)

    if not all_fg:
        print("ERROR: No FinnGen credible sets retrieved. Check network.")
        sys.exit(1)

    fg_df = pd.DataFrame(all_fg)
    print(f"\nTotal FinnGen CS variants: {len(fg_df)}")
    print(fg_df.groupby("phenocode")[["finngen_pip"]].agg(["count", "mean"]))

    # Step 2: Build position maps (only for target hg38 positions)
    print("\n--- Step 2: Building hg38->hg19 position maps ---")
    pos_maps = {}
    for phenocode in PHENO_MAP:
        # Collect target hg38 positions for this phenotype
        pheno_fg = fg_df[fg_df["phenocode"] == phenocode]
        targets = set(
            zip(pheno_fg["chromosome"].astype(int), pheno_fg["position_hg38"].astype(int))
        )
        pos_maps[phenocode] = build_position_map(phenocode, targets)

    # Step 3: Load our finemapping
    print("\n--- Step 3: Loading our finemapping results ---")
    our_fm = load_our_finemapping()
    if our_fm.empty:
        print("ERROR: No finemapping data loaded")
        sys.exit(1)

    # Step 4: Match variants
    print("\n--- Step 4: Matching variants ---")
    all_matched = []
    for phenocode, study_name in PHENO_MAP.items():
        pheno_fg = fg_df[fg_df["phenocode"] == phenocode]
        pos_map = pos_maps.get(phenocode, {})

        if pheno_fg.empty or not pos_map:
            print(f"  Skipping {phenocode}: no CS or no position map")
            continue

        matched = compare_pips(pheno_fg, pos_map, our_fm, study_name)
        if not matched.empty:
            all_matched.append(matched)
            print(
                f"  {phenocode} -> {study_name}: {len(matched)} variants matched"
            )
        else:
            print(f"  {phenocode} -> {study_name}: 0 variants matched")

    if not all_matched:
        print("\nERROR: No variants matched across any phenotype.")
        print("This may indicate a coordinate or allele mismatch issue.")
        sys.exit(1)

    matched_df = pd.concat(all_matched, ignore_index=True)
    print(f"\nTotal matched variants: {len(matched_df)}")

    # Step 5: Compute summaries
    print("\n--- Step 5: Computing concordance metrics ---")
    summary_df = compute_summary(matched_df)

    # Print key results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)

    if not summary_df.empty:
        good = summary_df[summary_df["good_cs"]]
        print(f"\nLoci with good CS: {len(good)}")
        print(f"Overall PIP correlation (SuSiE):  {matched_df[matched_df['good_cs']]['finngen_pip'].corr(matched_df[matched_df['good_cs']]['our_susie_pip']):.3f}")
        print(f"Overall PIP correlation (CARMA):  {matched_df[matched_df['good_cs']]['finngen_pip'].corr(matched_df[matched_df['good_cs']]['our_carma_pip']):.3f}")
        print(f"Overall PIP correlation (rec):    {matched_df[matched_df['good_cs']]['finngen_pip'].corr(matched_df[matched_df['good_cs']]['our_recommended_pip']):.3f}")

        if not good.empty:
            print(f"\nMean CS Jaccard:           {good['cs_jaccard'].mean():.3f}")
            print(f"Top variant agree (SuSiE): {good['top_variant_agree_susie'].sum()}/{len(good)}")
            print(f"Top variant agree (CARMA): {good['top_variant_agree_carma'].sum()}/{len(good)}")

        print("\nPer-locus summary:")
        cols = [
            "phenocode", "gene", "n_matched_variants", "good_cs",
            "pip_rho_susie", "pip_rho_carma", "cs_jaccard",
            "top_variant_agree_susie", "top_variant_agree_carma",
        ]
        print(summary_df[cols].to_string(index=False))

    # Step 6: Save outputs
    out_matched = os.path.join(RESULTS_DIR, "finngen_credset_comparison.csv")
    out_summary = os.path.join(RESULTS_DIR, "finngen_credset_summary.csv")

    matched_df.to_csv(out_matched, index=False)
    summary_df.to_csv(out_summary, index=False)
    print(f"\nSaved: {out_matched}")
    print(f"Saved: {out_summary}")

    # Step 7: Plot
    print("\n--- Step 6: Generating figure ---")
    fig_path = os.path.join(FIG_DIR, "finngen_credset_comparison.pdf")
    make_plots(matched_df, summary_df, fig_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
