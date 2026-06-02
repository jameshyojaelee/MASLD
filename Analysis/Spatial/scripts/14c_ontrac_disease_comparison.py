#!/usr/bin/env python3
"""
14c_ontrac_disease_comparison.py — Disease comparison along the DPT
diffusion-pseudotime niche trajectory (ONTraC fallback; ONTraC CLI unavailable).

PROVENANCE NOTE (F090/F091/F093/F094): the "niche trajectory" (NT) score is a DPT
diffusion-pseudotime rooted at the highest-hepatocyte spot (computed in 14b),
NOT a native ONTraC trajectory (the ONTraC CLI is unavailable, so 14b ran its DPT
fallback; see parameter_sweep.csv method=diffusion_pseudotime). Comments/columns that previously asserted
"ONTraC" are labeled as DPT diffusion-pseudotime for honest provenance. The NT
axis is therefore largely a hepatocyte-content gradient in multi-cellular Visium
spots (see F093 caveat below), not a validated disease-progression program.

Compares healthy vs MASLD NT (DPT) score distributions at the DONOR level
(F091: spot-level KS over thousands of autocorrelated spots was
pseudoreplication on an effective n=~5 donors), computes cell-type composition
along trajectory bins, and identifies trajectory-associated genes (TAGs) via
Spearman correlation on a log1p-CPM layer (F094: raw counts were depth-confounded)
with BH FDR correction.

Cross-references TAGs with dream DEGs, Conserved, and drug targets.

SLURM: --partition=cpu --cpus=16 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import ks_2samp, spearmanr
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_deconvolved_adata,
    load_dream_degs, load_conserved, load_drug_targets,
    save_csv, print_header, print_step,
)
from spatial_stats import ensure_lognorm

# Donor column for the GSE192741 spatial cohort (do NOT trust hardcoded
# inverted condition maps elsewhere). Healthy=JBO018/JBO022,
# Steatotic=JBO014/015/019.
DONOR_COL = "sample_id"


def load_nt_scores():
    """Load niche trajectory scores from 14b."""
    path = RESULTS_DIR / "ontrac" / "niche_trajectory_scores.csv"
    if not path.exists():
        print("  ERROR: Run 14b_ontrac_run.py first")
        sys.exit(1)
    df = pd.read_csv(path, index_col=0)
    print(f"  Loaded NT scores: {len(df)} spots")
    return df


def disease_trajectory_comparison(adata, nt_df):
    """Compare NT (DPT) score between conditions at the DONOR level (F091).

    Spot-level KS over thousands of spatially-autocorrelated, donor-nested spots
    is pseudoreplication (the old ks_pval=4.7e-63 was an n=spots artifact; the
    true n is the ~5 donors). Here we aggregate NT to one mean per donor, then
    run a donor-level Mann-Whitney across the per-donor means. We KEEP the column
    names (``ks_statistic``, ``ks_pval``) for downstream/back-compat but they now
    hold a DONOR-LEVEL Mann-Whitney U statistic and p-value, not a spot-level KS.
    A descriptive spot-level KS is also retained under ``spot_ks_*`` for context.
    """
    from scipy.stats import mannwhitneyu

    # Merge NT scores into adata.obs
    shared = adata.obs_names.intersection(nt_df.index)
    adata_sub = adata[shared].copy()
    adata_sub.obs["nt_score"] = nt_df.loc[shared, "nt_score"].values
    adata_sub.obs["niche_cluster"] = nt_df.loc[shared, "niche_cluster"].values

    conditions = adata_sub.obs["condition"].unique().tolist()
    cond_healthy = [c for c in conditions if "healthy" in c.lower() or c == "Healthy"]
    cond_disease = [c for c in conditions if c not in cond_healthy]

    if not cond_healthy or not cond_disease:
        print("  WARNING: Need both healthy and disease conditions for comparison")
        return pd.DataFrame(), adata_sub

    if DONOR_COL not in adata_sub.obs.columns:
        print(f"  WARNING: donor column '{DONOR_COL}' missing — cannot run "
              f"donor-level test; falling back to descriptive spot summary only.")

    # Per-donor mean NT score (the unit of evidence).
    donor_means = (adata_sub.obs.dropna(subset=["nt_score"])
                   .groupby([DONOR_COL, "condition"], observed=True)["nt_score"]
                   .mean().reset_index()) if DONOR_COL in adata_sub.obs.columns \
        else pd.DataFrame()

    results = []
    for disease in cond_disease:
        nt_healthy = adata_sub.obs.loc[
            adata_sub.obs["condition"].isin(cond_healthy), "nt_score"
        ].dropna()
        nt_disease = adata_sub.obs.loc[
            adata_sub.obs["condition"] == disease, "nt_score"
        ].dropna()

        # Descriptive spot-level KS (kept for context, NOT the headline).
        spot_ks_stat, spot_ks_pval = ks_2samp(nt_healthy, nt_disease)

        # Donor-level Mann-Whitney on per-donor means (the honest test).
        donor_stat, donor_pval, n_donor_h, n_donor_d = np.nan, np.nan, 0, 0
        if len(donor_means):
            mh = donor_means.loc[
                donor_means["condition"].isin(cond_healthy), "nt_score"].values
            md = donor_means.loc[
                donor_means["condition"] == disease, "nt_score"].values
            n_donor_h, n_donor_d = len(mh), len(md)
            if n_donor_h >= 1 and n_donor_d >= 1:
                try:
                    donor_stat, donor_pval = mannwhitneyu(
                        mh, md, alternative="two-sided")
                except ValueError:
                    donor_stat, donor_pval = np.nan, np.nan

        results.append({
            "comparison": f"Healthy_vs_{disease}",
            "n_donors_healthy": int(n_donor_h),
            "n_donors_disease": int(n_donor_d),
            "n_spots_healthy": len(nt_healthy),
            "n_spots_disease": len(nt_disease),
            "healthy_mean_nt": nt_healthy.mean(),
            "disease_mean_nt": nt_disease.mean(),
            "healthy_median_nt": nt_healthy.median(),
            "disease_median_nt": nt_disease.median(),
            # NOTE (F091): ks_* columns now hold a DONOR-LEVEL Mann-Whitney U,
            # not a spot-level KS. Spot KS retained separately for context.
            "ks_statistic": donor_stat,
            "ks_pval": donor_pval,
            "test": "donor_level_mannwhitney",
            "spot_ks_statistic": spot_ks_stat,
            "spot_ks_pval": spot_ks_pval,
        })
        print(f"    Healthy vs {disease}: donor-level MWU U={donor_stat}, "
              f"p={donor_pval:.3g} (n={n_donor_h} vs {n_donor_d} donors)")
        print(f"      [descriptive spot KS={spot_ks_stat:.3f}, "
              f"p={spot_ks_pval:.2e} — pseudoreplicated, NOT the headline]")
        print(f"      Healthy NT: {nt_healthy.mean():.3f} +/- {nt_healthy.std():.3f}")
        print(f"      Disease NT: {nt_disease.mean():.3f} +/- {nt_disease.std():.3f}")

    return pd.DataFrame(results), adata_sub


def trajectory_composition(adata_sub, n_bins=10):
    """Compute cell-type composition along NT score quantiles."""
    if "nt_score" not in adata_sub.obs.columns:
        return pd.DataFrame()

    nt = adata_sub.obs["nt_score"].dropna()
    if len(nt) == 0:
        return pd.DataFrame()

    adata_sub.obs["nt_bin"] = pd.qcut(
        adata_sub.obs["nt_score"], q=n_bins, labels=False, duplicates="drop"
    )

    # Get dominant cell type (from metadata or compute)
    ct_col = None
    for col in ["Cell_Type", "cell_type_dominant", "dominant_cell_type"]:
        if col in adata_sub.obs.columns:
            ct_col = col
            break

    if ct_col is None:
        print("  WARNING: No cell type column found for composition analysis")
        return pd.DataFrame()

    comp = pd.crosstab(
        adata_sub.obs["nt_bin"],
        adata_sub.obs[ct_col],
        normalize="index",
    )

    # Add bin statistics
    bin_stats = adata_sub.obs.groupby("nt_bin")["nt_score"].agg(["mean", "count"])
    comp = comp.join(bin_stats)

    return comp


def compute_trajectory_associated_genes(adata_sub, config):
    """Compute Spearman correlation of each gene with NT (DPT) score.

    F094: the deconvolved h5ad .X is RAW integer counts (depth-confounded), and
    NT is composition-derived so it tracks per-spot depth — correlating raw
    counts manufactures spurious associations. We correlate against a log1p-CPM
    layer instead (built once via ensure_lognorm). Spearman is rank-based, but
    normalization still removes the depth-driven rank inflation.

    Returns DataFrame with gene, rho, pval, padj, direction.
    """
    tag_padj = config["tag_padj_threshold"]
    tag_rho = config["tag_rho_threshold"]

    nt_scores = adata_sub.obs["nt_score"].values
    valid_mask = ~np.isnan(nt_scores)

    if valid_mask.sum() < 50:
        print("  WARNING: Too few spots with valid NT scores for TAG analysis")
        return pd.DataFrame()

    nt_valid = nt_scores[valid_mask]
    # F094: read a normalized (log1p-CPM) layer, not raw .X.
    lognorm_layer = ensure_lognorm(adata_sub)
    X = adata_sub.layers[lognorm_layer][valid_mask]

    n_genes = adata_sub.n_vars
    results = []

    print(f"  Computing Spearman correlations for {n_genes} genes...")
    for i in range(n_genes):
        if issparse(X):
            expr = np.asarray(X[:, i].todense()).flatten()
        else:
            expr = X[:, i].flatten()

        # Skip genes with zero variance
        if expr.std() == 0:
            continue

        rho, pval = spearmanr(nt_valid, expr)
        results.append({
            "gene": adata_sub.var_names[i],
            "spearman_rho": rho,
            "pval": pval,
        })

        if (i + 1) % 5000 == 0:
            print(f"    Processed {i + 1}/{n_genes} genes...")

    tag_df = pd.DataFrame(results)
    if len(tag_df) == 0:
        return tag_df

    # BH FDR correction
    pvals = tag_df["pval"].fillna(1.0).values
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    tag_df["padj"] = padj

    # Filter significant TAGs
    tag_df["is_tag"] = (tag_df["padj"] < tag_padj) & (tag_df["spearman_rho"].abs() > tag_rho)
    tag_df["direction"] = np.where(tag_df["spearman_rho"] > 0, "positive", "negative")

    # Sort by absolute rho
    tag_df = tag_df.sort_values("spearman_rho", ascending=False, key=abs)

    n_tag = tag_df["is_tag"].sum()
    n_pos = ((tag_df["is_tag"]) & (tag_df["direction"] == "positive")).sum()
    n_neg = ((tag_df["is_tag"]) & (tag_df["direction"] == "negative")).sum()
    print(f"  TAGs: {n_tag} (padj<{tag_padj}, |rho|>{tag_rho})")
    print(f"    Positive (increase with NT): {n_pos}")
    print(f"    Negative (decrease with NT): {n_neg}")
    if n_tag > 0 and n_neg / max(n_tag, 1) > 0.9:
        # F093 caveat: the NT (DPT) axis is rooted at the highest-hepatocyte spot,
        # so a near-total negative skew means TAGs are largely hepatocyte-identity
        # genes diluted by non-parenchymal infiltration (a spot-composition
        # gradient), NOT a within-cell-type disease program. Recommended
        # follow-up: correlate WITHIN hepatocyte-dominant spots or regress out
        # hepatocyte proportion before correlating.
        print(f"    CAVEAT (F093): {n_neg}/{n_tag} TAGs are NEGATIVE — this axis "
              f"is a hepatocyte-content/composition gradient, not a within-state "
              f"disease program. Interpret TAGs as composition-confounded.")

    return tag_df


def cross_reference_tags(tag_df):
    """Cross-reference TAGs with dream DEGs, Conserved, drug targets."""
    tags = set(tag_df[tag_df["is_tag"]]["gene"])
    if len(tags) == 0:
        print("  No TAGs to cross-reference")
        return pd.DataFrame()

    records = []

    # Dream DEGs
    try:
        dream = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
        gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
        deg_set = set(dream[gene_col].dropna())
        overlap_deg = tags & deg_set
        records.append({
            "reference": "Dream_DEGs",
            "n_tags": len(tags),
            "n_reference": len(deg_set),
            "n_overlap": len(overlap_deg),
            "pct_overlap": len(overlap_deg) / max(len(tags), 1) * 100,
        })
        print(f"  TAGs in dream DEGs: {len(overlap_deg)}/{len(tags)} "
              f"({len(overlap_deg) / max(len(tags), 1) * 100:.1f}%)")
    except Exception as e:
        print(f"  WARNING: Could not load dream DEGs: {e}")

    # Conserved
    cc = load_conserved()
    if cc:
        cc_set = set(cc)
        overlap_cc = tags & cc_set
        records.append({
            "reference": "Conserved",
            "n_tags": len(tags),
            "n_reference": len(cc_set),
            "n_overlap": len(overlap_cc),
            "pct_overlap": len(overlap_cc) / max(len(tags), 1) * 100,
        })
        print(f"  TAGs in Conserved: {len(overlap_cc)}/{len(tags)}")

    # Drug targets
    drug_df = load_drug_targets()
    if len(drug_df) > 0:
        dt_col = next(
            (c for c in drug_df.columns if "gene" in c.lower() or "target" in c.lower()),
            None,
        )
        if dt_col:
            dt_set = set(drug_df[dt_col].dropna())
            overlap_dt = tags & dt_set
            records.append({
                "reference": "Drug_targets",
                "n_tags": len(tags),
                "n_reference": len(dt_set),
                "n_overlap": len(overlap_dt),
                "pct_overlap": len(overlap_dt) / max(len(tags), 1) * 100,
            })
            print(f"  TAGs in drug targets: {len(overlap_dt)}/{len(tags)}")

    return pd.DataFrame(records)


def identify_novel_candidates(tag_df):
    """Identify TAGs that are NOT dream DEGs — potential novel spatial candidates."""
    tags = tag_df[tag_df["is_tag"]].copy()
    if len(tags) == 0:
        return pd.DataFrame()

    try:
        dream = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
        gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
        deg_set = set(dream[gene_col].dropna())
    except Exception:
        deg_set = set()

    tags["is_dream_deg"] = tags["gene"].isin(deg_set)
    novel = tags[~tags["is_dream_deg"]].copy()

    print(f"\n  Novel trajectory candidates (TAG but NOT dream DEG): {len(novel)}")
    if len(novel) > 0:
        top = novel.head(20)
        for _, row in top.iterrows():
            print(f"    {row['gene']}: rho={row['spearman_rho']:.3f}, padj={row['padj']:.2e}")

    return novel


def main():
    print_header("14c: Disease Comparison Along DPT Niche Trajectory "
                 "(ONTraC fallback; CLI unavailable)")

    config = load_config()
    ontrac_config = config["ontrac"]
    output_dir = RESULTS_DIR / "ontrac"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load NT scores
    print_step("Loading niche trajectory scores")
    nt_df = load_nt_scores()

    # Load spatial adata for expression data
    print_step("Loading deconvolved spatial data")
    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

    # Merge cell type info from ontrac metadata
    meta_path = RESULTS_DIR / "ontrac" / "input" / "ontrac_metadata.csv"
    if meta_path.exists():
        ontrac_meta = pd.read_csv(meta_path, index_col=0)
        shared = adata.obs_names.intersection(ontrac_meta.index)
        adata.obs.loc[shared, "Cell_Type"] = ontrac_meta.loc[shared, "Cell_Type"]

    # Disease trajectory comparison
    print_step("Comparing healthy vs disease trajectory distributions")
    disease_df, adata_sub = disease_trajectory_comparison(adata, nt_df)
    if len(disease_df) > 0:
        save_csv(disease_df, "disease_trajectory_comparison.csv", subdir="ontrac")

    # Cell-type composition along trajectory
    print_step("Computing cell-type composition along trajectory")
    comp_df = trajectory_composition(adata_sub)
    if len(comp_df) > 0:
        save_csv(comp_df, "trajectory_composition.csv", subdir="ontrac")
        print(f"  Composition bins: {len(comp_df)}")

    # Trajectory-associated genes
    print_step("Computing trajectory-associated genes (TAGs)")
    tag_df = compute_trajectory_associated_genes(adata_sub, ontrac_config)
    if len(tag_df) > 0:
        save_csv(tag_df, "trajectory_associated_genes.csv", subdir="ontrac")

        # Cross-reference
        print_step("Cross-referencing TAGs with bulk evidence")
        xref_df = cross_reference_tags(tag_df)
        if len(xref_df) > 0:
            save_csv(xref_df, "tag_cross_reference.csv", subdir="ontrac")

        # Novel candidates
        print_step("Identifying novel trajectory candidates")
        novel_df = identify_novel_candidates(tag_df)
        if len(novel_df) > 0:
            save_csv(novel_df, "novel_trajectory_candidates.csv", subdir="ontrac")

    print_header("14c: Complete (DPT diffusion-pseudotime; ONTraC fallback)")


if __name__ == "__main__":
    main()
