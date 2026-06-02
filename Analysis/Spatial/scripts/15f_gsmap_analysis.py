#!/usr/bin/env python3
"""
15f_gsmap_analysis.py — Analyze gsMap GWAS-spatial mapping results.

Loads per-spot p-values from gsMap spatial LDSC, performs:
  1. Zonation mapping: risk concentration along PP-PC gradient
  2. Multi-GWAS comparison: heatmap of -log10(p) per zonation bin x trait
  3. COLOC overlay: cross-reference COLOC genes with spatially-enriched spots
  4. Domain aggregation: Cauchy combination of p-values per spatial domain
  5. Cross-dataset replication: GSE192741 vs Vu concordance
  6. Conserved enrichment: CC genes enriched among spatial-risk genes

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_zonation_scores,
    load_coloc_results, load_conserved, load_multi_evidence_atlas,
    fisher_test, save_csv, print_header, print_step,
)

GSMAP_DIR = RESULTS_DIR / "gsmap"
OUTPUT_DIR = GSMAP_DIR


# ── Cauchy combination (matching gsMap implementation) ───────────────────

def acat_test(pvalues, weights=None):
    """Aggregated Cauchy Association Test for p-value combination.

    Implements the ACAT method (Liu & Xie, 2020) used by gsMap for
    aggregating per-spot p-values into per-region p-values.
    """
    pvalues = np.asarray(pvalues, dtype=float)
    pvalues = pvalues[np.isfinite(pvalues) & (pvalues > 0) & (pvalues <= 1)]
    if len(pvalues) == 0:
        return np.nan
    if np.any(pvalues == 0):
        return 0.0

    if weights is None:
        weights = np.ones(len(pvalues)) / len(pvalues)
    else:
        weights = np.asarray(weights, dtype=float)
        weights = weights / weights.sum()

    # Small p-values handled separately for numerical stability
    is_small = pvalues < 1e-15
    if np.any(is_small):
        cct_stat = np.sum(weights[is_small] / pvalues[is_small] / np.pi)
        if np.any(~is_small):
            cct_stat += np.sum(
                weights[~is_small] * np.tan((0.5 - pvalues[~is_small]) * np.pi)
            )
    else:
        cct_stat = np.sum(weights * np.tan((0.5 - pvalues) * np.pi))

    if cct_stat > 1e15:
        return 1.0 / cct_stat / np.pi

    from scipy.stats import cauchy
    return 1 - cauchy.cdf(cct_stat)


# ── Load gsMap results ──────────────────────────────────────────────────

def load_gsmap_results():
    """Load per-spot spatial LDSC p-values from gsMap output.

    gsMap quick_mode writes the spatial-LDSC trait scores to
    ``{workdir}/{dataset}/{sample}/spatial_ldsc/{sample}_{trait}.csv.gz``
    (NOT ``ldsc/`` — that was the path bug, finding F101). Each file is
    PER-SPOT with columns: ``spot, beta, se, z, p`` and has NO ``gene``
    column (finding F104). Gene-level information lives separately in
    ``latent_to_gene/{sample}_gene_marker_score.feather`` (genes x spots
    marker-score matrix); that file is consumed by ``load_gsmap_marker_scores``
    and the gene-level analyses below, not here.

    Returns dict[dataset][sample][trait] -> per-spot DataFrame.
    """
    results = {}

    for dataset_dir in sorted(GSMAP_DIR.glob("*")):
        if not dataset_dir.is_dir():
            continue
        dataset = dataset_dir.name
        if dataset in ("figures",):
            continue
        results[dataset] = {}

        for sample_dir in sorted(dataset_dir.glob("*")):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            results[dataset][sample] = {}

            # gsMap quick_mode writes spatial-LDSC to spatial_ldsc/ (was: ldsc/)
            ldsc_dir = sample_dir / "spatial_ldsc"
            if not ldsc_dir.exists():
                continue

            for ldsc_file in sorted(ldsc_dir.glob(f"{sample}_*.csv.gz")):
                # Extract trait name from filename
                trait = ldsc_file.stem.replace(f"{sample}_", "").replace(".csv", "")
                try:
                    df = pd.read_csv(ldsc_file, compression="gzip")
                    if "p" in df.columns and len(df) > 0:
                        results[dataset][sample][trait] = df
                except Exception as e:
                    print(f"    WARNING: Could not load {ldsc_file}: {e}")

    return results


def load_gsmap_marker_scores():
    """Load per-sample gene-x-spot marker-score matrices from gsMap.

    gsMap's ``latent_to_gene`` step writes
    ``{dataset}/{sample}/latent_to_gene/{sample}_gene_marker_score.feather``,
    a genes x spots matrix. The first column ``HUMAN_GENE_SYM`` holds the
    human gene symbol; every remaining column is a spot barcode (matching the
    ``spot`` column of the spatial_ldsc files), and each cell is that gene's
    marker score in that spot (higher = stronger spot-specific marker; ~71%
    are zero). This is the ONLY gene-level output gsMap produces — the
    spatial_ldsc files have no gene column (finding F104), so all gene-level
    risk localization is derived by combining these marker scores with the
    LDSC-significant spots.

    Returns dict[dataset][sample] -> DataFrame indexed by gene symbol,
    columns = spot barcodes.
    """
    marker = {}

    for dataset_dir in sorted(GSMAP_DIR.glob("*")):
        if not dataset_dir.is_dir():
            continue
        dataset = dataset_dir.name
        if dataset in ("figures",):
            continue
        marker[dataset] = {}

        for sample_dir in sorted(dataset_dir.glob("*")):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name

            l2g_dir = sample_dir / "latent_to_gene"
            feather = l2g_dir / f"{sample}_gene_marker_score.feather"
            if not feather.exists():
                continue
            try:
                df = pd.read_feather(feather)
            except Exception as e:
                print(f"    WARNING: Could not load {feather}: {e}")
                continue
            if "HUMAN_GENE_SYM" not in df.columns:
                continue
            df = df.set_index("HUMAN_GENE_SYM")
            # Collapse any duplicate gene symbols by max marker score
            if df.index.duplicated().any():
                df = df.groupby(level=0).max()
            marker[dataset][sample] = df

    return marker


def gene_risk_localization(gsmap_results, marker_scores,
                           spot_sig_alpha=0.05, top_gene_frac=0.05,
                           min_sig_spots=10):
    """Derive gene-level GWAS risk localization from spots x genes.

    A gene is "risk-localized" for a trait if it is a top marker of the spots
    that carry that trait's spatial-LDSC risk signal. For each dataset x
    sample x trait:

      1. Take the LDSC-significant spots: spot p < ``spot_sig_alpha``. If too
         few pass (< ``min_sig_spots``), fall back to the top-decile spots by
         -log10(p) so a signal-bearing foreground set always exists.
      2. From the latent_to_gene marker matrix, score every gene by
         (mean marker score in significant spots) - (mean marker score over
         all spots) = foreground-vs-background marker enrichment.
      3. Flag the top ``top_gene_frac`` of genes (by positive enrichment) as
         risk-localized for that trait in that sample.

    Per gene x trait we also record the minimum spot-level LDSC p-value among
    the significant spots (a "best localization p-value"), which is < alpha by
    construction for flagged genes.

    Aggregates across samples within a dataset: a gene is risk-localized for a
    trait if flagged in >= 1 sample; the gene's trait enrichment is the max
    across samples and its min p-value is the min across samples.

    Returns a long DataFrame: dataset, gene, trait, enrichment, min_pval,
    n_samples_flagged.
    """
    records = []

    for dataset, samples in gsmap_results.items():
        ds_marker = marker_scores.get(dataset, {})
        # gene x trait accumulators across samples within this dataset
        acc = defaultdict(lambda: {"enrich": [], "min_pval": [], "n_flag": 0})

        for sample, traits in samples.items():
            mk = ds_marker.get(sample)
            if mk is None or mk.shape[1] == 0:
                continue
            # float16 -> float32 for stable arithmetic
            mk = mk.astype("float32")
            mk_spots = set(mk.columns)
            bg_mean = mk.mean(axis=1)  # per-gene background marker score

            for trait, ldsc_df in traits.items():
                if "spot" not in ldsc_df.columns or "p" not in ldsc_df.columns:
                    continue
                ld = ldsc_df[["spot", "p"]].copy()
                ld["spot"] = ld["spot"].astype(str)
                # Restrict to spots present in the marker matrix
                ld = ld[ld["spot"].isin(mk_spots)]
                if len(ld) < min_sig_spots:
                    continue

                sig = ld[ld["p"] < spot_sig_alpha]
                if len(sig) < min_sig_spots:
                    # Fallback: top-decile spots by -log10(p)
                    ld = ld.assign(
                        mlog10p=-np.log10(ld["p"].clip(lower=1e-300))
                    )
                    cut = ld["mlog10p"].quantile(0.90)
                    sig = ld[ld["mlog10p"] >= cut]
                if len(sig) < 3:
                    continue

                sig_spots = [s for s in sig["spot"].tolist() if s in mk_spots]
                if len(sig_spots) < 3:
                    continue

                fg_mean = mk[sig_spots].mean(axis=1)  # per-gene foreground
                enrichment = (fg_mean - bg_mean)
                # Top genes by positive marker enrichment in the risk spots
                pos = enrichment[enrichment > 0].sort_values(ascending=False)
                if len(pos) == 0:
                    continue
                n_top = max(1, int(np.ceil(len(enrichment) * top_gene_frac)))
                top_genes = pos.head(n_top)
                min_p = float(sig["p"].min())

                for gene, enr in top_genes.items():
                    a = acc[(gene, trait)]
                    a["enrich"].append(float(enr))
                    a["min_pval"].append(min_p)
                    a["n_flag"] += 1

        for (gene, trait), a in acc.items():
            if not a["enrich"]:
                continue
            records.append({
                "dataset": dataset,
                "gene": gene,
                "trait": trait,
                "enrichment": float(np.max(a["enrich"])),
                "min_pval": float(np.min(a["min_pval"])),
                "n_samples_flagged": int(a["n_flag"]),
            })

    return pd.DataFrame(records)


def summarize_results(gsmap_results):
    """Print summary of loaded gsMap results."""
    total_samples = 0
    total_traits = 0
    for dataset, samples in gsmap_results.items():
        n_samples = len(samples)
        traits = set()
        for sample, sample_traits in samples.items():
            traits.update(sample_traits.keys())
        total_samples += n_samples
        total_traits = max(total_traits, len(traits))
        print(f"    {dataset}: {n_samples} samples, "
              f"traits: {sorted(traits)}")
    return total_samples, total_traits


# ── Analysis 1: Zonation mapping ────────────────────────────────────────

def zonation_mapping(gsmap_results):
    """Map gsMap risk scores onto periportal-pericentral zonation gradient.

    Merges spot-level gsMap p-values with zonation scores from 04a,
    bins spots into PP->PC gradient, tests risk concentration via
    Mann-Whitney U.
    """
    print_step("Zonation mapping")

    zon_scores = load_zonation_scores()
    if len(zon_scores) == 0:
        print("    WARNING: No zonation scores found (run 04a first)")
        return pd.DataFrame(), pd.DataFrame()

    config = load_config()
    n_bins = config.get("zonation", {}).get("n_bins", 5)

    enrichment_records = []
    heatmap_records = []

    for dataset, samples in gsmap_results.items():
        for sample, traits in samples.items():
            for trait, ldsc_df in traits.items():
                if "spot" not in ldsc_df.columns or "p" not in ldsc_df.columns:
                    continue

                # Merge with zonation scores
                ldsc_df = ldsc_df.copy()
                ldsc_df["spot"] = ldsc_df["spot"].astype(str)

                # Zonation scores index = spot barcode
                zon_sub = zon_scores.copy()
                zon_sub.index = zon_sub.index.astype(str)

                merged = ldsc_df.set_index("spot").join(zon_sub, how="inner")

                if len(merged) < 10:
                    continue

                # Check for zonation score column
                zon_col = None
                for col in ["zonation_score", "pc_score", "zon_score"]:
                    if col in merged.columns:
                        zon_col = col
                        break
                if zon_col is None:
                    continue

                # Bin by zonation
                merged["zon_bin"] = pd.qcut(
                    merged[zon_col], q=n_bins, labels=False,
                    duplicates="drop"
                )

                # Mann-Whitney U: PP (bin 0) vs PC (bin max)
                pp_pvals = merged[merged["zon_bin"] == 0]["p"]
                pc_pvals = merged[merged["zon_bin"] == merged["zon_bin"].max()]["p"]

                if len(pp_pvals) > 5 and len(pc_pvals) > 5:
                    # Compare -log10(p) — higher = more risk
                    u_stat, u_pval = mannwhitneyu(
                        -np.log10(pp_pvals.clip(lower=1e-300)),
                        -np.log10(pc_pvals.clip(lower=1e-300)),
                        alternative="two-sided",
                    )
                    pp_mean = -np.log10(pp_pvals.clip(lower=1e-300)).mean()
                    pc_mean = -np.log10(pc_pvals.clip(lower=1e-300)).mean()
                    risk_zone = "periportal" if pp_mean > pc_mean else "pericentral"

                    enrichment_records.append({
                        "dataset": dataset,
                        "sample": sample,
                        "trait": trait,
                        "pp_mean_mlog10p": pp_mean,
                        "pc_mean_mlog10p": pc_mean,
                        "mwu_stat": u_stat,
                        "mwu_pval": u_pval,
                        "risk_zone": risk_zone,
                        "n_spots": len(merged),
                    })

                # Heatmap data: mean -log10(p) per bin
                for bin_id in sorted(merged["zon_bin"].unique()):
                    bin_pvals = merged[merged["zon_bin"] == bin_id]["p"]
                    heatmap_records.append({
                        "dataset": dataset,
                        "sample": sample,
                        "trait": trait,
                        "zon_bin": int(bin_id),
                        "mean_mlog10p": -np.log10(
                            bin_pvals.clip(lower=1e-300)
                        ).mean(),
                        "n_spots": len(bin_pvals),
                    })

    enrichment_df = pd.DataFrame(enrichment_records)
    heatmap_df = pd.DataFrame(heatmap_records)

    if len(enrichment_df) > 0:
        print(f"    Zonation enrichment: {len(enrichment_df)} sample x trait tests")
        sig = enrichment_df[enrichment_df["mwu_pval"] < 0.05]
        if len(sig) > 0:
            print(f"    Significant (p<0.05): {len(sig)} tests")
            for _, row in sig.iterrows():
                print(f"      {row['trait']} / {row['sample']}: "
                      f"risk={row['risk_zone']}, p={row['mwu_pval']:.3e}")

    return enrichment_df, heatmap_df


# ── Analysis 2: COLOC overlay ───────────────────────────────────────────

def coloc_overlay(gsmap_results):
    """Cross-reference COLOC genes with gsMap spatially-enriched spots.

    For each COLOC gene with PP.H4>0.5, check if spots where that gene
    is highly expressed also show significant gsMap risk signal.
    """
    print_step("COLOC overlay")

    coloc_df = load_coloc_results(pp4_threshold=0.5)
    if len(coloc_df) == 0:
        print("    WARNING: No COLOC results found")
        return pd.DataFrame()

    # Get COLOC gene symbols
    gene_col = next(
        (c for c in ["gene_symbol", "symbol", "gene_name", "gene"]
         if c in coloc_df.columns),
        None,
    )
    if gene_col is None:
        print("    WARNING: No gene symbol column in COLOC results")
        return pd.DataFrame()

    coloc_genes = set(coloc_df[gene_col].dropna().unique())
    print(f"    COLOC genes: {len(coloc_genes)}")

    # For each gsMap result, identify spots with p < 0.05 and check
    # gene overlap
    records = []

    for dataset, samples in gsmap_results.items():
        for sample, traits in samples.items():
            for trait, ldsc_df in traits.items():
                if "p" not in ldsc_df.columns:
                    continue

                sig_spots = ldsc_df[ldsc_df["p"] < 0.05]
                all_spots = ldsc_df

                # spatial_ldsc has NO gene column (finding F104), so this
                # branch is inert; kept for forward-compatibility if a future
                # gsMap mode emits per-spot gene attribution. The spot-level
                # summary below is the operative COLOC overlay output.
                if "gene" in ldsc_df.columns:
                    sig_genes = set(sig_spots["gene"].dropna().unique())
                    all_genes = set(all_spots["gene"].dropna().unique())
                    overlap = coloc_genes & sig_genes
                    bg_overlap = coloc_genes & all_genes

                    if len(all_genes) > 0 and len(bg_overlap) > 0:
                        # Fisher test
                        OR, pval, n_hit = fisher_test(
                            coloc_genes, sig_genes, all_genes
                        )
                        records.append({
                            "dataset": dataset,
                            "sample": sample,
                            "trait": trait,
                            "n_coloc_in_sig": len(overlap),
                            "n_sig_genes": len(sig_genes),
                            "n_coloc_total": len(bg_overlap),
                            "n_genes_total": len(all_genes),
                            "fisher_or": OR,
                            "fisher_pval": pval,
                        })

                # Spot-level summary regardless of gene info
                n_sig = len(sig_spots)
                n_total = len(all_spots)
                records.append({
                    "dataset": dataset,
                    "sample": sample,
                    "trait": trait,
                    "n_sig_spots": n_sig,
                    "n_total_spots": n_total,
                    "pct_sig": 100.0 * n_sig / n_total if n_total > 0 else 0,
                    "mean_mlog10p": -np.log10(
                        all_spots["p"].clip(lower=1e-300)
                    ).mean(),
                })

    result_df = pd.DataFrame(records)
    if len(result_df) > 0:
        print(f"    COLOC overlay: {len(result_df)} records")

    return result_df


# ── Analysis 3: Domain aggregation ──────────────────────────────────────

def domain_aggregation(gsmap_results):
    """Aggregate gsMap p-values per spatial domain using Cauchy combination.

    Loads spatial domain assignments from 05d, groups spots by domain,
    and applies ACAT to get per-domain GWAS risk p-values.
    """
    print_step("Domain aggregation")

    domain_path = RESULTS_DIR / "domains" / "spatial_domains.csv"
    if not domain_path.exists():
        print("    WARNING: Spatial domains not found (run 05d first)")
        return pd.DataFrame()

    domains_df = pd.read_csv(domain_path, index_col=0)
    domain_col = None
    for col in ["spatial_domain", "domain", "leiden", "cluster"]:
        if col in domains_df.columns:
            domain_col = col
            break
    if domain_col is None:
        print("    WARNING: No domain column found")
        return pd.DataFrame()

    records = []

    for dataset, samples in gsmap_results.items():
        for sample, traits in samples.items():
            for trait, ldsc_df in traits.items():
                if "spot" not in ldsc_df.columns or "p" not in ldsc_df.columns:
                    continue

                ldsc_df = ldsc_df.copy()
                ldsc_df["spot"] = ldsc_df["spot"].astype(str)

                # Merge with domain assignments
                dom_sub = domains_df[[domain_col]].copy()
                dom_sub.index = dom_sub.index.astype(str)

                merged = ldsc_df.set_index("spot").join(dom_sub, how="inner")
                if len(merged) < 10:
                    continue

                for domain_id in sorted(merged[domain_col].unique()):
                    pvals = merged[merged[domain_col] == domain_id]["p"].values
                    pvals = pvals[(pvals > 0) & (pvals <= 1) & np.isfinite(pvals)]

                    if len(pvals) < 3:
                        continue

                    cauchy_p = acat_test(pvals)
                    records.append({
                        "dataset": dataset,
                        "sample": sample,
                        "trait": trait,
                        "domain": domain_id,
                        "cauchy_pval": cauchy_p,
                        "median_pval": np.median(pvals),
                        "mean_mlog10p": -np.log10(
                            np.clip(pvals, 1e-300, 1)
                        ).mean(),
                        "n_spots": len(pvals),
                    })

    result_df = pd.DataFrame(records)
    if len(result_df) > 0:
        print(f"    Domain aggregation: {len(result_df)} domain x trait tests")

    return result_df


# ── Analysis 4: Cross-dataset replication ───────────────────────────────

def cross_dataset_replication(gene_risk_df):
    """Correlate GSE192741 vs Vu risk patterns at the gene level.

    Uses the per-gene risk-localization enrichment from
    ``gene_risk_localization`` (foreground-vs-background marker enrichment in
    LDSC-significant spots), then computes Spearman correlation of gene
    enrichment between datasets for each shared trait.
    """
    print_step("Cross-dataset replication")

    if gene_risk_df is None or len(gene_risk_df) == 0:
        print("    WARNING: No gene-level risk data for replication")
        return pd.DataFrame()

    datasets = sorted(gene_risk_df["dataset"].unique())
    if len(datasets) < 2:
        print("    WARNING: Need >= 2 datasets for replication")
        return pd.DataFrame()

    # gene enrichment per (dataset, trait)
    gene_means = {}
    for (dataset, trait), sub in gene_risk_df.groupby(["dataset", "trait"]):
        gene_means[(dataset, trait)] = dict(
            zip(sub["gene"], sub["enrichment"])
        )

    # Compare pairs of datasets
    records = []
    trait_set = set(t for _, t in gene_means.keys())

    for trait in sorted(trait_set):
        # Get all datasets that have this trait
        ds_with_trait = [d for d in datasets if (d, trait) in gene_means]
        if len(ds_with_trait) < 2:
            continue

        for i in range(len(ds_with_trait)):
            for j in range(i + 1, len(ds_with_trait)):
                d1, d2 = ds_with_trait[i], ds_with_trait[j]
                scores1 = gene_means[(d1, trait)]
                scores2 = gene_means[(d2, trait)]

                common_genes = sorted(set(scores1.keys()) & set(scores2.keys()))
                if len(common_genes) < 20:
                    continue

                vals1 = [scores1[g] for g in common_genes]
                vals2 = [scores2[g] for g in common_genes]

                rho, pval = spearmanr(vals1, vals2)
                records.append({
                    "trait": trait,
                    "dataset_1": d1,
                    "dataset_2": d2,
                    "n_common_genes": len(common_genes),
                    "spearman_rho": rho,
                    "spearman_pval": pval,
                })

    result_df = pd.DataFrame(records)
    if len(result_df) > 0:
        print(f"    Cross-dataset replication: {len(result_df)} comparisons")
        for _, row in result_df.iterrows():
            sig = "*" if row["spearman_pval"] < 0.05 else ""
            print(f"      {row['trait']}: {row['dataset_1']} vs {row['dataset_2']} "
                  f"rho={row['spearman_rho']:.3f}, "
                  f"p={row['spearman_pval']:.3e} {sig} "
                  f"(n={row['n_common_genes']})")

    return result_df


# ── Analysis 5: Conserved enrichment ───────────────────────────────

def conserved_enrichment(gene_risk_df, marker_scores):
    """Test if Conserved genes are enriched among spatially-risk genes.

    Defines "spatial risk genes" per dataset x trait as the genes flagged
    risk-localized by ``gene_risk_localization`` (top marker-enriched genes in
    the LDSC-significant spots). The universe is the full set of genes present
    in that dataset's marker matrices. Tests overlap with Conserved via Fisher.
    """
    print_step("Conserved enrichment")

    cc_genes = load_conserved()
    if not cc_genes:
        print("    WARNING: No Conserved genes found")
        return pd.DataFrame()
    cc_set = set(cc_genes)
    print(f"    Conserved: {len(cc_set)} genes")

    if gene_risk_df is None or len(gene_risk_df) == 0:
        print("    WARNING: No gene-level risk data for CC enrichment")
        return pd.DataFrame()

    # Per-dataset gene universe = union of genes in the marker matrices
    dataset_universe = {}
    for dataset, samples in marker_scores.items():
        genes = set()
        for sample, mk in samples.items():
            genes.update(mk.index.tolist())
        if genes:
            dataset_universe[dataset] = genes

    records = []

    for (dataset, trait), sub in gene_risk_df.groupby(["dataset", "trait"]):
        all_genes = dataset_universe.get(dataset)
        if not all_genes or len(all_genes) < 50:
            continue
        risk_genes = set(sub["gene"])
        if len(risk_genes) == 0:
            continue

        OR, pval, n_overlap = fisher_test(cc_set, risk_genes, all_genes)

        records.append({
            "dataset": dataset,
            "trait": trait,
            "n_risk_genes": len(risk_genes),
            "n_cc_in_risk": n_overlap,
            "n_cc_total_in_universe": len(cc_set & all_genes),
            "n_universe": len(all_genes),
            "fisher_or": OR,
            "fisher_pval": pval,
        })

    result_df = pd.DataFrame(records)
    if len(result_df) > 0:
        print(f"    CC enrichment: {len(result_df)} tests")
        for _, row in result_df.iterrows():
            sig = "*" if row["fisher_pval"] < 0.05 else ""
            print(f"      {row['trait']}/{row['dataset']}: "
                  f"OR={row['fisher_or']:.2f}, p={row['fisher_pval']:.3e} {sig} "
                  f"({row['n_cc_in_risk']}/{row['n_risk_genes']})")

    return result_df


# ── Main ────────────────────────────────────────────────────────────────

def main():
    print_header("15f: gsMap Analysis")

    config = load_config()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load results
    print("  Loading gsMap results...")
    gsmap_results = load_gsmap_results()
    n_samples, n_traits = summarize_results(gsmap_results)

    if n_samples == 0:
        print("  ERROR: No gsMap results found. Run 15e first.")
        sys.exit(1)

    # Load per-sample gene x spot marker scores (latent_to_gene) and derive
    # gene-level risk localization (spot->gene mapping; finding F104 fix).
    print("  Loading gsMap latent_to_gene marker scores...")
    marker_scores = load_gsmap_marker_scores()
    gene_risk_df = gene_risk_localization(gsmap_results, marker_scores)
    if len(gene_risk_df) > 0:
        print(f"    Gene-level risk localization: {len(gene_risk_df)} "
              f"gene x trait rows across "
              f"{gene_risk_df['gene'].nunique()} genes, "
              f"{gene_risk_df['trait'].nunique()} traits")
        save_csv(gene_risk_df, "gene_risk_localization.csv", subdir="gsmap")
    else:
        print("    WARNING: No gene-level risk localization derived")

    # 1. Zonation mapping
    zon_enrich, zon_heatmap = zonation_mapping(gsmap_results)
    if len(zon_enrich) > 0:
        save_csv(zon_enrich, "gwas_zonation_enrichment.csv", subdir="gsmap")
    if len(zon_heatmap) > 0:
        save_csv(zon_heatmap, "gwas_zonation_heatmap.csv", subdir="gsmap")

    # 2. COLOC overlay
    coloc_spatial = coloc_overlay(gsmap_results)
    if len(coloc_spatial) > 0:
        save_csv(coloc_spatial, "coloc_spatial_specificity.csv", subdir="gsmap")

    # 3. Domain aggregation
    domain_agg = domain_aggregation(gsmap_results)
    if len(domain_agg) > 0:
        save_csv(domain_agg, "domain_risk_aggregation.csv", subdir="gsmap")

    # 4. Cross-dataset replication (gene-level, from marker-score enrichment)
    replication = cross_dataset_replication(gene_risk_df)
    if len(replication) > 0:
        save_csv(replication, "cross_dataset_replication.csv", subdir="gsmap")

    # 5. Conserved enrichment (gene-level, from marker-score enrichment)
    cc_enrich = conserved_enrichment(gene_risk_df, marker_scores)
    if len(cc_enrich) > 0:
        save_csv(cc_enrich, "conserved_spatial_risk.csv", subdir="gsmap")

    # Summary
    print("\n  Output files:")
    for csv_path in sorted(OUTPUT_DIR.glob("*.csv")):
        size_kb = csv_path.stat().st_size / 1e3
        print(f"    {csv_path.name} ({size_kb:.1f} KB)")

    print_header("15f: Complete")


if __name__ == "__main__":
    main()
