#!/usr/bin/env python3
"""
13c_multisp_communication.py — Multi-omic communication analysis with MultiSP domains.

Runs squidpy L-R analysis using MultiSP bimodal domains as the cluster key,
projects SCENIC+ regulon activity onto spots via cell-type-weighted
deconvolution, and builds multi-omic communication cards combining
L-R, regulon, and chromatin evidence.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=6:00:00
"""

import pathlib
import sys
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, C2L_PREFIX,
    load_config, load_deconvolved_adata,
    save_csv, print_header, print_step,
    harmonize_cell_types, strip_c2l_prefix,
)


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_multisp_adata():
    """Load spatial AnnData with MultiSP domains from 13b."""
    path = RESULTS_DIR / "multisp" / "spatial_with_multisp_domains.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Run 13b first: {path}")
    print_step(f"Loading MultiSP adata: {path}")
    return sc.read_h5ad(path)


def load_scenic_regulons():
    """Load SCENIC+ regulon data if available.

    Returns a dict: {regulon_name: {target_genes: [...], tf: str, ...}}
    and a DataFrame of regulon-cell-type activity if available.
    """
    scenic_dir = PROJECT_ROOT / "Analysis/ATAC/Human_Multiome/results/scenic_plus"
    regulon_data = {}
    regulon_activity = pd.DataFrame()

    if not scenic_dir.exists():
        print_step("WARNING: SCENIC+ results directory not found")
        return regulon_data, regulon_activity

    # Try loading regulon definitions
    for candidate in ["regulons.csv", "regulon_definitions.csv", "grn_output.csv"]:
        path = scenic_dir / candidate
        if path.exists():
            df = pd.read_csv(path)
            print_step(f"Loaded regulons from {candidate}: {len(df)} rows")
            # Parse into dict structure
            if "regulon" in df.columns and "gene" in df.columns:
                for reg_name, group in df.groupby("regulon"):
                    tf = reg_name.split("_")[0] if "_" in reg_name else reg_name
                    regulon_data[reg_name] = {
                        "tf": tf,
                        "target_genes": group["gene"].tolist(),
                        "n_targets": len(group),
                    }
            break

    # Try loading cell-type-level regulon activity
    for candidate in ["regulon_activity_per_celltype.csv", "auc_mtx.csv"]:
        path = scenic_dir / candidate
        if path.exists():
            regulon_activity = pd.read_csv(path, index_col=0)
            print_step(f"Loaded regulon activity: {regulon_activity.shape}")
            break

    return regulon_data, regulon_activity


def run_lr_with_multisp_domains(adata, config):
    """Run squidpy L-R analysis using MultiSP domains as cluster key."""
    comm_config = config["communication"]

    # Ensure multisp_domain is categorical
    if "multisp_domain" not in adata.obs.columns:
        raise ValueError("multisp_domain not found; run 13b first")
    adata.obs["multisp_domain"] = adata.obs["multisp_domain"].astype("category")

    # Build spatial neighbors
    sq.gr.spatial_neighbors(adata, coord_type="generic",
                            n_neighs=comm_config["n_neighbors"])

    # Run L-R analysis
    print_step(f"Running squidpy L-R (n_perms={comm_config['n_perms']})")
    try:
        sq.gr.ligrec(
            adata, n_perms=comm_config["n_perms"],
            cluster_key="multisp_domain",
            use_raw=False,
            transmitter_params={"categories": "ligand"},
            receiver_params={"categories": "receptor"},
        )
    except Exception as e:
        print(f"    WARNING: L-R analysis failed: {e}")
        return pd.DataFrame()

    # Extract results
    if "multisp_domain_ligrec" not in adata.uns:
        return pd.DataFrame()

    pvals = adata.uns["multisp_domain_ligrec"]["pvalues"]
    means = adata.uns["multisp_domain_ligrec"]["means"]

    all_results = []
    for (source, target) in pvals.keys():
        lr_pvals = pvals[(source, target)]
        lr_means = means[(source, target)]
        for lr_pair in lr_pvals.index:
            p = lr_pvals.loc[lr_pair]
            m = lr_means.loc[lr_pair]
            if np.isfinite(p):
                all_results.append({
                    "source_domain": source, "target_domain": target,
                    "lr_pair": lr_pair, "pvalue": p, "mean_expr": m,
                })

    if not all_results:
        return pd.DataFrame()

    results_df = pd.DataFrame(all_results)

    # BH FDR correction
    _, padj, _, _ = multipletests(results_df["pvalue"].values, method="fdr_bh")
    results_df["padj_bh"] = padj
    sig = results_df[results_df["padj_bh"] < comm_config["p_threshold"]].copy()

    print(f"    Total L-R tests: {len(results_df)}")
    print(f"    Significant (padj<{comm_config['p_threshold']}): {len(sig)}")

    return sig


def compare_lr_with_05b(multisp_lr):
    """Compare MultiSP-domain L-R results with existing 05b cell-type-based results."""
    lr_05b_path = RESULTS_DIR / "communication" / "differential_lr_pairs.csv"
    if not lr_05b_path.exists():
        print_step("WARNING: 05b differential L-R not found, skipping comparison")
        return pd.DataFrame()

    lr_05b = pd.read_csv(lr_05b_path)
    lr_pairs_05b = set(lr_05b["lr_pair"].dropna().unique())
    lr_pairs_multisp = set(multisp_lr["lr_pair"].dropna().unique()) if len(multisp_lr) > 0 else set()

    shared = lr_pairs_05b & lr_pairs_multisp
    only_05b = lr_pairs_05b - lr_pairs_multisp
    only_multisp = lr_pairs_multisp - lr_pairs_05b

    comparison = pd.DataFrame([{
        "n_05b": len(lr_pairs_05b),
        "n_multisp": len(lr_pairs_multisp),
        "n_shared": len(shared),
        "n_only_05b": len(only_05b),
        "n_only_multisp": len(only_multisp),
        "jaccard": len(shared) / max(len(shared | only_05b | only_multisp), 1),
    }])
    print(f"    05b vs MultiSP L-R: shared={len(shared)}, "
          f"only_05b={len(only_05b)}, only_multisp={len(only_multisp)}")
    return comparison


def project_regulon_activity(adata, regulon_activity, config):
    """Project regulon activity onto Visium spots via c2l-weighted pseudobulk.

    For each spot: regulon_score = sum(c2l_proportion_ct * regulon_activity_ct).
    """
    if len(regulon_activity) == 0:
        print_step("WARNING: No regulon activity data, skipping projection")
        return pd.DataFrame()

    # Get c2l proportions
    c2l_cols = [c for c in adata.obs.columns if c.startswith(C2L_PREFIX)]
    if not c2l_cols:
        print_step("WARNING: No c2l proportions found")
        return pd.DataFrame()

    c2l_names = [strip_c2l_prefix(c) for c in c2l_cols]

    # Determine which cell types match between regulon_activity and c2l
    reg_cts = regulon_activity.columns.tolist() if regulon_activity.shape[1] > regulon_activity.shape[0] else regulon_activity.index.tolist()

    # Try harmonization
    reg_cts_harmonized = harmonize_cell_types(reg_cts, source="atac")
    ct_map = {}
    for orig, harm in zip(reg_cts, reg_cts_harmonized):
        if harm in c2l_names:
            ct_map[orig] = c2l_cols[c2l_names.index(harm)]

    if len(ct_map) == 0:
        print_step("WARNING: No cell types matched between regulon activity and c2l")
        return pd.DataFrame()

    print_step(f"Matched {len(ct_map)} cell types for regulon projection")

    # Ensure regulon_activity has regulons as rows, cell types as columns
    if set(reg_cts) & set(ct_map.keys()):
        # Cell types are columns
        reg_mat = regulon_activity[list(ct_map.keys())].values  # (n_regulons, n_matched_cts)
        regulon_names = regulon_activity.index.tolist()
    else:
        # Cell types are rows — transpose
        reg_mat = regulon_activity.loc[list(ct_map.keys())].values.T
        regulon_names = regulon_activity.columns.tolist()

    # Get matched c2l proportions
    prop_cols = [ct_map[ct] for ct in ct_map.keys()]
    props = adata.obs[prop_cols].values.copy()
    prop_sums = props.sum(axis=1, keepdims=True)
    prop_sums[prop_sums == 0] = 1.0
    props = props / prop_sums

    # Weighted sum: (n_spots, n_cts) @ (n_cts, n_regulons) -> (n_spots, n_regulons)
    projected = props @ reg_mat.T if reg_mat.shape[0] != props.shape[1] else props @ reg_mat
    if projected.shape[1] != len(regulon_names):
        projected = props @ reg_mat

    spot_regulon_df = pd.DataFrame(projected, index=adata.obs_names, columns=regulon_names[:projected.shape[1]])
    print_step(f"Projected regulon activity: {spot_regulon_df.shape}")
    return spot_regulon_df


def correlate_regulons_with_domains(spot_regulon_df, adata):
    """Identify domain-specific regulons via Mann-Whitney U tests."""
    if len(spot_regulon_df) == 0:
        return pd.DataFrame()

    domains = sorted(adata.obs["multisp_domain"].unique())
    records = []

    for reg in spot_regulon_df.columns:
        values = spot_regulon_df[reg].values
        if np.std(values) < 1e-10:
            continue

        for domain in domains:
            in_domain = adata.obs["multisp_domain"] == domain
            vals_in = values[in_domain.values if hasattr(in_domain, 'values') else in_domain]
            vals_out = values[~(in_domain.values if hasattr(in_domain, 'values') else in_domain)]

            if len(vals_in) < 5 or len(vals_out) < 5:
                continue

            stat, pval = mannwhitneyu(vals_in, vals_out, alternative="greater")
            mean_in = np.mean(vals_in)
            mean_out = np.mean(vals_out)
            fc = mean_in / (mean_out + 1e-10)

            records.append({
                "regulon": reg, "domain": domain,
                "mean_in_domain": mean_in, "mean_outside": mean_out,
                "fold_change": fc, "mwu_pval": pval,
            })

    if not records:
        return pd.DataFrame()

    result = pd.DataFrame(records)
    _, padj, _, _ = multipletests(result["mwu_pval"].values, method="fdr_bh")
    result["padj_bh"] = padj

    sig = result[result["padj_bh"] < 0.05]
    print_step(f"Domain-specific regulons: {sig['regulon'].nunique()} regulons "
               f"across {sig['domain'].nunique()} domains (padj<0.05)")
    return result


def build_communication_cards(multisp_lr, regulon_domain_df, adata):
    """Build multi-omic communication cards per domain pair.

    Combines L-R evidence, regulon activity, and chromatin status.
    """
    if len(multisp_lr) == 0:
        return pd.DataFrame()

    cards = []
    for (src, tgt), lr_group in multisp_lr.groupby(["source_domain", "target_domain"]):
        card = {
            "source_domain": src,
            "target_domain": tgt,
            "n_lr_pairs": len(lr_group),
            "top_lr_pairs": "; ".join(lr_group.nlargest(5, "mean_expr")["lr_pair"].astype(str).tolist()),
            "mean_lr_strength": lr_group["mean_expr"].mean(),
        }

        # Add regulon evidence for source and target domains
        if len(regulon_domain_df) > 0:
            src_regs = regulon_domain_df[
                (regulon_domain_df["domain"] == src) & (regulon_domain_df["padj_bh"] < 0.05)
            ]
            tgt_regs = regulon_domain_df[
                (regulon_domain_df["domain"] == tgt) & (regulon_domain_df["padj_bh"] < 0.05)
            ]
            card["n_source_regulons"] = len(src_regs)
            card["n_target_regulons"] = len(tgt_regs)
            card["top_source_regulons"] = "; ".join(
                src_regs.nlargest(3, "fold_change")["regulon"].tolist()
            ) if len(src_regs) > 0 else ""
            card["top_target_regulons"] = "; ".join(
                tgt_regs.nlargest(3, "fold_change")["regulon"].tolist()
            ) if len(tgt_regs) > 0 else ""
        else:
            card["n_source_regulons"] = 0
            card["n_target_regulons"] = 0
            card["top_source_regulons"] = ""
            card["top_target_regulons"] = ""

        # Add chromatin evidence (if available in obsm)
        if "chromatin" in adata.obsm:
            src_mask = adata.obs["multisp_domain"] == src
            tgt_mask = adata.obs["multisp_domain"] == tgt
            src_chrom = np.mean(adata.obsm["chromatin"][src_mask.values if hasattr(src_mask, 'values') else src_mask], axis=0)
            tgt_chrom = np.mean(adata.obsm["chromatin"][tgt_mask.values if hasattr(tgt_mask, 'values') else tgt_mask], axis=0)
            rho, _ = spearmanr(src_chrom, tgt_chrom)
            card["chromatin_correlation"] = rho
        else:
            card["chromatin_correlation"] = np.nan

        # Evidence layers count
        n_layers = 1  # L-R always present
        if card["n_source_regulons"] > 0 or card["n_target_regulons"] > 0:
            n_layers += 1
        if not np.isnan(card["chromatin_correlation"]):
            n_layers += 1
        card["n_evidence_layers"] = n_layers

        cards.append(card)

    return pd.DataFrame(cards)


def main():
    print_header("13c: Multi-Omic Communication Analysis")

    config = load_config()
    output_dir = RESULTS_DIR / "multisp"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load data ────────────────────────────────────────────────────
    print_step("Step 1: Loading spatial data with MultiSP domains")
    adata = load_multisp_adata()
    print(f"    {adata.n_obs} spots, {adata.obs['multisp_domain'].nunique()} domains")

    # ── Step 2: L-R analysis with MultiSP domains ────────────────────────────
    print_step("Step 2: Running L-R analysis with MultiSP domain clusters")
    multisp_lr = run_lr_with_multisp_domains(adata, config)
    if len(multisp_lr) > 0:
        save_csv(multisp_lr, "multisp_differential_lr.csv", subdir="multisp")

        # Top interactions
        print(f"\n    Top 10 L-R interactions:")
        for _, row in multisp_lr.nsmallest(10, "padj_bh").iterrows():
            print(f"      Domain {row['source_domain']} -> {row['target_domain']}: "
                  f"{row['lr_pair']} (padj={row['padj_bh']:.4f})")

    # ── Step 3: Compare with 05b results ─────────────────────────────────────
    print_step("Step 3: Comparing with 05b cell-type-based L-R results")
    lr_comparison = compare_lr_with_05b(multisp_lr)
    if len(lr_comparison) > 0:
        save_csv(lr_comparison, "multisp_vs_05b_lr_comparison.csv", subdir="multisp")

    # ── Step 4: Project SCENIC+ regulon activity ─────────────────────────────
    print_step("Step 4: Loading and projecting SCENIC+ regulon activity")
    regulon_data, regulon_activity = load_scenic_regulons()
    print(f"    Regulons loaded: {len(regulon_data)} definitions, "
          f"activity matrix: {regulon_activity.shape if len(regulon_activity) > 0 else 'N/A'}")

    spot_regulon_df = project_regulon_activity(adata, regulon_activity, config)
    if len(spot_regulon_df) > 0:
        save_csv(spot_regulon_df, "spot_regulon_activity.csv", subdir="multisp")

    # ── Step 5: Correlate regulons with domains ──────────────────────────────
    print_step("Step 5: Identifying domain-specific regulons")
    regulon_domain_df = correlate_regulons_with_domains(spot_regulon_df, adata)
    if len(regulon_domain_df) > 0:
        save_csv(regulon_domain_df, "regulon_domain_activity.csv", subdir="multisp")

    # ── Step 6: Build multi-omic communication cards ─────────────────────────
    print_step("Step 6: Building multi-omic communication cards")
    cards = build_communication_cards(multisp_lr, regulon_domain_df, adata)
    if len(cards) > 0:
        save_csv(cards, "multi_omic_communication_cards.csv", subdir="multisp")
        print(f"\n    Communication cards: {len(cards)} domain pairs")
        multi_layer = cards[cards["n_evidence_layers"] >= 2]
        print(f"    Multi-layer cards (>=2 evidence layers): {len(multi_layer)}")

    # Summary
    print(f"\n  Summary:")
    print(f"    Significant L-R pairs: {len(multisp_lr)}")
    print(f"    Regulons projected: {len(spot_regulon_df.columns) if len(spot_regulon_df) > 0 else 0}")
    n_sig_regs = regulon_domain_df["regulon"].nunique() if len(regulon_domain_df) > 0 and "padj_bh" in regulon_domain_df.columns else 0
    print(f"    Domain-specific regulons: {n_sig_regs}")
    print(f"    Communication cards: {len(cards)}")

    print_header("13c: Complete")


if __name__ == "__main__":
    main()
