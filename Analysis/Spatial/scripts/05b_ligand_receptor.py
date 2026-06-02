#!/usr/bin/env python3
"""
05b_ligand_receptor.py — Spatial ligand-receptor interaction analysis.

Identifies spatially-resolved L-R interactions and builds cell-cell
communication networks comparing healthy vs. MASLD.

Donor-aware (2026-06-01 rigor pass):
  * F038 — ligrec runs on a log1p-CPM layer (`ensure_lognorm`), NOT the raw
    integer counts in `.X` (squidpy's per-cluster mean statistic on raw counts
    is depth-confounded).
  * F036 — the squidpy permutation test is run PER DONOR (sample_id) and the
    per-(source,target,LR) results are aggregated across donors, instead of
    pooling ~10^3 spatially-autocorrelated spots within a condition as if they
    were independent. The biological n is ~2 donors/condition (GSE192741:
    Healthy=JBO018/022, Steatotic=JBO014/015/019), so per-condition results are
    a donor-aggregated summary, not a spot-level permutation p-value. The
    `pvalue`/`padj_bh` columns are now the donor-median squidpy p (still named
    the same for downstream consumers, but their statistical meaning is the
    donor-aggregated value — see column note below).
  * F039 — `seed=42` passed to ligrec (and np.random.seed at module start) so
    the permutation numbers are reproducible.
  * F040 — the LR database is logged from config; squidpy's default resource is
    OmniPath/Intercell (not CellPhoneDB) — see warning in main().

SLURM: --partition=cpu --cpus=8 --mem=64G --time=4:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from statsmodels.stats.multitest import multipletests

np.random.seed(42)  # F039: lock the spatial_neighbors / permutation RNG

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, C2L_PREFIX, load_config, load_deconvolved_adata,
    save_csv, print_header,
)
from spatial_stats import ensure_lognorm

# Donor (biological replicate) column for GSE192741 — DO NOT pool across these.
DONOR_COL = "sample_id"
# Minimum spots per donor slice to attempt a ligrec permutation test.
MIN_SPOTS_PER_DONOR = 100


def _run_ligrec_one_slice(adata_slice, lognorm_layer, n_perms):
    """Run squidpy ligrec on a single donor slice; return per-(src,tgt,LR) rows.

    Reads expression from `lognorm_layer` (log1p-CPM) so the per-cluster mean
    statistic is not depth-confounded (F038). Returns a tidy DataFrame with one
    row per (source, target, lr_pair): pvalue + mean_expr for this donor.
    """
    # Point .X at the log-normalized layer for the duration of this test (F038).
    adata_slice = adata_slice.copy()
    adata_slice.X = adata_slice.layers[lognorm_layer]

    if "cell_type_dominant" in adata_slice.obs.columns:
        adata_slice.obs["cell_type_dominant"] = (
            adata_slice.obs["cell_type_dominant"]
            .str.replace(C2L_PREFIX, "", regex=False)
            .astype("category")
            .cat.remove_unused_categories()
        )

    sq.gr.spatial_neighbors(adata_slice, coord_type="generic", n_neighs=6)
    sq.gr.ligrec(
        adata_slice, n_perms=n_perms,
        cluster_key="cell_type_dominant",
        use_raw=False,
        seed=42,  # F039: deterministic permutation p-values
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
    )
    if "cell_type_dominant_ligrec" not in adata_slice.uns:
        return pd.DataFrame()

    pvals = adata_slice.uns["cell_type_dominant_ligrec"]["pvalues"]
    means = adata_slice.uns["cell_type_dominant_ligrec"]["means"]
    rows = []
    for (source, target) in pvals.keys():
        lr_pvals = pvals[(source, target)]
        lr_means = means[(source, target)]
        for lr_pair in lr_pvals.index:
            p = lr_pvals.loc[lr_pair]
            m = lr_means.loc[lr_pair]
            if np.isfinite(p):
                rows.append({
                    "source": source, "target": target,
                    "lr_pair": lr_pair, "pvalue": p, "mean_expr": m,
                })
    return pd.DataFrame(rows)


def run_ligrec(adata, condition, lognorm_layer, n_perms=10000):
    """Run squidpy L-R analysis for one condition, PER DONOR then aggregate (F036).

    Each donor slice is tested independently (no cross-donor / cross-slice
    pooling). Per-(source,target,lr_pair) results are aggregated across donors:
    `pvalue` = median donor p (donor-aggregated, NOT a spot-level p — name kept
    for downstream compatibility), `mean_expr` = mean of per-donor means,
    `n_donors` / `frac_donors_sig` expose the honest donor support.
    """
    adata_cond = adata[adata.obs["condition"] == condition]
    donors = adata_cond.obs[DONOR_COL].astype(str)
    per_donor = []
    n_used = 0
    for donor in sorted(donors.unique()):
        sl = adata_cond[donors.values == donor]
        if sl.n_obs < MIN_SPOTS_PER_DONOR:
            print(f"    {donor}: {sl.n_obs} spots < {MIN_SPOTS_PER_DONOR}, skipping")
            continue
        try:
            df_d = _run_ligrec_one_slice(sl, lognorm_layer, n_perms)
        except Exception as e:
            print(f"    WARNING: ligrec failed for {donor}: {e}")
            continue
        if len(df_d):
            df_d[DONOR_COL] = donor
            per_donor.append(df_d)
            n_used += 1
    if not per_donor:
        print(f"    WARNING: no donor with >= {MIN_SPOTS_PER_DONOR} spots for {condition}")
        return None

    alld = pd.concat(per_donor, ignore_index=True)
    # Aggregate across donors per (source, target, lr_pair). The condition has
    # only ~2 donors, so this is a donor-level summary, not a powered test.
    grp = alld.groupby(["source", "target", "lr_pair"], observed=True)
    agg = grp.agg(
        pvalue=("pvalue", "median"),          # donor-median p (F036 — see docstring)
        mean_expr=("mean_expr", "mean"),
        n_donors=("pvalue", "size"),
    ).reset_index()
    agg["frac_donors_sig"] = grp["pvalue"].apply(lambda s: float((s < 0.05).mean())).values
    agg["n_donors_total"] = n_used
    return agg


def extract_significant_lr(agg_df, p_threshold=0.01):
    """Extract significant L-R pairs from donor-aggregated ligrec results.

    `agg_df` is the per-(source,target,lr_pair) donor-aggregated table from
    run_ligrec. BH-FDR is applied to the donor-median p (F036): padj_bh is a
    donor-aggregated FDR, not a spot-level FDR (name kept for consumers).
    """
    if agg_df is None or len(agg_df) == 0:
        return pd.DataFrame()

    results_df = agg_df.copy()
    # Tag autocrine vs paracrine
    results_df["pair_type"] = np.where(
        results_df["source"] == results_df["target"], "autocrine", "paracrine")

    # Apply BH FDR correction across all donor-median p-values
    _, padj, _, _ = multipletests(results_df["pvalue"].values, method="fdr_bh")
    results_df["padj_bh"] = padj

    # Filter by FDR-corrected threshold
    sig = results_df[results_df["padj_bh"] < p_threshold].copy()
    return sig


def build_communication_network(sig_lr, p_threshold=0.01):
    """Build weighted cell-type communication network from L-R results."""
    import networkx as nx

    G = nx.DiGraph()
    if len(sig_lr) == 0:
        return G

    for (source, target), group in sig_lr.groupby(["source", "target"]):
        n_sig = len(group)
        mean_strength = group["mean_expr"].mean()
        G.add_edge(source, target, weight=n_sig, mean_strength=mean_strength)

    return G


def main():
    print_header("05b: Ligand-Receptor Interaction Analysis")

    config = load_config()
    comm_config = config["communication"]
    output_dir = RESULTS_DIR / "communication"
    output_dir.mkdir(parents=True, exist_ok=True)

    # F040: the config declares lr_database, but squidpy.gr.ligrec is called with
    # interactions=None and therefore uses its default OmniPath/Intercell
    # resource, NOT the declared database. Surface the mismatch instead of
    # silently disagreeing with the config/methods.
    declared_db = comm_config.get("lr_database", "OmniPath")
    print(f"  Configured lr_database: {declared_db}")
    if str(declared_db).lower() not in ("omnipath", "intercell"):
        print(f"    NOTE (F040): squidpy.gr.ligrec uses its default OmniPath/Intercell "
              f"resource; '{declared_db}' is NOT wired in. Report the LR database as "
              f"'OmniPath via squidpy.gr.ligrec default' in methods, or supply "
              f"interactions= built from {declared_db}.")

    adata = load_deconvolved_adata()
    print(f"  Loaded: {adata.n_obs} spots, {adata.obs[DONOR_COL].nunique()} donors")

    # F038: build a log1p-CPM layer once; ligrec reads this, never raw .X counts.
    lognorm_layer = ensure_lognorm(adata)
    print(f"  Using normalized layer for ligrec: '{lognorm_layer}'")

    conditions = adata.obs["condition"].unique().tolist()
    all_sig = {}

    for condition in conditions:
        cond_donors = adata.obs.loc[adata.obs["condition"] == condition, DONOR_COL].nunique()
        print(f"\n  Running L-R analysis: {condition} (n={cond_donors} donors)...")
        agg = run_ligrec(adata, condition, lognorm_layer,
                         n_perms=comm_config["n_perms"])
        if agg is None:
            continue

        sig_lr = extract_significant_lr(agg, p_threshold=comm_config["p_threshold"])
        all_sig[condition] = sig_lr
        save_csv(sig_lr, f"ligrec_{condition}.csv", subdir="communication")

        print(f"    Significant L-R pairs (donor-aggregated padj_bh<{comm_config['p_threshold']}): {len(sig_lr)}")
        if len(sig_lr) > 0:
            # Summary by pair type
            for pt, count in sig_lr["pair_type"].value_counts().items():
                print(f"      {pt}: {count}")
            print(f"    Top 10:")
            for _, row in sig_lr.sort_values("padj_bh").head(10).iterrows():
                print(f"      {row['source']} → {row['target']}: "
                      f"{row['lr_pair']} (padj={row['padj_bh']:.4f}, {row['pair_type']})")

        # Highlight MASLD-relevant pathways
        for axis_name, axis_pairs in [
            ("fibrosis", comm_config.get("fibrosis_lr", [])),
            ("inflammation", comm_config.get("inflammation_lr", [])),
            ("sinusoidal", comm_config.get("sinusoidal_lr", [])),
        ]:
            if axis_pairs and len(sig_lr) > 0:
                axis_hits = sig_lr[sig_lr["lr_pair"].isin(axis_pairs)]
                if len(axis_hits) > 0:
                    print(f"    {axis_name.upper()} axis hits: {len(axis_hits)}")

    # Differential L-R analysis with effect sizes
    if len(all_sig) >= 2:
        cond_list = list(all_sig.keys())
        cond_h = [c for c in cond_list if "healthy" in c.lower() or c == "Healthy"]
        cond_m = [c for c in cond_list if c not in cond_h]
        if cond_h and cond_m:
            sig_h = all_sig[cond_h[0]]
            sig_m = all_sig[cond_m[0]]
            lr_h = set(sig_h["lr_pair"]) if len(sig_h) > 0 else set()
            lr_m = set(sig_m["lr_pair"]) if len(sig_m) > 0 else set()
            emergent = lr_m - lr_h
            lost = lr_h - lr_m
            print(f"\n  Differential L-R:")
            print(f"    Disease-emergent: {len(emergent)}")
            print(f"    Disease-lost: {len(lost)}")

            # Build lookup dicts for fast access (lr_pair may be tuples)
            lr_to_row_m = {str(row["lr_pair"]): row for _, row in sig_m.iterrows()} if len(sig_m) > 0 else {}
            lr_to_row_h = {str(row["lr_pair"]): row for _, row in sig_h.iterrows()} if len(sig_h) > 0 else {}

            # Build differential table with effect sizes
            diff_rows = []
            for lr in emergent:
                lr_key = str(lr)
                row_m = lr_to_row_m.get(lr_key)
                mean_m = float(row_m["mean_expr"]) if row_m is not None else np.nan
                row_h_match = lr_to_row_h.get(lr_key)
                mean_h = float(row_h_match["mean_expr"]) if row_h_match is not None else 0.0
                diff_rows.append({
                    "lr_pair": lr, "category": "disease_emergent",
                    "mean_expr_disease": mean_m, "mean_expr_healthy": mean_h,
                    "delta_mean_expr": mean_m - mean_h,
                    "source": row_m["source"] if row_m is not None else np.nan,
                    "target": row_m["target"] if row_m is not None else np.nan,
                    "pair_type": row_m["pair_type"] if row_m is not None else np.nan,
                })
            for lr in lost:
                lr_key = str(lr)
                row_h = lr_to_row_h.get(lr_key)
                mean_h = float(row_h["mean_expr"]) if row_h is not None else np.nan
                row_m_match = lr_to_row_m.get(lr_key)
                mean_m = float(row_m_match["mean_expr"]) if row_m_match is not None else 0.0
                diff_rows.append({
                    "lr_pair": lr, "category": "disease_lost",
                    "mean_expr_disease": mean_m, "mean_expr_healthy": mean_h,
                    "delta_mean_expr": mean_m - mean_h,
                    "source": row_h["source"] if row_h is not None else np.nan,
                    "target": row_h["target"] if row_h is not None else np.nan,
                    "pair_type": row_h["pair_type"] if row_h is not None else np.nan,
                })
            diff_df = pd.DataFrame(diff_rows)
            save_csv(diff_df, "differential_lr_pairs.csv", subdir="communication")

    # Build communication networks
    import networkx as nx
    for condition, sig_lr in all_sig.items():
        G = build_communication_network(sig_lr)
        print(f"\n  {condition} network: {G.number_of_nodes()} nodes, "
              f"{G.number_of_edges()} edges, density={nx.density(G):.3f}")
        # Save as edge list
        edges = []
        for u, v, data in G.edges(data=True):
            edges.append({"source": u, "target": v, **data})
        save_csv(pd.DataFrame(edges), f"network_{condition}.csv",
                 subdir="communication")

    print_header("05b: Complete")


if __name__ == "__main__":
    main()
