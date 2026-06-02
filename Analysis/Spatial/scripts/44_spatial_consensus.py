#!/usr/bin/env python3
"""
44_spatial_consensus.py — Multi-cohort spatial consensus columns.

Aggregates per-cohort SVG / Moran's I results into consensus columns
(merged downstream by 06_integration.py). Currently integrates:

  - GSE192741 Visium (conditions: Healthy + Steatotic) — squidpy Moran's I.
    These are two CONDITIONS of ONE physical dataset (GSE192741), not two
    independent cohorts (F072), and are collapsed to one dataset for both the
    n_cohorts count and the Moran's-I mean.
  - Govaere2026 GeoMx SH-vs-PT — limma-style regional DE. |logFC| is a
    fold-change magnitude on a DIFFERENT, unbounded scale than Moran's I, so
    it is kept in its OWN column and never averaged with Moran's I (F074/F181);
    padj_bh < 0.1 binarizes the "svg" flag.

Outputs (Analysis/Spatial/results/integration/spatial_consensus.csv):
  - spatial_consensus_svg_n_cohorts        (int)  — number of independent
        DATASETS (not condition-level cohorts) with the gene SVG-significant.
  - spatial_consensus_morans_i_mean        (float) — mean of TRUE squidpy
        Moran's I across Visium datasets only, bounded ~[-1,1]. Does NOT
        include GeoMx |logFC|.
  - spatial_consensus_geomx_abs_logfc_mean (float) — mean |logFC| across GeoMx
        cohorts; a DE magnitude, NOT an autocorrelation statistic.
  - spatial_consensus_direction            (str)
        concordant_up | concordant_down | discordant | concordant_unsigned |
        single_cohort_only

NOTE: "direction" is only defined for cohorts with sign information
(Govaere logFC). The Visium Moran's I cohorts are autocorrelation
magnitudes (no sign); they count toward n_cohorts and morans_i_mean but
not toward direction. If only Visium cohort(s) contribute and there is
no signed cohort, direction = "single_cohort_only". With the current data
the consensus is dominated by one Visium dataset (GSE192741) and the
signed-direction axis has no positive instances (GeoMx never co-flags a
Visium SVG), so n_cohorts==2 means both GSE192741 conditions agree.
"""

import pathlib
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import RESULTS_DIR, print_header

SVG_PADJ_THRESH = 0.1
GOVAERE_LFC_FLOOR = 0.0  # |logFC| is used directly as a proxy magnitude


def _load_visium_svg_cohort(path: pathlib.Path, cohort_label: str,
                            dataset_label: str) -> pd.DataFrame:
    """Load a squidpy SVG CSV; returns long-format (gene, morans_i, svg, cohort).

    ``dataset_label`` is the underlying physical dataset (e.g. GSE192741), so
    that two conditions of the SAME Visium dataset are not counted as two
    independent cohorts downstream (F072). ``stat_type='morans_i'`` marks that
    the magnitude is a true squidpy Moran's I (bounded ~[-1,1]).
    """
    if not path.exists():
        print(f"  [skip] {path.name} not found")
        return pd.DataFrame()
    df = pd.read_csv(path, index_col=0)
    if "I" not in df.columns or "svg" not in df.columns:
        print(f"  [skip] {path.name} missing required cols (I, svg)")
        return pd.DataFrame()
    out = pd.DataFrame({
        "gene": df.index.astype(str),
        "magnitude": df["I"].astype(float),   # true squidpy Moran's I
        "svg": df["svg"].astype(bool),
        "logfc": np.nan,            # Visium SVG has no signed direction
        "cohort": cohort_label,
        "dataset": dataset_label,
        "stat_type": "morans_i",
    })
    print(f"  [{cohort_label}] {len(out)} genes, {out['svg'].sum()} SVG")
    return out


def _load_geomx_cohort(path: pathlib.Path, cohort_label: str,
                       dataset_label: str) -> pd.DataFrame:
    """Load a GeoMx regional DE CSV; returns long-format with signed logfc.

    ``stat_type='abs_logfc'`` marks that the magnitude is |logFC| (a fold-
    change magnitude, unbounded), NOT a Moran's I. Downstream this is kept in
    a SEPARATE consensus column so the two incommensurable scales are never
    averaged together (F074/F181).
    """
    if not path.exists():
        print(f"  [skip] {path.name} not found")
        return pd.DataFrame()
    df = pd.read_csv(path)
    needed = {"gene_symbol", "logFC", "padj_bh"}
    if not needed.issubset(df.columns):
        print(f"  [skip] {path.name} missing required cols ({needed - set(df.columns)})")
        return pd.DataFrame()
    df = df.dropna(subset=["gene_symbol", "logFC", "padj_bh"])
    # |logFC| is a regional-variability magnitude on a DIFFERENT scale than
    # Moran's I — stored under `magnitude` and tagged stat_type='abs_logfc'.
    out = pd.DataFrame({
        "gene": df["gene_symbol"].astype(str),
        "magnitude": df["logFC"].abs().astype(float),
        "svg": (df["padj_bh"].astype(float) < SVG_PADJ_THRESH) &
               (df["logFC"].abs() > GOVAERE_LFC_FLOOR),
        "logfc": df["logFC"].astype(float),
        "cohort": cohort_label,
        "dataset": dataset_label,
        "stat_type": "abs_logfc",
    })
    print(f"  [{cohort_label}] {len(out)} genes, {out['svg'].sum()} SVG (padj<{SVG_PADJ_THRESH})")
    return out


def build_consensus(long_df: pd.DataFrame) -> pd.DataFrame:
    """Collapse the long-format multi-cohort table to per-gene consensus.

    spatial_consensus_svg_n_cohorts : number of independent DATASETS (not
        condition-level cohorts) in which the gene is SVG-significant. The two
        GSE192741 Visium conditions (Healthy + Steatotic) collapse to a single
        dataset so they are not double-counted as two independent cohorts
        (F072). With current data the only second source (Govaere GeoMx) never
        co-flags a Visium SVG, so n_cohorts==2 means the gene was significant in
        both GSE192741 conditions, i.e. one physical Visium dataset.
    spatial_consensus_morans_i_mean : mean of TRUE squidpy Moran's I across the
        Visium cohorts only (stat_type=='morans_i'); bounded ~[-1,1]. GeoMx
        |logFC| is NOT mixed in here — that lives in
        spatial_consensus_geomx_abs_logfc_mean (different, unbounded scale).
        Same-dataset conditions are averaged to one value first so a gene in
        both GSE192741 conditions is not double-weighted (F074/F181).
    spatial_consensus_geomx_abs_logfc_mean : mean of |logFC| across GeoMx
        cohorts (stat_type=='abs_logfc'); a fold-change magnitude, NOT a
        spatial-autocorrelation statistic. Kept separate from Moran's I.
    spatial_consensus_direction     : derived from signed cohorts only.
    """
    if long_df.empty:
        return pd.DataFrame(columns=[
            "gene",
            "spatial_consensus_svg_n_cohorts",
            "spatial_consensus_morans_i_mean",
            "spatial_consensus_geomx_abs_logfc_mean",
            "spatial_consensus_direction",
        ])

    # n_cohorts = number of independent DATASETS where the gene is SVG-significant
    # (F072: collapse same-dataset conditions before counting). For each
    # (gene, dataset) take whether ANY cohort of that dataset flagged it SVG,
    # then count datasets with a True.
    svg_by_dataset = (long_df.groupby(["gene", "dataset"])["svg"].any()
                             .reset_index())
    n_cohorts = (svg_by_dataset[svg_by_dataset["svg"]]
                 .groupby("gene")["dataset"].nunique())

    # Moran's I mean: Visium (true Moran's I) ONLY (F074/F181). Average within
    # dataset first (collapse the two GSE192741 conditions), then across
    # datasets, so a gene present in both conditions is not double-weighted.
    visium = long_df[long_df["stat_type"] == "morans_i"]
    if not visium.empty:
        morans_mean = (visium.groupby(["gene", "dataset"])["magnitude"].mean()
                             .groupby("gene").mean())
    else:
        morans_mean = pd.Series(dtype=float)

    # GeoMx |logFC| mean kept in a SEPARATE column (different scale).
    geomx = long_df[long_df["stat_type"] == "abs_logfc"]
    if not geomx.empty:
        geomx_mean = (geomx.groupby(["gene", "dataset"])["magnitude"].mean()
                           .groupby("gene").mean())
    else:
        geomx_mean = pd.Series(dtype=float)

    # direction — only signed cohorts (logfc not NaN) participate
    signed = long_df.dropna(subset=["logfc"]).copy()
    # Restrict to cohorts where SVG flag is True for the direction call;
    # this avoids letting noisy non-significant signs decide direction.
    signed_sig = signed[signed["svg"]]
    if not signed_sig.empty:
        sign_per_cohort = (signed_sig.assign(sign=np.sign(signed_sig["logfc"]))
                                       .groupby(["gene", "cohort"])["sign"].mean()
                                       .reset_index())
        # per-gene: set of unique signs and number of signed cohorts
        per_gene = sign_per_cohort.groupby("gene")["sign"].agg(["nunique", "first", "count"])
        per_gene = per_gene.rename(columns={"nunique": "n_signs",
                                            "first": "sign_first",
                                            "count": "n_signed_cohorts"})
    else:
        per_gene = pd.DataFrame(columns=["n_signs", "sign_first", "n_signed_cohorts"])

    # Index over ALL genes seen in any cohort (a gene may appear only in GeoMx,
    # so morans_mean.index alone would drop it).
    all_genes = pd.Index(sorted(long_df["gene"].unique()))
    out = pd.DataFrame({"gene": all_genes})
    out["spatial_consensus_svg_n_cohorts"] = out["gene"].map(n_cohorts).fillna(0).astype(int)
    out["spatial_consensus_morans_i_mean"] = out["gene"].map(morans_mean)
    out["spatial_consensus_geomx_abs_logfc_mean"] = out["gene"].map(geomx_mean)

    # Direction labels:
    #   concordant_up   — >=1 SVG-significant signed cohort, all signs > 0
    #                     (or >= 2 signed cohorts agreeing on +)
    #   concordant_down — analogous for negative
    #   discordant      — >= 2 signed cohorts with conflicting signs
    #   concordant_unsigned — >= 2 SVG-significant cohorts but no signed
    #                         cohort flagged the gene (e.g. only the two
    #                         Visium Moran's-I conditions agree)
    #   single_cohort_only — fewer than 2 SVG-significant cohorts overall
    def _direction(g):
        n_c = int(n_cohorts.get(g, 0))
        if g in per_gene.index:
            row = per_gene.loc[g]
            n_signed = int(row["n_signed_cohorts"])
            n_signs = int(row["n_signs"])
            if n_signed >= 2 and n_signs > 1:
                return "discordant"
            if n_signed >= 2 and n_signs == 1:
                return "concordant_up" if row["sign_first"] > 0 else "concordant_down"
            # exactly one signed cohort with svg=True; let n_cohorts settle it
            if n_c >= 2:
                # combine 1 signed + at least 1 unsigned (Visium) → call by sign
                return "concordant_up" if row["sign_first"] > 0 else "concordant_down"
            return "concordant_up" if row["sign_first"] > 0 else "concordant_down"
        # no signed cohort flagged this gene
        if n_c >= 2:
            return "concordant_unsigned"
        return "single_cohort_only"

    out["spatial_consensus_direction"] = out["gene"].map(_direction)
    return out


def main():
    print_header("44: Spatial Consensus (multi-cohort SVG / Moran's I)")

    # ---------- Dataset 1: GSE192741 Visium (squidpy SVG, two CONDITIONS) ----------
    # Healthy + Steatotic are two conditions of ONE physical dataset (GSE192741),
    # tagged with the same dataset label so they collapse to one cohort (F072).
    visium_dfs = []
    healthy = _load_visium_svg_cohort(
        RESULTS_DIR / "svg" / "svgs_Healthy.csv",
        cohort_label="GSE192741_Healthy",
        dataset_label="GSE192741",
    )
    steatotic = _load_visium_svg_cohort(
        RESULTS_DIR / "svg" / "svgs_Steatotic.csv",
        cohort_label="GSE192741_Steatotic",
        dataset_label="GSE192741",
    )
    if not healthy.empty:
        visium_dfs.append(healthy)
    if not steatotic.empty:
        visium_dfs.append(steatotic)

    # ---------- Dataset 2: Govaere2026 GeoMx SH-vs-PT ----------
    geomx_pt = _load_geomx_cohort(
        RESULTS_DIR / "govaere2026" / "geomx_de_sh_vs_pt.csv",
        cohort_label="Govaere2026_GeoMx_SHvsPT",
        dataset_label="Govaere2026_GeoMx",
    )

    all_dfs = visium_dfs + ([geomx_pt] if not geomx_pt.empty else [])
    if not all_dfs:
        print("ERROR: no per-cohort spatial inputs available.")
        sys.exit(1)
    long_df = pd.concat(all_dfs, ignore_index=True)
    print(f"\n  Combined long-format rows: {len(long_df):,}")
    print(f"  Unique cohorts (conditions): {long_df['cohort'].nunique()} "
          f"({sorted(long_df['cohort'].unique().tolist())})")
    print(f"  Unique datasets: {long_df['dataset'].nunique()} "
          f"({sorted(long_df['dataset'].unique().tolist())})")
    print(f"  Unique genes: {long_df['gene'].nunique():,}")

    # ---------- Consensus ----------
    consensus = build_consensus(long_df)
    print(f"\n  Consensus genes: {len(consensus):,}")
    print(f"  n_cohorts dist: {consensus['spatial_consensus_svg_n_cohorts'].value_counts().sort_index().to_dict()}")
    print(f"  direction dist: {consensus['spatial_consensus_direction'].value_counts().to_dict()}")

    # ---------- Write ----------
    out_path = RESULTS_DIR / "integration" / "spatial_consensus.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    consensus.to_csv(out_path, index=False)
    print(f"\n  Saved: {out_path}")

    # ---------- Sanity ----------
    for g in ["GPNMB", "LPL", "FABP5"]:
        row = consensus[consensus["gene"] == g]
        if not row.empty:
            r = row.iloc[0]
            mi = r['spatial_consensus_morans_i_mean']
            gx = r['spatial_consensus_geomx_abs_logfc_mean']
            mi_s = f"{mi:.3f}" if pd.notna(mi) else "NA"
            gx_s = f"{gx:.3f}" if pd.notna(gx) else "NA"
            print(f"    {g}: n_cohorts={r['spatial_consensus_svg_n_cohorts']}, "
                  f"morans_i_mean={mi_s}, geomx_abs_logfc_mean={gx_s}, "
                  f"direction={r['spatial_consensus_direction']}")

    print_header("44: Complete")


if __name__ == "__main__":
    main()
