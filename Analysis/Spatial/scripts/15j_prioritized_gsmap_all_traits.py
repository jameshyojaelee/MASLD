#!/usr/bin/env python
"""15j — Prioritized-gene recovery among spatial GWAS-risk genes.

This is the main Fig. 4f data-prep path. It keeps the gsMap gene-risk
localization definition from 15f, then summarizes only the Fig. 4 prioritized
gene set across every gsMap-compatible GWAS trait in the current config.

Outputs:
  Analysis/Spatial/results/gsmap/gene_risk_localization_18trait.csv
  Analysis/Spatial/results/gsmap/prioritized_spatial_risk_18trait.csv

Run after 15e has produced spatial_ldsc outputs for all traits in
Analysis/Spatial/data/gsmap_gwas/gwas_config.yaml.
"""
from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

pd.set_option("compute.use_numexpr", False)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SCRIPTS = BASE / "Analysis/Spatial/scripts"
GS = BASE / "Analysis/Spatial/results/gsmap"
UDIR = BASE / "Analysis/Spatial/results/universe_validation"
GWAS_CONFIG = BASE / "Analysis/Spatial/data/gsmap_gwas/gwas_config.yaml"

TRAIT_META = {
    "ukbb_alt": ("ALT (UKB)", "Liver enzymes", 1),
    "mvp_alt": ("ALT (MVP EUR)", "Liver enzymes", 2),
    "ukbb_ast": ("AST (UKB)", "Liver enzymes", 3),
    "mvp_ast": ("AST (MVP EUR)", "Liver enzymes", 4),
    "ukbb_ggt": ("GGT (UKB)", "Liver enzymes", 5),
    "finngen_nafld": ("NAFLD (FinnGen)", "NAFLD / NASH", 10),
    "ghodsian_nafld": ("NAFLD (Ghodsian)", "NAFLD / NASH", 11),
    "mvp_nafld": ("NAFLD (MVP EUR)", "NAFLD / NASH", 12),
    "decode_nafld": ("NAFLD (deCODE)", "NAFLD / NASH", 13),
    "ukbb2023_nafld": ("NAFLD (UKB 2023)", "NAFLD / NASH", 14),
    "intermtn_nafld": ("NAFLD (Intermountain)", "NAFLD / NASH", 15),
    "anstee2020_nafld": ("NAFLD (Anstee 2020)", "NAFLD / NASH", 16),
    "nafld_2019": ("NAFLD (Namjou 2019)", "NAFLD / NASH", 17),
    "finngen_nash": ("NASH (FinnGen)", "NAFLD / NASH", 18),
    "pdff": ("PDFF (Pazoki 2022)", "Imaging liver fat", 30),
    "pdff_2021a": ("PDFF (2021a)", "Imaging liver fat", 31),
    "pdff_2021b": ("PDFF (2021b)", "Imaging liver fat", 32),
    "pdff_2022": ("PDFF (2022)", "Imaging liver fat", 33),
}

COHORT_META = {
    "gse192741": ("GSE192741", 1),
    "vu": ("Vu et al. 2025", 2),
}


def read_simple_yaml_keys(path: Path) -> list[str]:
    """Read top-level keys from the simple trait:path gsMap YAML config."""
    keys = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ":" not in line:
            continue
        key = line.split(":", 1)[0].strip()
        if key:
            keys.append(key)
    return keys


def read_gene_set(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def load_gsmap_results():
    """Load per-spot spatial LDSC p-values from gsMap output."""
    results = {}

    for dataset_dir in sorted(GS.glob("*")):
        if not dataset_dir.is_dir() or dataset_dir.name == "figures":
            continue
        dataset = dataset_dir.name
        results[dataset] = {}

        for sample_dir in sorted(dataset_dir.glob("*")):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            ldsc_dir = sample_dir / "spatial_ldsc"
            if not ldsc_dir.exists():
                continue

            results[dataset][sample] = {}
            for ldsc_file in sorted(ldsc_dir.glob(f"{sample}_*.csv.gz")):
                trait = ldsc_file.stem.replace(f"{sample}_", "").replace(".csv", "")
                try:
                    df = pd.read_csv(ldsc_file, compression="gzip")
                except Exception as e:
                    print(f"    WARNING: could not load {ldsc_file}: {e}")
                    continue
                if {"spot", "p"}.issubset(df.columns) and len(df) > 0:
                    results[dataset][sample][trait] = df

    return results


def load_gsmap_marker_scores():
    """Load per-sample gene-x-spot marker-score matrices from gsMap."""
    marker = {}

    for dataset_dir in sorted(GS.glob("*")):
        if not dataset_dir.is_dir() or dataset_dir.name == "figures":
            continue
        dataset = dataset_dir.name
        marker[dataset] = {}

        for sample_dir in sorted(dataset_dir.glob("*")):
            if not sample_dir.is_dir():
                continue
            sample = sample_dir.name
            feather = sample_dir / "latent_to_gene" / f"{sample}_gene_marker_score.feather"
            if not feather.exists():
                continue
            try:
                df = pd.read_feather(feather)
            except Exception as e:
                print(f"    WARNING: could not load {feather}: {e}")
                continue
            if "HUMAN_GENE_SYM" not in df.columns:
                continue
            df = df.set_index("HUMAN_GENE_SYM")
            if df.index.duplicated().any():
                df = df.groupby(level=0).max()
            marker[dataset][sample] = df

    return marker


def gene_risk_localization(gsmap_results, marker_scores,
                           spot_sig_alpha=0.05, top_gene_frac=0.05,
                           min_sig_spots=10):
    """Derive gene-level spatial GWAS-risk localization from spot p-values.

    This mirrors 15f_gsmap_analysis.py without importing scanpy-heavy utilities:
    select LDSC-significant spots per sample/trait, score genes by marker-score
    enrichment in those spots, and flag the top positive-enrichment genes.
    """
    records = []

    for dataset, samples in gsmap_results.items():
        ds_marker = marker_scores.get(dataset, {})
        acc = defaultdict(lambda: {"enrich": [], "min_pval": [], "n_flag": 0})

        for sample, traits in samples.items():
            mk = ds_marker.get(sample)
            if mk is None or mk.shape[1] == 0:
                continue
            mk = mk.astype("float32")
            mk_spots = set(mk.columns)
            bg_mean = mk.mean(axis=1)

            for trait, ldsc_df in traits.items():
                if "spot" not in ldsc_df.columns or "p" not in ldsc_df.columns:
                    continue
                ld = ldsc_df[["spot", "p"]].copy()
                ld["spot"] = ld["spot"].astype(str)
                ld = ld[ld["spot"].isin(mk_spots)]
                if len(ld) < min_sig_spots:
                    continue

                sig = ld[ld["p"] < spot_sig_alpha]
                if len(sig) < min_sig_spots:
                    ld = ld.assign(mlog10p=-np.log10(ld["p"].clip(lower=1e-300)))
                    cut = ld["mlog10p"].quantile(0.90)
                    sig = ld[ld["mlog10p"] >= cut]
                if len(sig) < 3:
                    continue

                sig_spots = [s for s in sig["spot"].tolist() if s in mk_spots]
                if len(sig_spots) < 3:
                    continue

                fg_mean = mk[sig_spots].mean(axis=1)
                enrichment = fg_mean - bg_mean
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

        for (gene, trait), vals in acc.items():
            records.append({
                "dataset": dataset,
                "gene": gene,
                "trait": trait,
                "enrichment": max(vals["enrich"]),
                "min_pval": min(vals["min_pval"]),
                "n_samples_flagged": vals["n_flag"],
            })

    return pd.DataFrame(records)


def fisher_test(foreground: set[str], risk: set[str], universe: set[str]):
    from scipy.stats import fisher_exact

    fg = foreground & universe
    rk = risk & universe
    a = len(fg & rk)
    b = len(fg - rk)
    c = len(rk - fg)
    d = len(universe - fg - rk)
    if min(a + b, c + d, a + c, b + d) == 0:
        return np.nan, np.nan, (a, b, c, d)
    odds, pval = fisher_exact([[a, b], [c, d]])
    return odds, pval, (a, b, c, d)


def wald_ci(a: int, b: int, c: int, d: int, z: float = 1.959964):
    if min(a, b, c, d) == 0:
        return np.nan, np.nan
    log_or = math.log((a * d) / (b * c))
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return math.exp(log_or - z * se), math.exp(log_or + z * se)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Write outputs for available traits instead of failing on missing traits.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=GS,
        help="Output directory; defaults to Analysis/Spatial/results/gsmap.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    config_traits = read_simple_yaml_keys(GWAS_CONFIG)
    expected_traits = [t for t in config_traits if t in TRAIT_META]
    unknown = [t for t in config_traits if t not in TRAIT_META]
    if unknown:
        print(f"WARNING: traits in config without display metadata: {unknown}")
    if not expected_traits:
        raise SystemExit(f"No recognized gsMap traits found in {GWAS_CONFIG}")

    print("=" * 70)
    print("  15j: Prioritized spatial GWAS-risk genes across all gsMap traits")
    print("=" * 70)
    print(f"  expected traits: {len(expected_traits)}")

    print("\n  Loading gsMap spatial_ldsc results and marker scores...")
    gsmap_results = load_gsmap_results()
    marker_scores = load_gsmap_marker_scores()
    loc = gene_risk_localization(gsmap_results, marker_scores)

    observed_pairs = set(zip(loc["dataset"], loc["trait"]))
    expected_pairs = {
        (cohort, trait)
        for cohort in COHORT_META
        for trait in expected_traits
    }
    missing_pairs = sorted(expected_pairs - observed_pairs)
    if missing_pairs and not args.allow_missing:
        print("\nERROR: Missing gsMap risk-localization outputs:")
        for cohort, trait in missing_pairs:
            print(f"  {cohort:10} {trait}")
        print("\nRun 15e_run_gsmap.sh with the 18-trait gwas_config.yaml, then rerun 15j.")
        return 2

    if missing_pairs:
        print("\nWARNING: writing partial outputs; missing pairs:")
        for cohort, trait in missing_pairs:
            print(f"  {cohort:10} {trait}")

    out_loc = args.output_dir / "gene_risk_localization_18trait.csv"
    loc.to_csv(out_loc, index=False)
    print(f"\n  wrote {out_loc} ({len(loc):,} rows)")

    prioritized = read_gene_set(UDIR / "prioritized_universe_FINAL.txt")
    universe = {}
    for cohort, sample_map in marker_scores.items():
        genes = set()
        for marker in sample_map.values():
            genes.update(str(g) for g in marker.index)
        if genes:
            universe[cohort] = genes

    risk = {
        (cohort, trait): set(sub["gene"])
        for (cohort, trait), sub in loc.groupby(["dataset", "trait"])
    }

    rows = []
    for trait in expected_traits:
        trait_label, family, trait_order = TRAIT_META[trait]
        for cohort, (cohort_label, cohort_order) in COHORT_META.items():
            if cohort not in universe or (cohort, trait) not in risk:
                if args.allow_missing:
                    continue
                raise SystemExit(f"Missing universe or risk set for {cohort}/{trait}")

            U = universe[cohort]
            R = risk[(cohort, trait)] & U
            OR, pval, (a, b, c, d) = fisher_test(prioritized, R, U)
            lo, hi = wald_ci(a, b, c, d)
            n_prior_universe = len(prioritized & U)
            rows.append({
                "trait": trait,
                "trait_label": trait_label,
                "trait_family": family,
                "trait_order": trait_order,
                "cohort": cohort,
                "cohort_label": cohort_label,
                "cohort_order": cohort_order,
                "n_risk_genes": len(R),
                "n_prior_in_risk": a,
                "n_prior_total": len(prioritized),
                "n_prior_in_universe": n_prior_universe,
                "n_universe": len(U),
                "frac_risk_prioritized": a / len(R) if R else np.nan,
                "frac_prioritized_recovered": a / n_prior_universe if n_prior_universe else np.nan,
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "fisher_or": OR,
                "fisher_pval": pval,
                "ci_low": lo,
                "ci_high": hi,
                "neglog10p": -math.log10(max(pval, 1e-300)) if pval == pval else np.nan,
            })

    summary = pd.DataFrame(rows).sort_values(["trait_order", "cohort_order"])
    out_summary = args.output_dir / "prioritized_spatial_risk_18trait.csv"
    summary.to_csv(out_summary, index=False)
    print(f"  wrote {out_summary} ({len(summary):,} rows)")

    print("\n=== prioritized spatial GWAS-risk genes ===")
    cols = [
        "trait_label", "cohort_label", "n_risk_genes", "n_prior_in_risk",
        "frac_prioritized_recovered", "fisher_or", "fisher_pval",
    ]
    with pd.option_context("display.max_rows", 80, "display.width", 140):
        print(summary[cols].to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
