#!/usr/bin/env python
"""15k — Phenotype-level spatial genetic-risk enrichment for Fig 4f.

Clumps the 18 gsMap-compatible EUR GWAS into 6 liver phenotypes (ALT, AST, GGT,
NAFLD, NASH, PDFF) and summarizes, per phenotype x spatial cohort, how strongly
the Fig4 prioritized MASLD gene set is enriched among spatial GWAS-risk genes.

Design (approved 2026-07-09):
  * Per-study Fisher enrichment ORs already exist in prioritized_spatial_risk_18trait.csv
    (one OR per trait x cohort). We do NOT meta-analyze them into a pooled CI —
    studies within a phenotype share samples (UKB / MVP recur), so a parametric
    pooled CI would fake precision.
  * Each phenotype x cohort is summarized by the MEDIAN per-study OR; the whisker
    is the across-study RANGE (min-max). This is #study-invariant and shows
    consistency, which is the panel's message.
  * Single-study phenotypes (GGT, NASH) show that study's OR + parametric 95% CI
    (flagged n_studies == 1).
  * Driver genes: prioritized genes recovered as spatial-risk genes in the most
    (trait, cohort) pairs of a phenotype (cross-study robust, both cohorts).

Inputs (produced by 15j / 15h2 after 15e):
  Analysis/Spatial/results/gsmap/prioritized_spatial_risk_18trait.csv
  Analysis/Spatial/results/gsmap/gene_risk_localization_18trait.csv
  Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt

Outputs (isolated; canonical files untouched):
  Analysis/Spatial/results/gsmap/phenotype_spatial_risk.csv    (6 pheno x 2 cohort)
  Analysis/Spatial/results/gsmap/phenotype_driver_genes.csv

Env: rnaseq (or spatial). pandas/numpy only.
"""
from pathlib import Path

import numpy as np
import pandas as pd
pd.set_option("compute.use_numexpr", False)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GS = BASE / "Analysis/Spatial/results/gsmap"
UDIR = BASE / "Analysis/Spatial/results/universe_validation"

# ── phenotype clumping (trait -> phenotype) ──────────────────────────────────
PHENO = {
    "ukbb_alt": "ALT", "mvp_alt": "ALT",
    "ukbb_ast": "AST", "mvp_ast": "AST",
    "ukbb_ggt": "GGT",
    "finngen_nafld": "NAFLD", "ghodsian_nafld": "NAFLD", "mvp_nafld": "NAFLD",
    "decode_nafld": "NAFLD", "ukbb2023_nafld": "NAFLD", "intermtn_nafld": "NAFLD",
    "anstee2020_nafld": "NAFLD", "nafld_2019": "NAFLD",
    "finngen_nash": "NASH",
    "pdff": "PDFF", "pdff_2021a": "PDFF", "pdff_2021b": "PDFF", "pdff_2022": "PDFF",
}
# display order + family
PHENO_ORDER = ["ALT", "AST", "GGT", "NAFLD", "NASH", "PDFF"]
PHENO_FAMILY = {"ALT": "Liver enzymes", "AST": "Liver enzymes", "GGT": "Liver enzymes",
                "NAFLD": "Disease diagnosis", "NASH": "Disease diagnosis",
                "PDFF": "Imaging liver fat"}
COHORT_LABEL = {"gse192741": "GSE192741", "vu": "Vu et al. 2025"}
N_DRIVERS = 2


def main():
    print("=" * 70)
    print("  15k: phenotype-level spatial genetic-risk enrichment (Fig 4f)")
    print("=" * 70)

    prior = set(l.strip() for l in open(UDIR / "prioritized_universe_FINAL.txt") if l.strip())
    print(f"  prioritized universe: {len(prior):,} genes")

    # ── 1. per-study ORs -> phenotype summary (median + range) ──
    per = pd.read_csv(GS / "prioritized_spatial_risk_18trait.csv")
    per["phenotype"] = per["trait"].map(PHENO)
    unmapped = per[per["phenotype"].isna()]["trait"].unique()
    assert len(unmapped) == 0, f"unmapped traits: {unmapped}"

    rows = []
    for (pheno, cohort), sub in per.groupby(["phenotype", "cohort"]):
        ors = sub["fisher_or"].to_numpy(dtype=float)
        n = len(sub)
        med = float(np.median(ors))
        if n == 1:
            lo = float(sub["ci_low"].iloc[0])
            hi = float(sub["ci_high"].iloc[0])
            itype = "ci"
        else:
            lo, hi = float(np.min(ors)), float(np.max(ors))
            itype = "range"
        rows.append(dict(
            phenotype=pheno, cohort=cohort, cohort_label=COHORT_LABEL[cohort],
            family=PHENO_FAMILY[pheno], n_studies=n,
            or_summary=med, or_lo=lo, or_hi=hi, interval_type=itype,
            n_below1=int((ors < 1).sum()),
            or_min=float(np.min(ors)), or_max=float(np.max(ors)),
            median_prior_recovered=int(np.median(sub["n_prior_in_risk"])),
        ))
    ph = pd.DataFrame(rows)
    ph["phenotype"] = pd.Categorical(ph["phenotype"], categories=PHENO_ORDER, ordered=True)
    ph = ph.sort_values(["phenotype", "cohort_label"]).reset_index(drop=True)
    ph.to_csv(GS / "phenotype_spatial_risk.csv", index=False)
    print(f"\n  wrote phenotype_spatial_risk.csv ({len(ph)} rows)")

    # ── 2. driver genes: prioritized genes recovered across the most studies ──
    loc = pd.read_csv(GS / "gene_risk_localization_18trait.csv")
    loc = loc[loc["gene"].isin(prior)].copy()
    loc["phenotype"] = loc["trait"].map(PHENO)
    loc = loc[loc["phenotype"].notna()]

    drv_rows = []
    for pheno in PHENO_ORDER:
        sub = loc[loc["phenotype"] == pheno]
        # support = # distinct (trait, cohort) pairs where gene is a spatial-risk gene
        support = (sub.groupby("gene")
                      .agg(n_pairs=("trait", lambda s: len(s)),
                           n_cohorts=("dataset", lambda s: s.nunique()),
                           n_traits=("trait", lambda s: s.nunique()),
                           mean_enr=("enrichment", "mean"))
                      .reset_index())
        # robust drivers: present in BOTH cohorts, then most (trait,cohort) support
        support = support.sort_values(
            ["n_cohorts", "n_pairs", "mean_enr"], ascending=[False, False, False])
        top = support.head(N_DRIVERS)
        for _, r in top.iterrows():
            drv_rows.append(dict(phenotype=pheno, gene=r["gene"],
                                 n_pairs=int(r["n_pairs"]),
                                 n_cohorts=int(r["n_cohorts"]),
                                 n_traits=int(r["n_traits"]),
                                 mean_enrichment=round(float(r["mean_enr"]), 3)))
    drv = pd.DataFrame(drv_rows)
    drv["phenotype"] = pd.Categorical(drv["phenotype"], categories=PHENO_ORDER, ordered=True)
    drv = drv.sort_values(["phenotype", "n_pairs"], ascending=[True, False]).reset_index(drop=True)
    drv.to_csv(GS / "phenotype_driver_genes.csv", index=False)
    print(f"  wrote phenotype_driver_genes.csv ({len(drv)} rows)")

    # ── report ──
    print("\n=== phenotype x cohort enrichment (median OR [interval], n studies) ===")
    print(f"{'phenotype':10} {'cohort':14} {'nS':>3} {'medOR':>6} "
          f"{'lo':>6} {'hi':>6} {'type':>6} {'<1':>3}")
    for r in ph.itertuples():
        print(f"{str(r.phenotype):10} {r.cohort_label:14} {r.n_studies:3d} "
              f"{r.or_summary:6.2f} {r.or_lo:6.2f} {r.or_hi:6.2f} "
              f"{r.interval_type:>6} {r.n_below1:3d}")

    print("\n=== driver genes (top prioritized recovered per phenotype) ===")
    for pheno in PHENO_ORDER:
        g = drv[drv["phenotype"] == pheno]
        lab = ", ".join(f"{x.gene}(pairs={x.n_pairs},coh={x.n_cohorts})"
                        for x in g.itertuples())
        print(f"  {pheno:10} {lab}")


if __name__ == "__main__":
    main()
