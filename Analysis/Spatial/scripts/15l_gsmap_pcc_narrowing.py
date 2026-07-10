#!/usr/bin/env python
"""15l — gsMap PCC gene-narrowing table for Fig 4f (spatial-narrowing panel).

gsMap's native per-gene "narrow risk genes via spatial context" statistic is
PCC (gsMap/diagnosis.py): the Pearson correlation, across spots, between a gene's
spatial specificity (GSS) and the trait's spatial heritability (-log10 P of
spatial LDSC). High PCC = the gene's tissue-expression pattern tracks where the
GWAS risk concentrates. gsMap writes it per sample x trait as
  report/{trait}/{sample}_{trait}_Gene_Diagnostic_Info.csv  [Gene, Annotation, Median_GSS, PCC]

This script aggregates PCC per (gene, trait, cohort) across the cohort's samples
and flags which genes are in the Fig4 prioritized set and the genetic/COLOC set,
so the figure can show that spatial concordance surfaces the causal genes.

NOTE: the run used --annotation condition (single-valued per sample), so the
`Annotation`/`Median_GSS` fields are global and uninformative — we use PCC only,
which is annotation-independent.

Output (isolated; canonical untouched):
  Analysis/Spatial/results/gsmap/gsmap_pcc_by_trait.csv

Env: rnaseq (or spatial). pandas/numpy/scipy.
"""
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
pd.set_option("compute.use_numexpr", False)

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GS = BASE / "Analysis/Spatial/results/gsmap"
UDIR = BASE / "Analysis/Spatial/results/universe_validation"

COHORT_LABEL = {"gse192741": "GSE192741", "vu": "Vu et al. 2025"}
DIAG_GLOB = str(GS / "*/*/report/*/*_Gene_Diagnostic_Info.csv")
# path: .../gsmap/{cohort}/{sample}/report/{trait}/{sample}_{trait}_Gene_Diagnostic_Info.csv
PATH_RE = re.compile(r"/gsmap/(?P<cohort>[^/]+)/(?P<sample>[^/]+)/report/(?P<trait>[^/]+)/")


def rd(fname):
    return set(l.strip() for l in open(UDIR / fname) if l.strip())


def main():
    print("=" * 70)
    print("  15l: gsMap PCC gene-narrowing aggregation (Fig 4f)")
    print("=" * 70)

    prioritized = rd("prioritized_universe_FINAL.txt")
    genetic = rd("universe_genetic.txt")   # COLOC/genetic arm (same set 15h2/15k use)
    print(f"  prioritized={len(prioritized):,}  genetic/COLOC={len(genetic):,}")

    files = sorted(glob.glob(DIAG_GLOB))
    print(f"  Gene_Diagnostic_Info files: {len(files)} (expect 270 = 15 samples x 18 traits)")

    frames = []
    for f in files:
        m = PATH_RE.search(f)
        if not m:
            print(f"  WARN: unparsed path {f}")
            continue
        df = pd.read_csv(f, usecols=["Gene", "PCC"])
        df["cohort"] = m.group("cohort")
        df["trait"] = m.group("trait")
        df["sample"] = m.group("sample")
        frames.append(df)
    allg = pd.concat(frames, ignore_index=True)
    print(f"  rows loaded: {len(allg):,}")

    # ── aggregate per (gene, trait, cohort) across samples ──
    agg = (allg.groupby(["trait", "cohort", "Gene"])
                .agg(pcc_mean=("PCC", "mean"),
                     pcc_median=("PCC", "median"),
                     n_samples=("PCC", "size"))
                .reset_index()
                .rename(columns={"Gene": "gene"}))
    agg["cohort_label"] = agg["cohort"].map(COHORT_LABEL)
    agg["is_prioritized"] = agg["gene"].isin(prioritized)
    agg["is_genetic"] = agg["gene"].isin(genetic)          # COLOC/GWAS-genetic
    agg["is_convergent"] = agg["is_prioritized"] & agg["is_genetic"]
    agg["rank"] = agg.groupby(["trait", "cohort"])["pcc_mean"].rank(ascending=False, method="min")
    agg = agg.sort_values(["trait", "cohort", "rank"]).reset_index(drop=True)

    out = GS / "gsmap_pcc_by_trait.csv"
    agg.to_csv(out, index=False)
    print(f"\n  wrote {out.name}  ({len(agg):,} rows, "
          f"{agg['trait'].nunique()} traits x {agg['cohort'].nunique()} cohorts)")

    # ── VERIFICATION: does high PCC surface genetic/prioritized genes? ──
    from scipy.stats import mannwhitneyu
    print("\n=== does spatial concordance (PCC) discriminate genetic-risk genes? ===")
    print(f"{'trait':16} {'cohort':10} {'PCC_genetic':>11} {'PCC_bg':>8} {'MWU_p':>10} {'top15_genetic':>13}")
    focus = ["mvp_nafld", "finngen_nafld", "ghodsian_nafld", "ukbb_alt", "mvp_alt", "pdff"]
    for trait in focus:
        for cohort in ["gse192741", "vu"]:
            sub = agg[(agg.trait == trait) & (agg.cohort == cohort)]
            if sub.empty:
                continue
            gen = sub.loc[sub.is_genetic, "pcc_mean"]
            bg = sub.loc[~sub.is_genetic, "pcc_mean"]
            if len(gen) < 5 or len(bg) < 5:
                continue
            p = mannwhitneyu(gen, bg, alternative="greater").pvalue
            top15 = sub.nsmallest(15, "rank")
            n_gen_top = int(top15.is_genetic.sum())
            print(f"{trait:16} {cohort:10} {gen.median():11.3f} {bg.median():8.3f} "
                  f"{p:10.2e} {n_gen_top:6d}/15")

    # ── show top narrowed genes for the likely primary trait ──
    print("\n=== top-15 PCC-narrowed genes: mvp_nafld / GSE192741 ===")
    t = agg[(agg.trait == "mvp_nafld") & (agg.cohort == "gse192741")].nsmallest(15, "rank")
    for r in t.itertuples():
        flags = []
        if r.is_prioritized: flags.append("prioritized")
        if r.is_genetic: flags.append("COLOC")
        print(f"  {int(r.rank):2d}. {r.gene:12} PCC={r.pcc_mean:.3f} "
              f"(n={r.n_samples}) {', '.join(flags)}")


if __name__ == "__main__":
    main()
