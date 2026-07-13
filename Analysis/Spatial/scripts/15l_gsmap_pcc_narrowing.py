#!/usr/bin/env python
"""15l — gsMap PCC gene-narrowing table for Fig 4f (spatial-narrowing panel).

gsMap's native per-gene "narrow risk genes via spatial context" statistic is
PCC (gsMap/diagnosis.py): the Pearson correlation, across spots, between a gene's
spatial specificity (GSS) and the trait's spatial heritability (-log10 P of
spatial LDSC). High PCC = the gene's tissue-expression pattern tracks where the
GWAS risk concentrates. gsMap writes it per sample x trait as
  report/{trait}/{sample}_{trait}_Gene_Diagnostic_Info.csv  [Gene, Annotation, Median_GSS, PCC]

This script aggregates PCC per (gene, trait, cohort) across samples and attaches
the frozen manuscript evidence classes. The primary classes are defined only by
canonical TREAT DE and Tier-1/2 SuSiE-COLOC, so the spatial analysis is not part
of its own gene-set definition. Broad historical universe flags are retained as
explicit sensitivity columns only.

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
RELEASE_ID = os.environ.get("MANUSCRIPT_RELEASE_ID", "2026-07-10-r1")
CLASS_FILE = (BASE / "RNA-seq/results/manuscript_release" / RELEASE_ID /
              "evidence_class_table.tsv")

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

    prioritized_broad = rd("prioritized_universe_FINAL.txt")
    genetic_broad = rd("universe_genetic.txt")
    if not CLASS_FILE.exists():
        raise FileNotFoundError(
            f"Missing frozen evidence classes: {CLASS_FILE}. "
            "Run scripts/manuscript/build_evidence_class_release.R first."
        )
    classes = pd.read_csv(CLASS_FILE, sep="\t", low_memory=False)
    classes = classes[[
        "symbol", "analysis_release_id", "primary_evidence_class",
        "sensitivity_evidence_class", "genetic_trait_scope",
        "genetic_confidence", "bulk_AveExpr", "gene_biotype",
    ]].drop_duplicates("symbol")
    print(
        f"  release={RELEASE_ID}  primary SuSiE={sum(classes.genetic_confidence == 'susie'):,}  "
        f"convergent={sum(classes.primary_evidence_class == 'convergent'):,}"
    )
    print(
        f"  broad sensitivity: prioritized={len(prioritized_broad):,}  "
        f"genetic={len(genetic_broad):,}"
    )

    files = sorted(glob.glob(DIAG_GLOB))
    print(f"  Gene_Diagnostic_Info files: {len(files)} (expect 270 = 15 samples x 18 traits)")

    frames = []
    for f in files:
        m = PATH_RE.search(f)
        if not m:
            print(f"  WARN: unparsed path {f}")
            continue
        df = pd.read_csv(f, usecols=["Gene", "Median_GSS", "PCC"])
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
                     median_gss_mean=("Median_GSS", "mean"),
                     n_samples=("PCC", "size"))
                .reset_index()
                .rename(columns={"Gene": "gene"}))
    agg["cohort_label"] = agg["cohort"].map(COHORT_LABEL)
    agg = agg.merge(classes, left_on="gene", right_on="symbol", how="left")
    agg["analysis_release_id"] = agg["analysis_release_id"].fillna(RELEASE_ID)
    agg["is_susie_genetic"] = agg["primary_evidence_class"].isin(
        ["genetic_only", "convergent"]
    )
    agg["is_treat_disease_state"] = agg["primary_evidence_class"].isin(
        ["disease_state_only", "convergent"]
    )
    agg["is_convergent"] = agg["primary_evidence_class"].eq("convergent")
    agg["is_genetic_only"] = agg["primary_evidence_class"].eq("genetic_only")
    agg["is_disease_state_only"] = agg["primary_evidence_class"].eq("disease_state_only")
    agg["is_prioritized_broad_sensitivity"] = agg["gene"].isin(prioritized_broad)
    agg["is_genetic_broad_sensitivity"] = agg["gene"].isin(genetic_broad)
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
            gen = sub.loc[sub.is_susie_genetic, "pcc_mean"]
            bg = sub.loc[~sub.is_susie_genetic, "pcc_mean"]
            if len(gen) < 5 or len(bg) < 5:
                continue
            p = mannwhitneyu(gen, bg, alternative="greater").pvalue
            top15 = sub.nsmallest(15, "rank")
            n_gen_top = int(top15.is_susie_genetic.sum())
            print(f"{trait:16} {cohort:10} {gen.median():11.3f} {bg.median():8.3f} "
                  f"{p:10.2e} {n_gen_top:6d}/15")

    # ── show top narrowed genes for the likely primary trait ──
    print("\n=== top-15 PCC-narrowed genes: mvp_nafld / GSE192741 ===")
    t = agg[(agg.trait == "mvp_nafld") & (agg.cohort == "gse192741")].nsmallest(15, "rank")
    for r in t.itertuples():
        flags = []
        if r.is_convergent: flags.append("convergent")
        elif r.is_genetic_only: flags.append("genetic-only")
        elif r.is_disease_state_only: flags.append("disease-state-only")
        print(f"  {int(r.rank):2d}. {r.gene:12} PCC={r.pcc_mean:.3f} "
              f"(n={r.n_samples}) {', '.join(flags)}")


if __name__ == "__main__":
    main()
