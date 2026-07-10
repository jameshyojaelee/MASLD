#!/usr/bin/env python
"""15h2 — EUR-expansion gsMap spatial-risk enrichment (11 new traits).

Companion to 15h_arm_decomposition_gsmap.py, run AFTER 15e computes spatial_ldsc
for the 11 added EUR traits. Reuses 15f's exact localization functions (imported,
not copied) so the risk-gene definition is identical to the canonical 7-trait run.

Steps:
  1. Regenerate gene-level risk localization across ALL traits now on disk
     (7 existing + 11 new) via 15f.gene_risk_localization. Written to an ISOLATED
     file (gene_risk_localization_18trait.csv) — the canonical
     gene_risk_localization.csv (consumed by fig4g) is left untouched.
  2. Validate: the 7 existing traits' per-(dataset,trait) risk-gene counts must
     reproduce the on-disk canonical file (deterministic check).
  3. Fisher enrichment (identical to 15h: sample OR, Wald CI) for the 11 NEW
     traits x {gse192741, vu} x 4 partitions {Convergent, Prioritized union,
     Transcriptomic-only, Genetic-only}. Written to eur_expansion_spatial_risk.csv.

Env: rnaseq (or spatial). Report-only; does NOT touch fig4g_gsmap_risk_in_tissue.R.
"""
import importlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
pd.set_option("compute.use_numexpr", False)  # silence numexpr/numpy ABI warning

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SCRIPTS = BASE / "Analysis/Spatial/scripts"
GS = BASE / "Analysis/Spatial/results/gsmap"
UDIR = BASE / "Analysis/Spatial/results/universe_validation"
sys.path.insert(0, str(SCRIPTS))
g15f = importlib.import_module("15f_gsmap_analysis")

NEW_TRAITS = [
    "mvp_nafld", "mvp_alt", "mvp_ast", "decode_nafld", "ukbb2023_nafld",
    "intermtn_nafld", "anstee2020_nafld", "nafld_2019",
    "pdff_2021a", "pdff_2021b", "pdff_2022",
]
EXISTING = ["finngen_nafld", "finngen_nash", "ghodsian_nafld", "pdff",
            "ukbb_alt", "ukbb_ast", "ukbb_ggt"]


def rd(f):
    return set(l.strip() for l in open(UDIR / f) if l.strip())


def fisher_test(set_a, set_b, universe):
    from scipy.stats import fisher_exact
    a = len(set_a & set_b); b = len(set_a - set_b)
    c = len(set_b - set_a); d = len(universe - set_a - set_b)
    if min(a + b, c + d, a + c, b + d) == 0:
        return np.nan, np.nan, (a, b, c, d)
    odds, pval = fisher_exact([[a, b], [c, d]])
    return odds, pval, (a, b, c, d)


def wald_ci(a, b, c, d, z=1.959964):
    if min(a, b, c, d) == 0:
        return np.nan, np.nan
    lo = math.log((a * d) / (b * c)); se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return math.exp(lo - z * se), math.exp(lo + z * se)


def main():
    print("=" * 70)
    print("  15h2: EUR-expansion gsMap spatial-risk enrichment")
    print("=" * 70)

    # ── 1. regenerate localization across all traits on disk ──
    print("\n  Loading spatial_ldsc results + marker scores ...")
    res = g15f.load_gsmap_results()
    mk = g15f.load_gsmap_marker_scores()
    loc = g15f.gene_risk_localization(res, mk)
    out_loc = GS / "gene_risk_localization_18trait.csv"
    loc.to_csv(out_loc, index=False)
    traits_found = sorted(loc["trait"].unique())
    print(f"  regenerated localization: {len(loc):,} rows, "
          f"{len(traits_found)} traits -> {out_loc.name}")
    print(f"  traits: {traits_found}")

    # ── 2. validate reproduction of the 7 existing traits ──
    canon = pd.read_csv(GS / "gene_risk_localization.csv")
    print("\n  Validation (existing 7 traits reproduce canonical counts):")
    ok_all = True
    for (ds, tr), sub in canon.groupby(["dataset", "trait"]):
        n_canon = len(sub)
        n_new = len(loc[(loc.dataset == ds) & (loc.trait == tr)])
        flag = "OK" if n_canon == n_new else "DIFF"
        if n_canon != n_new:
            ok_all = False
        print(f"    {ds:10} {tr:15} canon={n_canon:5d} regen={n_new:5d} {flag}")
    print(f"  reproduction: {'ALL MATCH' if ok_all else 'MISMATCH — investigate'}")

    # ── 3. partitions + per-dataset universe (identical to 15h) ──
    rna_sig = rd("universe_transcriptomic.txt")
    coloc = rd("universe_genetic.txt")
    FINAL = rd("prioritized_universe_FINAL.txt")
    GROUPS = [
        ("Convergent",          rna_sig & coloc),
        ("Prioritized union",   FINAL),
        ("Transcriptomic-only", rna_sig - coloc),
        ("Genetic-only",        coloc - rna_sig),
    ]
    print(f"\n  partitions: Convergent={len(rna_sig & coloc)} "
          f"FINAL={len(FINAL)} T-only={len(rna_sig - coloc)} "
          f"G-only={len(coloc - rna_sig)}")

    import pyarrow.feather as pf
    universe = {}
    for ds_dir in sorted(GS.glob("*")):
        if not ds_dir.is_dir() or ds_dir.name in ("figures",):
            continue
        genes = set()
        for s_dir in sorted(ds_dir.glob("*")):
            fe = s_dir / "latent_to_gene" / f"{s_dir.name}_gene_marker_score.feather"
            if fe.exists():
                genes.update(pf.read_table(
                    fe, columns=["HUMAN_GENE_SYM"]
                ).column("HUMAN_GENE_SYM").to_pylist())
        if genes:
            universe[ds_dir.name] = genes

    risk = {(ds, tr): set(sub["gene"])
            for (ds, tr), sub in loc.groupby(["dataset", "trait"])}

    # ── 4. Fisher enrichment for the 11 new traits ──
    rows = []
    for trait in NEW_TRAITS:
        for ds in ["gse192741", "vu"]:
            if (ds, trait) not in risk:
                print(f"    WARNING: no risk set for {ds}/{trait} (spatial_ldsc missing?)")
                continue
            for label, S in GROUPS:
                OR, p, (a, b, c, d) = fisher_test(S, risk[(ds, trait)], universe[ds])
                lo, hi = wald_ci(a, b, c, d)
                rows.append(dict(trait=trait, cohort=ds, partition=label,
                                 n_set=len(S), a=a, fisher_or=OR,
                                 ci_low=lo, ci_high=hi, fisher_pval=p,
                                 neglog10p=-math.log10(max(p, 1e-300))
                                 if p == p else np.nan))

    res_df = pd.DataFrame(rows)
    out = GS / "eur_expansion_spatial_risk.csv"
    res_df.to_csv(out, index=False)
    print(f"\n  saved {out}  ({len(res_df)} rows)")

    # ── compact report table ──
    print("\n=== EUR-EXPANSION spatial-risk enrichment (OR [95% CI], p) ===")
    print(f"{'trait':16} {'cohort':10} {'partition':20} {'OR':>6} "
          f"{'CI_low':>7} {'CI_high':>7} {'pval':>11} {'a':>5}")
    for r in res_df.itertuples():
        orr = f"{r.fisher_or:.2f}" if r.fisher_or == r.fisher_or else "NA"
        lo = f"{r.ci_low:.2f}" if r.ci_low == r.ci_low else "NA"
        hi = f"{r.ci_high:.2f}" if r.ci_high == r.ci_high else "NA"
        print(f"{r.trait:16} {r.cohort:10} {r.partition:20} {orr:>6} "
              f"{lo:>7} {hi:>7} {r.fisher_pval:11.2e} {r.a:5d}")

    # ── Convergent-partition replication verdict (both cohorts) ──
    print("\n=== CONVERGENT-partition replication (both cohorts OR>1, p<0.05) ===")
    conv = res_df[res_df.partition == "Convergent"]
    for trait in NEW_TRAITS:
        g = conv[(conv.trait == trait) & (conv.cohort == "gse192741")]
        v = conv[(conv.trait == trait) & (conv.cohort == "vu")]
        if len(g) and len(v):
            gr, gp = g.fisher_or.iloc[0], g.fisher_pval.iloc[0]
            vr, vp = v.fisher_or.iloc[0], v.fisher_pval.iloc[0]
            both = (gr > 1 and gp < 0.05 and vr > 1 and vp < 0.05)
            print(f"  {trait:16} gse192741 OR={gr:.2f} p={gp:.1e} | "
                  f"vu OR={vr:.2f} p={vp:.1e} -> "
                  f"{'REPLICATES' if both else 'no'}")


if __name__ == "__main__":
    main()
