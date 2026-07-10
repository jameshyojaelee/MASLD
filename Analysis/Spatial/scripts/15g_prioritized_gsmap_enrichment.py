#!/usr/bin/env python
"""Recompute gsMap spatial-risk enrichment with set_a = prioritized MASLD universe
(Fig4 target list), replicating 15f_gsmap_analysis.py `conserved_enrichment` exactly
minus the conserved-core -> prioritized swap. Validates the pipeline by first
reproducing the on-disk conserved-core fisher_or, then emits the prioritized version.
"""
import sys, math
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow.feather as pf
from scipy.stats import fisher_exact

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
sys.path.insert(0, str(BASE / "Analysis/Spatial/scripts"))
from spatial_utils import load_conserved  # noqa: E402  (identical set as 15f)

GSMAP = BASE / "Analysis/Spatial/results/gsmap"
PRIOR_TXT = BASE / "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt"

LIVER = ["ukbb_alt", "ukbb_ast", "ukbb_ggt", "finngen_nafld"]  # headline replicating traits

def fisher_test(set_a, set_b, universe):
    """Exact replica of spatial_utils.fisher_test (sample OR from scipy)."""
    a = len(set_a & set_b)
    b = len(set_a - set_b)
    c = len(set_b - set_a)
    d = len(universe - set_a - set_b)
    if min(a + b, c + d, a + c, b + d) == 0:
        return np.nan, np.nan, 0, (a, b, c, d)
    odds, pval = fisher_exact([[a, b], [c, d]])
    return odds, pval, a, (a, b, c, d)

def wald_ci(a, b, c, d, z=1.959964):
    """Wald CI on log(sample OR); all cells large here so this is well-behaved."""
    if min(a, b, c, d) == 0:
        return np.nan, np.nan
    logor = math.log((a * d) / (b * c))
    se = math.sqrt(1/a + 1/b + 1/c + 1/d)
    return math.exp(logor - z * se), math.exp(logor + z * se)

# ── universe per dataset: union of HUMAN_GENE_SYM across sample feathers ──
def build_universe():
    uni = {}
    for ds_dir in sorted(GSMAP.glob("*")):
        if not ds_dir.is_dir() or ds_dir.name in ("figures",):
            continue
        genes = set()
        for s_dir in sorted(ds_dir.glob("*")):
            fe = s_dir / "latent_to_gene" / f"{s_dir.name}_gene_marker_score.feather"
            if fe.exists():
                col = pf.read_table(fe, columns=["HUMAN_GENE_SYM"]).column("HUMAN_GENE_SYM").to_pylist()
                genes.update(col)
        if genes:
            uni[ds_dir.name] = genes
    return uni

universe = build_universe()
for d, g in universe.items():
    print(f"[universe] {d}: {len(g)} genes")

# ── sets ──
prioritized = set(l.strip() for l in open(PRIOR_TXT) if l.strip())
conserved = set(load_conserved())
print(f"[sets] prioritized={len(prioritized)}  conserved(validation)={len(conserved)}")

# ── risk genes per (dataset, trait) from gene_risk_localization.csv ──
grl = pd.read_csv(GSMAP / "gene_risk_localization.csv")
ondisk = pd.read_csv(GSMAP / "conserved_spatial_risk.csv")
ondisk_map = {(r.dataset, r.trait): (r.fisher_or, r.fisher_pval) for r in ondisk.itertuples()}

rows = []
print("\n=== VALIDATION: reproduce on-disk conserved-core fisher_or ===")
print(f"{'dataset':10} {'trait':15} {'repro_OR':>9} {'ondisk_OR':>9} {'match':>6}")
for (ds, trait), sub in grl.groupby(["dataset", "trait"]):
    uni = universe.get(ds)
    if not uni:
        continue
    risk = set(sub["gene"])
    # validation: conserved
    cor, cp, ca, _ = fisher_test(conserved, risk, uni)
    od = ondisk_map.get((ds, trait), (np.nan, np.nan))
    match = "OK" if (not np.isnan(cor) and not np.isnan(od[0]) and abs(cor - od[0]) < 1e-6) else "DIFF"
    print(f"{ds:10} {trait:15} {cor:9.4f} {od[0]:9.4f} {match:>6}")
    # target: prioritized
    por, pp, pa_, (a, b, c, d) = fisher_test(prioritized, risk, uni)
    lo, hi = wald_ci(a, b, c, d)
    rows.append(dict(dataset=ds, trait=trait, n_risk_genes=len(risk),
                     n_prior_in_risk=a, n_prior_total=len(prioritized),
                     n_prior_in_universe=len(prioritized & uni), n_universe=len(uni),
                     a=a, b=b, c=c, d=d,
                     fisher_or=por, fisher_pval=pp, ci_low=lo, ci_high=hi))

res = pd.DataFrame(rows)
res["neglog10p"] = -np.log10(res["fisher_pval"].clip(lower=1e-300))
res["is_liver"] = res["trait"].isin(LIVER)

print("\n=== PRIORITIZED-TARGET enrichment (set_a = prioritized_universe_FINAL) ===")
print(f"{'dataset':10} {'trait':15} {'OR':>6} {'CI_low':>7} {'CI_high':>7} {'pval':>11} {'a':>5} {'liver':>6}")
for r in res.sort_values(["is_liver", "trait", "dataset"], ascending=[False, True, True]).itertuples():
    print(f"{r.dataset:10} {r.trait:15} {r.fisher_or:6.2f} {r.ci_low:7.2f} {r.ci_high:7.2f} "
          f"{r.fisher_pval:11.2e} {r.a:5d} {'YES' if r.is_liver else '':>6}")

# ── replication verdict for liver traits ──
print("\n=== REPLICATION VERDICT (liver-enzyme + FinnGen-NAFLD, both cohorts) ===")
for trait in LIVER:
    sub = res[res.trait == trait]
    g = sub[sub.dataset == "gse192741"]
    v = sub[sub.dataset == "vu"]
    if len(g) and len(v):
        gr, gp = g.fisher_or.iloc[0], g.fisher_pval.iloc[0]
        vr, vp = v.fisher_or.iloc[0], v.fisher_pval.iloc[0]
        both = (gr > 1 and gp < 0.05 and vr > 1 and vp < 0.05)
        print(f"  {trait:15} gse192741 OR={gr:.2f} p={gp:.1e} | vu OR={vr:.2f} p={vp:.1e}  "
              f"-> {'REPLICATES' if both else 'DOES NOT REPLICATE'}")

out = GSMAP / "prioritized_spatial_risk.csv"
res.to_csv(out, index=False)
print(f"\n[saved] {out}")
