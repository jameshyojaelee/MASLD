#!/usr/bin/env python
"""F1 follow-up: 1-Mb block bootstrap intervals for the F1 fractions, and the per-pair detail
tables the assignment asks for (arm plausibility, TOP-LD coverage strata, double-het donors).

The resampling unit is the 1-Mb block (chrN~Mb), fixed by IMPLEMENTATION_SPEC.md section 2.
Seed 20260914, 10,000 draws. Every statistic is recomputed end to end inside each draw from the
resampled rows, not from a stored per-pair summary.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1]).resolve()
TAB = OUT / "tables"
SEED = 20260914
NBOOT = 10_000

m = pd.read_csv(TAB / "pair_phase_substrate.tsv", sep="\t")
mp = m[~m.is_gnmt].copy()
mp["block_1mb"] = mp.chrom + "~" + (mp.pos1 // 1_000_000).astype(int).astype(str)

# The TOP-LD lookup keys on position only. Check that the TOP-LD entry at each matched position
# carries the same two alleles as the P5 SNV, so a co-located indel cannot supply a pair's LD.
inf = pd.read_csv(TAB / "topld_variant_info.tsv", sep="\t")
alle = inf.set_index(["chrom", "pos38"])[["topld_ref", "topld_alt"]].to_dict("index")


def allele_conc(c, p, r, a):
    e = alle.get((c, p))
    if e is None:
        return "absent"
    return "same_alleles" if {e["topld_ref"], e["topld_alt"]} == {r, a} else "different_alleles"


mp["v1_topld_allele_state"] = [allele_conc(r.chrom, r.pos1, r.ref1, r.alt1)
                               for r in mp.itertuples()]
mp["v2_topld_allele_state"] = [allele_conc(r.chrom, r.pos2, r.ref2, r.alt2)
                               for r in mp.itertuples()]
mp["topld_allele_concordant"] = ((mp.v1_topld_allele_state == "same_alleles")
                                 & (mp.v2_topld_allele_state == "same_alleles"))
mp[["signal_uid", "gene", "v1", "v2", "v1_topld_allele_state", "v2_topld_allele_state",
    "topld_allele_concordant", "topld_pair_state", "topld_dprime", "topld_r2",
    "panel_Dprime", "panel_r2"]].to_csv(TAB / "topld_allele_concordance.tsv",
                                        sep="\t", index=False)
print("TOP-LD allele concordance x coverage:")
print(pd.crosstab(mp.topld_allele_concordant, mp.topld_pair_state))
print("covered pairs with an allele-discordant TOP-LD entry: "
      f"{int(((mp.topld_pair_state == 'covered') & ~mp.topld_allele_concordant).sum())}")


def boot(df: pd.DataFrame, fn, label: str) -> dict:
    """Block bootstrap: resample 1-Mb blocks with replacement, recompute fn on the drawn rows."""
    rng = np.random.default_rng(SEED)
    blocks = df.block_1mb.to_numpy()
    uniq, inv = np.unique(blocks, return_inverse=True)
    idx_by_block = [np.flatnonzero(inv == i) for i in range(len(uniq))]
    obs = fn(df)
    draws = np.empty(NBOOT)
    for b in range(NBOOT):
        pick = rng.integers(0, len(uniq), len(uniq))
        rows = np.concatenate([idx_by_block[p] for p in pick])
        draws[b] = fn(df.iloc[rows])
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return dict(statistic=label, observed=float(obs), ci95_lo=float(lo), ci95_hi=float(hi),
                n_rows=int(len(df)), n_blocks_1mb=int(len(uniq)), n_boot=NBOOT, seed=SEED)


rows = []
rows.append(boot(mp, lambda d: (d.dprime_best > 0.9).mean(),
                 "share of P5 pairs with D' > 0.9 (best available source, missing counts as not)"))
rows.append(boot(mp[mp.topld_pair_state == "covered"], lambda d: (d.topld_dprime > 0.9).mean(),
                 "share with TOP-LD D' > 0.9, TOP-LD-covered pairs only"))
rows.append(boot(mp[mp.panel_state == "ok"], lambda d: (d.panel_Dprime > 0.9).mean(),
                 "share with 1000G-EUR-panel D' > 0.9, panel-computable pairs only"))
rows.append(boot(mp[mp.panel_state == "ok"], lambda d: (d.arm_min_freq < 0.01).mean(),
                 "share whose rarer single-variant arm has EUR haplotype frequency < 0.01"))
rows.append(boot(mp[mp.panel_state == "ok"], lambda d: (d.arm_min_freq < 0.005).mean(),
                 "share whose rarer single-variant arm has EUR haplotype frequency < 0.005"))
rows.append(boot(mp[mp.panel_state == "ok"],
                 lambda d: ((d.panel_n_alt1_ref2 == 0) | (d.panel_n_ref1_alt2 == 0)).mean(),
                 "share with a single-variant arm seen in 0 of 758 EUR haplotypes"))
rows.append(boot(mp[mp.panel_state == "ok"], lambda d: d.arm_min_freq.median(),
                 "median EUR haplotype frequency of the rarer single-variant arm"))
rows.append(boot(mp[mp.topld_pair_state == "covered"], lambda d: d.topld_r2.median(),
                 "median TOP-LD r2, covered pairs"))
ci = pd.DataFrame(rows)
ci.to_csv(TAB / "f1_block_bootstrap_intervals.tsv", sep="\t", index=False)
print(ci.to_string(index=False))

# ---- strata tables the assignment asks for, reported not summarised away
strat = []
for name, sub in [("all_P5_pairs", mp),
                  ("sep_le_500bp", mp[mp.separation_bp <= 500]),
                  ("sep_501_5000bp", mp[(mp.separation_bp > 500) & (mp.separation_bp <= 5000)]),
                  ("sep_gt_5000bp", mp[mp.separation_bp > 5000]),
                  ("universe_A_direct", mp[mp.universe == "A_direct"]),
                  ("universe_B_direct", mp[mp.universe == "B_direct"]),
                  ("universe_C_enzyme", mp[mp.universe == "C_enzyme"]),
                  ("p5_scored_only", mp[mp.p5_scored])]:
    ok = sub[sub.panel_state == "ok"]
    cov = sub[sub.topld_pair_state == "covered"]
    strat.append(dict(
        stratum=name, n_pairs=len(sub), n_blocks_1mb=sub.block_1mb.nunique(),
        n_topld_covered=len(cov), n_panel_ok=len(ok),
        frac_dprime_gt_09_all=(sub.dprime_best > 0.9).mean(),
        frac_dprime_gt_09_topld_covered=(cov.topld_dprime > 0.9).mean() if len(cov) else np.nan,
        median_topld_r2=cov.topld_r2.median() if len(cov) else np.nan,
        median_panel_r2=ok.panel_r2.median() if len(ok) else np.nan,
        frac_arm_min_lt_001=(ok.arm_min_freq < 0.01).mean() if len(ok) else np.nan,
        frac_arm_zero=((ok.panel_n_alt1_ref2 == 0) | (ok.panel_n_ref1_alt2 == 0)).mean()
        if len(ok) else np.nan,
        median_arm_min_freq=ok.arm_min_freq.median() if len(ok) else np.nan,
        frac_coupling_positive=(cov.topld_corr_sign == "+").mean() if len(cov) else np.nan,
    ))
st = pd.DataFrame(strat)
st.to_csv(TAB / "f1_strata.tsv", sep="\t", index=False)
print("\n" + st.to_string(index=False))

# ---- which of the two arms is the rare one
ok = mp[mp.panel_state == "ok"].copy()
ok["rarer_arm"] = np.where(ok.arm_v1_freq <= ok.arm_v2_freq, "v1_alone", "v2_alone")
print("\nrarer arm:\n" + ok.rarer_arm.value_counts().to_string())
print("\nboth arms < 0.01:", int(((ok.arm_v1_freq < 0.01) & (ok.arm_v2_freq < 0.01)).sum()),
      " exactly one arm < 0.01:",
      int(((ok.arm_v1_freq < 0.01) ^ (ok.arm_v2_freq < 0.01)).sum()),
      " neither arm < 0.01:", int(((ok.arm_v1_freq >= 0.01) & (ok.arm_v2_freq >= 0.01)).sum()))

# ---- double-het donors, the only read-backed phase candidates
ase = pd.read_csv(TAB / "ase_pair_het_donors.tsv", sep="\t")
# One variant pair can appear under several signals; the unit here is the unordered SITE PAIR,
# not the signal-by-pair row, and donors are never pooled across the two deposits.
ase["site_pair"] = [" | ".join(sorted([r.v1, r.v2])) for r in ase.itertuples()]
sites = (ase.groupby(["dataset", "site_pair"])
         .agg(separation_bp=("separation_bp", "first"),
              both_sites_measured=("both_sites_measured", "max"),
              n_donors_het_at_both=("n_donors_het_at_both", "max"),
              donors_het_at_both=("donors_het_at_both", "first"),
              genes=("gene", lambda s: ",".join(sorted({str(x) for x in s}))),
              n_signals=("signal_uid", "nunique"))
         .reset_index().sort_values(["dataset", "separation_bp"]))
sites.to_csv(TAB / "f1_read_backed_phase_candidates.tsv", sep="\t", index=False)
summ = (sites.groupby("dataset")
        .apply(lambda d: pd.Series(dict(
            n_site_pairs_both_measured=int(d.both_sites_measured.sum()),
            n_site_pairs_with_double_het_donor=int((d.n_donors_het_at_both > 0).sum()),
            n_site_pairs_le_500bp_both_measured=int(
                ((d.separation_bp <= 500) & d.both_sites_measured).sum()),
            n_site_pairs_le_500bp_with_double_het=int(
                ((d.separation_bp <= 500) & (d.n_donors_het_at_both > 0)).sum()),
            max_donors_at_one_site_pair=int(d.n_donors_het_at_both.max()),
            max_donors_le_500bp=int(d.loc[d.separation_bp <= 500, "n_donors_het_at_both"].max()
                                    if (d.separation_bp <= 500).any() else 0))),
               include_groups=False)
        .reset_index())
summ.to_csv(TAB / "f1_read_backed_phase_summary.tsv", sep="\t", index=False)
print("\nread-backed candidates, unordered site pairs:\n" + summ.to_string(index=False))
print("\nsite pairs <= 500 bp with at least one double-het donor:\n" +
      sites[(sites.separation_bp <= 500) & (sites.n_donors_het_at_both > 0)][
          ["dataset", "genes", "site_pair", "separation_bp", "n_donors_het_at_both",
           "donors_het_at_both"]].to_string(index=False))
print("\nall site pairs with at least one double-het donor, any separation: "
      + str(int((sites.n_donors_het_at_both > 0).sum())) + " dataset-rows, "
      + str(sites.loc[sites.n_donors_het_at_both > 0, "site_pair"].nunique())
      + " distinct site pairs")

json.dump(dict(
    rarer_arm=ok.rarer_arm.value_counts().to_dict(),
    both_arms_lt_001=int(((ok.arm_v1_freq < 0.01) & (ok.arm_v2_freq < 0.01)).sum()),
    exactly_one_arm_lt_001=int(((ok.arm_v1_freq < 0.01) ^ (ok.arm_v2_freq < 0.01)).sum()),
    neither_arm_lt_001=int(((ok.arm_v1_freq >= 0.01) & (ok.arm_v2_freq >= 0.01)).sum()),
    n_blocks_1mb_all=int(mp.block_1mb.nunique()),
), open(TAB / "f1_arm_detail.json", "w"), indent=1)
print("\ndone")
