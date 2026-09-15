#!/usr/bin/env python3
"""B3 step 3: print the numbers that go into RESULTS.md, and the two prediction verdicts."""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

PRIMARY_SCORER = "CHIP_HISTONE"
PRIMARY_STAT = "mean_abs"


def fmt(r):
    return (f"{r.observed:+.4g} [{r.ci_lo:+.4g}, {r.ci_hi:+.4g}] "
            f"p={'<=2.0e-04' if r.p_boot <= r.p_floor + 1e-12 else f'{r.p_boot:.2g}'} q={r.q_bh:.2g}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tab = os.path.join(a.out, "tables")
    pd.set_option("display.width", 250)
    res = pd.read_csv(os.path.join(tab, "b3_contrasts.tsv"), sep="\t")
    meta = json.load(open(os.path.join(tab, "b3_boot_meta.json")))
    print("BOOT META:", meta)

    prim = res[(res.family == "primary_native")]
    print("\n=== PRIMARY FAMILY (adult-liver / AVI, native scale), %d tests ===" % len(prim))
    for stat in ["mean_abs", "p95_abs", "conc_top5"]:
        print(f"\n--- {stat} ---")
        p = prim[prim.stat == stat].pivot_table(index="scorer", columns="contrast",
                                                values=["observed", "ci_lo", "ci_hi", "q_bh"])
        for sc in sorted(prim.scorer.unique()):
            line = [f"{sc:13s}"]
            for c in sorted(prim.contrast.unique()):
                r = prim[(prim.scorer == sc) & (prim.stat == stat) & (prim.contrast == c)]
                if len(r) == 0:
                    continue
                r = r.iloc[0]
                line.append(f"{c[:2]} {fmt(r)}")
            print("  " + " | ".join(line))

    print("\n=== B3.1 verdict channel: %s adult_liver %s ===" % (PRIMARY_SCORER, PRIMARY_STAT))
    key = prim[(prim.scorer == PRIMARY_SCORER) & (prim.stat == PRIMARY_STAT)].set_index("contrast")
    if "T1_skill_Q5_minus_Q1" not in key.index:
        print("  verdict channel absent from this table; nothing to score")
        return
    t1 = key.loc["T1_skill_Q5_minus_Q1"]
    b31 = bool(t1.observed > 0 and t1.ci_lo > 0)
    print("  T1 Q5-Q1 =", fmt(t1), "-> B3.1", "MET" if b31 else "NOT MET")
    t3, t4, t2, t5 = (key.loc["T3_promoter_std"], key.loc["T4_log10_dist_tss_std"],
                      key.loc["T2_skill_slope_std"], key.loc["T5_skill_slope_std_adjusted"])
    b32 = bool(t3.observed > 0 and t3.ci_lo > 0 and t4.observed < 0 and t4.ci_hi < 0)
    print("  T3 promoter     =", fmt(t3))
    print("  T4 log10distTSS =", fmt(t4))
    print("  T2 skill slope  =", fmt(t2))
    print("  T5 skill adj    =", fmt(t5), "-> B3.2", "MET" if b32 else "NOT MET")

    print("\n=== direction agreement across the six scorers (native, primary group) ===")
    for c in ["T1_skill_Q5_minus_Q1", "T2_skill_slope_std", "T3_promoter_std",
              "T4_log10_dist_tss_std", "T5_skill_slope_std_adjusted"]:
        s = prim[(prim.contrast == c) & (prim.stat == PRIMARY_STAT)]
        pos = int((s.observed > 0).sum())
        sig = int((s.q_bh < 0.05).sum())
        excl = int(((s.ci_lo > 0) | (s.ci_hi < 0)).sum())
        print(f"  {c:32s} positive {pos}/{len(s)}  CI excludes 0 {excl}/{len(s)}  BH q<0.05 {sig}/{len(s)}")

    print("\n=== quintile profile (native, primary group) ===")
    qp = pd.read_csv(os.path.join(tab, "b3_quintile_profile.tsv"), sep="\t")
    print(qp.pivot_table(index=["scorer", "stat"], columns="quintile", values="mean").to_string())

    print("\n=== HepG2 secondary family, %s, T1 ===" % PRIMARY_STAT)
    h = res[(res.family == "hepg2_native") & (res.stat == PRIMARY_STAT) &
            (res.contrast == "T1_skill_Q5_minus_Q1")]
    for _, r in h.iterrows():
        print(f"  {r.scorer:13s} {fmt(r)}")

    print("\n=== quantile-scale family (secondary), %s, T1 ===" % PRIMARY_STAT)
    q = res[(res.family == "primary_quantile") & (res.stat == PRIMARY_STAT) &
            (res.contrast == "T1_skill_Q5_minus_Q1")]
    for _, r in q.iterrows():
        print(f"  {r.scorer:13s} {fmt(r)}")

    print("\n=== BH survivors per family ===")
    print(res.groupby("family").agg(tests=("q_bh", "size"),
                                    bh_lt_05=("q_bh", lambda s: int((s < 0.05).sum())),
                                    at_p_floor=("p_boot", lambda s: int((s <= 2.0 / 10001 + 1e-12).sum()))
                                    ).to_string())


if __name__ == "__main__":
    main()
