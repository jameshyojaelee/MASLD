#!/usr/bin/env python
"""
rbh_fp_resolution_report.py — Quantify how RBH filters paralogs out of L4_blast.

For a set of canonical multi-partner lncRNAs (NEAT1, MEG3, HOTAIR, XIST, H19,
MIR122, MALAT1) and 50 random human lncRNAs with >= 5 m->h partners, report:
    n_partners_before  partners in m->h L4_blast.tsv
    n_partners_after   partners surviving RBH

Also reports the aggregate before/after counts at the table level.

Usage:
    python rbh_fp_resolution_report.py \\
        --m2h data/external/orthologs/layers/L4_blast.tsv \\
        --rbh data/external/orthologs/layers/L4_blast_rbh.tsv \\
        --out Cas13_Library_Design/reviews/rbh_fp_resolution.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


CANONICAL = ["NEAT1", "MEG3", "HOTAIR", "XIST", "H19", "MIR122", "MALAT1"]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--m2h", required=True, type=Path,
                   help="L4_blast.tsv (gene-level m->h)")
    p.add_argument("--rbh", required=True, type=Path,
                   help="L4_blast_rbh.tsv (RBH-passing pairs)")
    p.add_argument("--out", required=True, type=Path,
                   help="Output TSV for the per-gene resolution table")
    p.add_argument("--n-random", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    m2h = pd.read_csv(args.m2h, sep="\t", low_memory=False)
    rbh = pd.read_csv(args.rbh, sep="\t", low_memory=False)

    n_m2h_pairs = len(m2h)
    n_rbh_pairs = len(rbh)
    print(f"[fp] m->h pairs: {n_m2h_pairs}", file=sys.stderr)
    print(f"[fp] RBH pairs : {n_rbh_pairs}", file=sys.stderr)
    print(f"[fp] RBH retention: {n_rbh_pairs / max(1, n_m2h_pairs):.1%}",
          file=sys.stderr)

    # Build sets for membership tests
    m2h_set = set(zip(m2h["mouse_ensembl"], m2h["human_ensembl"]))
    rbh_set = set(zip(rbh["mouse_ensembl"], rbh["human_ensembl"]))
    if not rbh_set.issubset(m2h_set):
        # RBH should always be a subset of m->h (since it's defined by being
        # m->h best AND h->m best). If not, sanity-flag it.
        extra = len(rbh_set - m2h_set)
        print(f"[fp] WARNING: {extra} RBH pairs are not in m->h L4_blast.tsv",
              file=sys.stderr)

    # Per-human-gene partner counts
    m2h["sym_upper"] = m2h["human_symbol"].astype(str).str.upper()
    rbh["sym_upper"] = rbh["human_symbol"].astype(str).str.upper()
    m2h_by_hgene = m2h.groupby("human_ensembl")
    rbh_by_hgene = rbh.groupby("human_ensembl")

    rows = []

    # Canonical cases
    for sym in CANONICAL:
        sub_m = m2h[m2h["sym_upper"] == sym]
        if sub_m.empty:
            rows.append({
                "kind": "canonical",
                "human_symbol": sym,
                "human_ensembl": "",
                "n_before": 0,
                "n_after": 0,
                "kept_mouse_symbols": "",
                "dropped_mouse_symbols": "",
            })
            continue
        for hg, before in sub_m.groupby("human_ensembl"):
            after = rbh[rbh["human_ensembl"] == hg]
            kept = sorted(after["mouse_symbol"].dropna().astype(str).unique())
            dropped = sorted(set(before["mouse_symbol"].dropna().astype(str))
                             - set(kept))
            rows.append({
                "kind": "canonical",
                "human_symbol": sym,
                "human_ensembl": hg,
                "n_before": len(before),
                "n_after": len(after),
                "kept_mouse_symbols": ",".join(kept),
                "dropped_mouse_symbols": ",".join(dropped),
            })

    # Random sample of human genes with >= 5 m->h partners
    counts = m2h.groupby("human_ensembl").size().rename("n")
    multi_partner_genes = counts[counts >= 5].index.tolist()
    rng = np.random.default_rng(args.seed)
    sample = rng.choice(multi_partner_genes,
                        size=min(args.n_random, len(multi_partner_genes)),
                        replace=False)
    for hg in sample:
        before = m2h[m2h["human_ensembl"] == hg]
        after = rbh[rbh["human_ensembl"] == hg]
        kept = sorted(after["mouse_symbol"].dropna().astype(str).unique())
        dropped = sorted(set(before["mouse_symbol"].dropna().astype(str))
                         - set(kept))
        rows.append({
            "kind": "random_multi_partner",
            "human_symbol": before["human_symbol"].iloc[0] if not before.empty else "",
            "human_ensembl": hg,
            "n_before": len(before),
            "n_after": len(after),
            "kept_mouse_symbols": ",".join(kept),
            "dropped_mouse_symbols": ",".join(dropped),
        })

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[fp] wrote {len(df)} rows to {args.out}", file=sys.stderr)

    # Console summary
    print("\nCanonical multi-partner lncRNAs (before -> after RBH):",
          file=sys.stderr)
    canon = df[df["kind"] == "canonical"]
    for _, r in canon.iterrows():
        print(f"  {r['human_symbol']:8s}  {r['n_before']:>2d} -> {r['n_after']:>2d}  "
              f"kept={r['kept_mouse_symbols'] or '-'}  "
              f"dropped={r['dropped_mouse_symbols'] or '-'}", file=sys.stderr)

    rand = df[df["kind"] == "random_multi_partner"]
    if not rand.empty:
        before_sum = rand["n_before"].sum()
        after_sum = rand["n_after"].sum()
        print(f"\nRandom 50-gene panel (>=5 partners each):", file=sys.stderr)
        print(f"  total partners before RBH: {before_sum}", file=sys.stderr)
        print(f"  total partners after  RBH: {after_sum}", file=sys.stderr)
        print(f"  reduction: {1 - after_sum / max(1, before_sum):.1%}",
              file=sys.stderr)

    print(f"\nAggregate L4_blast.tsv:  {n_m2h_pairs} pairs", file=sys.stderr)
    print(f"Aggregate L4_blast_rbh.tsv: {n_rbh_pairs} pairs "
          f"({n_rbh_pairs / max(1, n_m2h_pairs):.1%} retained)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
