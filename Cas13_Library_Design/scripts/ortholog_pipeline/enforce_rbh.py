#!/usr/bin/env python
"""
enforce_rbh.py — Reciprocal Best Hit (Wolf & Koonin 2012, Genome Biol).

A pair (M, H) is an RBH iff:
    argmax_h bitscore(M -> h)  ==  H     (M's best human partner is H)
  AND
    argmax_m bitscore(H -> m)  ==  M     (H's best mouse partner is M)

This is the most basic ortholog-inference test. Without it, paralogs (e.g. the
4 false NEAT1 mouse partners Frmd8os/9430037G07Rik/9330159M07Rik/Gm45470)
get inflated to Tier M alongside the true ortholog (Neat1).

Inputs:
    --m2h  L4_blast.tsv (gene-level, blastn mouse query vs human db)
    --h2m  L4_blast_h2m_gene.tsv (gene-level, blastn human query vs mouse db)
    --out  output TSV path

Output schema (matches L4_blast.tsv plus h2m columns + is_rbh flag).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def best_partner(df: pd.DataFrame, query_col: str, target_col: str,
                 score_col: str) -> pd.DataFrame:
    """Return one row per query_col: the target with max score_col, with ties
    broken deterministically (by target_col ascending) for reproducibility."""
    # Sort: score DESC, target ASC (stable tiebreak)
    df_sorted = df.sort_values([score_col, target_col],
                               ascending=[False, True], kind="mergesort")
    return df_sorted.drop_duplicates(subset=[query_col], keep="first").copy()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--m2h", required=True, type=Path,
                   help="Gene-level mouse->human BLAST pairs (L4_blast.tsv)")
    p.add_argument("--h2m", required=True, type=Path,
                   help="Gene-level human->mouse BLAST pairs (L4_blast_h2m_gene.tsv)")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    print(f"[rbh] loading m->h pairs from {args.m2h}", file=sys.stderr)
    m2h = pd.read_csv(args.m2h, sep="\t", low_memory=False)
    for c in ("blast_bitscore", "blast_pident", "blast_length", "blast_evalue"):
        if c in m2h.columns:
            m2h[c] = pd.to_numeric(m2h[c], errors="coerce")
    print(f"[rbh]   m->h gene pairs: {len(m2h)} "
          f"(unique mouse genes: {m2h['mouse_ensembl'].nunique()})", file=sys.stderr)

    print(f"[rbh] loading h->m pairs from {args.h2m}", file=sys.stderr)
    h2m = pd.read_csv(args.h2m, sep="\t", low_memory=False)
    for c in ("blast_bitscore_h2m", "blast_pident_h2m", "blast_length_h2m",
              "blast_evalue_h2m"):
        if c in h2m.columns:
            h2m[c] = pd.to_numeric(h2m[c], errors="coerce")
    print(f"[rbh]   h->m gene pairs: {len(h2m)} "
          f"(unique human genes: {h2m['human_ensembl'].nunique()})", file=sys.stderr)

    # Step 1: For each mouse gene, find its single best human partner.
    print(f"[rbh] finding m->h best partner per mouse gene...", file=sys.stderr)
    m2h_best = best_partner(m2h, "mouse_ensembl", "human_ensembl",
                            "blast_bitscore")
    print(f"[rbh]   mouse genes with a best human partner: {len(m2h_best)}",
          file=sys.stderr)

    # Step 2: For each human gene, find its single best mouse partner.
    print(f"[rbh] finding h->m best partner per human gene...", file=sys.stderr)
    h2m_best = best_partner(h2m, "human_ensembl", "mouse_ensembl",
                            "blast_bitscore_h2m")
    print(f"[rbh]   human genes with a best mouse partner: {len(h2m_best)}",
          file=sys.stderr)

    # Step 3: Inner-merge: the pair must appear in BOTH "best" tables.
    keys = ["mouse_ensembl", "human_ensembl"]
    m2h_keys = m2h_best[keys + ["blast_pident", "blast_length", "blast_evalue",
                                "blast_bitscore", "mouse_symbol", "human_symbol"]]
    h2m_keys = h2m_best[keys + ["blast_pident_h2m", "blast_length_h2m",
                                "blast_evalue_h2m", "blast_bitscore_h2m"]]

    rbh = m2h_keys.merge(h2m_keys, on=keys, how="inner")
    print(f"[rbh] reciprocal best hit pairs: {len(rbh)}", file=sys.stderr)

    # Diagnostics: how many m->h "best" pairs got pruned by RBH?
    pruned_m2h = len(m2h_best) - len(rbh)
    print(f"[rbh]   pruned non-reciprocal m->h best pairs: {pruned_m2h}",
          file=sys.stderr)
    print(f"[rbh]   RBH coverage of m->h best: {len(rbh) / max(1, len(m2h_best)):.1%}",
          file=sys.stderr)

    # Step 4: Emit common-schema output. RBH pairs are stronger evidence than
    # plain BLAST, so we mark them with their own tier flag and keep "M".
    out = pd.DataFrame({
        "mouse_ensembl": rbh["mouse_ensembl"],
        "mouse_symbol": rbh["mouse_symbol"],
        "human_ensembl": rbh["human_ensembl"],
        "human_symbol": rbh["human_symbol"],
        "blast_pident": rbh["blast_pident"],
        "blast_length": rbh["blast_length"],
        "blast_evalue": rbh["blast_evalue"],
        "blast_bitscore": rbh["blast_bitscore"],
        "blast_pident_h2m": rbh["blast_pident_h2m"],
        "blast_length_h2m": rbh["blast_length_h2m"],
        "blast_bitscore_h2m": rbh["blast_bitscore_h2m"],
        "is_rbh": 1,
        "tier_M_blast_rbh": 1,
        "confidence_tier": "M",
        "provenance_sources": "blast_rbh",
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[rbh] wrote {len(out)} RBH pairs to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
