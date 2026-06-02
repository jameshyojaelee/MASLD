#!/usr/bin/env python
"""
compare_blast_mmseqs2.py — Compare L4 BLAST-RBH vs L5 MMseqs2-RBH layers.

Reports overlap, unique pairs, pident correlation for shared pairs,
and canonical lncRNA spot-checks.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
LAYERS = PROJECT / "data/external/orthologs/layers"

CANONICAL_PAIRS = [
    ("Meg3", "MEG3"),
    ("Malat1", "MALAT1"),
    ("Neat1", "NEAT1"),
    ("Hotair", "HOTAIR"),
    ("Xist", "XIST"),
    ("H19", "H19"),
]


def main() -> int:
    blast_path = LAYERS / "L4_blast_rbh.tsv"
    mmseqs_path = LAYERS / "L5_mmseqs2_rbh.tsv"

    if not blast_path.exists():
        print(f"ERROR: {blast_path} not found", file=sys.stderr)
        return 1
    if not mmseqs_path.exists():
        print(f"ERROR: {mmseqs_path} not found", file=sys.stderr)
        return 1

    blast = pd.read_csv(blast_path, sep="\t", low_memory=False)
    mmseqs = pd.read_csv(mmseqs_path, sep="\t", low_memory=False)

    print("=" * 70)
    print("L4 BLAST-RBH vs L5 MMseqs2-RBH comparison")
    print("=" * 70)

    # ── Counts ────────────────────────────────────────────────────────────
    print(f"\nL4 BLAST-RBH pairs:  {len(blast):,}")
    print(f"L5 MMseqs2-RBH pairs: {len(mmseqs):,}")

    # ── Overlap by (mouse_ensembl, human_ensembl) ─────────────────────────
    blast_pairs = set(zip(blast["mouse_ensembl"], blast["human_ensembl"]))
    mmseqs_pairs = set(zip(mmseqs["mouse_ensembl"], mmseqs["human_ensembl"]))

    overlap = blast_pairs & mmseqs_pairs
    blast_only = blast_pairs - mmseqs_pairs
    mmseqs_only = mmseqs_pairs - blast_pairs

    print(f"\nOverlap (both):    {len(overlap):,}")
    print(f"BLAST-only:        {len(blast_only):,}")
    print(f"MMseqs2-only:      {len(mmseqs_only):,}")
    if len(blast_pairs) > 0:
        print(f"Jaccard index:     {len(overlap) / len(blast_pairs | mmseqs_pairs):.3f}")
        print(f"BLAST recall by MMseqs2: {len(overlap) / len(blast_pairs):.1%}")
        print(f"MMseqs2 recall by BLAST: {len(overlap) / len(mmseqs_pairs):.1%}")

    # ── pident correlation for overlapping pairs ──────────────────────────
    if len(overlap) > 10:
        print(f"\n--- pident correlation for {len(overlap):,} overlapping pairs ---")
        blast_idx = blast.set_index(["mouse_ensembl", "human_ensembl"])
        mmseqs_idx = mmseqs.set_index(["mouse_ensembl", "human_ensembl"])

        overlap_list = list(overlap)
        b_pident = blast_idx.loc[overlap_list, "blast_pident"].values
        m_pident = mmseqs_idx.loc[overlap_list, "mmseqs2_pident"].values

        r, p_r = stats.pearsonr(b_pident, m_pident)
        rho, p_rho = stats.spearmanr(b_pident, m_pident)
        print(f"  Pearson r  = {r:.4f} (p = {p_r:.2e})")
        print(f"  Spearman rho = {rho:.4f} (p = {p_rho:.2e})")
        diff = m_pident - b_pident
        print(f"  Mean pident diff (MMseqs2 - BLAST): {np.mean(diff):.2f}")
        print(f"  Median pident diff: {np.median(diff):.2f}")

    # ── Repeat risk distribution comparison ───────────────────────────────
    print("\n--- Repeat risk distribution ---")
    if "repeat_risk_flag" in blast.columns:
        print("\nBLAST-RBH (L4_blast.tsv, not L4_blast_rbh which lacks this col):")
        blast_all = pd.read_csv(LAYERS / "L4_blast.tsv", sep="\t", low_memory=False)
        if "repeat_risk_flag" in blast_all.columns:
            print(blast_all["repeat_risk_flag"].value_counts().to_string())
    print("\nMMseqs2-RBH (L5):")
    print(mmseqs["repeat_risk_flag"].value_counts().to_string())

    # ── Canonical pair spot-check ─────────────────────────────────────────
    print("\n--- Canonical lncRNA pair spot-check ---")
    print(f"{'Mouse':<12} {'Human':<12} {'BLAST-RBH':<12} {'MMseqs2-RBH':<12}")
    print("-" * 50)
    for m_sym, h_sym in CANONICAL_PAIRS:
        in_blast = any(
            (str(r.get("mouse_symbol", "")).upper() == m_sym.upper()
             and str(r.get("human_symbol", "")).upper() == h_sym.upper())
            for _, r in blast.iterrows()
        ) if len(blast) < 5000 else (
            blast[
                (blast["mouse_symbol"].str.upper() == m_sym.upper())
                & (blast["human_symbol"].str.upper() == h_sym.upper())
            ].shape[0] > 0
        )
        in_mmseqs = mmseqs[
            (mmseqs["mouse_symbol"].str.upper() == m_sym.upper())
            & (mmseqs["human_symbol"].str.upper() == h_sym.upper())
        ].shape[0] > 0

        b_mark = "YES" if in_blast else "no"
        m_mark = "YES" if in_mmseqs else "no"
        print(f"{m_sym:<12} {h_sym:<12} {b_mark:<12} {m_mark:<12}")

    # ── MMseqs2-only sample: top 10 by bitscore ──────────────────────────
    if mmseqs_only:
        print(f"\n--- Top 10 MMseqs2-only pairs by bitscore ---")
        mmseqs_only_df = mmseqs[
            mmseqs.apply(
                lambda r: (r["mouse_ensembl"], r["human_ensembl"]) in mmseqs_only,
                axis=1,
            )
        ].nlargest(10, "mmseqs2_bitscore")
        print(mmseqs_only_df[[
            "mouse_symbol", "human_symbol",
            "mmseqs2_pident", "mmseqs2_length", "mmseqs2_bitscore",
            "repeat_risk_flag",
        ]].to_string(index=False))

    # ── BLAST-only sample: top 10 by bitscore ────────────────────────────
    if blast_only:
        print(f"\n--- Top 10 BLAST-only pairs by bitscore ---")
        blast_only_df = blast[
            blast.apply(
                lambda r: (r["mouse_ensembl"], r["human_ensembl"]) in blast_only,
                axis=1,
            )
        ].nlargest(10, "blast_bitscore")
        cols_to_show = ["mouse_symbol", "human_symbol",
                        "blast_pident", "blast_length", "blast_bitscore"]
        cols_to_show = [c for c in cols_to_show if c in blast_only_df.columns]
        print(blast_only_df[cols_to_show].to_string(index=False))

    print("\n" + "=" * 70)
    print("Comparison complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
