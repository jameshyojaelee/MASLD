#!/usr/bin/env python
"""
parse_mmseqs2.py — Parse MMseqs2 easy-rbh output to gene-level L5 ortholog table.

MMseqs2 easy-rbh produces transcript-level reciprocal best hits. This script:
  1. Maps transcript IDs to gene IDs using the same side-TSV files from the
     BLAST pipeline (*_headers.tsv).
  2. Aggregates to gene level (best bitscore per mouse-human gene pair).
  3. Applies repeat_risk_flag logic (same as blast_to_orthotable.py).
  4. Emits common-schema L5 layer.

Usage:
  python parse_mmseqs2.py \
    --mmseqs2 mmseqs2_work/mmseqs2_rbh_result \
    --mouse-headers blast_work/mouse_vM37_lncRNA_headers.tsv \
    --human-headers blast_work/human_v47_lncRNA_headers.tsv \
    --out data/external/orthologs/layers/L5_mmseqs2_rbh.tsv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_headers(tsv_path: Path) -> dict[str, tuple[str, str]]:
    """Load tx_id -> (gene_id_base, gene_name) from side-TSV."""
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    df["gene_id_base"] = df["gene_id"].str.split(".").str[0]
    return dict(zip(df["tx_id"], zip(df["gene_id_base"], df["gene_name"])))


def repeat_risk(aln_len: float, pident: float) -> str:
    """Classify repeat risk per Shah et al. (2019) / Smit et al. (2020).

    At ~90 Myr mouse-human divergence, short high-identity fragments are
    classic LINE/SINE repeat signatures rather than true orthology.
    """
    if aln_len < 100:
        return "high_risk" if pident > 85 else "medium_risk"
    return "low_risk"


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mmseqs2", required=True, type=Path,
                   help="MMseqs2 easy-rbh output (BLAST-like tab)")
    p.add_argument("--mouse-headers", required=True, type=Path,
                   help="Mouse transcript header TSV")
    p.add_argument("--human-headers", required=True, type=Path,
                   help="Human transcript header TSV")
    p.add_argument("--out", required=True, type=Path,
                   help="Output L5 layer TSV")
    p.add_argument("--min-pident", type=float, default=50.0,
                   help="Minimum percent identity (default: 50)")
    p.add_argument("--min-length", type=int, default=50,
                   help="Minimum alignment length in bp (default: 50)")
    p.add_argument("--max-evalue", type=float, default=1e-5,
                   help="Maximum E-value (default: 1e-5)")
    args = p.parse_args()

    # ── Load header mappings ──────────────────────────────────────────────
    print("[mmseqs2] loading transcript -> gene mappings...", file=sys.stderr)
    mouse_lookup = load_headers(args.mouse_headers)
    human_lookup = load_headers(args.human_headers)
    print(f"[mmseqs2]   mouse transcripts: {len(mouse_lookup)}", file=sys.stderr)
    print(f"[mmseqs2]   human transcripts: {len(human_lookup)}", file=sys.stderr)

    # ── Read MMseqs2 output ───────────────────────────────────────────────
    print(f"[mmseqs2] reading {args.mmseqs2}...", file=sys.stderr)
    cols = ["qseqid", "sseqid", "pident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df = pd.read_csv(args.mmseqs2, sep="\t", names=cols, header=None)
    print(f"[mmseqs2]   raw RBH hits: {len(df)}", file=sys.stderr)

    # ── Filter ────────────────────────────────────────────────────────────
    before = len(df)
    df = df[
        (df["pident"] >= args.min_pident)
        & (df["alnlen"] >= args.min_length)
        & (df["evalue"] <= args.max_evalue)
    ].copy()
    print(f"[mmseqs2]   after filters (pident>={args.min_pident}, "
          f"len>={args.min_length}, evalue<={args.max_evalue}): "
          f"{len(df)} hits (dropped {before - len(df)})", file=sys.stderr)

    # ── Map tx -> gene ────────────────────────────────────────────────────
    df["mouse_ensembl"] = df["qseqid"].map(
        lambda t: mouse_lookup.get(t, (None, None))[0])
    df["mouse_symbol"] = df["qseqid"].map(
        lambda t: mouse_lookup.get(t, (None, None))[1])
    df["human_ensembl"] = df["sseqid"].map(
        lambda t: human_lookup.get(t, (None, None))[0])
    df["human_symbol"] = df["sseqid"].map(
        lambda t: human_lookup.get(t, (None, None))[1])

    unmapped_q = df["mouse_ensembl"].isna().sum()
    unmapped_s = df["human_ensembl"].isna().sum()
    if unmapped_q > 0 or unmapped_s > 0:
        print(f"[mmseqs2]   WARNING: unmapped query={unmapped_q}, "
              f"target={unmapped_s}", file=sys.stderr)
    df = df.dropna(subset=["mouse_ensembl", "human_ensembl"])
    print(f"[mmseqs2]   after gene-mapping dropna: {len(df)}", file=sys.stderr)

    # ── Aggregate to gene level (best bitscore per pair) ──────────────────
    df = (df.sort_values("bitscore", ascending=False)
            .drop_duplicates(subset=["mouse_ensembl", "human_ensembl"],
                             keep="first"))
    print(f"[mmseqs2]   gene-level pairs: {len(df)}", file=sys.stderr)

    # ── Repeat risk flag ──────────────────────────────────────────────────
    df["repeat_risk_flag"] = [
        repeat_risk(row["alnlen"], row["pident"])
        for _, row in df.iterrows()
    ]
    risk_counts = df["repeat_risk_flag"].value_counts()
    print("[mmseqs2]   repeat_risk_flag distribution:", file=sys.stderr)
    for k, v in risk_counts.items():
        print(f"[mmseqs2]     {k}: {v} ({v / len(df):.1%})", file=sys.stderr)

    # ── Emit common-schema L5 table ───────────────────────────────────────
    out = pd.DataFrame({
        "mouse_ensembl": df["mouse_ensembl"].values,
        "mouse_symbol": df["mouse_symbol"].values,
        "mouse_biotype": "lncRNA",
        "human_ensembl": df["human_ensembl"].values,
        "human_symbol": df["human_symbol"].values,
        "human_biotype": "lncRNA",
        "tier_M_mmseqs2_rbh": 1,
        "mmseqs2_pident": df["pident"].values,
        "mmseqs2_length": df["alnlen"].values,
        "mmseqs2_evalue": df["evalue"].values,
        "mmseqs2_bitscore": df["bitscore"].values,
        "repeat_risk_flag": df["repeat_risk_flag"].values,
        "confidence_tier": "M",
        "provenance_sources": "mmseqs2_easy-rbh",
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"\n[mmseqs2] wrote {len(out)} gene-level RBH pairs to {args.out}",
          file=sys.stderr)

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"[mmseqs2] unique mouse_ensembl: {out['mouse_ensembl'].nunique()}",
          file=sys.stderr)
    print(f"[mmseqs2] unique human_ensembl: {out['human_ensembl'].nunique()}",
          file=sys.stderr)
    print(f"[mmseqs2] pident summary:", file=sys.stderr)
    print(out["mmseqs2_pident"].describe().to_string(), file=sys.stderr)
    print(f"[mmseqs2] bitscore summary:", file=sys.stderr)
    print(out["mmseqs2_bitscore"].describe().to_string(), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
