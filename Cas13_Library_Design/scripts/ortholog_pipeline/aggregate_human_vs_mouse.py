#!/usr/bin/env python
"""
aggregate_human_vs_mouse.py — Collapse human->mouse blastn transcript hits to
gene-level (human_gene, mouse_gene) pairs, keeping the best alignment per pair.

This is the reverse-direction analogue of blast_to_orthotable.py. We do NOT
write tier flags here -- the output is consumed only by enforce_rbh.py, which
intersects forward and reverse BLAST to produce L4_blast_rbh.tsv.

Filters mirror the forward pipeline: pident >= 50, length >= 50, evalue <= 1e-5.

Output columns:
    human_ensembl, human_symbol, mouse_ensembl, mouse_symbol,
    blast_pident_h2m, blast_length_h2m, blast_evalue_h2m, blast_bitscore_h2m
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_fasta_headers(fa_path: Path) -> pd.DataFrame:
    """Load side-emitted *_headers.tsv (tx_id, gene_id, gene_name, length)."""
    tsv_path = Path(str(fa_path).replace(".fa", "_headers.tsv"))
    if tsv_path.exists():
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        df["gene_id_base"] = df["gene_id"].astype(str).str.split(".").str[0]
        return df
    rows = []
    with open(fa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            line = line[1:].strip()
            parts = line.split("|")
            tx_id = parts[0]
            gene_id = parts[1] if len(parts) > 1 else ""
            gene_name = parts[5] if len(parts) > 5 else ""
            rows.append((tx_id, gene_id, gene_name))
    df = pd.DataFrame(rows, columns=["tx_id", "gene_id", "gene_name"])
    df["gene_id_base"] = df["gene_id"].str.split(".").str[0]
    return df


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--blast", required=True, type=Path,
                   help="blast outfmt 6 (qseqid=human tx, sseqid=mouse tx)")
    p.add_argument("--mouse-fa", required=True, type=Path)
    p.add_argument("--human-fa", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--min-pident", type=float, default=50.0)
    p.add_argument("--min-length", type=int, default=50)
    p.add_argument("--max-evalue", type=float, default=1e-5)
    args = p.parse_args()

    print(f"[h2m] parsing FASTA headers...", file=sys.stderr)
    mouse_tx = parse_fasta_headers(args.mouse_fa)
    human_tx = parse_fasta_headers(args.human_fa)
    print(f"[h2m]   mouse transcripts: {len(mouse_tx)}", file=sys.stderr)
    print(f"[h2m]   human transcripts: {len(human_tx)}", file=sys.stderr)

    print(f"[h2m] reading blast output...", file=sys.stderr)
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    bl = pd.read_csv(args.blast, sep="\t", names=cols)
    print(f"[h2m]   raw blast hits: {len(bl)}", file=sys.stderr)

    bl = bl[
        (bl["pident"] >= args.min_pident)
        & (bl["length"] >= args.min_length)
        & (bl["evalue"] <= args.max_evalue)
    ].copy()
    print(f"[h2m]   after pident>={args.min_pident} length>={args.min_length} "
          f"evalue<={args.max_evalue}: {len(bl)} hits", file=sys.stderr)

    # qseqid=human, sseqid=mouse
    human_lookup = dict(zip(human_tx["tx_id"], zip(human_tx["gene_id_base"],
                                                   human_tx["gene_name"])))
    mouse_lookup = dict(zip(mouse_tx["tx_id"], zip(mouse_tx["gene_id_base"],
                                                   mouse_tx["gene_name"])))

    bl["human_ensembl"] = bl["qseqid"].map(lambda t: human_lookup.get(t, (None, None))[0])
    bl["human_symbol"] = bl["qseqid"].map(lambda t: human_lookup.get(t, (None, None))[1])
    bl["mouse_ensembl"] = bl["sseqid"].map(lambda t: mouse_lookup.get(t, (None, None))[0])
    bl["mouse_symbol"] = bl["sseqid"].map(lambda t: mouse_lookup.get(t, (None, None))[1])

    bl = bl.dropna(subset=["human_ensembl", "mouse_ensembl"])
    print(f"[h2m]   after gene-mapping dropna: {len(bl)} hits", file=sys.stderr)

    # Best alignment per (human_gene, mouse_gene) pair by bitscore
    bl = bl.sort_values("bitscore", ascending=False).drop_duplicates(
        subset=["human_ensembl", "mouse_ensembl"], keep="first"
    )
    print(f"[h2m]   gene-level pairs (best per pair): {len(bl)}", file=sys.stderr)

    out = pd.DataFrame({
        "human_ensembl": bl["human_ensembl"],
        "human_symbol": bl["human_symbol"],
        "mouse_ensembl": bl["mouse_ensembl"],
        "mouse_symbol": bl["mouse_symbol"],
        "blast_pident_h2m": bl["pident"],
        "blast_length_h2m": bl["length"],
        "blast_evalue_h2m": bl["evalue"],
        "blast_bitscore_h2m": bl["bitscore"],
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[h2m] wrote {len(out)} human->mouse lncRNA pairs to {args.out}",
          file=sys.stderr)
    print(f"[h2m] unique human_ensembl: {out['human_ensembl'].nunique()}",
          file=sys.stderr)
    print(f"[h2m] unique mouse_ensembl: {out['mouse_ensembl'].nunique()}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
