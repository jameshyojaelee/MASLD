#!/usr/bin/env python
"""
blast_to_orthotable.py — Aggregate blastn transcript-level hits to gene-level
mouse↔human lncRNA pairs in the common ortholog-table schema.

Inputs:
  --blast   blast outfmt 6 (qseqid sseqid pident length ... evalue bitscore)
  --mouse-fa  mouse lncRNA FASTA (for transcript_id → gene_id mapping)
  --human-fa  human lncRNA FASTA (same)
  --out     output TSV path

Strategy:
  1. Parse FASTA headers (GENCODE format: tx|gene|...|gene_name|length|biotype|)
     to build transcript→(gene_id, gene_name) maps.
  2. For each blast hit, lift to (mouse_gene, human_gene) pair and keep the
     best alignment per pair (max bitscore).
  3. Filter at e-value < 1e-5, pident ≥ 50%, length ≥ 50 bp (sensible defaults).
  4. Assign repeat_risk_flag per pair (P0-6 fix): at 90 Myr divergence, short
     high-identity fragments are likely LINE/SINE repeat-derived, not truly
     orthologous. RepeatMasker is unavailable on this cluster, so we flag based
     on alignment characteristics (Shah et al. 2019 signature):
       high_risk   : alignment_length < 100 AND pident > 85%
       medium_risk : alignment_length < 100 AND pident 70–85%
       low_risk    : alignment_length >= 100
  5. Emit common-schema TSV with tier_M_blast=1, confidence_tier=M.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_fasta_headers(fa_path: Path) -> pd.DataFrame:
    """Read the side-TSV emitted by the sbatch step (tx_id, gene_id, gene_name, length).

    The sbatch step writes simplified FASTAs (header = tx_id only) and a side
    TSV with the full mapping. This loader prefers the TSV but falls back to
    parsing GENCODE pipe-delimited headers if the TSV is missing.
    """
    tsv_path = Path(str(fa_path).replace(".fa", "_headers.tsv"))
    if tsv_path.exists():
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        df["gene_id_base"] = df["gene_id"].astype(str).str.split(".").str[0]
        return df
    # Fallback: parse FASTA headers directly (GENCODE pipe-delimited)
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
    p.add_argument("--blast", required=True, type=Path)
    p.add_argument("--mouse-fa", required=True, type=Path)
    p.add_argument("--human-fa", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--min-pident", type=float, default=50.0)
    p.add_argument("--min-length", type=int, default=50)
    p.add_argument("--max-evalue", type=float, default=1e-5)
    args = p.parse_args()

    print(f"[blast] parsing FASTA headers...", file=sys.stderr)
    mouse_tx = parse_fasta_headers(args.mouse_fa)
    human_tx = parse_fasta_headers(args.human_fa)
    print(f"[blast]  mouse transcripts: {len(mouse_tx)}", file=sys.stderr)
    print(f"[blast]  human transcripts: {len(human_tx)}", file=sys.stderr)

    print(f"[blast] reading blast output...", file=sys.stderr)
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    bl = pd.read_csv(args.blast, sep="\t", names=cols)
    print(f"[blast]  raw blast hits: {len(bl)}", file=sys.stderr)

    # Filter at thresholds
    bl = bl[
        (bl["pident"] >= args.min_pident)
        & (bl["length"] >= args.min_length)
        & (bl["evalue"] <= args.max_evalue)
    ].copy()
    print(f"[blast]  after pident>={args.min_pident} length>={args.min_length} "
          f"evalue<={args.max_evalue}: {len(bl)} hits", file=sys.stderr)

    # Map transcript → gene (use unversioned IDs)
    mouse_lookup = dict(zip(mouse_tx["tx_id"], zip(mouse_tx["gene_id_base"],
                                                   mouse_tx["gene_name"])))
    human_lookup = dict(zip(human_tx["tx_id"], zip(human_tx["gene_id_base"],
                                                   human_tx["gene_name"])))

    bl["mouse_ensembl"] = bl["qseqid"].map(lambda t: mouse_lookup.get(t, (None, None))[0])
    bl["mouse_symbol"] = bl["qseqid"].map(lambda t: mouse_lookup.get(t, (None, None))[1])
    bl["human_ensembl"] = bl["sseqid"].map(lambda t: human_lookup.get(t, (None, None))[0])
    bl["human_symbol"] = bl["sseqid"].map(lambda t: human_lookup.get(t, (None, None))[1])

    bl = bl.dropna(subset=["mouse_ensembl", "human_ensembl"])
    print(f"[blast]  after gene-mapping dropna: {len(bl)} hits", file=sys.stderr)

    # Per (mouse_gene, human_gene), keep best alignment by bitscore
    bl = bl.sort_values("bitscore", ascending=False).drop_duplicates(
        subset=["mouse_ensembl", "human_ensembl"], keep="first"
    )
    print(f"[blast]  gene-level pairs (best per pair): {len(bl)}", file=sys.stderr)

    # ---- P0-6: Repeat risk flagging ----------------------------------------
    # At ~90 Myr mouse-human divergence, short high-identity BLAST fragments
    # are classic LINE/SINE repeat signatures rather than true orthology.
    # RepeatMasker is not available on this cluster, so we flag based on
    # alignment characteristics per Shah et al. (2019) / Smit et al. (2020).
    #   high_risk   : length < 100 AND pident > 85  (repeat hallmark)
    #   medium_risk : length < 100 AND pident [70, 85]
    #   low_risk    : length >= 100
    def _repeat_risk(row) -> str:
        aln_len = float(row.get("length", 0) or 0)
        pident = float(row.get("pident", 0) or 0)
        if aln_len < 100:
            return "high_risk" if pident > 85 else "medium_risk"
        return "low_risk"

    bl["repeat_risk_flag"] = bl.apply(_repeat_risk, axis=1)
    risk_counts = bl["repeat_risk_flag"].value_counts()
    print(f"[blast]  repeat_risk_flag distribution:", file=sys.stderr)
    for k, v in risk_counts.items():
        print(f"[blast]    {k}: {v} ({v / len(bl):.1%})", file=sys.stderr)

    # Common-schema output
    out = pd.DataFrame({
        "mouse_ensembl": bl["mouse_ensembl"],
        "mouse_symbol": bl["mouse_symbol"],
        "mouse_biotype": "lncRNA",
        "human_ensembl": bl["human_ensembl"],
        "human_symbol": bl["human_symbol"],
        "human_biotype": "lncRNA",
        "tier_H_biomart": 0,
        "tier_H_mirbase": 0,
        "tier_H_mirgenedb": 0,
        "tier_M_phasej_synteny": 0,
        "tier_M_lncbook": 0,
        "tier_M_noncode": 0,
        "tier_M_hezroni2015": 0,
        "tier_M_snodb": 0,
        "tier_M_blast": 1,
        "blast_pident": bl["pident"],
        "blast_length": bl["length"],
        "blast_evalue": bl["evalue"],
        "blast_bitscore": bl["bitscore"],
        "repeat_risk_flag": bl["repeat_risk_flag"],
        "tier_M_ortho2align": 0,
        "tier_M_ortho2align_pvalue": pd.NA,
        "tier_M_ortho2align_identity": pd.NA,
        "tier_M_pseudogene_parent": 0,
        "tier_L_liftover": 0,
        "tier_L_liftover_overlap_bp": pd.NA,
        "tier_L_reciprocal_liftover": 0,
        "confidence_tier": "M",
        "evidence_count": 1,
        "provenance_sources": "blastn_dc-megablast",
        "notes": pd.NA,
    })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[blast] wrote {len(out)} mouse↔human lncRNA pairs to {args.out}",
          file=sys.stderr)

    # Summary stats
    print(f"[blast] unique mouse_ensembl: {out['mouse_ensembl'].nunique()}",
          file=sys.stderr)
    print(f"[blast] unique human_ensembl: {out['human_ensembl'].nunique()}",
          file=sys.stderr)
    print(f"[blast] pident summary:", file=sys.stderr)
    print(out["blast_pident"].describe().to_string(), file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
