#!/usr/bin/env python
"""
prep_orthofinder_fastas.py — Prepare protein FASTAs for OrthoFinder.

OrthoFinder expects one FASTA per species with unique headers. GENCODE
protein FASTAs contain multiple isoforms per gene. This script:
  1. Parses GENCODE protein FASTAs (human v49, mouse vM38)
  2. Selects the longest isoform per gene (canonical OrthoFinder practice)
  3. Writes clean FASTAs with headers: >ENSG_base|gene_symbol
  4. Writes an ID mapping TSV for downstream re-annotation

GENCODE header format:
  >ENSP...|ENST...|ENSG...|...|...|transcript_name|gene_symbol|length

Usage:
  python prep_orthofinder_fastas.py
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HUMAN_FA = PROJECT / "data/external/orthologs/gencode.v49.pc_translations.fa.gz"
MOUSE_FA = PROJECT / "data/external/orthologs/gencode.vM38.pc_translations.fa.gz"
WORK = (PROJECT / "Cas13_Library_Design/scripts/ortholog_pipeline"
        / "orthofinder_work")
OUT_DIR = WORK / "proteomes"
MAP_DIR = WORK / "id_maps"


def strip_version(eid: str) -> str:
    """ENSG00000186092.7 -> ENSG00000186092"""
    return eid.split(".")[0]


def parse_gencode_protein_fasta(fasta_path: Path) -> dict:
    """Parse GENCODE protein FASTA, return {ensg_base: {ensp, symbol, length, seq}}
    keeping only the longest isoform per gene."""
    genes: dict[str, dict] = {}

    with gzip.open(fasta_path, "rt") as f:
        current_header = None
        current_seq_parts = []

        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                # Process previous record
                if current_header is not None:
                    _store_if_longest(genes, current_header,
                                      "".join(current_seq_parts))
                current_header = line[1:]  # strip ">"
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        # Last record
        if current_header is not None:
            _store_if_longest(genes, current_header,
                              "".join(current_seq_parts))

    return genes


def _store_if_longest(genes: dict, header: str, seq: str):
    """Store protein if it's longer than what we have for this gene."""
    fields = header.split("|")
    if len(fields) < 7:
        return  # malformed header

    ensp = fields[0]
    enst = fields[1]
    ensg = strip_version(fields[2])
    symbol = fields[6] if len(fields) > 6 else ensg
    prot_len = len(seq)

    if ensg not in genes or prot_len > genes[ensg]["length"]:
        genes[ensg] = {
            "ensp": strip_version(ensp),
            "enst": strip_version(enst),
            "symbol": symbol,
            "length": prot_len,
            "seq": seq,
        }


def write_fasta(genes: dict, out_path: Path):
    """Write one-protein-per-gene FASTA with header >ENSG_base|symbol."""
    with open(out_path, "w") as f:
        for ensg in sorted(genes.keys()):
            info = genes[ensg]
            f.write(f">{ensg}|{info['symbol']}\n")
            # Write sequence in 60-char lines
            seq = info["seq"]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
    print(f"  Wrote {len(genes)} proteins to {out_path}", file=sys.stderr)


def write_id_map(genes: dict, out_path: Path, species: str):
    """Write ID mapping TSV: ensg_base, symbol, ensp_base, enst_base, prot_length."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("ensembl_gene\tgene_symbol\tensembl_protein\t"
                "ensembl_transcript\tprotein_length\tspecies\n")
        for ensg in sorted(genes.keys()):
            info = genes[ensg]
            f.write(f"{ensg}\t{info['symbol']}\t{info['ensp']}\t"
                    f"{info['enst']}\t{info['length']}\t{species}\n")
    print(f"  Wrote ID map to {out_path}", file=sys.stderr)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)

    print("Parsing human protein FASTA...", file=sys.stderr)
    human = parse_gencode_protein_fasta(HUMAN_FA)
    print(f"  {len(human)} human genes (longest isoform per gene)",
          file=sys.stderr)
    write_fasta(human, OUT_DIR / "Human.fa")
    write_id_map(human, MAP_DIR / "human_id_map.tsv", "human")

    print("Parsing mouse protein FASTA...", file=sys.stderr)
    mouse = parse_gencode_protein_fasta(MOUSE_FA)
    print(f"  {len(mouse)} mouse genes (longest isoform per gene)",
          file=sys.stderr)
    write_fasta(mouse, OUT_DIR / "Mouse.fa")
    write_id_map(mouse, MAP_DIR / "mouse_id_map.tsv", "mouse")

    # Summary
    print(f"\nSummary:", file=sys.stderr)
    print(f"  Human: {len(human):,} genes  "
          f"(from {HUMAN_FA.name})", file=sys.stderr)
    print(f"  Mouse: {len(mouse):,} genes  "
          f"(from {MOUSE_FA.name})", file=sys.stderr)
    print(f"  Total sequences for OrthoFinder: "
          f"{len(human) + len(mouse):,}", file=sys.stderr)
    print(f"\nProteome dir ready: {OUT_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
