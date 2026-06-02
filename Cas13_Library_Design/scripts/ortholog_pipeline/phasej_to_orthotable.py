#!/usr/bin/env python
"""
phasej_to_orthotable.py — Convert Phase J ncRNA synteny output into the
master-ortholog-table common schema.

Phase J pipeline (RNA-seq/55_ncrna_conservation.R) computed flanking-PCG
synteny for all 12,671 human lncRNAs and recorded mouse lncRNAs in the
syntenic interval as `mouse_lnc_between` (semicolon-separated symbols).

This converter:
  1. Reads lncrna_synteny_conservation.csv
  2. Filters to synteny_status == "Synteny_conserved" (4,674 rows)
  3. Explodes mouse_lnc_between on ";"
  4. Resolves mouse_symbol → ENSMUSG via GENCODE vM38 metadata
  5. Resolves human_symbol → ENSG via GENCODE v49 metadata
  6. Writes L1.5_phasej_synteny.tsv in common schema
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEFAULT_PHASEJ = PROJECT / "RNA-seq/results/ncrna/lncrna_synteny_conservation.csv"
DEFAULT_MOUSE_META = PROJECT / "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"
DEFAULT_HUMAN_META = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
DEFAULT_OUT = PROJECT / "data/external/orthologs/layers/L1.5_phasej_synteny.tsv"


def load_mouse_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={
        "mouse_ensembl_base": "mouse_ensembl",
        "mouse_symbol_gtf": "mouse_symbol",
        "mouse_biotype": "mouse_biotype",
    })
    # Some mouse symbols appear with multiple Ensembl IDs; keep first.
    df = df.drop_duplicates(subset=["mouse_symbol"], keep="first")
    return df[["mouse_symbol", "mouse_ensembl", "mouse_biotype"]]


def load_human_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={
        "ensembl_base": "human_ensembl",
        "gene_name": "human_symbol",
        "gene_biotype": "human_biotype",
    })
    df = df.drop_duplicates(subset=["human_symbol"], keep="first")
    return df[["human_symbol", "human_ensembl", "human_biotype"]]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phasej", type=Path, default=DEFAULT_PHASEJ)
    p.add_argument("--mouse-meta", type=Path, default=DEFAULT_MOUSE_META)
    p.add_argument("--human-meta", type=Path, default=DEFAULT_HUMAN_META)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    print(f"[phasej] reading {args.phasej}", file=sys.stderr)
    pj = pd.read_csv(args.phasej)
    n_total = len(pj)
    pj = pj[pj["synteny_status"] == "Synteny_conserved"].copy()
    print(f"[phasej] {len(pj)}/{n_total} rows with Synteny_conserved", file=sys.stderr)

    pj["mouse_symbols_list"] = pj["mouse_lnc_between"].fillna("").str.split(";")
    exploded = pj.explode("mouse_symbols_list")
    exploded = exploded[exploded["mouse_symbols_list"].astype(str).str.len() > 0]
    exploded = exploded.rename(columns={
        "lncrna": "human_symbol",
        "mouse_symbols_list": "mouse_symbol",
        "upstream_pcg": "upstream_pcg_human",
        "downstream_pcg": "downstream_pcg_human",
        "mouse_upstream": "mouse_upstream_pcg",
        "mouse_downstream": "mouse_downstream_pcg",
    })
    exploded["mouse_symbol"] = exploded["mouse_symbol"].astype(str).str.strip()
    print(f"[phasej] exploded to {len(exploded)} (mouse, human) pairs", file=sys.stderr)

    mouse_meta = load_mouse_metadata(args.mouse_meta)
    human_meta = load_human_metadata(args.human_meta)
    print(f"[phasej] mouse meta {len(mouse_meta)}, human meta {len(human_meta)}",
          file=sys.stderr)

    merged = exploded.merge(mouse_meta, on="mouse_symbol", how="left")
    merged = merged.merge(human_meta, on="human_symbol", how="left")

    n_mouse_missing = merged["mouse_ensembl"].isna().sum()
    n_human_missing = merged["human_ensembl"].isna().sum()
    print(f"[phasej] resolution: mouse ENSMUSG missing {n_mouse_missing}, "
          f"human ENSG missing {n_human_missing}", file=sys.stderr)

    # Common schema columns
    out = pd.DataFrame({
        "mouse_ensembl": merged["mouse_ensembl"],
        "mouse_symbol": merged["mouse_symbol"],
        "mouse_biotype": merged["mouse_biotype"],
        "human_ensembl": merged["human_ensembl"],
        "human_symbol": merged["human_symbol"],
        "human_biotype": merged["human_biotype"],
        "tier_H_biomart": 0,
        "tier_H_mirbase": 0,
        "tier_H_mirgenedb": 0,
        "tier_M_phasej_synteny": 1,
        "tier_M_lncbook": 0,
        "tier_M_noncode": 0,
        "tier_M_hezroni2015": 0,
        "tier_M_snodb": 0,
        "tier_M_ortho2align": 0,
        "tier_M_ortho2align_pvalue": pd.NA,
        "tier_M_ortho2align_identity": pd.NA,
        "tier_M_pseudogene_parent": 0,
        "tier_L_liftover": 0,
        "tier_L_liftover_overlap_bp": pd.NA,
        "tier_L_reciprocal_liftover": 0,
        "confidence_tier": "M",
        "evidence_count": 1,
        "provenance_sources": "phasej_synteny_55",
        "notes": (
            "synteny_anchors_human:"
            + merged["upstream_pcg_human"].fillna("").astype(str)
            + "/"
            + merged["downstream_pcg_human"].fillna("").astype(str)
            + ";mouse:"
            + merged["mouse_upstream_pcg"].fillna("").astype(str)
            + "/"
            + merged["mouse_downstream_pcg"].fillna("").astype(str)
        ),
    })

    # Drop unresolvable pairs (no ENSMUSG and no ENSG)
    fully_unresolved = out["mouse_ensembl"].isna() & out["human_ensembl"].isna()
    out = out[~fully_unresolved].copy()
    print(f"[phasej] writing {len(out)} pairs to {args.out}", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    # Summary stats
    print(f"[phasej] biotype distribution (mouse_biotype):", file=sys.stderr)
    print(out["mouse_biotype"].value_counts(dropna=False).to_string(), file=sys.stderr)
    print(f"[phasej] unique mouse symbols: {out['mouse_symbol'].nunique()}",
          file=sys.stderr)
    print(f"[phasej] unique human symbols: {out['human_symbol'].nunique()}",
          file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
