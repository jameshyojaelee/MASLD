#!/usr/bin/env python
"""
pull_guides.py — dynamic Cas13 guide puller for ANY mouse gene target(s).

Use this to add guides for new candidate targets without rebuilding the whole
library. Accepts mouse symbols, mouse Ensembl gene IDs, or human symbols (mapped
to the mouse ortholog via the project ortholog table). Queries the cached parquet
index and applies the same constitutive-first, transcript-spread selection as the
full library build.

Examples:
    # three mouse genes by symbol, 10 guides each
    python pull_guides.py --genes Pnpla3,Col1a1,Meg3 --id-type mouse_symbol --out /tmp/g.csv

    # human symbols (orthology-mapped), 6 guides each
    python pull_guides.py --genes THRB,FASN --id-type human_symbol --n 6

    # from a file (one gene per line), mouse Ensembl IDs
    python pull_guides.py --gene-file new_targets.txt --id-type ensmusg --out new_guides.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guides_config as cfg
import guide_selection as gs


def _read_genes(args) -> list[str]:
    genes: list[str] = []
    if args.genes:
        genes += [g.strip() for g in args.genes.split(",") if g.strip()]
    if args.gene_file:
        with open(args.gene_file) as fh:
            genes += [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not genes:
        sys.exit("ERROR: provide --genes and/or --gene-file")
    # de-dup preserving order
    seen, uniq = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq


def resolve_targets(lazy, genes: list[str], id_type: str) -> pd.DataFrame:
    """Map user queries -> DataFrame[query, gene_id_base] (+ mapping notes)."""
    if id_type == "ensmusg":
        return pd.DataFrame({
            "query": genes,
            "gene_id_base": [g.split(".")[0] for g in genes],
        })
    if id_type == "mouse_symbol":
        m = gs.symbols_to_ids(lazy, genes)
        return pd.DataFrame({
            "query": genes,
            "gene_id_base": [m.get(g.upper()) for g in genes],
        })
    if id_type == "human_symbol":
        mapped = gs.human_to_mouse_ids(genes)
        return mapped[["query", "gene_id_base"]]
    sys.exit(f"ERROR: unknown --id-type {id_type}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", help="comma-separated gene list")
    ap.add_argument("--gene-file", help="file with one gene per line")
    ap.add_argument("--id-type", default="mouse_symbol",
                    choices=["mouse_symbol", "ensmusg", "human_symbol"])
    ap.add_argument("--n", type=int, default=cfg.N_GUIDES_DEFAULT,
                    help=f"guides per gene (default {cfg.N_GUIDES_DEFAULT})")
    ap.add_argument("--release", default=cfg.DEFAULT_RELEASE,
                    choices=list(cfg.RELEASES.keys()))
    ap.add_argument("--out", help="output CSV path (default: stdout)")
    ap.add_argument("--summary", help="optional per-gene coverage CSV path")
    args = ap.parse_args()

    genes = _read_genes(args)
    lazy = gs.load_index(args.release)
    targets = resolve_targets(lazy, genes, args.id_type)

    n_unmapped = targets["gene_id_base"].isna().sum()
    if n_unmapped:
        miss = targets.loc[targets["gene_id_base"].isna(), "query"].tolist()
        print(f"[pull] WARNING: {n_unmapped} unmapped query(ies): {miss}", file=sys.stderr)

    guides_df, summary_df = gs.select_guides(lazy, targets, n=args.n)

    short = summary_df[summary_df["short_flag"] & ~summary_df["missing_flag"]]
    missing = summary_df[summary_df["missing_flag"]]
    print(f"[pull] {len(summary_df)} gene(s): "
          f"{len(guides_df)} guides | "
          f"{len(short)} short (<{args.n}) | {len(missing)} with no guides",
          file=sys.stderr)

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        guides_df.to_csv(args.out, index=False)
        print(f"[pull] wrote {args.out}", file=sys.stderr)
    else:
        guides_df.to_csv(sys.stdout, index=False)
    if args.summary:
        summary_df.to_csv(args.summary, index=False)
        print(f"[pull] wrote {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
