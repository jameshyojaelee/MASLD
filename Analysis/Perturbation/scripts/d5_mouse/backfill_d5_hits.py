"""Backfill the empty mouse_gene (symbol) column in d5_mouse_orthologs.csv
using the new symbol-based ortholog table at
data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz.

Match priority per row:
  1. (human_gene, mouse_ensembl) join — strictest, preserves the original pair
  2. mouse_ensembl join — recovers pairs where the human symbol drifted
  3. human_gene join — last resort

Subagent: d5-ortholog-mapper.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HITS_CSV = PROJECT_ROOT / "Analysis/Perturbation/data/hits/d5_mouse_orthologs.csv"
ORTHO_TSV = PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"


def main() -> None:
    hits = pd.read_csv(HITS_CSV)
    ortho = pd.read_csv(ORTHO_TSV, sep="\t")

    n = len(hits)
    print(f"hits rows: {n:,}")
    print(f"ortho rows: {len(ortho):,}")
    print(f"hits with empty mouse_gene BEFORE: {(hits['mouse_gene'].isna() | (hits['mouse_gene']=='')).sum()}")

    # Build lookup maps
    pair_map = (
        ortho.dropna(subset=["mouse_gene_symbol"])
             .set_index(["human_gene_symbol", "mouse_gene_ensembl"])["mouse_gene_symbol"]
             .to_dict()
    )
    mouse_map = (
        ortho.dropna(subset=["mouse_gene_symbol"])
             .drop_duplicates("mouse_gene_ensembl")
             .set_index("mouse_gene_ensembl")["mouse_gene_symbol"]
             .to_dict()
    )
    human_map = (
        ortho.dropna(subset=["mouse_gene_symbol"])
             .drop_duplicates("human_gene_symbol")
             .set_index("human_gene_symbol")["mouse_gene_symbol"]
             .to_dict()
    )

    def lookup(row) -> tuple[str | None, str]:
        hg, me = row["human_gene"], row["mouse_ensembl"]
        if pd.notna(hg) and pd.notna(me):
            v = pair_map.get((hg, me))
            if v:
                return v, "pair"
        if pd.notna(me):
            v = mouse_map.get(me)
            if v:
                return v, "mouse_ensembl"
        if pd.notna(hg):
            v = human_map.get(hg)
            if v:
                return v, "human_symbol"
        return None, "unmapped"

    mapped, sources = [], []
    for _, row in hits.iterrows():
        v, src = lookup(row)
        mapped.append(v)
        sources.append(src)

    hits["mouse_gene"] = mapped
    hits["mouse_gene_source"] = sources

    print(f"hits with empty mouse_gene AFTER: {hits['mouse_gene'].isna().sum()}")
    print("source breakdown:")
    print(pd.Series(sources).value_counts().to_string())

    hits.to_csv(HITS_CSV, index=False)
    print(f"Wrote: {HITS_CSV}")


if __name__ == "__main__":
    main()
