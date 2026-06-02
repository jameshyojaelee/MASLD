"""
ortholog_lookup.py — High-level Python API for the master ortholog table.

Use this for any downstream analysis where you need to map mouse↔human genes.

Example:
    from ortholog_pipeline.ortholog_lookup import (
        get_mouse_ortholog, get_human_ortholog, OrthologLookup
    )

    # Single-gene lookup (returns best ortholog at min_tier)
    get_mouse_ortholog("MEG3", min_tier="M")
    # → {"mouse_symbol": "Meg3", "mouse_ensembl": "ENSMUSG00000021268",
    #    "confidence_tier": "M", "provenance_sources": "blast"}

    get_human_ortholog("Mir122", min_tier="H")
    # → {"human_symbol": "MIR122", ..., "confidence_tier": "H",
    #    "provenance_sources": "mirbase"}

    # Batch lookup
    lookup = OrthologLookup()
    df = lookup.map_genes(["MEG3", "HOTAIR", "MIR122", "ALB"], "human_to_mouse")
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Iterable, Literal

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MASTER = PROJECT / "data/external/orthologs/master_ortholog_table.tsv.gz"

TIER_RANK = {"H": 0, "M": 1, "L": 2}
TierLevel = Literal["H", "M", "L"]

# miRNA single-form -> paralog fallback. HGNC has no bare MIRLET7A gene
# (only MIRLET7A1/2/3); MGI has no bare Mirlet7a (only Mirlet7a-1, -2). When
# a user queries the single-form symbol, the L2_mirbase layer maps the
# miRBase mature family to the paralog-specific GENCODE rows. To make the
# single-form query resolve, on exact-match miss we try matching the same
# stem with an appended paralog suffix (\d+ for HGNC concat, -\d+ for MGI
# hyphenated). Returns the FIRST matching row (smallest paralog id ranks
# first because lexicographic on "MIRLET7A1" < "MIRLET7A2").
_MIRNA_QUERY_PREFIX_RE = re.compile(r"^(MIR|Mir|mir)")


def _strip_version(s):
    if pd.isna(s):
        return s
    return str(s).split(".")[0]


class OrthologLookup:
    """Cached, vectorized lookup against the master ortholog table.

    Instantiate once, then call .map_genes() or .get_one() for any number
    of lookups. The table is loaded lazily on first use.
    """

    def __init__(self, master_path: Path = MASTER,
                 exclude_v49_absent: bool = True):
        self._path = master_path
        self._exclude_v49_absent = exclude_v49_absent
        self._df: Optional[pd.DataFrame] = None

    @property
    def df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self._path, sep="\t", low_memory=False,
                                   dtype=str)
            # Strip versions defensively
            for c in ("mouse_ensembl", "human_ensembl"):
                if c in self._df.columns:
                    self._df[c] = self._df[c].map(_strip_version)
            # Exclude v49-absent human genes by default (P1-13 fix).
            # 27 human Ensembl IDs are absent from GENCODE v49 metadata
            # (retired/merged IDs). Rows are flagged, not deleted, for
            # provenance. Set exclude_v49_absent=False to include them.
            if self._exclude_v49_absent and "v49_absent" in self._df.columns:
                n_before = len(self._df)
                self._df = self._df[
                    pd.to_numeric(self._df["v49_absent"],
                                  errors="coerce").fillna(0) == 0
                ].reset_index(drop=True)
                n_dropped = n_before - len(self._df)
                if n_dropped > 0:
                    import sys
                    print(f"[ortholog_lookup] Excluded {n_dropped} rows with "
                          f"v49_absent=1 (GENCODE v49-absent human genes)",
                          file=sys.stderr)
            # Tier rank for fast filtering
            self._df["_tier_rank"] = self._df["confidence_tier"].map(TIER_RANK)
        return self._df

    def _filter_tier(self, df: pd.DataFrame, min_tier: TierLevel) -> pd.DataFrame:
        return df[df["_tier_rank"] <= TIER_RANK[min_tier]]

    def get_one(self, query: str, direction: str = "human_to_mouse",
                min_tier: TierLevel = "M",
                by: str = "symbol") -> Optional[dict]:
        """Return best ortholog dict for a single query gene.

        Args:
          query: gene symbol or Ensembl ID
          direction: "human_to_mouse" or "mouse_to_human"
          min_tier: H (strict), M (default), or L (most permissive)
          by: "symbol" (case-insensitive) or "ensembl"

        Returns:
          dict with keys mouse_*, human_*, confidence_tier, provenance_sources,
          or None if no match at requested tier.
        """
        d = self.df
        d = self._filter_tier(d, min_tier)
        if direction == "human_to_mouse":
            src_sym, src_id = "human_symbol", "human_ensembl"
        else:
            src_sym, src_id = "mouse_symbol", "mouse_ensembl"
        if by == "symbol":
            q = query.upper()
            hits = d[d[src_sym].astype(str).str.upper() == q]
            # miRNA single-form fallback: MIRLET7A -> match MIRLET7A1 /
            # MIRLET7A-1 / etc. Only triggered on exact-match miss for queries
            # that look like miRNA symbols (start with MIR / Mir / mir).
            if hits.empty and _MIRNA_QUERY_PREFIX_RE.match(query):
                # Pattern: query optionally followed by a hyphen and digits,
                # or directly followed by digits. Anchored at start and end.
                pat = re.compile(rf"^{re.escape(q)}-?\d+$")
                mask = d[src_sym].astype(str).str.upper().str.match(pat)
                hits = d[mask.fillna(False)]
        elif by == "ensembl":
            q = _strip_version(query)
            hits = d[d[src_id] == q]
        else:
            raise ValueError(f"by must be 'symbol' or 'ensembl', got {by}")
        if hits.empty:
            return None
        # Sort by tier rank (best tier first), then by BLAST bitscore (higher
        # = more sequence similarity), then by evidence count
        if "blast_bitscore" in hits.columns:
            hits = hits.copy()
            hits["_bitscore_num"] = pd.to_numeric(hits["blast_bitscore"],
                                                   errors="coerce").fillna(0)
        else:
            hits = hits.copy()
            hits["_bitscore_num"] = 0
        # Use whichever evidence_count column exists (renamed in v4 rebuild)
        _ev_col = ("independent_evidence_count"
                    if "independent_evidence_count" in hits.columns
                    else "raw_evidence_count"
                    if "raw_evidence_count" in hits.columns
                    else "evidence_count"
                    if "evidence_count" in hits.columns
                    else None)
        if _ev_col:
            sort_cols = ["_tier_rank", "_bitscore_num", _ev_col, src_sym]
        else:
            sort_cols = ["_tier_rank", "_bitscore_num", src_sym]
        hits = hits.sort_values(
            sort_cols,
            ascending=[True, False] + [False] * (len(sort_cols) - 3) + [True],
        )
        r = hits.iloc[0]
        ev = int(r.get("independent_evidence_count")
                 or r.get("raw_evidence_count")
                 or r.get("evidence_count") or 0)
        return {
            "mouse_symbol": r.get("mouse_symbol"),
            "mouse_ensembl": r.get("mouse_ensembl"),
            "human_symbol": r.get("human_symbol"),
            "human_ensembl": r.get("human_ensembl"),
            "confidence_tier": r.get("confidence_tier"),
            "provenance_sources": r.get("provenance_sources"),
            "evidence_count": ev,
        }

    def map_genes(self, genes: Iterable[str],
                  direction: str = "human_to_mouse",
                  min_tier: TierLevel = "M",
                  by: str = "symbol",
                  keep_all_orthologs: bool = False) -> pd.DataFrame:
        """Batch map a list of query genes.

        Args:
          genes: iterable of gene symbols / Ensembl IDs
          direction, min_tier, by: as in get_one()
          keep_all_orthologs: if True, returns ALL ortholog rows per query
            (useful for many-many cases like miRNA families). If False
            (default), returns only the best ortholog per query.

        Returns:
          DataFrame with columns query, mouse_symbol, mouse_ensembl,
          human_symbol, human_ensembl, confidence_tier, provenance_sources,
          evidence_count. Rows for queries with no match show NaN.
        """
        d = self.df
        d = self._filter_tier(d, min_tier)
        if direction == "human_to_mouse":
            src_sym, src_id = "human_symbol", "human_ensembl"
        else:
            src_sym, src_id = "mouse_symbol", "mouse_ensembl"
        genes_list = list(genes)
        if by == "symbol":
            queries = [g.upper() for g in genes_list]
            d_search = d.copy()
            d_search["_match_key"] = d_search[src_sym].astype(str).str.upper()
        else:
            queries = [_strip_version(g) for g in genes_list]
            d_search = d.copy()
            d_search["_match_key"] = d_search[src_id]
        qdf = pd.DataFrame({"query": genes_list, "_match_key": queries})
        merged = qdf.merge(d_search, on="_match_key", how="left")
        # miRNA single-form fallback: for symbol queries that returned no
        # match AND look like miRNA family symbols, append paralog suffix
        # variants and try the merge again. Triggers only on misses.
        if by == "symbol":
            miss_mask = merged[src_sym].isna()
            miss_queries = (merged.loc[miss_mask, "query"].dropna()
                            .drop_duplicates().tolist())
            mirna_misses = [
                q for q in miss_queries
                if isinstance(q, str) and _MIRNA_QUERY_PREFIX_RE.match(q)
            ]
            if mirna_misses:
                # Pre-uppercase symbol column once
                sym_upper = d_search["_match_key"]
                extra_rows = []
                for q in mirna_misses:
                    pat = re.compile(rf"^{re.escape(q.upper())}-?\d+$")
                    hit_mask = sym_upper.astype(str).str.match(pat).fillna(False)
                    hits = d_search[hit_mask]
                    if not hits.empty:
                        # Sort by tier rank then symbol (smallest paralog first)
                        hits = hits.copy()
                        hits["query"] = q
                        # Place row IDs first deterministically by sorting on the
                        # uppercase symbol (alpha-numeric -> "MIRLET7A1" < "...A2")
                        hits = hits.sort_values(
                            ["_tier_rank", "_match_key"],
                            ascending=[True, True],
                        )
                        extra_rows.append(hits)
                if extra_rows:
                    extra = pd.concat(extra_rows, ignore_index=True)
                    # Replace the miss rows with the fallback hits.
                    merged = pd.concat(
                        [merged[~merged["query"].isin(mirna_misses)], extra],
                        ignore_index=True,
                        sort=False,
                    )
        if "blast_bitscore" in merged.columns:
            merged["_bitscore_num"] = pd.to_numeric(
                merged["blast_bitscore"], errors="coerce").fillna(0)
        else:
            merged["_bitscore_num"] = 0
        # Use whichever evidence_count column exists (renamed in v4 rebuild)
        _ev_col = ("independent_evidence_count"
                    if "independent_evidence_count" in merged.columns
                    else "raw_evidence_count"
                    if "raw_evidence_count" in merged.columns
                    else "evidence_count"
                    if "evidence_count" in merged.columns
                    else None)
        if not keep_all_orthologs:
            sort_cols = ["query", "_tier_rank", "_bitscore_num"]
            asc = [True, True, False]
            if _ev_col:
                sort_cols.append(_ev_col)
                asc.append(False)
            merged = merged.sort_values(
                sort_cols, ascending=asc
            ).drop_duplicates(subset=["query"], keep="first")
        # Normalize evidence_count column name for output
        if _ev_col and _ev_col != "evidence_count":
            merged = merged.rename(columns={_ev_col: "evidence_count"})
        elif _ev_col is None:
            merged["evidence_count"] = 0
        keep = ["query", "mouse_symbol", "mouse_ensembl",
                "human_symbol", "human_ensembl",
                "confidence_tier", "provenance_sources",
                "evidence_count"]
        return merged[keep].reset_index(drop=True)

    def summary(self) -> dict:
        """Stats on the master table."""
        d = self.df
        return {
            "total_pairs": len(d),
            "tier_H": (d["confidence_tier"] == "H").sum(),
            "tier_M": (d["confidence_tier"] == "M").sum(),
            "tier_L": (d["confidence_tier"] == "L").sum(),
            "unique_mouse_ensembl": d["mouse_ensembl"].nunique(),
            "unique_human_ensembl": d["human_ensembl"].nunique(),
        }


# Module-level convenience functions
_DEFAULT_LOOKUP: Optional[OrthologLookup] = None


def _default() -> OrthologLookup:
    global _DEFAULT_LOOKUP
    if _DEFAULT_LOOKUP is None:
        _DEFAULT_LOOKUP = OrthologLookup()
    return _DEFAULT_LOOKUP


def get_mouse_ortholog(human_gene: str, min_tier: TierLevel = "M",
                       by: str = "symbol") -> Optional[dict]:
    """Convenience: find mouse ortholog of a human gene."""
    return _default().get_one(human_gene, "human_to_mouse", min_tier, by)


def get_human_ortholog(mouse_gene: str, min_tier: TierLevel = "M",
                       by: str = "symbol") -> Optional[dict]:
    """Convenience: find human ortholog of a mouse gene."""
    return _default().get_one(mouse_gene, "mouse_to_human", min_tier, by)


def map_human_to_mouse(genes: Iterable[str], min_tier: TierLevel = "M",
                       **kwargs) -> pd.DataFrame:
    """Convenience: batch map human → mouse."""
    return _default().map_genes(genes, "human_to_mouse", min_tier, **kwargs)


def map_mouse_to_human(genes: Iterable[str], min_tier: TierLevel = "M",
                       **kwargs) -> pd.DataFrame:
    """Convenience: batch map mouse → human."""
    return _default().map_genes(genes, "mouse_to_human", min_tier, **kwargs)


if __name__ == "__main__":
    # CLI smoke test
    lookup = OrthologLookup()
    print("=== Master ortholog table summary ===")
    for k, v in lookup.summary().items():
        print(f"  {k}: {v:,}" if isinstance(v, int) else f"  {k}: {v}")
    print()
    print("=== Spot-check: human → mouse @ tier M ===")
    for g in ["MEG3", "HOTAIR", "MIR122", "ALB", "TP53", "BRCA1",
              "PNPLA3", "HSD17B13", "MALAT1", "NEAT1", "XIST", "H19"]:
        r = get_mouse_ortholog(g, min_tier="M")
        if r:
            print(f"  {g:10s} → {r['mouse_symbol']!s:20s} {r['mouse_ensembl']!s:20s} "
                  f"[{r['confidence_tier']}] {r['provenance_sources']}")
        else:
            r_l = get_mouse_ortholog(g, min_tier="L")
            if r_l:
                print(f"  {g:10s} → {r_l['mouse_symbol']!s:20s} (only available at Tier L)")
            else:
                print(f"  {g:10s} → no ortholog")
