#!/usr/bin/env python
"""
guide_selection.py — core query + per-gene guide selection for the Cas13 library.

Shared by build_library_guides.py (full roster) and pull_guides.py (ad-hoc
genes). Reads the cached parquet index (see build_guide_index.py), never the
3.8 GB CSV directly.

Selection objective = gene-level knockdown maximising isoform coverage:
  1. Among a gene's eligible guides, prefer those hitting the MOST isoforms
     (constitutive); descend isoform-coverage tiers only until >= N candidates.
  2. Within that candidate set, spread guides along the transcript (bin by
     position, take the best combined_score per bin) so guides don't cluster.
  3. Tie-break / fill by combined_score (TIGER/Cas13Design composite).

The upstream pool is already post-basic-criteria; we additionally re-assert the
two cheap guards (no homopolymer run; TIGER>=0.75 OR Cas13Design>=0.75) so the
output is guaranteed clean and the pass-rate is logged. Off-target cleanliness is
trusted from the source (mismatch counts the on-target in the 0-mm slot).
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
import guides_config as cfg

# ortholog_lookup lives one level up under scripts/
sys.path.insert(0, str(cfg.CAS13LIB / "scripts"))

_POS_RE = re.compile(r":(\d+)\s*$")


# ---------------------------------------------------------------------------
# Index access
# ---------------------------------------------------------------------------
def load_index(release: str = cfg.DEFAULT_RELEASE) -> pl.LazyFrame:
    path = cfg.index_path(release)
    if not path.exists():
        sys.exit(
            f"ERROR: guide index not found: {path}\n"
            f"Build it first:  python build_guide_index.py --release {release}"
        )
    return pl.scan_parquet(path)


def symbols_to_ids(lazy: pl.LazyFrame, symbols: Iterable[str]) -> dict:
    """Resolve mouse gene symbols -> gene_id_base via one indexed scan.

    Returns {SYMBOL_UPPER: gene_id_base}. If a symbol maps to >1 gene_id (rare),
    the first by gene_id_base is kept.
    """
    syms = [str(s).upper() for s in symbols]
    hit = (
        lazy.with_columns(pl.col("symbol").str.to_uppercase().alias("_su"))
        .filter(pl.col("_su").is_in(syms))
        .select(["_su", "gene_id_base"])
        .unique()
        .sort(["_su", "gene_id_base"])
        .collect()
    )
    out: dict[str, str] = {}
    for su, gid in hit.iter_rows():
        out.setdefault(su, gid)
    return out


def human_to_mouse_ids(human_symbols: Iterable[str],
                       min_tier: str = cfg.ORTHOLOG_MIN_TIER) -> pd.DataFrame:
    """Map human symbols -> mouse gene_id_base via the project ortholog table.

    Returns DataFrame[query, gene_id_base, mouse_symbol, mapped].
    """
    from ortholog_pipeline.ortholog_lookup import map_human_to_mouse  # noqa
    m = map_human_to_mouse(list(human_symbols), min_tier=min_tier)
    m = m.rename(columns={"mouse_ensembl": "gene_id_base"})
    m["mapped"] = m["gene_id_base"].notna()
    return m[["query", "gene_id_base", "mouse_symbol", "mapped"]]


# ---------------------------------------------------------------------------
# Row helpers
# ---------------------------------------------------------------------------
def _n_isoforms(tx_id_set: Optional[str]) -> int:
    if not tx_id_set:
        return 0
    return tx_id_set.count("|") + 1


def _position(tx_id_pos: Optional[str]) -> Optional[int]:
    """Representative transcript position = trailing int of the first segment."""
    if not tx_id_pos:
        return None
    first = tx_id_pos.split("|", 1)[0]
    m = _POS_RE.search(first)
    return int(m.group(1)) if m else None


def _has_homopolymer(seq: str) -> bool:
    return any(h in seq for h in cfg.HOMOPOLYMERS)


def _passes_score(tiger: Optional[float], cas13: Optional[float]) -> bool:
    t = tiger if tiger is not None else 0.0
    c = cas13 if cas13 is not None else 0.0
    return (t >= cfg.TIGER_MIN) or (c >= cfg.CAS13_MIN)


# ---------------------------------------------------------------------------
# Per-gene selection
# ---------------------------------------------------------------------------
def _spread_select(candidate: list[dict], n: int) -> list[dict]:
    """Pick N guides spread along the transcript, best combined_score per region."""
    if len(candidate) <= n:
        return list(candidate)
    have_pos = [r["position"] for r in candidate if r["position"] is not None]
    if len(set(have_pos)) < n:
        return sorted(candidate, key=lambda r: -r["combined_score"])[:n]
    pmin, pmax = min(have_pos), max(have_pos)
    width = max((pmax - pmin) / n, 1e-9)
    bins: dict[int, list[dict]] = defaultdict(list)
    for r in candidate:
        p = r["position"]
        b = n if p is None else min(int((p - pmin) / width), n - 1)  # n = no-pos bucket
        bins[b].append(r)
    selected: list[dict] = []
    for b in range(n):
        if bins.get(b):
            best = max(bins[b], key=lambda r: r["combined_score"])
            selected.append(best)
            bins[b].remove(best)
    if len(selected) < n:
        leftover = [r for b in bins for r in bins[b]]
        leftover.sort(key=lambda r: -r["combined_score"])
        selected += leftover[: (n - len(selected))]
    return selected[:n]


def _select_for_gene(rows: list[dict], n: int) -> list[dict]:
    """Constitutive-first candidate set, then transcript-spread top-N."""
    if not rows:
        return []
    iso_tiers = sorted({r["n_isoforms_targeted"] for r in rows}, reverse=True)
    candidate: list[dict] = []
    for iso in iso_tiers:
        candidate += [r for r in rows if r["n_isoforms_targeted"] == iso]
        if len(candidate) >= n:
            break
    selected = _spread_select(candidate, n)
    selected.sort(
        key=lambda r: (r["position"] if r["position"] is not None else 1 << 60,
                       -r["combined_score"])
    )
    return selected


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def select_guides(lazy: pl.LazyFrame,
                  targets: pd.DataFrame,
                  n: int = cfg.N_GUIDES_DEFAULT) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select up to N guides per target gene.

    Args:
      lazy:    parquet index LazyFrame (load_index()).
      targets: DataFrame with columns [query, gene_id_base] (one row per gene).
               gene_id_base may be NA (unmapped) -> reported as missing.
      n:       guides per gene.

    Returns:
      guides_df  — one row per selected guide (schema below).
      summary_df — one row per target gene (coverage / shortfall report).
    """
    targets = targets.dropna(subset=["gene_id_base"]).drop_duplicates("gene_id_base")
    want_ids = targets["gene_id_base"].tolist()
    q_by_id = dict(zip(targets["gene_id_base"], targets["query"]))

    # One predicate-pushdown scan for all wanted genes.
    pool = (
        lazy.filter(pl.col("gene_id_base").is_in(want_ids))
        .collect()
    )

    # Group rows per gene.
    groups: dict[str, list[dict]] = defaultdict(list)
    n_pool_pass = 0
    n_pool_total = 0
    for r in pool.iter_rows(named=True):
        n_pool_total += 1
        seq = r["guide_seq"] or ""
        # Re-assert the two cheap basic-criteria guards.
        if _has_homopolymer(seq) or not _passes_score(r["tiger_score"], r["cas13_score"]):
            continue
        n_pool_pass += 1
        r["n_isoforms_targeted"] = _n_isoforms(r["tx_id_set"])
        r["position"] = _position(r["tx_id_pos"])
        groups[r["gene_id_base"]].append(r)

    if n_pool_total:
        print(f"[select] guard pass-rate: {n_pool_pass:,}/{n_pool_total:,} "
              f"({100*n_pool_pass/n_pool_total:.1f}%) of fetched pool rows")

    guide_records: list[dict] = []
    summary_records: list[dict] = []
    for gid in want_ids:
        rows = groups.get(gid, [])
        query = q_by_id.get(gid, gid)
        n_avail = len(rows)
        selected = _select_for_gene(rows, n)
        symbol = selected[0]["symbol"] if selected else (rows[0]["symbol"] if rows else None)
        biotype = selected[0]["biotype"] if selected else (rows[0]["biotype"] if rows else None)
        for rank, r in enumerate(selected, start=1):
            guide_records.append({
                "guide_id": f"{r['symbol']}_g{rank:02d}",
                "gene_id_mouse": gid,
                "gene_symbol_mouse": r["symbol"],
                "biotype": r["biotype"],
                "guide_seq": r["guide_seq"],
                "target_seq": r["target_seq"],
                "region": r["region"],
                "tiger_score": r["tiger_score"],
                "cas13_score": r["cas13_score"],
                "combined_score": r["combined_score"],
                "n_isoforms_targeted": r["n_isoforms_targeted"],
                "single_isoform": r["n_isoforms_targeted"] == 1,
                "tx_id_set": r["tx_id_set"],
                "tx_id_pos": r["tx_id_pos"],
                "position": r["position"],
                "n_target": r["n_target"],
                "any_indel": r["any_indel"],
                "mismatch": r["mismatch"],
                "rank_within_gene": rank,
                "n_available_pool": n_avail,
            })
        summary_records.append({
            "query": query,
            "gene_id_mouse": gid,
            "gene_symbol_mouse": symbol,
            "biotype": biotype,
            "n_available_pool": n_avail,
            "n_selected": len(selected),
            "short_flag": len(selected) < n,
            "missing_flag": n_avail == 0,
        })

    guides_df = pd.DataFrame(guide_records)
    summary_df = pd.DataFrame(summary_records)
    return guides_df, summary_df
