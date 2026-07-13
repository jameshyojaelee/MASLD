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
def _txpos(r: dict) -> dict:
    """Parse a guide's tx_id_pos ('TX:pos|TX:pos|...') into {transcript_id: pos},
    cached on the row. Used for transcript-AWARE overlap detection: two guides can
    only overlap if they target a COMMON transcript."""
    tp = r.get("_txpos")
    if tp is None:
        tp = {}
        for seg in (r.get("tx_id_pos") or "").split("|"):
            if ":" in seg:
                t, p = seg.rsplit(":", 1)
                try:
                    tp[t] = int(p)
                except ValueError:
                    pass
        r["_txpos"] = tp
    return tp


def _far_enough(r: dict, selected: list[dict], min_spacing: int) -> bool:
    """True iff guide r does NOT overlap any already-selected guide. Overlap is
    TRANSCRIPT-AWARE: two guides conflict only when they target a COMMON transcript
    and their positions ON THAT TRANSCRIPT are < min_spacing apart. Guides on
    disjoint isoforms cannot overlap (different sequence) and never conflict.

    This fixes the multi-isoform blind spot where the single `position` (the
    FIRST-listed transcript only) made two guides at the SAME site look far apart
    when they happened to list different transcripts first (e.g. Adora1 g01/g03,
    2 nt apart on the shared transcript but 779 'apart' by first-segment position).
    Falls back to the representative `position` axis only when a pair lacks
    transcript coordinates."""
    rt = _txpos(r)
    for s in selected:
        st = _txpos(s)
        if rt and st:
            if any(abs(rt[t] - st[t]) < min_spacing for t in (rt.keys() & st.keys())):
                return False
            # disjoint transcripts -> no possible overlap -> not a conflict
        else:                                   # missing tx coords -> safe fallback
            p, q = r.get("position"), s.get("position")
            if p is not None and q is not None and abs(p - q) < min_spacing:
                return False
    return True


def _spread_select(candidate: list[dict], n: int,
                   min_spacing: int = cfg.MIN_GUIDE_SPACING) -> list[dict]:
    """Pick up to N guides spread along the transcript, each >= min_spacing bp from
    every other selected guide (no overlap when min_spacing >= GUIDE_LEN).

    Transcript-wide spread is preserved by position binning; the min-distance
    constraint is enforced at EVERY step (both the per-bin pick and the score-ordered
    fill), which closes the prior overlap leak. If a gene has fewer than N
    non-overlapping guides, fewer are returned (clean) -- an overlapping guide is
    never added (non-overlap takes priority over hitting N)."""
    if not candidate:
        return []
    have_pos = [r["position"] for r in candidate if r["position"] is not None]
    selected: list[dict] = []
    # --- Binned pass: one non-overlapping guide per transcript bin (spread) ---------
    if len(set(have_pos)) >= n:
        pmin, pmax = min(have_pos), max(have_pos)
        width = max((pmax - pmin) / n, 1e-9)
        bins: dict[int, list[dict]] = defaultdict(list)
        for r in candidate:
            p = r["position"]
            b = n if p is None else min(int((p - pmin) / width), n - 1)  # n = no-pos bucket
            bins[b].append(r)
        for b in range(n):
            for r in sorted(bins.get(b, []), key=lambda r: -r["combined_score"]):
                if _far_enough(r, selected, min_spacing):   # skip if it overlaps a pick
                    selected.append(r)
                    break
            if len(selected) >= n:
                break
    # --- Greedy fill: best score first, non-overlapping only ------------------------
    # Covers small / low-position-diversity pools and tops up bins that were skipped.
    if len(selected) < n:
        sel_ids = {id(s) for s in selected}
        for r in sorted(candidate, key=lambda r: -r["combined_score"]):
            if id(r) in sel_ids:
                continue
            if _far_enough(r, selected, min_spacing):
                selected.append(r)
                sel_ids.add(id(r))
                if len(selected) >= n:
                    break
    return selected[:n]


def _has_cds(r: dict) -> bool:
    """Guide overlaps the CODING SEQUENCE (region 'CDS', '5'UTR|CDS', '3'UTR|CDS')."""
    return "CDS" in str(r.get("region") or "").upper()


def _constitutive_candidate(pool: list[dict], n: int) -> list[dict]:
    """Constitutive-first candidate set: take whole isoform-coverage tiers
    (most isoforms first) until >= n candidates accumulate."""
    if not pool:
        return []
    iso_tiers = sorted({r["n_isoforms_targeted"] for r in pool}, reverse=True)
    cand: list[dict] = []
    for iso in iso_tiers:
        cand += [r for r in pool if r["n_isoforms_targeted"] == iso]
        if len(cand) >= n:
            break
    return cand


def _select_for_gene(rows: list[dict], n: int) -> list[dict]:
    """CDS-first, constitutive-first candidate set, then transcript-spread top-N.

    Region preference (RfxCas13d targets the mature mRNA, but CDS-targeting is the
    validated knockdown design): prefer guides overlapping the CODING SEQUENCE;
    only fall back to UTR-only guides when a gene has < n CDS candidates. lncRNA
    genes have no CDS (region 'lncRNA') so they pass straight through the fallback
    branch unchanged. Within the chosen pool: constitutive-first (max isoform
    coverage) -> transcript-spread -> combined_score.
    """
    if not rows:
        return []
    cds_rows = [r for r in rows if _has_cds(r)]
    utr_rows = [r for r in rows if not _has_cds(r)]
    # Progressive widening so the >=MIN_GUIDE_SPACING spacing constraint doesn't starve
    # genes that have plenty of guides at lower isoform tiers: try the constitutive-first
    # CDS set, then ALL CDS rows, then CDS + UTR -- stopping as soon as N non-overlapping
    # guides are selectable. Constitutive-first is preserved as the FIRST attempt.
    candidate = _constitutive_candidate(cds_rows, n)            # prefer CDS, constitutive
    selected = _spread_select(candidate, n)
    if len(selected) < n and len(cds_rows) > len(candidate):    # widen to all CDS rows
        selected = _spread_select(cds_rows, n)
    if len(selected) < n and utr_rows:                          # UTR fallback
        selected = _spread_select(cds_rows + utr_rows, n)
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
