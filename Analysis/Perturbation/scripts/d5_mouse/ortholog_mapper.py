"""D5 ortholog mapper — bridges human atlas hits ↔ mouse gene symbols.

Canonical table (preferred):
  data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz
  cols: human_gene_symbol, human_gene_ensembl,
        mouse_gene_symbol, mouse_gene_ensembl,
        ortholog_type, conservation_score

Schema-agnostic: also accepts the legacy Ensembl-only table
  (mouse_ensembl_gene_id, human_ensembl_gene_id, orthology_type) — symbol
  lookups will return [] in that mode but ortholog_type joins still work.

Provides:
  - human_to_mouse(symbol) -> list[str]
  - mouse_to_human(symbol) -> list[str]
  - map_hit_table(df, column) -> df with mouse_gene column joined

Subagent: d5-ortholog-mapper.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ORTHO_CANDIDATES = [
    PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz",
    PROJECT_ROOT / "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz",
    PROJECT_ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz",
    PROJECT_ROOT / "data/external/orthologs/mouse_human_orthologs.tsv.gz",
]


def _resolve_ortho_path() -> Path:
    for p in ORTHO_CANDIDATES:
        if p.exists():
            return p
    # Skeleton fallback so the runner can still produce valid JSON.
    return ORTHO_CANDIDATES[0]


ORTHO_PATH = _resolve_ortho_path()

_cache: dict[str, pd.DataFrame] = {}


def _load() -> pd.DataFrame:
    if "df" not in _cache:
        if not ORTHO_PATH.exists():
            _cache["df"] = pd.DataFrame(
                columns=["human_symbol", "mouse_symbol", "ortholog_type"]
            )
        else:
            _cache["df"] = pd.read_csv(ORTHO_PATH, sep="\t")
    return _cache["df"]


def _symbol_cols(df: pd.DataFrame) -> tuple[str | None, str | None]:
    """Return (human_col, mouse_col) symbol column names, or (None, None) if absent."""
    h = next((c for c in (
        "human_gene_symbol", "human_symbol", "hgnc_symbol", "external_gene_name_h"
    ) if c in df.columns), None)
    m = next((c for c in (
        "mouse_gene_symbol", "mouse_symbol", "mgi_symbol", "external_gene_name_m"
    ) if c in df.columns), None)
    return h, m


def _type_col(df: pd.DataFrame) -> str | None:
    return next((c for c in ("ortholog_type", "orthology_type") if c in df.columns), None)


def human_to_mouse(symbol: str) -> list[str]:
    df = _load()
    h, m = _symbol_cols(df)
    if not h or not m:
        return []
    sub = df[df[h].astype(str).str.upper() == symbol.upper()]
    return sorted(sub[m].dropna().unique().tolist())


def mouse_to_human(symbol: str) -> list[str]:
    df = _load()
    h, m = _symbol_cols(df)
    if not h or not m:
        return []
    sub = df[df[m] == symbol]
    return sorted(sub[h].dropna().unique().tolist())


def map_hit_table(hits: pd.DataFrame, *, column: str = "human_gene") -> pd.DataFrame:
    df = _load()
    h, m = _symbol_cols(df)
    if not h or not m:
        # Pass through unchanged when symbol columns are missing (skeleton mode)
        return hits.copy()
    tcol = _type_col(df) or "ortholog_type"
    cols = [h, m] + ([tcol] if tcol in df.columns else [])
    sub = df[cols].rename(columns={h: column, tcol: "ortholog_type"})
    merged = hits.merge(sub, on=column, how="left")
    if m in merged.columns:
        if "mouse_gene" in merged.columns and m != "mouse_gene":
            merged["mouse_gene"] = merged[m].fillna(merged["mouse_gene"])
        else:
            merged["mouse_gene"] = merged[m]
    return merged
