"""Cached data loaders for the MASLD Atlas Streamlit app.

All loaders use ``st.cache_data`` so cold-boot on Hugging Face Spaces stays
fast. Large per-gene JSON files use a small LRU cache so memory stays bounded
(we cannot hold all 34K gene profiles in RAM).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from . import paths


def _read_json(path: Path) -> Any:
    if not path.exists():
        raise FileNotFoundError(
            f"Required data file is missing: {path}. "
            "Set MASLD_DATA_DIR or run `make build-parquet` from the atlas root."
        )
    with path.open("r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Summary + featured
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def load_atlas_summary() -> dict:
    return _read_json(paths.summary_json())


@st.cache_data(ttl=3600, show_spinner=False)
def load_featured_genes() -> list[dict]:
    data = _read_json(paths.featured_genes())
    assert isinstance(data, list)
    return data


# ---------------------------------------------------------------------------
# Gene index (compact) — used by the Explorer page
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def load_gene_index() -> pd.DataFrame:
    """Flatten ``gene_index.json`` into a DataFrame.

    gene_index.json uses the M1-M7 aligned evidence key scheme; each row has
    an optional ``evidence`` dict. We materialise the eight evidence columns
    (``s1_human`` ... ``s8_proteomics``) as scalar float columns for fast
    filtering / sorting.
    """
    entries = _read_json(paths.gene_index())
    df = pd.DataFrame(entries)

    evidence_cols = [
        "s1_human",
        "s2_mouse",
        "s3_genetic",
        "s4_essential",
        "s5_epigenomic",
        "s6_spatial",
        "s7_singlecell",
        "s8_proteomics",
    ]

    ev = pd.json_normalize(df.get("evidence", pd.Series([{}] * len(df))).fillna({}))
    for col in evidence_cols:
        if col not in ev.columns:
            ev[col] = 0.0
        ev[col] = pd.to_numeric(ev[col], errors="coerce").fillna(0.0)

    df = pd.concat([df.drop(columns=["evidence"], errors="ignore"), ev[evidence_cols]], axis=1)

    # Derived convenience flags
    if "dream_padj" in df.columns and "dream_logfc" in df.columns:
        df["is_deg"] = (df["dream_padj"] < 0.05) & (df["dream_logfc"].abs() > 0.3)
    else:
        df["is_deg"] = False

    df["coloc_linked"] = df["s3_genetic"].astype(float) > 0.5
    return df


# ---------------------------------------------------------------------------
# Per-gene profile (on-demand; small LRU in the cache)
# ---------------------------------------------------------------------------

@st.cache_data(ttl=300, max_entries=500, show_spinner=False)
def load_gene_profile(symbol: str) -> dict | None:
    path = paths.gene_profile_path(symbol)
    if not path.exists():
        return None
    with path.open("r") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Convergence + Bayesian
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def load_convergence_matrix() -> pd.DataFrame:
    """Return convergence modalities as m1..m8 columns + ``gene``.

    Handles three shapes: dict-of-genes wrapper, bare list, and the current
    production form where ``modalities`` is a list of 8 floats.
    """
    raw = _read_json(paths.convergence_matrix())
    rows = raw["genes"] if isinstance(raw, dict) and "genes" in raw else raw
    df = pd.DataFrame(rows)

    modality_cols = [f"m{i}" for i in range(1, 9)]

    if "modalities" in df.columns:
        sample = df["modalities"].dropna().iloc[0] if len(df) else None
        if isinstance(sample, list):
            mod_df = pd.DataFrame(
                df["modalities"].tolist(), columns=modality_cols, index=df.index
            )
            df = pd.concat([df.drop(columns=["modalities"]), mod_df], axis=1)
        else:
            mods = pd.json_normalize(df["modalities"].fillna({}))
            for col in modality_cols:
                if col not in mods.columns:
                    mods[col] = 0
            df = pd.concat(
                [df.drop(columns=["modalities"]), mods[modality_cols]], axis=1
            )
    else:
        for col in modality_cols:
            if col not in df.columns:
                df[col] = 0

    # Normalise symbol column
    if "gene" not in df.columns and "symbol" in df.columns:
        df = df.rename(columns={"symbol": "gene"})

    for col in modality_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    if "count" not in df.columns:
        # Derive count if missing: #modalities with score > 0.2 (matches Next.js default)
        df["count"] = (df[modality_cols] > 0.2).sum(axis=1)
    return df


@st.cache_data(ttl=3600, show_spinner=False)
def load_bayesian_ranking() -> pd.DataFrame:
    """Canonical Bayesian ranking (Script 46d, 2026-04-23 onward).

    Replaces the legacy 46b archetype-similarity score (archived under
    archive/46b_retired_2026-04-23/). Each row carries posterior_prob, tier,
    concordance_state, dominant_stage_S1, druggability_tier, and the
    coloc_best_pp4 / coloc_best_gwas anchor.
    """
    raw = _read_json(paths.bayesian_ranking())
    rows = raw["genes"] if isinstance(raw, dict) and "genes" in raw else raw
    return pd.DataFrame(rows)


@st.cache_data(ttl=3600, show_spinner=False)
def load_proteomics_summary() -> dict:
    return _read_json(paths.proteomics_summary())


# ---------------------------------------------------------------------------
# Wide atlas mega-table (Parquet)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Network (portal_export_v2)
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, max_entries=256, show_spinner=False)
def load_gene_graph(symbol: str) -> dict | None:
    """Load portal_export_v2 per-gene graph (v4 schema).

    Returns ``None`` if the file is absent (some gene symbols have no network
    neighbors exported). Never raises — callers handle ``None`` as an empty
    neighborhood.
    """
    path = paths.gene_graph_path(symbol)
    if not path.exists():
        return None
    with path.open("r") as f:
        return json.load(f)


@st.cache_data(ttl=3600, show_spinner=False)
def load_community_labels() -> dict:
    path = paths.community_labels()
    if not path.exists():
        return {}
    with path.open("r") as f:
        return json.load(f)


@st.cache_data(ttl=3600, show_spinner=False)
def load_layer_metadata() -> list[dict]:
    """Return the portal_export_v2 layer metadata list (colors, descriptions)."""
    path = paths.data_dir() / "network" / "portal_export_v2" / "layer_metadata.json"
    if not path.exists():
        return []
    with path.open("r") as f:
        raw = json.load(f)
    return raw.get("layers", []) if isinstance(raw, dict) else raw


@st.cache_data(ttl=3600, show_spinner=False)
def load_filter_indexes() -> dict:
    """Load all filter_indexes/*.json into a dict keyed by basename."""
    out: dict[str, Any] = {}
    d = paths.filter_indexes_dir()
    if not d.exists():
        return out
    for p in d.glob("*.json"):
        try:
            with p.open("r") as f:
                out[p.stem] = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue
    return out


@st.cache_data(ttl=3600, show_spinner=False)
def load_query_bundles() -> dict:
    """Load Q2_*.json + Q3_*.json (and any other query_bundle files)."""
    out: dict[str, Any] = {}
    d = paths.query_bundles_dir()
    if not d.exists():
        return out
    for p in d.glob("*.json"):
        try:
            with p.open("r") as f:
                out[p.stem] = json.load(f)
        except (json.JSONDecodeError, OSError):
            continue
    return out


# ---------------------------------------------------------------------------
# Cross-species, drugs, knowledge graph
# ---------------------------------------------------------------------------

@st.cache_data(ttl=3600, show_spinner=False)
def load_cross_species_data() -> dict | None:
    """Return cross_species.json if present (else None; caller synthesises)."""
    path = paths.cross_species()
    if not path.exists():
        return None
    with path.open("r") as f:
        return json.load(f)


@st.cache_data(ttl=3600, show_spinner=False)
def load_drugs_data() -> dict | None:
    """Return drug_pipeline.json if present."""
    path = paths.drug_pipeline()
    if not path.exists():
        return None
    with path.open("r") as f:
        return json.load(f)


@st.cache_data(ttl=3600, show_spinner=False)
def load_clinical_drugs() -> pd.DataFrame:
    """Load clinical MASLD drug validation table.

    Prefers ``clinical_drugs.json`` in DATA_DIR (if ever exported); otherwise
    falls back to the CSV at
    ``RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv``.
    """
    json_path = paths.data_dir() / "clinical_drugs.json"
    if json_path.exists():
        with json_path.open("r") as f:
            raw = json.load(f)
        rows = raw["drugs"] if isinstance(raw, dict) and "drugs" in raw else raw
        return pd.DataFrame(rows)
    csv_path = paths.clinical_drugs_csv()
    if csv_path.exists():
        return pd.read_csv(csv_path)
    return pd.DataFrame()


@st.cache_data(ttl=3600, show_spinner=False)
def load_knowledge_graph() -> dict | None:
    path = paths.knowledge_graph()
    if not path.exists():
        return None
    with path.open("r") as f:
        return json.load(f)


@st.cache_data(ttl=3600, show_spinner=False)
def load_atlas_parquet(columns: tuple[str, ...] | None = None) -> pd.DataFrame:
    """Load the wide multi-evidence atlas as a DataFrame.

    ``columns`` accepts a tuple (required for hashing by ``st.cache_data``) and
    is forwarded to :func:`pandas.read_parquet`. Restricting columns keeps the
    Atlas / Explorer pages responsive on the free HF tier.
    """
    path = paths.atlas_parquet()
    if not path.exists():
        raise FileNotFoundError(
            f"atlas.parquet not found at {path}. "
            "Run `python -m data.build_parquet` from masld-atlas-portal/."
        )
    cols = list(columns) if columns else None
    return pd.read_parquet(path, columns=cols)
