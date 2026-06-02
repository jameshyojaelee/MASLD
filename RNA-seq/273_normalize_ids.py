#!/usr/bin/env python3
"""
273_normalize_ids.py
====================

Fixes a critical bug in the multi-evidence network pipeline: different edge
files use different gene identifiers (Ensembl IDs vs HGNC symbols), causing
downstream joins (262c confidence, 263 noisy-or composite, 265 multiplex
assembly, etc.) to silently drop all cross-layer matches.

What this does
--------------
1. Builds a canonical ENSG -> HGNC symbol map from network_nodes.csv.
   Version suffixes (e.g. ENSG00000198886.5) are stripped before building
   the map so that both versioned and un-versioned ENSG IDs resolve.
2. For each edge layer (ppi, coexpr, regulon, lr, genetic, pathway, cerna):
     - Detects the ID type of gene_a / gene_b from the first two rows
       (already_symbol / converted_from_ensembl / mixed).
     - Strips any ENSG version suffixes.
     - Maps ENSG IDs to HGNC symbols via the node map.
     - Drops edges whose endpoints have no symbol mapping (gene not in the
       network node set).
     - Preserves raw_score and metadata columns verbatim.
3. Backs up each original edge file to
   results/network/edges_backup_pre_id_norm/ before overwriting.
   Idempotent: skips backup if the backup already exists.
4. Writes results/network/id_normalization_summary.csv with per-layer counts.

SBATCH
------
Submit via logs/network_pipeline/scripts/273_idnorm.sbatch (written by
run_network_remediation.sh). Activates the `spatial` env.
"""

from __future__ import annotations

import shutil
import time
from pathlib import Path

import pandas as pd

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE      = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq")
NET_DIR   = BASE / "results" / "network"
NODES_CSV = NET_DIR / "network_nodes.csv"
BACKUP    = NET_DIR / "edges_backup_pre_id_norm"
SUMMARY   = NET_DIR / "id_normalization_summary.csv"

EDGE_FILES = [
    "edges_ppi.csv",
    "edges_coexpr.csv",
    "edges_regulon.csv",
    "edges_lr.csv",
    "edges_genetic.csv",
    "edges_pathway.csv",
    "edges_cerna.csv",
]

KEEP_COLS = ["gene_a", "gene_b", "raw_score", "metadata"]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def strip_version(x: str) -> str:
    """Strip trailing .N version from an Ensembl ID. Non-ENSG strings pass through."""
    if isinstance(x, str) and x.startswith("ENSG") and "." in x:
        return x.split(".", 1)[0]
    return x


def build_ensg_to_symbol_map(nodes_csv: Path) -> dict[str, str]:
    """Build an ENSG -> HGNC symbol map from the network node set."""
    nodes = pd.read_csv(nodes_csv, usecols=["human_symbol", "ensembl_id"])
    nodes = nodes.dropna(subset=["human_symbol", "ensembl_id"])
    nodes["ensembl_id"] = nodes["ensembl_id"].astype(str).map(strip_version)
    # If the same ENSG maps to multiple symbols, keep first occurrence.
    nodes = nodes.drop_duplicates(subset=["ensembl_id"], keep="first")
    return dict(zip(nodes["ensembl_id"], nodes["human_symbol"]))


def detect_mode(df: pd.DataFrame) -> str:
    """Inspect the first two rows to classify ID type.

    Returns one of: already_symbol | converted_from_ensembl | mixed | empty
    """
    if df.empty:
        return "empty"
    probe = df.head(2)
    a_ensg = probe["gene_a"].astype(str).str.startswith("ENSG")
    b_ensg = probe["gene_b"].astype(str).str.startswith("ENSG")
    ensg_frac = (a_ensg.sum() + b_ensg.sum()) / (2 * len(probe))
    if ensg_frac == 1.0:
        return "converted_from_ensembl"
    if ensg_frac == 0.0:
        return "already_symbol"
    return "mixed"


def normalize_layer(
    csv_path: Path,
    ensg2sym: dict[str, str],
    valid_symbols: set[str],
) -> tuple[pd.DataFrame, int, int, str]:
    """Load an edge file, map ENSG endpoints to symbols, drop unmapped edges.

    Returns (normalized_df, n_input, n_output, mode).
    """
    df = pd.read_csv(csv_path)
    n_in = len(df)

    # Safety: ensure required columns exist; preserve extras if present.
    for c in ["gene_a", "gene_b"]:
        if c not in df.columns:
            raise ValueError(f"{csv_path.name}: missing required column {c!r}")

    mode = detect_mode(df)

    # Strip version suffixes defensively on every row (cheap).
    df["gene_a"] = df["gene_a"].astype(str).map(strip_version)
    df["gene_b"] = df["gene_b"].astype(str).map(strip_version)

    def to_symbol(x: str) -> str | float:
        if x.startswith("ENSG"):
            return ensg2sym.get(x, pd.NA)
        # Already a symbol; keep only if present in the node symbol set.
        return x if x in valid_symbols else pd.NA

    df["gene_a"] = df["gene_a"].map(to_symbol)
    df["gene_b"] = df["gene_b"].map(to_symbol)

    df = df.dropna(subset=["gene_a", "gene_b"])
    n_out = len(df)

    # Preserve raw_score + metadata; also keep any other columns already present.
    ordered = [c for c in KEEP_COLS if c in df.columns]
    extras = [c for c in df.columns if c not in ordered]
    df = df[ordered + extras]

    return df, n_in, n_out, mode


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    print("=== 273_normalize_ids.py ===")
    print(f"Base dir: {BASE}")
    print(f"Network:  {NET_DIR}")

    if not NODES_CSV.exists():
        raise FileNotFoundError(f"Missing node file: {NODES_CSV}")

    # --- 1. Build canonical map ------------------------------------------------
    print(f"\n[1/3] Building ENSG -> symbol map from {NODES_CSV.name}")
    ensg2sym = build_ensg_to_symbol_map(NODES_CSV)
    valid_symbols = set(ensg2sym.values())
    print(f"      {len(ensg2sym):,} ENSG IDs -> {len(valid_symbols):,} unique symbols")

    # --- 2. Back up originals --------------------------------------------------
    BACKUP.mkdir(parents=True, exist_ok=True)
    print(f"\n[2/3] Backing up originals to {BACKUP}")
    for name in EDGE_FILES:
        src = NET_DIR / name
        dst = BACKUP / name
        if not src.exists():
            print(f"      SKIP (missing): {name}")
            continue
        if dst.exists():
            print(f"      SKIP (backup exists): {name}")
            continue
        shutil.copy2(src, dst)
        print(f"      backed up: {name}")

    # --- 3. Normalize each layer ----------------------------------------------
    print(f"\n[3/3] Normalizing edge files")
    rows = []
    for name in EDGE_FILES:
        src = NET_DIR / name
        if not src.exists():
            print(f"      SKIP (missing): {name}")
            continue
        t_layer = time.time()
        norm, n_in, n_out, mode = normalize_layer(src, ensg2sym, valid_symbols)
        norm.to_csv(src, index=False)
        dropped = n_in - n_out
        pct = 100.0 * dropped / n_in if n_in else 0.0
        layer = name.replace("edges_", "").replace(".csv", "")
        rows.append({
            "layer": layer,
            "n_edges_input": n_in,
            "n_edges_output": n_out,
            "n_edges_dropped": dropped,
            "pct_dropped": round(pct, 3),
            "mode": mode,
        })
        print(f"      {layer:<10} in={n_in:>9,}  out={n_out:>9,}  "
              f"drop={dropped:>8,} ({pct:5.2f}%)  mode={mode}  "
              f"[{time.time() - t_layer:5.1f}s]")

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY, index=False)
    print(f"\nWrote summary: {SUMMARY}")

    # --- Print summary table --------------------------------------------------
    print("\n=== ID Normalization Summary ===")
    print(f"{'Layer':<10} {'Input':>10} {'Output':>10} {'Dropped':>10} {'%Drop':>7}  Mode")
    for r in rows:
        print(f"{r['layer']:<10} {r['n_edges_input']:>10,} {r['n_edges_output']:>10,} "
              f"{r['n_edges_dropped']:>10,} {r['pct_dropped']:>6.2f}%  {r['mode']}")

    print(f"\nTotal elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
