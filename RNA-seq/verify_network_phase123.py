#!/usr/bin/env python3
"""Verification script for network Phase 1/2/3 rehab.

Reads the rebuilt portal gene_graphs JSONs for a handful of spot-check
genes and prints an edge-count breakdown by bucket. Intended to confirm:

- CYP3A4 (stable hepatocyte metabolic hub) gains `neighbors_d_stable`
  edges it was missing before.
- HNF4A (TF) shows up in `neighbors_d_regulon_out` (it's a SCENIC+ TF).
- COL1A1 (fibrosis collagen) retains F2-dynamic edges.
- RORA, THRB (drug targets) -- layer-rich baseline, should remain rich.
- Compares `edge_counts` section of each JSON to confirm the new layers
  appear.
"""
from __future__ import annotations

import json
from pathlib import Path

BASE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
GRAPHS_DIR = BASE / "RNA-seq" / "results" / "network" / "portal_export_v2" / "gene_graphs"
LAYER_COUNTS_CSV = BASE / "RNA-seq" / "results" / "network" / "layer_counts.csv"

SPOT_CHECK = ["CYP3A4", "HNF4A", "COL1A1", "RORA", "THRB", "APOB", "HMGCS2"]


def main() -> int:
    if LAYER_COUNTS_CSV.exists():
        print("=== Atlas-wide layer counts (from 290) ===")
        print(LAYER_COUNTS_CSV.read_text())

    if not GRAPHS_DIR.exists():
        print(f"[FAIL] gene graphs dir not found: {GRAPHS_DIR}")
        return 1

    print("\n=== Per-gene neighborhood breakdown ===")
    for gene in SPOT_CHECK:
        p = GRAPHS_DIR / f"{gene}.json"
        if not p.exists():
            print(f"\n[{gene}] (no JSON found)")
            continue
        doc = json.loads(p.read_text())
        neighbors = doc.get("neighbors", {})
        ec = doc.get("edge_counts", {})
        badges = doc.get("badges", {})
        print(f"\n[{gene}]")
        print(f"  edge_counts: {ec}")
        if badges:
            print(f"  badges (Phase-3 collapsed layers): {badges}")
        for bucket, rows in neighbors.items():
            if rows:
                print(f"  {bucket}: {len(rows)} edges")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
