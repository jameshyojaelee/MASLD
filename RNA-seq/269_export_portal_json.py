#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=net_269_export
#SBATCH --output=logs/net_269_export_%j.out
#SBATCH --error=logs/net_269_export_%j.err
# ===========================================================================
# Script 269: Export Portal JSON
# ===========================================================================
# Purpose: Stage 4a of the Bayesian multiplex gene network pipeline. Export
#          pre-computed network data as JSON files optimized for the Next.js
#          portal. Produces a self-contained directory with community hierarchy,
#          global layout, layer metadata, search index, and per-gene graphs.
#
# Input:
#   - RNA-seq/results/network/gene_graphs/*.json       (Script 267)
#   - RNA-seq/results/network/communities/communities.json (Script 268)
#   - RNA-seq/results/network/communities/community_labels.csv (Script 268)
#   - RNA-seq/results/network/global_layout.csv         (Script 265)
#   - RNA-seq/results/network/bayesian_diagnostics.csv  (Script 264)
#
# Output:
#   portal_export/
#   ├── communities.json
#   ├── global_layout.json
#   ├── layer_metadata.json
#   ├── search_index.json
#   └── gene_graphs/
#       ├── THRB.json
#       └── ...
#
# Environment: spatial (json, pandas, multiprocessing)
# ===========================================================================

import json
import os
import shutil
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
GENE_GRAPH_DIR = NET_DIR / "gene_graphs"
COMM_DIR = NET_DIR / "communities"
OUTDIR = NET_DIR / "portal_export"
OUTDIR_GRAPHS = OUTDIR / "gene_graphs"

# Layer color palette and descriptions
LAYER_DEFS = {
    "ppi":      {"color": "#e74c3c", "description": "Protein-protein interactions (STRING physical)"},
    "coexpr":   {"color": "#3498db", "description": "Co-expression (Pearson > threshold across 1,444 samples)"},
    "regulon":  {"color": "#2ecc71", "description": "Transcription factor regulon membership (SCENIC+)"},
    "lr":       {"color": "#f39c12", "description": "Ligand-receptor signaling (LIANA consensus)"},
    "genetic":  {"color": "#9b59b6", "description": "Genetic co-association (COLOC + TWAS)"},
    "pathway":  {"color": "#1abc9c", "description": "Shared pathway membership (MSigDB Hallmark + Reactome)"},
    "spatial":  {"color": "#e67e22", "description": "Spatial co-localization (Visium Moran's bivariate)"},
    "cosmos":   {"color": "#34495e", "description": "Mechanistic signaling paths (COSMOS)"},
    "cerna":    {"color": "#e91e63", "description": "ceRNA shared miRNA regulation"},
    "xspecies": {"color": "#607d8b", "description": "Cross-species concordance (human-mouse ortholog)"},
}

REQUIRED_GENE_GRAPH_FIELDS = {"gene", "community", "neighbors_composite", "halo"}

t0 = time.time()


def log(msg):
    elapsed = time.time() - t0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ===========================================================================
# 1. Validate a single gene graph JSON (used in parallel pool)
# ===========================================================================
def validate_gene_graph(path_str):
    """Validate a gene graph JSON. Returns (path, size, error_or_None)."""
    path = Path(path_str)
    try:
        with open(path) as f:
            data = json.load(f)
        missing = REQUIRED_GENE_GRAPH_FIELDS - set(data.keys())
        if missing:
            return (path_str, path.stat().st_size, f"missing fields: {missing}")
        return (path_str, path.stat().st_size, None)
    except json.JSONDecodeError as e:
        return (path_str, path.stat().st_size if path.exists() else 0, f"invalid JSON: {e}")
    except Exception as e:
        return (path_str, 0, str(e))


# ===========================================================================
# 2. Build layer_metadata.json
# ===========================================================================
def build_layer_metadata(diag_path):
    """Build layer metadata from bayesian_diagnostics.csv + hardcoded palette."""
    log("Building layer_metadata.json...")
    layers = []

    diag = None
    if diag_path.exists():
        diag = pd.read_csv(diag_path)
        log(f"  Loaded diagnostics: {len(diag)} rows")
    else:
        log(f"  WARNING: {diag_path} not found — pi0 and total_edges will be null")

    for name, meta in LAYER_DEFS.items():
        entry = {
            "name": name,
            "color": meta["color"],
            "description": meta["description"],
            "pi0": None,
            "total_edges": None,
        }
        if diag is not None:
            # Match by layer name column
            for col_candidate in ["layer", "layer_name", "name"]:
                if col_candidate in diag.columns:
                    row = diag[diag[col_candidate] == name]
                    if len(row) == 1:
                        if "pi0" in diag.columns:
                            val = row["pi0"].values[0]
                            entry["pi0"] = float(val) if pd.notna(val) else None
                        if "total_edges" in diag.columns:
                            val = row["total_edges"].values[0]
                            entry["total_edges"] = int(val) if pd.notna(val) else None
                        elif "n_edges" in diag.columns:
                            val = row["n_edges"].values[0]
                            entry["total_edges"] = int(val) if pd.notna(val) else None
                    break
        layers.append(entry)

    out_path = OUTDIR / "layer_metadata.json"
    with open(out_path, "w") as f:
        json.dump({"layers": layers}, f, indent=2)
    log(f"  Wrote {out_path} ({out_path.stat().st_size:,} bytes)")
    return layers


# ===========================================================================
# 3. Build communities.json
# ===========================================================================
def build_communities():
    """Merge community hierarchy and labels into a single JSON."""
    log("Building communities.json...")

    comm_json_path = COMM_DIR / "communities.json"
    labels_path = COMM_DIR / "community_labels.csv"

    if not comm_json_path.exists():
        log(f"  ERROR: {comm_json_path} not found — skipping communities.json")
        return

    with open(comm_json_path) as f:
        communities = json.load(f)
    log(f"  Loaded communities hierarchy: {len(communities)} top-level keys")

    # Merge labels if available
    if labels_path.exists():
        labels = pd.read_csv(labels_path)
        log(f"  Loaded community labels: {len(labels)} entries")
        label_map = {}
        for _, row in labels.iterrows():
            cid = str(row.get("community_id", row.get("community", "")))
            label = row.get("label", row.get("name", ""))
            enrichment = row.get("enrichment", row.get("top_pathway", ""))
            label_map[cid] = {"label": str(label), "enrichment": str(enrichment)}

        # Attach labels to community data if it is a dict of communities
        if isinstance(communities, dict) and "communities" in communities:
            for comm in communities.get("communities", []):
                cid = str(comm.get("id", comm.get("community_id", "")))
                if cid in label_map:
                    comm["label"] = label_map[cid]["label"]
                    comm["enrichment"] = label_map[cid]["enrichment"]
        elif isinstance(communities, list):
            for comm in communities:
                cid = str(comm.get("id", comm.get("community_id", "")))
                if cid in label_map:
                    comm["label"] = label_map[cid]["label"]
                    comm["enrichment"] = label_map[cid]["enrichment"]

        communities = {"communities": communities} if isinstance(communities, list) else communities
        communities["label_count"] = len(label_map)
    else:
        log("  WARNING: community_labels.csv not found — no enrichment labels")

    out_path = OUTDIR / "communities.json"
    with open(out_path, "w") as f:
        json.dump(communities, f, separators=(",", ":"))
    log(f"  Wrote {out_path} ({out_path.stat().st_size:,} bytes)")


# ===========================================================================
# 4. Build global_layout.json
# ===========================================================================
def build_global_layout(comm_assign_path):
    """Merge FA2 layout coordinates with community assignments. Normalize to [0,1]."""
    log("Building global_layout.json...")

    layout_path = NET_DIR / "global_layout.csv"
    if not layout_path.exists():
        log(f"  ERROR: {layout_path} not found — skipping global_layout.json")
        return

    layout = pd.read_csv(layout_path)
    log(f"  Loaded FA2 layout: {len(layout)} genes")

    # Identify coordinate columns
    x_col = next((c for c in layout.columns if c.lower() in ("x", "fa2_x", "pos_x")), None)
    y_col = next((c for c in layout.columns if c.lower() in ("y", "fa2_y", "pos_y")), None)
    gene_col = next((c for c in layout.columns if c.lower() in
                      ("gene", "human_symbol", "symbol", "gene_symbol", "node")), None)

    if not all([x_col, y_col, gene_col]):
        log(f"  ERROR: Could not identify x/y/gene columns in global_layout.csv "
            f"(found: {list(layout.columns[:10])})")
        return

    # Normalize coordinates to [0, 1]
    x_vals = layout[x_col].values.astype(float)
    y_vals = layout[y_col].values.astype(float)
    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()
    x_range = x_max - x_min if x_max != x_min else 1.0
    y_range = y_max - y_min if y_max != y_min else 1.0
    x_norm = (x_vals - x_min) / x_range
    y_norm = (y_vals - y_min) / y_range

    # Merge community assignments if available
    comm_macro = {}
    comm_meso = {}
    if comm_assign_path and comm_assign_path.exists():
        ca = pd.read_csv(comm_assign_path)
        ca_gene_col = next((c for c in ca.columns if c.lower() in
                            ("gene", "human_symbol", "symbol", "node")), None)
        if ca_gene_col:
            for _, row in ca.iterrows():
                g = str(row[ca_gene_col])
                if "community_macro" in ca.columns:
                    comm_macro[g] = row["community_macro"]
                elif "macro" in ca.columns:
                    comm_macro[g] = row["macro"]
                if "community_meso" in ca.columns:
                    comm_meso[g] = row["community_meso"]
                elif "meso" in ca.columns:
                    comm_meso[g] = row["meso"]
        log(f"  Merged community assignments: {len(comm_macro)} macro, {len(comm_meso)} meso")

    genes = []
    for i in range(len(layout)):
        symbol = str(layout[gene_col].iloc[i])
        entry = {
            "symbol": symbol,
            "x": round(float(x_norm[i]), 5),
            "y": round(float(y_norm[i]), 5),
        }
        if symbol in comm_macro:
            val = comm_macro[symbol]
            entry["community_macro"] = int(val) if pd.notna(val) else None
        if symbol in comm_meso:
            val = comm_meso[symbol]
            entry["community_meso"] = int(val) if pd.notna(val) else None
        genes.append(entry)

    out_path = OUTDIR / "global_layout.json"
    with open(out_path, "w") as f:
        json.dump({"genes": genes}, f, separators=(",", ":"))
    log(f"  Wrote {out_path} ({out_path.stat().st_size:,} bytes, {len(genes)} genes)")


# ===========================================================================
# 5. Build search_index.json
# ===========================================================================
def build_search_index(gene_graph_results):
    """Build compact search index from validated gene graphs.

    Format: {genes: [{s: symbol, d: degree, c: community_macro, n: top_neighbor}]}
    Target: < 300 KB.
    """
    log("Building search_index.json...")

    entries = []
    for path_str, size, error in gene_graph_results:
        if error is not None:
            continue
        path = Path(path_str)
        try:
            with open(path) as f:
                data = json.load(f)
            symbol = data.get("gene", path.stem)
            neighbors = data.get("neighbors_composite", [])
            degree = len(neighbors)
            community = data.get("community", {})
            comm_macro = community.get("macro", community.get("community_macro", None))

            # Top neighbor: highest composite weight
            top_neighbor = None
            if neighbors:
                # neighbors may be list of dicts or list of [gene, weight] pairs
                if isinstance(neighbors[0], dict):
                    sorted_n = sorted(neighbors, key=lambda x: x.get("weight", 0), reverse=True)
                    top_neighbor = sorted_n[0].get("gene", sorted_n[0].get("symbol", None))
                elif isinstance(neighbors[0], (list, tuple)) and len(neighbors[0]) >= 2:
                    sorted_n = sorted(neighbors, key=lambda x: x[1], reverse=True)
                    top_neighbor = sorted_n[0][0]

            entry = {"s": symbol, "d": degree}
            if comm_macro is not None:
                entry["c"] = int(comm_macro) if not isinstance(comm_macro, str) else comm_macro
            if top_neighbor is not None:
                entry["n"] = top_neighbor
            entries.append(entry)
        except Exception:
            continue

    # Sort by degree descending for quick lookup of high-connectivity genes
    entries.sort(key=lambda x: x.get("d", 0), reverse=True)

    out_path = OUTDIR / "search_index.json"
    with open(out_path, "w") as f:
        json.dump({"genes": entries}, f, separators=(",", ":"))

    size_kb = out_path.stat().st_size / 1024
    log(f"  Wrote {out_path} ({size_kb:.1f} KB, {len(entries)} genes)")
    if size_kb > 300:
        log(f"  WARNING: search_index.json exceeds 300 KB target ({size_kb:.1f} KB)")


# ===========================================================================
# 6. Validate and copy gene_graphs
# ===========================================================================
def validate_and_copy_gene_graphs():
    """Validate all gene graph JSONs in parallel, copy valid ones to portal_export."""
    log("Validating gene graphs...")

    if not GENE_GRAPH_DIR.exists():
        log(f"  ERROR: {GENE_GRAPH_DIR} not found — skipping gene graph export")
        return []

    json_files = sorted(GENE_GRAPH_DIR.glob("*.json"))
    log(f"  Found {len(json_files)} gene graph files")

    if len(json_files) == 0:
        return []

    # Parallel validation
    n_workers = min(cpu_count(), 8)
    with Pool(n_workers) as pool:
        results = pool.map(validate_gene_graph, [str(p) for p in json_files])

    n_valid = sum(1 for _, _, err in results if err is None)
    n_invalid = sum(1 for _, _, err in results if err is not None)
    log(f"  Validated: {n_valid} valid, {n_invalid} invalid")

    if n_invalid > 0:
        for path_str, _, err in results:
            if err is not None:
                log(f"    INVALID: {Path(path_str).name} — {err}")
                if n_invalid > 10:
                    log(f"    ... (showing first 10 of {n_invalid})")
                    break

    # Copy valid files to portal export
    OUTDIR_GRAPHS.mkdir(parents=True, exist_ok=True)
    n_copied = 0
    for path_str, _, err in results:
        if err is None:
            src = Path(path_str)
            dst = OUTDIR_GRAPHS / src.name
            shutil.copy2(src, dst)
            n_copied += 1

    log(f"  Copied {n_copied} gene graphs to {OUTDIR_GRAPHS}")
    return results


# ===========================================================================
# Main
# ===========================================================================
def main():
    log("=== Script 269: Export Portal JSON ===")
    log(f"Network dir: {NET_DIR}")

    # Create output directory
    OUTDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR_GRAPHS.mkdir(parents=True, exist_ok=True)

    # Community assignments path (used by global_layout)
    comm_assign_path = COMM_DIR / "community_assignments.csv"
    if not comm_assign_path.exists():
        # Fallback: check alternate location
        alt = NET_DIR / "community_assignments.csv"
        comm_assign_path = alt if alt.exists() else None

    # --- Step 1: Layer metadata ---
    diag_path = NET_DIR / "bayesian_diagnostics.csv"
    build_layer_metadata(diag_path)

    # --- Step 2: Communities ---
    build_communities()

    # --- Step 3: Global layout ---
    build_global_layout(comm_assign_path)

    # --- Step 4: Validate and copy gene graphs ---
    gene_graph_results = validate_and_copy_gene_graphs()

    # --- Step 5: Search index (depends on gene graph validation) ---
    build_search_index(gene_graph_results)

    # --- Step 6: Summary ---
    log("")
    log("=== Portal Export Summary ===")

    total_files = 0
    total_size = 0
    for p in OUTDIR.rglob("*.json"):
        total_files += 1
        total_size += p.stat().st_size

    log(f"Total files: {total_files}")
    log(f"Total size:  {total_size / 1024 / 1024:.2f} MB")

    # Gene graph size distribution
    graph_files = list(OUTDIR_GRAPHS.glob("*.json"))
    if graph_files:
        sizes = [p.stat().st_size for p in graph_files]
        largest = max(sizes)
        smallest = min(sizes)
        largest_name = graph_files[sizes.index(largest)].name
        smallest_name = graph_files[sizes.index(smallest)].name
        log(f"Gene graphs:  {len(graph_files)} files")
        log(f"  Largest:    {largest_name} ({largest / 1024:.1f} KB)")
        log(f"  Smallest:   {smallest_name} ({smallest / 1024:.1f} KB)")
        log(f"  Median:     {sorted(sizes)[len(sizes)//2] / 1024:.1f} KB")
    else:
        log("Gene graphs:  0 files")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
