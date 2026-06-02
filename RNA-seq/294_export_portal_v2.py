#!/usr/bin/env python3
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --job-name=net_294_portal
#SBATCH --output=logs/net_294_portal_%j.out
#SBATCH --error=logs/net_294_portal_%j.err
# ===========================================================================
# Script 294: Export Portal JSON v2 (MASLD Network v2 - Architecture C)
# ===========================================================================
# Purpose: Replaces v1 (269_export_portal_json.py). Reads the edge annotation
#          atlas directly (3.37M edges, 40+ cols) rather than composite scores.
#          Emits a per-gene JSON graph keyed to the new STRING + 5 D-type
#          (D-F2, D-COLOC, D-LR, D-ceRNA, D-XS) layer design, plus F0-F1 /
#          F3-F4 community JSON and transition table for the portal sidebar.
#
# Inputs:
#   RNA-seq/results/network/edge_annotation_atlas.parquet
#   RNA-seq/results/network/network_nodes.csv
#   RNA-seq/results/network/communities_f2/communities_F01.csv
#   RNA-seq/results/network/communities_f2/communities_F34.csv
#   RNA-seq/results/network/communities_f2/community_transition_table.csv
#   RNA-seq/results/network/communities_f2/community_enrichment_F01.csv
#   RNA-seq/results/network/communities_f2/community_enrichment_F34.csv
#   RNA-seq/results/multi_evidence/multi_evidence_atlas.csv   (subset of cols)
#   RNA-seq/results/network/edge_annotation_summary.csv       (optional)
#
# Output: RNA-seq/results/network/portal_export_v2/
#   manifest.json, search_index.json, layer_metadata.json,
#   communities_F01.json, communities_F34.json, community_transitions.json,
#   stage_manifest.json, gene_graphs/{SYMBOL}.json
#
# Environment: spatial (pandas, pyarrow, numpy)
# ===========================================================================

from __future__ import annotations

import json
import math
import os
import sys
import time
import warnings
from datetime import datetime, timezone
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
NET_DIR = BASE / "RNA-seq" / "results" / "network"
ATLAS_PARQUET = NET_DIR / "edge_annotation_atlas.parquet"
NODES_CSV = NET_DIR / "network_nodes.csv"
COMM_DIR = NET_DIR / "communities_f2"
COMM_F01_CSV = COMM_DIR / "communities_F01.csv"
COMM_F34_CSV = COMM_DIR / "communities_F34.csv"
COMM_TRANS_CSV = COMM_DIR / "community_transition_table.csv"
ENRICH_F01_CSV = COMM_DIR / "community_enrichment_F01.csv"
ENRICH_F34_CSV = COMM_DIR / "community_enrichment_F34.csv"
ATLAS_CSV = BASE / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"
EDGE_SUMMARY_CSV = NET_DIR / "edge_annotation_summary.csv"

OUTDIR = NET_DIR / "portal_export_v2"
OUTDIR_GRAPHS = OUTDIR / "gene_graphs"

# ---------------------------------------------------------------------------
# Layer palette (STRING + 5 D-types)
# ---------------------------------------------------------------------------

LAYER_DEFS: dict[str, dict[str, str]] = {
    "S":        {"color": "#34495e", "description": "STRING v12 high-confidence PPI baseline (score >= 0.7)"},
    "D-F2":     {"color": "#e74c3c", "description": "Stage-dynamic co-expression edges (F2-transition specific; continuous |delta_r|)"},
    "D-STABLE": {"color": "#1abc9c", "description": "Stable co-expression across all fibrosis stages (|r_mean| >= 0.6, stage-invariant)"},
    "D-REG":    {"color": "#d35400", "description": "Directed TF -> target edges from SCENIC+ regulons (hepatocyte + disease)"},
    "D-COLOC":  {"color": "#9b59b6", "description": "GWAS-eQTL colocalized gene pairs (PP4 >= 0.3 shared locus + eQTL sentinel sharing)"},
    "D-LR":     {"color": "#f39c12", "description": "Ligand-receptor signaling edges (LIANA consensus across cell-type pairs)"},
    "D-ceRNA":  {"color": "#e91e63", "description": "Competing endogenous RNA edges (shared miRNA regulation)"},
    "D-XS":     {"color": "#2ecc71", "description": "Cross-species conserved edges (WGCNA module preservation + cell-type concordance)"},
}
LAYER_ORDER = ["S", "D-F2", "D-STABLE", "D-REG", "D-COLOC", "D-LR", "D-ceRNA", "D-XS"]

# Gene-level columns to pull from the multi-evidence atlas (kept compact for <10KB JSONs)
ATLAS_USECOLS = [
    "human_symbol",
    "ensembl_id",
    "gene_biotype",
    "dream_logFC",
    "dream_padj",
    "is_deg",
    "coloc_susie_best_pp4",
    "is_conserved",
    "sex_class",
    "zonation_class",
    "ferroptosis_class",
    "dgidb_druggable",
    "attribution_class",
]

TOP_K_STRING = 20
TOP_K_PER_DTYPE = 20

t0 = time.time()


def log(msg: str) -> None:
    elapsed = time.time() - t0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _scalar(v: Any) -> Any:
    """Convert a pandas/NumPy scalar to a JSON-safe python primitive."""
    if v is None:
        return None
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    try:
        if pd.isna(v):
            return None
    except (ValueError, TypeError):
        pass
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        f = float(v)
        return None if math.isnan(f) else f
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, (np.ndarray,)):
        return [_scalar(x) for x in v.tolist()]
    return v


def _round(v: Any, nd: int = 4) -> Any:
    s = _scalar(v)
    if isinstance(s, float):
        return round(s, nd)
    return s


def _write_json(path: Path, obj: Any, compact: bool = True) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    sep = (",", ":") if compact else (", ", ": ")
    with open(path, "w") as f:
        json.dump(obj, f, separators=sep, default=_scalar)
    return path.stat().st_size


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_atlas() -> pd.DataFrame:
    log(f"Loading edge atlas from {ATLAS_PARQUET.name} ...")
    if not ATLAS_PARQUET.exists():
        log(f"  ERROR: atlas parquet missing: {ATLAS_PARQUET}")
        sys.exit(1)
    df = pd.read_parquet(ATLAS_PARQUET)
    log(f"  Loaded {len(df):,} edges, {df.shape[1]} columns")
    return df


def load_nodes() -> pd.DataFrame:
    log(f"Loading nodes from {NODES_CSV.name} ...")
    if not NODES_CSV.exists():
        log(f"  ERROR: nodes CSV missing: {NODES_CSV}")
        sys.exit(1)
    nodes = pd.read_csv(NODES_CSV, low_memory=False)
    log(f"  Loaded {len(nodes):,} nodes")
    return nodes


def load_gene_atlas() -> pd.DataFrame:
    """Load multi-evidence atlas, only required columns."""
    if not ATLAS_CSV.exists():
        log(f"  WARNING: atlas CSV missing: {ATLAS_CSV}")
        return pd.DataFrame(columns=ATLAS_USECOLS)
    # Determine which requested columns actually exist
    header = pd.read_csv(ATLAS_CSV, nrows=0).columns.tolist()
    cols = [c for c in ATLAS_USECOLS if c in header]
    missing = [c for c in ATLAS_USECOLS if c not in header]
    if missing:
        log(f"  WARNING: atlas missing columns: {missing}")
    log(f"Loading multi-evidence atlas ({len(cols)} columns) ...")
    df = pd.read_csv(ATLAS_CSV, usecols=cols, low_memory=False)
    log(f"  Loaded {len(df):,} atlas rows")
    return df


def load_communities() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    log("Loading community assignments + enrichment ...")
    f01 = pd.read_csv(COMM_F01_CSV) if COMM_F01_CSV.exists() else pd.DataFrame()
    f34 = pd.read_csv(COMM_F34_CSV) if COMM_F34_CSV.exists() else pd.DataFrame()
    trans = pd.read_csv(COMM_TRANS_CSV) if COMM_TRANS_CSV.exists() else pd.DataFrame()
    enr01 = pd.read_csv(ENRICH_F01_CSV) if ENRICH_F01_CSV.exists() else pd.DataFrame()
    enr34 = pd.read_csv(ENRICH_F34_CSV) if ENRICH_F34_CSV.exists() else pd.DataFrame()
    log(f"  F01: {len(f01)} genes, F34: {len(f34)} genes, "
        f"transitions: {len(trans)}, enr01: {len(enr01)}, enr34: {len(enr34)}")
    return f01, f34, trans, enr01, enr34


# ---------------------------------------------------------------------------
# Community JSONs
# ---------------------------------------------------------------------------

def top_hallmark_per_community(enr: pd.DataFrame) -> dict[int, dict[str, Any]]:
    """Pick the smallest-pvalue Hallmark gene_set per community_id."""
    if enr.empty or "community_id" not in enr.columns:
        return {}
    enr = enr.copy()
    if "pvalue" in enr.columns:
        enr = enr.sort_values("pvalue", ascending=True, kind="mergesort")
    out: dict[int, dict[str, Any]] = {}
    for cid, sub in enr.groupby("community_id", sort=False):
        row = sub.iloc[0]
        out[int(cid)] = {
            "top_hallmark": _scalar(row.get("gene_set")),
            "top_hallmark_p": _round(row.get("pvalue"), 6),
            "overlap": _scalar(row.get("overlap")),
            "set_size": _scalar(row.get("set_size")),
            "comm_size": _scalar(row.get("community_size") or row.get("comm_size_in_universe")),
        }
    return out


def build_communities_json(comm_df: pd.DataFrame, enr_df: pd.DataFrame, trans_df: pd.DataFrame,
                            side: str) -> dict[str, Any]:
    """Build communities_F01.json / communities_F34.json."""
    if comm_df.empty:
        return {"side": side, "communities": []}
    top_hit = top_hallmark_per_community(enr_df)

    # dominance_change = share of community members that persist to the other side
    dom: dict[int, dict[str, Any]] = {}
    if not trans_df.empty:
        key_col = "F01_id" if side == "F01" else "F34_id"
        other_col = "F34_id" if side == "F01" else "F01_id"
        n_col = "n_F01" if side == "F01" else "n_F34"
        if key_col in trans_df.columns:
            for cid, sub in trans_df.groupby(key_col):
                total = int(sub[n_col].iloc[0]) if n_col in sub.columns else 0
                shared_max = int(sub["n_shared_genes"].max()) if "n_shared_genes" in sub.columns else 0
                dom[int(cid)] = {
                    "n_genes": total,
                    "best_partner": _scalar(sub.sort_values("n_shared_genes", ascending=False).iloc[0][other_col]),
                    "best_shared": shared_max,
                    "dominance_change": round(shared_max / total, 4) if total else None,
                }

    out_communities = []
    gid_col = "community_id" if "community_id" in comm_df.columns else "macro_id"
    for cid, sub in comm_df.groupby(gid_col):
        cid_i = int(cid)
        enrich = top_hit.get(cid_i, {})
        dom_i = dom.get(cid_i, {})
        out_communities.append({
            "id": cid_i,
            "label": f"{side}-C{cid_i}",
            "n_genes": int(len(sub)),
            "top_hallmark": enrich.get("top_hallmark"),
            "top_hallmark_p": enrich.get("top_hallmark_p"),
            "overlap": enrich.get("overlap"),
            "set_size": enrich.get("set_size"),
            "dominance_change": dom_i.get("dominance_change"),
            "best_partner_other_side": dom_i.get("best_partner"),
            "best_shared": dom_i.get("best_shared"),
        })
    out_communities.sort(key=lambda c: -c["n_genes"])
    return {"side": side, "n_communities": len(out_communities), "communities": out_communities}


def build_transitions_json(trans_df: pd.DataFrame) -> dict[str, Any]:
    if trans_df.empty:
        return {"transitions": []}
    records = []
    for _, row in trans_df.iterrows():
        records.append({
            "f01_id": _scalar(row.get("F01_id")),
            "f34_id": _scalar(row.get("F34_id")),
            "n_shared": _scalar(row.get("n_shared_genes")),
            "n_F01": _scalar(row.get("n_F01")),
            "n_F34": _scalar(row.get("n_F34")),
            "jaccard": _round(row.get("jaccard"), 4),
            "F01_class": _scalar(row.get("F01_class")),
            "F34_class": _scalar(row.get("F34_class")),
            "classification": _scalar(row.get("classification")),
        })
    return {"n_transitions": len(records), "transitions": records}


# ---------------------------------------------------------------------------
# Per-gene edge packaging
# ---------------------------------------------------------------------------

EDGE_ANNOTATION_COLS = [
    "string_score", "in_string_ge700", "is_contested",
    "r_F01", "r_F2", "r_F34", "delta_max", "emergence_stage",
    "direction_at_F2", "fdr_f2", "loco_replication_fraction",
    # D-F2 deconvolution-adjusted annotation (281d):
    "f2_delta_max_raw", "f2_delta_max_resid",
    "f2_deconv_concordance", "f2_composition_driven",
    # D-COLOC (including expanded channel info from 287):
    "pp4_min", "gwas_list", "locus_id", "channel", "pp4_thresh",
    # D-LR consensus (288):
    "score_diff_lr", "cell_type_pairs", "emergence_stage_lr",
    "direction", "n_contexts",
    "differential_at_F2", "sex_specific", "metabolic_hep_mac",
    # D-REG (293):
    "reg_tf", "reg_target", "regulon_id", "regulon_source",
    "regulon_activity_diff", "regulon_activity_padj", "regulon_directed",
    # D-STABLE (286):
    "stable_r_mean", "stable_min_r", "stable_max_r",
    "stable_range_r", "stable_delta_max",
    # D-XS (297):
    "xs_channel", "xs_module_id",
    "xs_preservation_Z", "xs_celltype",
    # D-ceRNA expanded (289):
    "n_shared_mirnas", "shared_mirnas", "provenance",
    # Node annotations:
    "celltype_driver_a", "celltype_driver_b",
    "coloc_linked_gene_a", "coloc_linked_gene_b",
    "coloc_linked", "druggable_pair",
    "druggable_a", "druggable_b",
    "sex_class_a", "sex_class_b",
    "conserved_mouse_a", "conserved_mouse_b",
    "ferroptosis_a", "ferroptosis_b",
    "zonation_a", "zonation_b",
    "provenance_hash",
]


def _edge_record(row: pd.Series, gene: str) -> dict[str, Any]:
    """Build a compact edge record from an atlas row, oriented at `gene`."""
    a = row["gene_a"]; b = row["gene_b"]
    partner = b if a == gene else a
    rec: dict[str, Any] = {
        "partner": _scalar(partner),
        "type": _scalar(row.get("type")),
    }
    for col in EDGE_ANNOTATION_COLS:
        if col not in row.index:
            continue
        v = row[col]
        if isinstance(v, float):
            rec[col] = _round(v, 4)
        else:
            rec[col] = _scalar(v)
    # Drop null fields to keep payloads small
    return {k: v for k, v in rec.items() if v is not None and v != ""}


def _rank_key(rec: dict[str, Any]) -> float:
    """Ranking key for Top-K per type. Larger is better."""
    t = rec.get("type", "")
    if "S" in t.split(","):
        return float(rec.get("string_score") or 0.0)
    if "D-F2" in t:
        return float(rec.get("delta_max") or 0.0)
    if "D-STABLE" in t:
        return float(rec.get("stable_r_mean") or 0.0)
    if "D-REG" in t:
        # Strong regulon edges = low padj (= high 1 - padj)
        padj = rec.get("regulon_activity_padj")
        return float(1.0 - padj) if padj is not None else 0.0
    if "D-COLOC" in t:
        v = rec.get("pp4_min")
        return float(v) if v is not None else 0.0
    if "D-LR" in t:
        return float(rec.get("score_diff_lr") or 0.0)
    if "D-ceRNA" in t:
        return float(rec.get("n_shared_mirnas") or rec.get("delta_max") or 0.0)
    if "D-XS" in t:
        return float(rec.get("xs_preservation_Z") or 0.0)
    return 0.0


# ---------------------------------------------------------------------------
# Gene graph writer (runs in worker processes)
# ---------------------------------------------------------------------------

_WORKER_STATE: dict[str, Any] = {}


def _worker_init(edges_pkl: str, gene_meta: dict[str, dict[str, Any]],
                  community_meta: dict[str, dict[int, dict[str, Any]]],
                  outdir: str,
                  collapsed_layers: set[str]) -> None:
    # Load once per worker. The edges DataFrame is shared via pickle path.
    _WORKER_STATE["edges"] = pd.read_pickle(edges_pkl)
    _WORKER_STATE["by_a"] = _WORKER_STATE["edges"].groupby("gene_a", sort=False).indices
    _WORKER_STATE["by_b"] = _WORKER_STATE["edges"].groupby("gene_b", sort=False).indices
    _WORKER_STATE["gene_meta"] = gene_meta
    _WORKER_STATE["community_meta"] = community_meta
    _WORKER_STATE["outdir"] = Path(outdir)
    _WORKER_STATE["collapsed_layers"] = collapsed_layers


def _build_single_gene(gene: str) -> tuple[str, int, str | None]:
    try:
        edges = _WORKER_STATE["edges"]
        by_a = _WORKER_STATE["by_a"]
        by_b = _WORKER_STATE["by_b"]
        idx_a = by_a.get(gene, np.empty(0, dtype=int))
        idx_b = by_b.get(gene, np.empty(0, dtype=int))
        if len(idx_a) + len(idx_b) == 0:
            sub = edges.iloc[[]]
        else:
            idx = np.concatenate([idx_a, idx_b])
            sub = edges.iloc[idx]

        # Build records, bucket by type. `type` is a comma-joined tag list
        # like "S,D-F2" so an edge is added to every matching bucket.
        buckets: dict[str, list[dict[str, Any]]] = {t: [] for t in LAYER_ORDER}
        r_F01_vals: list[float] = []
        r_F2_vals: list[float] = []
        r_F34_vals: list[float] = []
        emergence_stages: list[str] = []
        contested_list: list[dict[str, Any]] = []
        druggable_pairs: list[dict[str, Any]] = []

        for _, row in sub.iterrows():
            rec = _edge_record(row, gene)
            t = rec.get("type") or ""
            tags = [s for s in t.split(",") if s]
            for tag in tags:
                if tag in buckets:
                    buckets[tag].append(rec)
            # Collect trajectory info
            for src, dst in (("r_F01", r_F01_vals), ("r_F2", r_F2_vals), ("r_F34", r_F34_vals)):
                v = rec.get(src)
                if isinstance(v, (int, float)):
                    dst.append(float(v))
            if rec.get("emergence_stage"):
                emergence_stages.append(rec["emergence_stage"])
            if rec.get("is_contested"):
                contested_list.append(rec)
            if rec.get("druggable_pair"):
                druggable_pairs.append(rec)

        neighbors: dict[str, list[dict[str, Any]]] = {}

        # Top STRING edges (K=20)
        s_edges = sorted(buckets.get("S", []), key=_rank_key, reverse=True)[:TOP_K_STRING]
        neighbors["neighbors_string"] = s_edges

        # D-F2: split by emergence. All six categories from the atlas are
        # emitted so the portal can filter by F-stage. Top-K per category.
        df2 = buckets.get("D-F2", [])
        for stage, key in [
            ("F2_emerging",     "neighbors_d_f2_emerging"),
            ("F2_dissolving",   "neighbors_d_f2_dissolving"),
            ("transient_F2",    "neighbors_d_f2_transient"),
            ("progressive_up",  "neighbors_d_f2_progressive_up"),
            ("progressive_down","neighbors_d_f2_progressive_down"),
            ("F01_specific",    "neighbors_d_f2_f01_specific"),
            ("F34_specific",    "neighbors_d_f2_f34_specific"),
        ]:
            subset = [e for e in df2 if e.get("emergence_stage") == stage]
            neighbors[key] = sorted(subset, key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]

        # D-STABLE: Top-K (potentially very many; cap per-gene).
        stable_edges = sorted(buckets.get("D-STABLE", []), key=_rank_key, reverse=True)
        neighbors["neighbors_d_stable"] = stable_edges[:TOP_K_PER_DTYPE]

        # D-REG: directed TF -> target. Split by whether this gene is the TF
        # (out-degree) or the target (in-degree).
        reg_edges = buckets.get("D-REG", [])
        reg_out = [e for e in reg_edges if e.get("reg_tf") == gene]
        reg_in = [e for e in reg_edges if e.get("reg_target") == gene]
        reg_other = [e for e in reg_edges if e.get("reg_tf") != gene and e.get("reg_target") != gene]
        neighbors["neighbors_d_regulon_out"] = sorted(reg_out, key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]
        neighbors["neighbors_d_regulon_in"] = sorted(reg_in, key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]
        if reg_other:
            neighbors["neighbors_d_regulon_other"] = sorted(reg_other, key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]

        # D-COLOC: all (should be few)
        neighbors["neighbors_d_coloc"] = sorted(buckets.get("D-COLOC", []), key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]

        # D-LR: all (now ~1-2K total across all genes, few per gene)
        neighbors["neighbors_d_lr"] = sorted(buckets.get("D-LR", []), key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]

        # D-ceRNA + D-XS: Top-K
        neighbors["top_d_cerna"] = sorted(buckets.get("D-ceRNA", []), key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]
        neighbors["top_d_xs"] = sorted(buckets.get("D-XS", []), key=_rank_key, reverse=True)[:TOP_K_PER_DTYPE]

        neighbors["contested"] = contested_list[:TOP_K_PER_DTYPE]
        neighbors["druggable_pairs"] = druggable_pairs[:TOP_K_PER_DTYPE]

        # Gene-level attributes
        meta = _WORKER_STATE["gene_meta"].get(gene, {})

        # Community membership (from community assignments)
        comm_meta = _WORKER_STATE["community_meta"]
        f01_entry = comm_meta.get("F01", {}).get(gene)
        f34_entry = comm_meta.get("F34", {}).get(gene)

        # Trajectory summary
        def _mean(xs: list[float]) -> float | None:
            return round(float(np.mean(xs)), 4) if xs else None

        # Most common emergence_stage among edges incident to this gene
        emerg_mode = None
        if emergence_stages:
            vals, counts = np.unique(emergence_stages, return_counts=True)
            emerg_mode = str(vals[int(np.argmax(counts))])

        # Phase-3 collapse-to-badge: layers below the floor surface as
        # node-level badges rather than as toggleable layer tabs. The
        # badge tells the reader this gene participates in the sparse
        # layer without wasting a UI slot on a near-empty tab.
        collapsed = _WORKER_STATE.get("collapsed_layers", set())
        badges: dict[str, Any] = {}
        for tag in collapsed:
            n_here = len(buckets.get(tag, []))
            if n_here > 0:
                badges[tag] = n_here

        doc = {
            "gene": gene,
            "ensembl_id": meta.get("ensembl_id"),
            "biotype": meta.get("gene_biotype"),
            "attributes": {
                "is_deg": meta.get("is_deg"),
                "dream_logFC": meta.get("dream_logFC"),
                "dream_padj": meta.get("dream_padj"),
                "coloc_susie_best_pp4": meta.get("coloc_susie_best_pp4"),
                "is_conserved": meta.get("is_conserved"),
                "sex_class": meta.get("sex_class"),
                "dgidb_druggable": meta.get("dgidb_druggable"),
                "ferroptosis_class": meta.get("ferroptosis_class"),
                "zonation_class": meta.get("zonation_class"),
                "attribution_class": meta.get("attribution_class"),
            },
            "community": {
                "f01": f01_entry,
                "f34": f34_entry,
            },
            "neighbors": neighbors,
            "edge_counts": {t: len(buckets.get(t, [])) for t in LAYER_ORDER},
            "badges": badges,
            "f_stage_trajectory": {
                "r_F01": _mean(r_F01_vals),
                "r_F2": _mean(r_F2_vals),
                "r_F34": _mean(r_F34_vals),
                "emergence_stage": emerg_mode,
            },
        }

        # Strip None from attributes to shrink payload
        doc["attributes"] = {k: v for k, v in doc["attributes"].items() if v is not None}

        outdir: Path = _WORKER_STATE["outdir"]
        out_path = outdir / f"{gene}.json"
        size = _write_json(out_path, doc, compact=True)
        return (gene, size, None)
    except Exception as exc:  # pragma: no cover - diagnostic
        return (gene, 0, f"{type(exc).__name__}: {exc}")


# ---------------------------------------------------------------------------
# Manifest + search_index + layer metadata
# ---------------------------------------------------------------------------

def build_layer_metadata(edges: pd.DataFrame,
                          summary_df: pd.DataFrame | None,
                          collapsed_layers: set[str] | None = None) -> list[dict[str, Any]]:
    log("Building layer_metadata.json ...")
    collapsed_layers = collapsed_layers or set()
    # Per-layer edge counts via the `type` tag list (multi-tag edges count
    # in every tagged bucket). Use str.contains with word boundaries —
    # 15M row `apply` loops are unusable, pandas regex is ~1000x faster.
    def count_layer(tag: str) -> int:
        col = edges["type"].fillna("")
        pat = r"(?:^|,)" + tag + r"(?:,|$)"
        return int(col.str.contains(pat, regex=True).sum())
    layers: list[dict[str, Any]] = []
    gold_map: dict[str, float] = {}
    if summary_df is not None and "type" in summary_df.columns:
        for _, row in summary_df.iterrows():
            t = str(row["type"])
            for col in ("gold_fold_enrichment", "fold_enrichment", "gold_enrichment"):
                if col in summary_df.columns and pd.notna(row.get(col)):
                    gold_map[t] = float(row[col])
                    break
    for t in LAYER_ORDER:
        meta = LAYER_DEFS[t]
        layers.append({
            "type": t,
            "color": meta["color"],
            "description": meta["description"],
            "n_edges": count_layer(t),
            "gold_fold_enrichment": gold_map.get(t),
            "collapsed_to_badge": t in collapsed_layers,
        })
    size = _write_json(OUTDIR / "layer_metadata.json", {"layers": layers}, compact=False)
    log(f"  Wrote layer_metadata.json ({size:,} bytes)")
    return layers


def build_search_index(nodes: pd.DataFrame, gene_meta: dict[str, dict[str, Any]],
                         degree_map: dict[str, int],
                         comm_meta: dict[str, dict[int, dict[str, Any]]]) -> None:
    log("Building search_index.json ...")
    entries: list[dict[str, Any]] = []
    f01 = comm_meta.get("F01", {})
    f34 = comm_meta.get("F34", {})
    for sym in nodes["human_symbol"].dropna().astype(str).tolist():
        meta = gene_meta.get(sym, {})
        e: dict[str, Any] = {
            "s": sym,
            "d": int(degree_map.get(sym, 0)),
        }
        f01e = f01.get(sym)
        f34e = f34.get(sym)
        if f01e and f01e.get("macro_id") is not None:
            e["c_macro_F01"] = f01e.get("macro_id")
        if f34e and f34e.get("macro_id") is not None:
            e["c_macro_F34"] = f34e.get("macro_id")
        if meta.get("is_deg"):
            e["is_deg"] = True
        pp4 = meta.get("coloc_susie_best_pp4")
        if isinstance(pp4, (int, float)) and pp4 >= 0.5:
            e["is_coloc"] = True
        if meta.get("dgidb_druggable"):
            e["is_drug"] = True
        entries.append(e)
    entries.sort(key=lambda x: -x["d"])
    size = _write_json(OUTDIR / "search_index.json", {"genes": entries}, compact=True)
    size_kb = size / 1024
    log(f"  Wrote search_index.json ({size_kb:.1f} KB, {len(entries)} genes)")


def build_manifest(edges: pd.DataFrame, layers: list[dict[str, Any]],
                    summary_df: pd.DataFrame | None, n_gene_graphs: int) -> None:
    log("Building manifest.json ...")
    criterion_status: dict[str, Any] = {}
    if summary_df is not None:
        for _, row in summary_df.iterrows():
            t = str(row.get("type", "unknown"))
            criterion_status[t] = {
                col: _scalar(row[col]) for col in summary_df.columns if col != "type"
            }
    manifest = {
        "version": "v2",
        "architecture": "C",
        "build_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "atlas_parquet": str(ATLAS_PARQUET.relative_to(BASE)),
        "n_edges_total": int(len(edges)),
        "n_edges_per_type": {t: int(c) for t, c in edges["type"].value_counts().items()},
        "n_gene_graphs": n_gene_graphs,
        "layers": [ly["type"] for ly in layers],
        "criterion_status": criterion_status,
    }
    size = _write_json(OUTDIR / "manifest.json", manifest, compact=False)
    log(f"  Wrote manifest.json ({size:,} bytes)")


def build_stage_manifest() -> None:
    """Emit per-F-stage sample counts if the upstream file exists; else a placeholder."""
    stage_path = NET_DIR / "stage" / "stage_sample_counts.csv"
    payload: dict[str, Any] = {"F01": None, "F2": None, "F34": None}
    if stage_path.exists():
        try:
            df = pd.read_csv(stage_path)
            by_stage: dict[str, dict[str, int]] = {}
            stage_col = next((c for c in ("f_stage", "stage") if c in df.columns), None)
            cohort_col = next((c for c in ("cohort", "dataset") if c in df.columns), None)
            n_col = next((c for c in ("n_samples", "n") if c in df.columns), None)
            if stage_col and cohort_col and n_col:
                for _, row in df.iterrows():
                    by_stage.setdefault(str(row[stage_col]), {})[str(row[cohort_col])] = int(row[n_col])
                payload = by_stage
        except Exception as exc:
            log(f"  WARNING: failed to parse stage_sample_counts.csv: {exc}")
    else:
        log(f"  NOTE: {stage_path} not present; writing placeholder stage_manifest.json")
    _write_json(OUTDIR / "stage_manifest.json", payload, compact=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    log("=== Script 294: Export Portal JSON v2 ===")
    log(f"Output dir: {OUTDIR}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    OUTDIR_GRAPHS.mkdir(parents=True, exist_ok=True)

    # Load data
    edges = load_atlas()
    nodes = load_nodes()
    gene_atlas = load_gene_atlas()
    f01_df, f34_df, trans_df, enr01_df, enr34_df = load_communities()

    summary_df = None
    if EDGE_SUMMARY_CSV.exists():
        try:
            summary_df = pd.read_csv(EDGE_SUMMARY_CSV)
            log(f"Loaded edge_annotation_summary.csv ({len(summary_df)} rows)")
        except Exception as exc:
            log(f"WARNING: failed to read edge summary: {exc}")

    # Build community JSONs
    log("Writing community JSONs ...")
    _write_json(OUTDIR / "communities_F01.json",
                build_communities_json(f01_df, enr01_df, trans_df, "F01"), compact=False)
    _write_json(OUTDIR / "communities_F34.json",
                build_communities_json(f34_df, enr34_df, trans_df, "F34"), compact=False)
    _write_json(OUTDIR / "community_transitions.json",
                build_transitions_json(trans_df), compact=False)

    # Prep gene -> community lookup
    def _comm_lookup(df: pd.DataFrame) -> dict[str, dict[str, Any]]:
        if df.empty or "gene" not in df.columns:
            return {}
        keep = [c for c in ("gene", "community_id", "macro_id", "meso_id", "micro_id") if c in df.columns]
        out: dict[str, dict[str, Any]] = {}
        for _, row in df[keep].iterrows():
            g = str(row["gene"])
            entry = {col: _scalar(row[col]) for col in keep if col != "gene"}
            out[g] = entry
        return out

    comm_meta = {"F01": _comm_lookup(f01_df), "F34": _comm_lookup(f34_df)}

    # Attach top-Hallmark per macro community so each gene JSON has its label
    top01 = top_hallmark_per_community(enr01_df)
    top34 = top_hallmark_per_community(enr34_df)
    for g, entry in comm_meta["F01"].items():
        cid = entry.get("macro_id") or entry.get("community_id")
        if cid is not None and int(cid) in top01:
            entry.update({"top_hallmark": top01[int(cid)]["top_hallmark"],
                          "top_hallmark_p": top01[int(cid)]["top_hallmark_p"]})
    for g, entry in comm_meta["F34"].items():
        cid = entry.get("macro_id") or entry.get("community_id")
        if cid is not None and int(cid) in top34:
            entry.update({"top_hallmark": top34[int(cid)]["top_hallmark"],
                          "top_hallmark_p": top34[int(cid)]["top_hallmark_p"]})

    # Build gene-level metadata lookup
    gene_meta: dict[str, dict[str, Any]] = {}
    if not gene_atlas.empty and "human_symbol" in gene_atlas.columns:
        for _, row in gene_atlas.iterrows():
            sym = str(row["human_symbol"])
            gene_meta[sym] = {c: _scalar(row[c]) for c in gene_atlas.columns if c != "human_symbol"}
    # Ensure all node symbols have at least ensembl/biotype from network_nodes.csv
    for _, row in nodes.iterrows():
        sym = str(row.get("human_symbol"))
        meta = gene_meta.setdefault(sym, {})
        for col in ("ensembl_id", "gene_biotype", "dream_logFC", "dream_padj",
                    "is_deg", "coloc_susie_best_pp4", "is_conserved", "sex_class",
                    "zonation_class", "ferroptosis_class", "dgidb_druggable",
                    "attribution_class"):
            if col in row.index and meta.get(col) in (None, ""):
                meta[col] = _scalar(row[col])

    # Layer metadata (uses layer_counts.csv from 290; collapsed layers
    # will be rendered as node badges by the portal, not as layer tabs).
    _layer_counts_tmp = NET_DIR / "layer_counts.csv"
    _collapsed_tmp: set[str] = set()
    if _layer_counts_tmp.exists():
        try:
            _lc_tmp = pd.read_csv(_layer_counts_tmp)
            _collapsed_tmp = set(
                _lc_tmp.loc[_lc_tmp["collapse_to_badge"].fillna(False).astype(bool), "layer"].tolist()
            )
        except Exception:
            pass
    layers = build_layer_metadata(edges, summary_df, _collapsed_tmp)

    # Degree map (unique partners across all edge types)
    log("Computing degree map ...")
    deg_a = edges.groupby("gene_a").size()
    deg_b = edges.groupby("gene_b").size()
    degree_series = deg_a.add(deg_b, fill_value=0).astype(int)
    degree_map = degree_series.to_dict()

    # Phase-3 collapse: load per-layer counts written by 290 and flag
    # any layer below the edge-count floor for node-level badge display.
    collapsed_layers: set[str] = set()
    layer_counts_csv = NET_DIR / "layer_counts.csv"
    if layer_counts_csv.exists():
        try:
            lc = pd.read_csv(layer_counts_csv)
            collapsed_layers = set(
                lc.loc[lc["collapse_to_badge"].fillna(False).astype(bool), "layer"].tolist()
            )
            log(f"Collapse-to-badge layers (from 290): {sorted(collapsed_layers)}")
        except Exception as exc:
            log(f"[warn] failed to read {layer_counts_csv.name}: {exc}")

    # Persist edges to a tmp pickle for worker processes (avoid duplicating in RAM * N workers
    # via fork-safe shared mmapped storage)
    tmp_pkl = OUTDIR / ".edges_tmp.pkl"
    log(f"Serializing edges to {tmp_pkl} for worker processes ...")
    edges.to_pickle(tmp_pkl)

    # Gene list to process
    gene_list = sorted(nodes["human_symbol"].dropna().astype(str).unique().tolist())
    log(f"Writing gene graphs for {len(gene_list):,} genes ...")

    n_workers = min(cpu_count(), 16)
    if n_workers > 1:
        with Pool(n_workers, initializer=_worker_init,
                  initargs=(str(tmp_pkl), gene_meta, comm_meta,
                            str(OUTDIR_GRAPHS), collapsed_layers)) as pool:
            results = pool.map(_build_single_gene, gene_list, chunksize=128)
    else:
        _worker_init(str(tmp_pkl), gene_meta, comm_meta,
                     str(OUTDIR_GRAPHS), collapsed_layers)
        results = [_build_single_gene(g) for g in gene_list]

    # Remove tmp pickle
    try:
        tmp_pkl.unlink()
    except OSError:
        pass

    n_ok = sum(1 for _, _, err in results if err is None)
    n_err = sum(1 for _, _, err in results if err is not None)
    sizes = [s for _, s, err in results if err is None]
    log(f"  Wrote {n_ok} gene graphs, {n_err} errors")
    if sizes:
        log(f"  Sizes KB  min={min(sizes)/1024:.1f}  "
            f"median={np.median(sizes)/1024:.1f}  max={max(sizes)/1024:.1f}")
    if n_err > 0:
        for gene, _, err in results:
            if err is not None:
                log(f"    ERROR {gene}: {err}")
                break  # show first

    # Search index
    build_search_index(nodes, gene_meta, degree_map, comm_meta)

    # Stage manifest (placeholder if missing)
    build_stage_manifest()

    # Manifest (last; depends on n_gene_graphs)
    build_manifest(edges, layers, summary_df, n_ok)

    # Summary
    total_bytes = sum(p.stat().st_size for p in OUTDIR.rglob("*.json"))
    log("")
    log("=== Portal Export v2 Summary ===")
    log(f"Files in {OUTDIR}: {sum(1 for _ in OUTDIR.rglob('*.json'))}")
    log(f"Total size:  {total_bytes / 1024 / 1024:.2f} MB")
    log(f"Done in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    main()
