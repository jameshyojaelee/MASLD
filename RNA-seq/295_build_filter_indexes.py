#!/usr/bin/env python3
"""
295_build_filter_indexes.py
----------------------------
Extend the portal's search_index.json with per-gene filter fields AND
build pre-built inverted filter indexes + pre-computed query bundles
for fast multi-filter queries in the Next.js portal.

Outputs (under results/network/portal_export_v2/):
  - search_index.json                     (rewritten, v3 schema, expanded)
  - filter_indexes/by_f_stage.json
  - filter_indexes/by_coloc.json
  - filter_indexes/by_druggable.json
  - filter_indexes/by_comm_F01.json
  - filter_indexes/by_comm_F34.json
  - filter_indexes/by_sex_class.json
  - filter_indexes/by_cohort_support.json
  - filter_indexes/by_community_pair.json
  - filter_indexes/index_manifest.json
  - query_bundles/Q2_coloc_f2_switch.json
  - query_bundles/Q3_fibrogenic_druggable.json
  - query_bundles/THRB_stratified.json

Usage (SLURM):
  sbatch RNA-seq/295_build_filter_indexes.py
"""
#SBATCH --job-name=net_295_indexes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/net_295_indexes_%j.out
#SBATCH --error=RNA-seq/logs/net_295_indexes_%j.err

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from collections import defaultdict, Counter

import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

PORTAL_DIR   = PROJECT_ROOT / "RNA-seq/results/network/portal_export_v2"
# Gene graphs live at results/network/gene_graphs/ (not under portal_export_v2).
GRAPH_DIR    = PROJECT_ROOT / "RNA-seq/results/network/gene_graphs"
SEARCH_IDX   = PORTAL_DIR / "search_index.json"

ATLAS_PARQUET = PROJECT_ROOT / "RNA-seq/results/network/edge_annotation_atlas.parquet"
ATLAS_CSV     = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
COMM_F01_CSV  = PROJECT_ROOT / "RNA-seq/results/network/communities_f2/communities_F01.csv"
COMM_F34_CSV  = PROJECT_ROOT / "RNA-seq/results/network/communities_f2/communities_F34.csv"
COMM_LABELS_JSON = PROJECT_ROOT / "RNA-seq/results/network/communities_f2/community_labels.json"

OUT_FILTERS = PORTAL_DIR / "filter_indexes"
OUT_QUERIES = PORTAL_DIR / "query_bundles"
OUT_FILTERS.mkdir(parents=True, exist_ok=True)
OUT_QUERIES.mkdir(parents=True, exist_ok=True)


SCHEMA_VERSION = 3


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _log(msg: str, t0: float | None = None) -> None:
    stamp = time.strftime("%H:%M:%S")
    if t0 is not None:
        print(f"[{stamp}] {msg} ({time.time() - t0:.1f}s)", flush=True)
    else:
        print(f"[{stamp}] {msg}", flush=True)


def _save_json(obj, path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, separators=(",", ":"))
    return path.stat().st_size


def _safe_str(x) -> str | None:
    if x is None:
        return None
    if isinstance(x, float) and np.isnan(x):
        return None
    s = str(x).strip()
    return s if s and s.lower() != "nan" else None


# -----------------------------------------------------------------------------
# Step 1. Load existing search index + community tables + atlas
# -----------------------------------------------------------------------------
def load_search_index() -> list[dict]:
    t0 = time.time()
    with open(SEARCH_IDX) as fh:
        obj = json.load(fh)
    genes = obj["genes"] if isinstance(obj, dict) and "genes" in obj else obj
    # dedup (existing file has dup rows e.g. CAPG)
    seen = set()
    uniq = []
    for g in genes:
        s = g.get("s")
        if not s or s in seen:
            continue
        seen.add(s)
        uniq.append(g)
    _log(f"loaded search_index: {len(genes)} rows -> {len(uniq)} unique genes", t0)
    return uniq


def load_community_tables() -> tuple[dict, dict]:
    t0 = time.time()
    f01 = pd.read_csv(COMM_F01_CSV, usecols=["gene", "macro_id"])
    f34 = pd.read_csv(COMM_F34_CSV, usecols=["gene", "macro_id"])
    m01 = dict(zip(f01["gene"].astype(str), f01["macro_id"].astype("Int64")))
    m34 = dict(zip(f34["gene"].astype(str), f34["macro_id"].astype("Int64")))
    _log(f"loaded community tables F01={len(m01)} F34={len(m34)}", t0)
    return m01, m34


def load_atlas_gene_level() -> pd.DataFrame:
    """Load multi_evidence_atlas and extract the filter-relevant columns only."""
    t0 = time.time()
    want = [
        "human_symbol",
        "gene_biotype",
        "bulk_padj",
        "bulk_logFC",
        "is_conserved",
        "coloc_susie_best_pp4",
        "coloc_n_gwas_h4_08",
        "n_coloc_sources",
        "cross_ancestry_coloc_replication",
        "dgidb_druggable",
        "opentargets_drug",
        "sex_class",
        "zonation_class",
        "ferroptosis_class",
        "progression_gene_class",
        "sources_active",
        "f2_inflection_logFC",
        "f2_inflection_padj",
        "spatial_n_datasets",
    ]
    head = pd.read_csv(ATLAS_CSV, nrows=0)
    use = [c for c in want if c in head.columns]
    df = pd.read_csv(ATLAS_CSV, usecols=use, low_memory=False)
    df = df.rename(columns={"human_symbol": "gene"}).set_index("gene")
    _log(f"loaded multi_evidence_atlas: {df.shape}", t0)
    return df


# -----------------------------------------------------------------------------
# Step 2. Derive per-gene filter fields from gene_graphs + atlas
# -----------------------------------------------------------------------------
def _f_stage_from_trajectory(traj: dict | None, g_attrs: dict | None) -> str | None:
    """Return one of: F01_specific, F2_emerging, F2_dissolving, F34_specific,
    stable, NA. Uses edge-level trajectory ratios stored in gene_graph."""
    if not traj:
        return None
    r01 = traj.get("r_F01")
    r2  = traj.get("r_F2")
    r34 = traj.get("r_F34")
    emer = traj.get("emergence_stage")
    if emer and isinstance(emer, str):
        emer_l = emer.lower()
        if "transient_f2" in emer_l or "f2_emerg" in emer_l:
            return "F2_emerging"
        if "f2_dissolv" in emer_l:
            return "F2_dissolving"
        if emer_l in ("f01_specific", "f34_specific", "stable"):
            return emer_l.replace("f01", "F01").replace("f34", "F34")
    # Fallback: infer from ratios
    try:
        vals = {"F01": float(r01 or 0), "F2": float(r2 or 0), "F34": float(r34 or 0)}
    except (TypeError, ValueError):
        return None
    if not any(vals.values()):
        return None
    top = max(vals, key=vals.get)
    return {"F01": "F01_specific", "F2": "F2_emerging", "F34": "F34_specific"}[top]


def derive_f_stage_from_edges(min_edges: int = 3) -> dict[str, dict]:
    """Roll up per-gene f_stage from edge_annotation_atlas.parquet.

    For each gene, find all D-F2 edges it participates in (type == 'D-F2'),
    count emergence_stage categories, and assign modal category as f_stage.
    Genes with fewer than min_edges D-F2 edges -> f_stage = 'undefined'.

    Returns: {gene: {'f_stage': str, 'n_d_f2_edges': int,
                     'n_F2_emerging': int, 'n_F2_dissolving': int, ...}}
    """
    t0 = time.time()
    df = pd.read_parquet(ATLAS_PARQUET, columns=["gene_a", "gene_b", "type", "emergence_stage"])
    df = df[df["type"] == "D-F2"].copy()
    _log(f"  D-F2 edges: {len(df):,}")

    # Stack into long form (one row per (gene, edge_idx))
    a = df[["gene_a", "emergence_stage"]].rename(columns={"gene_a": "gene"})
    b = df[["gene_b", "emergence_stage"]].rename(columns={"gene_b": "gene"})
    long = pd.concat([a, b], ignore_index=True)
    long = long.dropna(subset=["gene"])

    # Pivot: counts per gene per emergence_stage
    counts = long.groupby(["gene", "emergence_stage"]).size().unstack(fill_value=0)
    counts["n_d_f2_edges"] = counts.sum(axis=1)

    stage_cols = [c for c in counts.columns if c != "n_d_f2_edges"]
    def _pick(row):
        if row["n_d_f2_edges"] < min_edges:
            return "undefined"
        sub = row[stage_cols]
        if sub.sum() == 0:
            return "undefined"
        return str(sub.idxmax())
    counts["f_stage"] = counts.apply(_pick, axis=1)

    # Build output dict with renamed columns
    rename = {
        "F2_emerging": "n_F2_emerging",
        "F2_dissolving": "n_F2_dissolving",
        "transient_F2": "n_transient_F2",
        "progressive_up": "n_progressive_up",
        "progressive_down": "n_progressive_down",
        "F34_specific": "n_F34_specific",
    }
    out = {}
    for gene, row in counts.iterrows():
        rec = {"f_stage": row["f_stage"], "n_d_f2_edges": int(row["n_d_f2_edges"])}
        for src, dst in rename.items():
            if src in counts.columns:
                rec[dst] = int(row[src])
        out[str(gene)] = rec
    # Report distribution
    dist = Counter(v["f_stage"] for v in out.values())
    _log(f"derive_f_stage_from_edges: {len(out)} genes; distribution={dict(dist)}", t0)
    return out


def load_gene_graph(symbol: str) -> dict | None:
    p = GRAPH_DIR / f"{symbol}.json"
    if not p.exists():
        return None
    try:
        with open(p) as fh:
            return json.load(fh)
    except Exception:
        return None


def derive_per_gene_fields(
    idx: list[dict],
    atlas: pd.DataFrame,
    m01: dict,
    m34: dict,
    f_stage_map: dict[str, dict] | None = None,
) -> list[dict]:
    """Expand each search_index record with filter fields."""
    t0 = time.time()
    atlas_genes = set(atlas.index.astype(str))

    # COLOC threshold (PP.H4 >= 0.8, or any source_h4_08 hit)
    def _is_coloc(row) -> bool:
        v = row.get("coloc_susie_best_pp4")
        if pd.notna(v) and float(v) >= 0.8:
            return True
        v = row.get("coloc_n_gwas_h4_08")
        if pd.notna(v) and float(v) >= 1:
            return True
        return False

    def _is_druggable(row) -> bool:
        for c in ("dgidb_druggable", "opentargets_drug"):
            v = row.get(c)
            if pd.notna(v):
                if isinstance(v, (int, float)) and float(v) >= 1:
                    return True
                s = str(v).strip().lower()
                if s and s not in ("0", "false", "nan", "none", ""):
                    return True
        return False

    out: list[dict] = []
    n_graphs = 0
    for i, g in enumerate(idx):
        sym = g["s"]
        rec = {
            "s": sym,
            "d": g.get("d", 0),
            "is_deg": bool(g.get("is_deg", False)),
        }
        # atlas-derived fields
        if sym in atlas_genes:
            row = atlas.loc[sym]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            row_d = row.to_dict()
            rec["is_coloc"]     = _is_coloc(row_d)
            rec["is_druggable"] = _is_druggable(row_d)
            sx = _safe_str(row_d.get("sex_class"))
            if sx: rec["sex_class"] = sx
            zn = _safe_str(row_d.get("zonation_class"))
            if zn: rec["zonation"] = zn
            pc = _safe_str(row_d.get("progression_gene_class"))
            if pc: rec["progression_class"] = pc
            fc = _safe_str(row_d.get("ferroptosis_class"))
            if fc: rec["ferroptosis_class"] = fc
            bt = _safe_str(row_d.get("gene_biotype"))
            if bt: rec["biotype"] = bt
            sa = row_d.get("sources_active")
            if pd.notna(sa):
                rec["sources_active"] = int(sa)
            sd = row_d.get("spatial_n_datasets")
            if pd.notna(sd):
                rec["spatial_n_datasets"] = int(sd)
            cc = row_d.get("is_conserved")
            if pd.notna(cc):
                rec["is_conserved"] = bool(cc)
        else:
            rec["is_coloc"]     = False
            rec["is_druggable"] = False

        # community from table
        if sym in m01:
            v = m01[sym]
            if pd.notna(v):
                rec["comm_F01"] = int(v)
        if sym in m34:
            v = m34[sym]
            if pd.notna(v):
                rec["comm_F34"] = int(v)
        # fall back to graph community if not in table
        # (deferred to graph-pass below to save memory if needed)

        # graph-derived fields
        gg = load_gene_graph(sym)
        if gg is not None:
            n_graphs += 1
            ec = gg.get("edge_counts") or {}
            neigh = gg.get("neighbors") or {}
            traj = gg.get("f_stage_trajectory") or {}
            comm = gg.get("community") or {}

            if "comm_F01" not in rec and comm.get("f01") is not None:
                try: rec["comm_F01"] = int(comm["f01"])
                except Exception: pass
            if "comm_F34" not in rec and comm.get("f34") is not None:
                try: rec["comm_F34"] = int(comm["f34"])
                except Exception: pass

            rec["n_d_f2_emerging"]   = len(neigh.get("top_d_f2_emerging", []) or [])
            rec["n_d_f2_dissolving"] = len(neigh.get("top_d_f2_dissolving", []) or [])
            rec["n_d_coloc"]         = len(neigh.get("top_d_coloc", []) or [])
            rec["n_d_lr"]            = len(neigh.get("top_d_lr", []) or [])
            rec["n_contested"]       = len(neigh.get("all_contested", []) or [])
            rec["n_druggable_pairs"] = len(neigh.get("top_druggable_pairs", []) or [])

            for k, v in ec.items():
                try: rec[f"ec_{k}"] = int(v)
                except Exception: pass

            # NOTE: legacy edge-level emergence_stage in gene_graph is NOT a
            # gene-level label; ignore it and use the parquet-derived rollup below.

            # attributes may override is_deg if authoritative
            attrs = gg.get("attributes") or {}
            if "is_deg" in attrs:
                rec["is_deg"] = bool(attrs["is_deg"])

        # Gene-level f_stage rollup from edge parquet (authoritative)
        if f_stage_map and sym in f_stage_map:
            fs = f_stage_map[sym]
            rec["f_stage"] = fs["f_stage"]
            rec["n_d_f2_edges"] = fs["n_d_f2_edges"]
            for k in ("n_F2_emerging", "n_F2_dissolving", "n_transient_F2",
                      "n_progressive_up", "n_progressive_down", "n_F34_specific"):
                if k in fs:
                    rec[k] = fs[k]

        if (i + 1) % 2000 == 0:
            _log(f"  expanded {i + 1}/{len(idx)} genes (graphs={n_graphs})")

        out.append(rec)

    _log(f"expanded {len(out)} genes ({n_graphs} gene_graphs loaded)", t0)
    return out


# -----------------------------------------------------------------------------
# Step 3. Build inverted indexes
# -----------------------------------------------------------------------------
def build_inverted_indexes(expanded: list[dict]) -> dict[str, Path]:
    t0 = time.time()
    by_f_stage = defaultdict(list)
    by_coloc    = []
    by_drug     = []
    by_f01      = defaultdict(list)
    by_f34      = defaultdict(list)
    by_sex      = defaultdict(list)
    by_cohort   = {"ge_5": [], "ge_7": [], "ge_9": []}
    by_pair     = defaultdict(list)   # "F01:X->F34:Y" -> [gene,...]

    for r in expanded:
        s = r["s"]
        if r.get("f_stage"):
            by_f_stage[r["f_stage"]].append(s)
        if r.get("is_coloc"):
            by_coloc.append(s)
        if r.get("is_druggable"):
            by_drug.append(s)
        if "comm_F01" in r:
            by_f01[str(r["comm_F01"])].append(s)
        if "comm_F34" in r:
            by_f34[str(r["comm_F34"])].append(s)
        if r.get("sex_class"):
            by_sex[r["sex_class"]].append(s)
        # Cohort support proxy = sources_active
        sa = r.get("sources_active")
        if sa is not None:
            if sa >= 5: by_cohort["ge_5"].append(s)
            if sa >= 7: by_cohort["ge_7"].append(s)
            if sa >= 9: by_cohort["ge_9"].append(s)
        if "comm_F01" in r and "comm_F34" in r:
            by_pair[f"F01:{r['comm_F01']}->F34:{r['comm_F34']}"].append(s)

    out = {}
    out["by_f_stage.json"]        = _save_json(dict(by_f_stage),   OUT_FILTERS / "by_f_stage.json")
    out["by_coloc.json"]          = _save_json(by_coloc,            OUT_FILTERS / "by_coloc.json")
    out["by_druggable.json"]      = _save_json(by_drug,             OUT_FILTERS / "by_druggable.json")
    out["by_comm_F01.json"]       = _save_json(dict(by_f01),        OUT_FILTERS / "by_comm_F01.json")
    out["by_comm_F34.json"]       = _save_json(dict(by_f34),        OUT_FILTERS / "by_comm_F34.json")
    out["by_sex_class.json"]      = _save_json(dict(by_sex),        OUT_FILTERS / "by_sex_class.json")
    out["by_cohort_support.json"] = _save_json(by_cohort,           OUT_FILTERS / "by_cohort_support.json")
    out["by_community_pair.json"] = _save_json(dict(by_pair),       OUT_FILTERS / "by_community_pair.json")

    _log(
        "inverted indexes: "
        f"f_stage={ {k: len(v) for k,v in by_f_stage.items()} } "
        f"coloc={len(by_coloc)} drug={len(by_drug)} "
        f"F01={len(by_f01)} F34={len(by_f34)} sex={len(by_sex)} "
        f"pairs={len(by_pair)}",
        t0,
    )
    return out


# -----------------------------------------------------------------------------
# Step 4. Pre-computed query bundles
# -----------------------------------------------------------------------------
FIBROGENIC_HINTS = (
    "fibro", "ecm", "collagen", "stellate", "myofib", "wound", "emt",
)


def build_query_bundles(expanded: list[dict]) -> dict[str, Path]:
    t0 = time.time()
    by_sym = {r["s"]: r for r in expanded}

    # Q2: COLOC + F2 switch (emerging OR dissolving)
    q2_strict = [
        r["s"] for r in expanded
        if r.get("is_coloc") and r.get("f_stage") in ("F2_emerging", "F2_dissolving")
    ]
    q2 = q2_strict
    q2_mode = "strict_modal_F2_switch"
    # Fallback: any gene with >=1 F2_emerging or F2_dissolving edge + is_coloc
    if len(q2) < 5:
        q2 = [
            r["s"] for r in expanded
            if r.get("is_coloc")
            and (r.get("n_F2_emerging", 0) + r.get("n_F2_dissolving", 0)) >= 1
        ]
        q2_mode = "any_F2_switch_edge (fallback from strict modal)"
    q2_obj = {
        "query": "COLOC PP4>=0.8 AND F2-switch (emerging|dissolving)",
        "mode": q2_mode,
        "n": len(q2),
        "genes": q2,
    }
    _save_json(q2_obj, OUT_QUERIES / "Q2_coloc_f2_switch.json")

    # Q3: fibrogenic F34 communities ∩ druggable.
    # Primary: use community_labels.json F34 labels containing 'Fibrogenic'.
    # Fallback 1: label top_hallmarks contain EMT / MYOGENESIS / APICAL_JUNCTION.
    # Fallback 2: drop fibrogenic constraint -> druggable in any 'expanding' F34 community.
    fib_comms_from_labels: set[int] = set()
    fallback_emerging_comms: set[int] = set()
    try:
        with open(COMM_LABELS_JSON) as fh:
            comm_labels = json.load(fh)
        f34_labels = comm_labels.get("F34", {})
        hallmark_hits = ("EPITHELIAL_MESENCHYMAL_TRANSITION", "MYOGENESIS", "APICAL_JUNCTION")
        for cid, info in f34_labels.items():
            try:
                cid_int = int(cid)
            except (TypeError, ValueError):
                continue
            label = (info.get("label") or "").lower()
            tops = info.get("top_hallmarks") or []
            if "fibrogenic" in label:
                fib_comms_from_labels.add(cid_int)
            elif any(h in tops for h in hallmark_hits):
                fib_comms_from_labels.add(cid_int)
            shift = (info.get("dominant_shift") or "").lower()
            tclass = (info.get("transition_class") or "").lower()
            if shift == "expanding" or tclass == "expanding":
                fallback_emerging_comms.add(cid_int)
    except Exception as e:
        _log(f"  WARN could not load community_labels.json: {e}")
        comm_labels = {}
    # Identify fibrogenic F34 communities by looking for genes in each community
    # that carry a fibrogenic/ECM signal. We use a simple heuristic: a community
    # is "fibrogenic" if any of its genes has zonation_class != pericentral AND
    # progression_class containing "fibrosis" OR ferroptosis_class set OR
    # community genes include known ECM markers (COL1A1/2, ACTA2, TIMP1, PDGFRB).
    # Ensure we always have some fibrogenic community set; back-fill with ECM markers
    ecm_markers = {"COL1A1", "COL1A2", "COL3A1", "ACTA2", "TIMP1", "PDGFRB", "LOX",
                   "LOXL1", "LOXL2", "TGFB1", "TGFBR1", "PDGFA", "FN1"}
    comm_f34_members = defaultdict(set)
    for r in expanded:
        if "comm_F34" in r:
            comm_f34_members[r["comm_F34"]].add(r["s"])
    fib_comms: set[int] = set(fib_comms_from_labels)
    if not fib_comms:
        for cid, members in comm_f34_members.items():
            if members & ecm_markers:
                fib_comms.add(int(cid))

    q3 = [
        r["s"] for r in expanded
        if r.get("is_druggable") and r.get("comm_F34") in fib_comms
    ]
    q3_mode = "fibrogenic_F34"
    # Fallback: druggable in any 'expanding' F34 community
    if len(q3) < 5 and fallback_emerging_comms:
        q3 = [
            r["s"] for r in expanded
            if r.get("is_druggable") and r.get("comm_F34") in fallback_emerging_comms
        ]
        q3_mode = "druggable_in_expanding_F34 (fibrogenic fallback)"
    # Last resort: just druggable F34 genes with an emerging f_stage
    if len(q3) < 5:
        q3 = [
            r["s"] for r in expanded
            if r.get("is_druggable") and r.get("f_stage") in ("F2_emerging", "F34_specific")
        ]
        q3_mode = "druggable_with_emerging_fstage (loose fallback)"

    q3_obj = {
        "query": "F34 fibrogenic community AND druggable",
        "mode": q3_mode,
        "n_fibrogenic_communities": sorted(int(c) for c in fib_comms),
        "n": len(q3),
        "genes": q3,
    }
    _save_json(q3_obj, OUT_QUERIES / "Q3_fibrogenic_druggable.json")

    # THRB stratified: attach its graph + split neighbor lists by f-stage bucket
    thrb_graph = load_gene_graph("THRB")
    thrb_bundle = {
        "symbol": "THRB",
        "record": by_sym.get("THRB", {}),
        "tabs": {"F01": [], "F2": [], "F34": []},
    }
    if thrb_graph:
        neigh = thrb_graph.get("neighbors") or {}
        # edges of each type already carry stage-specific lists in the graph
        for tab, keys in {
            "F01": ["top_string"],
            "F2":  ["top_d_f2_emerging", "top_d_f2_dissolving"],
            "F34": ["top_d_coloc", "top_d_lr", "top_d_cerna", "top_d_xs",
                    "top_druggable_pairs"],
        }.items():
            merged = []
            for k in keys:
                for item in (neigh.get(k) or []):
                    if isinstance(item, dict):
                        merged.append(item)
                    else:
                        merged.append({"partner": item})
            thrb_bundle["tabs"][tab] = merged
        thrb_bundle["graph_attributes"]  = thrb_graph.get("attributes")
        thrb_bundle["graph_edge_counts"] = thrb_graph.get("edge_counts")
        thrb_bundle["f_stage_trajectory"] = thrb_graph.get("f_stage_trajectory")
    _save_json(thrb_bundle, OUT_QUERIES / "THRB_stratified.json")

    _log(
        f"query bundles: Q2={len(q2)} Q3={len(q3)} THRB=ok (fib_communities={len(fib_comms)})",
        t0,
    )
    return {}


# -----------------------------------------------------------------------------
# Step 5. Write expanded search_index.json (v3)
# -----------------------------------------------------------------------------
def write_expanded_search_index(expanded: list[dict]) -> int:
    t0 = time.time()
    obj = {
        "schema_version": SCHEMA_VERSION,
        "n": len(expanded),
        "genes": expanded,
    }
    size = _save_json(obj, SEARCH_IDX)
    _log(f"wrote search_index.json (v{SCHEMA_VERSION}, {size/1e6:.2f} MB)", t0)
    return size


# -----------------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------------
def write_manifest(expanded: list[dict], filter_sizes: dict[str, int]) -> None:
    total = sum(
        p.stat().st_size for p in OUT_FILTERS.glob("*.json")
    ) + sum(p.stat().st_size for p in OUT_QUERIES.glob("*.json"))
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "n_genes": len(expanded),
        "filter_indexes_bytes": {k: v for k, v in filter_sizes.items()},
        "filter_indexes_total_bytes": total,
        "fields_added": [
            "is_coloc", "is_druggable", "f_stage",
            "comm_F01", "comm_F34", "sex_class", "zonation",
            "progression_class", "ferroptosis_class", "biotype",
            "sources_active", "spatial_n_datasets", "is_conserved",
            "n_d_f2_emerging", "n_d_f2_dissolving", "n_d_coloc",
            "n_d_lr", "n_contested", "n_druggable_pairs",
        ],
        "inverted_indexes": sorted(p.name for p in OUT_FILTERS.glob("*.json")),
        "query_bundles": sorted(p.name for p in OUT_QUERIES.glob("*.json")),
    }
    _save_json(manifest, OUT_FILTERS / "index_manifest.json")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    t_all = time.time()
    _log(f"PORTAL_DIR = {PORTAL_DIR}")

    idx    = load_search_index()
    m01, m34 = load_community_tables()
    atlas  = load_atlas_gene_level()

    f_stage_map = derive_f_stage_from_edges(min_edges=3)

    expanded = derive_per_gene_fields(idx, atlas, m01, m34, f_stage_map=f_stage_map)

    filter_sizes = build_inverted_indexes(expanded)
    build_query_bundles(expanded)
    write_expanded_search_index(expanded)
    write_manifest(expanded, filter_sizes)

    _log(f"ALL DONE in {time.time() - t_all:.1f}s")


if __name__ == "__main__":
    main()
