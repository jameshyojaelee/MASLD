#!/usr/bin/env python3
"""
290_edge_annotation_atlas.py — MASLD Network v2, Architecture C, Phase 3.

Builds the unified edge annotation atlas (publishable deliverable).

Combines Type-S STRING backbone (combined_score/1000 >= 0.7) with Type-D
disease-native edges (F2, COLOC, LR, ceRNA) into a single per-edge table
annotated with per-edge provenance + per-gene flags propagated from the
multi-evidence node atlas + v1 community reference assignments.

Architecture C has NO composite score. Value = annotation richness.

Outputs:
    results/network/edge_annotation_atlas.parquet
    results/network/edge_annotation_atlas.csv
    results/network/edge_annotation_summary.csv

SBATCH:
    #!/bin/bash
    #SBATCH --job-name=net_290_annot
    #SBATCH --partition=cpu
    #SBATCH --mem=64G
    #SBATCH --cpus-per-task=8
    #SBATCH --time=48:00:00
    #SBATCH --output=logs/net_290_annot_%j.out
    #SBATCH --error=logs/net_290_annot_%j.err
    source ~/.bashrc
    micromamba activate spatial
    python 290_edge_annotation_atlas.py
"""

from pathlib import Path
import pandas as pd
import numpy as np
import hashlib
import time
import sys

SCRIPT_VERSION = "2026-04-23"

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
NETWORK_DIR = PROJECT_ROOT / "RNA-seq" / "results" / "network"
ATLAS_NODE_PATH = PROJECT_ROOT / "RNA-seq" / "results" / "multi_evidence" / "multi_evidence_atlas.csv"
NODES_PATH = NETWORK_DIR / "network_nodes.csv"
COMMUNITY_PATH = NETWORK_DIR / "communities" / "community_assignments.csv"

STRING_PATH = NETWORK_DIR / "edges_ppi.csv"
D_F2_LOCO_PATH = NETWORK_DIR / "edges_d_f2_loco.csv"
D_F2_PATH = NETWORK_DIR / "edges_d_f2.csv"
D_COLOC_PATH = NETWORK_DIR / "edges_d_coloc.csv"
D_LR_PATH = NETWORK_DIR / "edges_d_lr.csv"
D_CERNA_PATH = NETWORK_DIR / "edges_d_cerna.csv"

# New Phase 1/2/3 layers (optional; absent files are skipped).
D_F2_CONTINUOUS_PATH   = NETWORK_DIR / "edges_d_f2_continuous.csv"
D_F2_DECONV_FLAG_PATH  = NETWORK_DIR / "edges_d_f2_deconv_flag.csv"
D_REGULON_PATH         = NETWORK_DIR / "edges_d_regulon.csv"
D_STABLE_PATH          = NETWORK_DIR / "edges_d_stable.csv"
D_COLOC_EXPANDED_PATH  = NETWORK_DIR / "edges_d_coloc_expanded.csv"
D_LR_CONSENSUS_PATH    = NETWORK_DIR / "edges_d_lr_consensus.csv"
D_CERNA_TRIPLET_PATH   = NETWORK_DIR / "edges_d_cerna_triplet.csv"
D_XS_PATH              = NETWORK_DIR / "edges_d_xs_conserved.csv"

# When an "expanded" file is present we prefer it over the legacy strict
# variant; the portal atlas propagates both source flags so layer semantics
# remain traceable.
PREFER_EXPANDED = True

# Floor below which a type collapses from a toggleable layer to a
# node-level badge (Phase 3 collapse rule). 294 reads these via the atlas
# summary and writes per-gene badges instead of neighbor buckets.
LAYER_COUNT_FLOOR = 500

OUT_PARQUET = NETWORK_DIR / "edge_annotation_atlas.parquet"
OUT_CSV = NETWORK_DIR / "edge_annotation_atlas.csv"
OUT_SUMMARY = NETWORK_DIR / "edge_annotation_summary.csv"

STRING_THRESHOLD = 0.7


def timed_load(path, loader_msg, loader=pd.read_csv, **kwargs):
    if not Path(path).exists():
        print(f"[WARN] missing: {path} -- skipping {loader_msg}")
        return None
    t0 = time.time()
    df = loader(path, **kwargs)
    print(f"[load] {loader_msg}: {len(df):,} rows from {path.name} in {time.time()-t0:.1f}s")
    return df


def canonicalize(df, a_col="gene_a", b_col="gene_b"):
    df = df.copy()
    df[a_col] = df[a_col].astype(str)
    df[b_col] = df[b_col].astype(str)
    swap = df[a_col] > df[b_col]
    a_new = np.where(swap, df[b_col], df[a_col])
    b_new = np.where(swap, df[a_col], df[b_col])
    df[a_col] = a_new
    df[b_col] = b_new
    df = df[df[a_col] != df[b_col]]
    return df


def detect_pair_cols(df):
    candidates = [
        ("gene_a", "gene_b"),
        ("geneA", "geneB"),
        ("source", "target"),
        ("from", "to"),
        ("protein1", "protein2"),
        ("symbol_a", "symbol_b"),
        ("gene1", "gene2"),
        ("ligand", "receptor"),
        ("mrna_a", "mrna_b"),
    ]
    cols_lower = {c.lower(): c for c in df.columns}
    for a, b in candidates:
        if a.lower() in cols_lower and b.lower() in cols_lower:
            return cols_lower[a.lower()], cols_lower[b.lower()]
    raise ValueError(f"Could not detect gene pair columns among: {list(df.columns)}")


def load_string():
    df = timed_load(STRING_PATH, "Type-S STRING")
    if df is None:
        return None
    a, b = detect_pair_cols(df)
    score_col = None
    for c in ["raw_score", "combined_score", "score", "string_score"]:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        print("[WARN] STRING: no score column detected")
        df["string_score"] = np.nan
    else:
        vals = df[score_col].astype(float)
        if vals.max() > 1.5:
            vals = vals / 1000.0
        df["string_score"] = vals
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    df = df[["gene_a", "gene_b", "string_score"]]
    df = df[df["string_score"] >= STRING_THRESHOLD]
    df = canonicalize(df)
    df = df.groupby(["gene_a", "gene_b"], as_index=False)["string_score"].max()
    df["in_string_ge700"] = True
    print(f"[filt] STRING edges >= {STRING_THRESHOLD}: {len(df):,}")
    return df


def load_d_f2():
    # Phase-2 continuous variant preferred when present; falls back to the
    # gated LOCO file, then the pre-LOCO file. The schema is a superset so
    # the downstream merge is unchanged.
    if PREFER_EXPANDED and D_F2_CONTINUOUS_PATH.exists():
        path = D_F2_CONTINUOUS_PATH
    elif D_F2_LOCO_PATH.exists():
        path = D_F2_LOCO_PATH
    else:
        path = D_F2_PATH
    df = timed_load(path, f"Type-D F2 ({path.name})")
    if df is None:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    keep_cols = ["gene_a", "gene_b"]
    rename_map = {}
    for cand, out in [
        ("r_F01", "r_F01"), ("r_f01", "r_F01"),
        ("r_F2", "r_F2"), ("r_f2", "r_F2"),
        ("r_F34", "r_F34"), ("r_f34", "r_F34"),
        ("delta_max", "delta_max"),
        ("emergence_stage", "emergence_stage"),
        ("fdr_f2", "fdr_f2"), ("fdr_F2", "fdr_f2"), ("padj_f2", "fdr_f2"),
        ("loco_replication_fraction", "loco_replication_fraction"),
        ("loco_rep_frac", "loco_replication_fraction"),
        ("replicates_loco_fraction", "loco_replication_fraction"),
    ]:
        if cand in df.columns and out not in rename_map.values():
            rename_map[cand] = out
    df = df.rename(columns=rename_map)
    for want in ["r_F01", "r_F2", "r_F34", "delta_max", "emergence_stage",
                 "fdr_f2", "loco_replication_fraction"]:
        if want not in df.columns:
            df[want] = np.nan
        keep_cols.append(want)
    df = df[keep_cols]
    df = canonicalize(df)
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_d_coloc():
    # Phase-3 expanded variant (PP4 >= 0.3 shared-locus + eQTL sentinel)
    # preferred when present. Adds `channel` / `pp4_thresh` columns.
    path = D_COLOC_EXPANDED_PATH if (PREFER_EXPANDED and
                                     D_COLOC_EXPANDED_PATH.exists()) else D_COLOC_PATH
    df = timed_load(path, f"Type-D COLOC ({path.name})")
    if df is None:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    rename_map = {}
    for cand, out in [("pp4_min", "pp4_min"), ("min_pp4", "pp4_min"), ("pp4", "pp4_min"),
                      ("gwas_list", "gwas_list"), ("gwas", "gwas_list"),
                      ("locus_id", "locus_id"), ("locus", "locus_id")]:
        if cand in df.columns and out not in rename_map.values():
            rename_map[cand] = out
    df = df.rename(columns=rename_map)
    for want in ["pp4_min", "gwas_list", "locus_id"]:
        if want not in df.columns:
            df[want] = np.nan
    df = df[["gene_a", "gene_b", "pp4_min", "gwas_list", "locus_id"]]
    df = canonicalize(df)
    agg = df.groupby(["gene_a", "gene_b"], as_index=False).agg({
        "pp4_min": "max",
        "gwas_list": lambda s: ",".join(sorted(set(str(x) for x in s if pd.notna(x)))),
        "locus_id": lambda s: ",".join(sorted(set(str(x) for x in s if pd.notna(x)))),
    })
    return agg


def load_d_lr():
    path = (D_LR_CONSENSUS_PATH if (PREFER_EXPANDED and D_LR_CONSENSUS_PATH.exists())
            else D_LR_PATH)
    df = timed_load(path, f"Type-D LR ({path.name})")
    if df is None:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    rename_map = {}
    for cand, out in [("score_diff_lr", "score_diff_lr"), ("score_diff", "score_diff_lr"),
                      ("lr_score_diff", "score_diff_lr"),
                      ("cell_type_pairs", "cell_type_pairs"), ("celltype_pairs", "cell_type_pairs"),
                      ("emergence_stage_lr", "emergence_stage_lr"),
                      ("emergence_stage", "emergence_stage_lr")]:
        if cand in df.columns and out not in rename_map.values():
            rename_map[cand] = out
    df = df.rename(columns=rename_map)
    for want in ["score_diff_lr", "cell_type_pairs", "emergence_stage_lr"]:
        if want not in df.columns:
            df[want] = np.nan
    df = df[["gene_a", "gene_b", "score_diff_lr", "cell_type_pairs", "emergence_stage_lr"]]
    df = canonicalize(df)
    agg = df.groupby(["gene_a", "gene_b"], as_index=False).agg({
        "score_diff_lr": "max",
        "cell_type_pairs": lambda s: ",".join(sorted(set(str(x) for x in s if pd.notna(x)))),
        "emergence_stage_lr": lambda s: ",".join(sorted(set(str(x) for x in s if pd.notna(x)))),
    })
    return agg


def load_d_cerna():
    path = (D_CERNA_TRIPLET_PATH if (PREFER_EXPANDED and D_CERNA_TRIPLET_PATH.exists())
            else D_CERNA_PATH)
    df = timed_load(path, f"Type-D ceRNA ({path.name})")
    if df is None or len(df) == 0:
        print("[info] ceRNA: 0 edges")
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    rename_map = {}
    for cand, out in [("n_shared_mirnas", "n_shared_mirnas"), ("n_shared", "n_shared_mirnas"),
                      ("shared_mirnas", "shared_mirnas"), ("shared_miRNAs", "shared_mirnas")]:
        if cand in df.columns and out not in rename_map.values():
            rename_map[cand] = out
    df = df.rename(columns=rename_map)
    for want in ["n_shared_mirnas", "shared_mirnas"]:
        if want not in df.columns:
            df[want] = np.nan
    df = df[["gene_a", "gene_b", "n_shared_mirnas", "shared_mirnas"]]
    df = canonicalize(df)
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_d_regulon():
    """D-REG: directed TF -> target edges from SCENIC+ regulons (Script 293)."""
    df = timed_load(D_REGULON_PATH, "Type-D REG")
    if df is None or len(df) == 0:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    keep = ["gene_a", "gene_b"]
    for cand, out in [("tf", "reg_tf"), ("target", "reg_target"),
                      ("regulon_id", "regulon_id"),
                      ("regulon_source", "regulon_source"),
                      ("regulon_activity_diff", "regulon_activity_diff"),
                      ("activity_padj", "regulon_activity_padj"),
                      ("directed", "regulon_directed")]:
        if cand in df.columns:
            df = df.rename(columns={cand: out})
            keep.append(out)
    df = df[keep]
    # D-REG is DIRECTED (TF -> target). Do NOT canonicalize gene order —
    # but the atlas merge is undirected-keyed; to survive that, we emit
    # the canonical pair and keep tf/target as metadata columns.
    df["regulon_directed"] = True
    swap = df["gene_a"] > df["gene_b"]
    ga_canon = np.where(swap, df["gene_b"], df["gene_a"])
    gb_canon = np.where(swap, df["gene_a"], df["gene_b"])
    df["gene_a"] = ga_canon
    df["gene_b"] = gb_canon
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_d_stable():
    """D-STABLE: stable co-expression across F01/F2/F34 (Script 286)."""
    df = timed_load(D_STABLE_PATH, "Type-D STABLE")
    if df is None or len(df) == 0:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    keep = ["gene_a", "gene_b"]
    rename_map = {
        "r_F01": "stable_r_F01", "r_F2": "stable_r_F2", "r_F34": "stable_r_F34",
        "r_mean": "stable_r_mean", "min_r": "stable_min_r", "max_r": "stable_max_r",
        "range_r": "stable_range_r", "delta_max": "stable_delta_max",
    }
    for src, dst in rename_map.items():
        if src in df.columns:
            df = df.rename(columns={src: dst})
            keep.append(dst)
    df = df[keep]
    df = canonicalize(df)
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_d_xs():
    """D-XS: cross-species conserved edges (Script 297)."""
    df = timed_load(D_XS_PATH, "Type-D XS")
    if df is None or len(df) == 0:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    keep = ["gene_a", "gene_b"]
    for cand, out in [("channel", "xs_channel"),
                      ("module_id", "xs_module_id"),
                      ("preservation_Z", "xs_preservation_Z"),
                      ("preservation_Z_mean", "xs_preservation_Z_mean"),
                      ("celltype", "xs_celltype")]:
        if cand in df.columns:
            df = df.rename(columns={cand: out})
            keep.append(out)
    df = df[keep]
    df = canonicalize(df)
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_d_f2_deconv_flag():
    """Annotation-only: merges deconv_concordance + composition_driven onto
    existing D-F2 edges. Returns a canonicalized table with (gene_a,
    gene_b) + flag columns. When this table is merged, rows without a
    match just get NaNs — no new edges are created."""
    df = timed_load(D_F2_DECONV_FLAG_PATH, "D-F2 deconv flag")
    if df is None or len(df) == 0:
        return None
    a, b = detect_pair_cols(df)
    df = df.rename(columns={a: "gene_a", b: "gene_b"})
    keep = ["gene_a", "gene_b"]
    for cand, out in [("delta_max", "f2_delta_max_raw"),
                      ("delta_max_resid", "f2_delta_max_resid"),
                      ("deconv_concordance", "f2_deconv_concordance"),
                      ("composition_driven", "f2_composition_driven")]:
        if cand in df.columns:
            df = df.rename(columns={cand: out})
            keep.append(out)
    df = df[keep]
    df = canonicalize(df)
    df = df.drop_duplicates(["gene_a", "gene_b"])
    return df


def load_node_atlas():
    t0 = time.time()
    if not ATLAS_NODE_PATH.exists():
        print(f"[WARN] node atlas missing: {ATLAS_NODE_PATH}")
        return None
    wanted = ["symbol", "gene_symbol", "gene", "hgnc_symbol", "human_symbol",
              "dgidb_druggable", "coloc_susie_best_pp4", "sex_class",
              "is_conserved", "ferroptosis_class", "zonation_class",
              "attribution_class"]
    header = pd.read_csv(ATLAS_NODE_PATH, nrows=0).columns.tolist()
    usecols = [c for c in wanted if c in header]
    df = pd.read_csv(ATLAS_NODE_PATH, usecols=usecols, low_memory=False)
    sym_col = None
    for c in ["symbol", "gene_symbol", "gene", "hgnc_symbol", "human_symbol"]:
        if c in df.columns:
            sym_col = c
            break
    if sym_col is None:
        print("[WARN] no symbol column in node atlas")
        return None
    df = df.rename(columns={sym_col: "symbol"})
    df["symbol"] = df["symbol"].astype(str)
    df = df.drop_duplicates("symbol")
    print(f"[load] node atlas: {len(df):,} genes, {len(df.columns)} cols in {time.time()-t0:.1f}s")
    return df


def load_communities():
    if not COMMUNITY_PATH.exists():
        print(f"[WARN] community assignments missing: {COMMUNITY_PATH}")
        return None
    df = pd.read_csv(COMMUNITY_PATH)
    sym_col = None
    for c in ["symbol", "gene", "gene_symbol", "node"]:
        if c in df.columns:
            sym_col = c
            break
    if sym_col is None:
        print("[WARN] no symbol column in community file")
        return None
    df = df.rename(columns={sym_col: "symbol"})
    keep = ["symbol"]
    for c in ["community_macro", "community_meso", "community_micro",
              "macro", "meso", "micro"]:
        if c in df.columns:
            keep.append(c)
    df = df[keep].drop_duplicates("symbol")
    df = df.rename(columns={"macro": "community_macro",
                            "meso": "community_meso",
                            "micro": "community_micro"})
    print(f"[load] communities: {len(df):,} nodes, cols={list(df.columns)}")
    return df


def merge_edges(base, incoming, source_tag):
    if incoming is None or len(incoming) == 0:
        return base
    incoming = incoming.copy()
    incoming["__src_" + source_tag] = True
    if base is None:
        return incoming
    merged = pd.merge(base, incoming, on=["gene_a", "gene_b"], how="outer")
    return merged


def propagate_gene_annotations(edges, atlas):
    if atlas is None:
        for col in ["celltype_driver", "coloc_linked_gene", "druggable",
                    "sex_class", "conserved_mouse", "ferroptosis", "zonation"]:
            edges[f"{col}_a"] = np.nan
            edges[f"{col}_b"] = np.nan
        edges["coloc_linked"] = False
        edges["druggable_pair"] = False
        return edges

    amap = atlas.set_index("symbol")

    def gene_col(series, atlas_col, default=np.nan):
        if atlas_col not in amap.columns:
            return pd.Series([default] * len(series), index=series.index)
        return series.map(amap[atlas_col])

    edges["celltype_driver_a"] = gene_col(edges["gene_a"], "attribution_class", "unknown")
    edges["celltype_driver_b"] = gene_col(edges["gene_b"], "attribution_class", "unknown")
    edges["celltype_driver_a"] = edges["celltype_driver_a"].fillna("unknown")
    edges["celltype_driver_b"] = edges["celltype_driver_b"].fillna("unknown")

    pp4_a = gene_col(edges["gene_a"], "coloc_susie_best_pp4")
    pp4_b = gene_col(edges["gene_b"], "coloc_susie_best_pp4")
    edges["coloc_pp4_a"] = pd.to_numeric(pp4_a, errors="coerce")
    edges["coloc_pp4_b"] = pd.to_numeric(pp4_b, errors="coerce")
    edges["coloc_linked"] = ((edges["coloc_pp4_a"] > 0.5) | (edges["coloc_pp4_b"] > 0.5)).fillna(False)

    drug_a = gene_col(edges["gene_a"], "dgidb_druggable").fillna(False).astype(bool)
    drug_b = gene_col(edges["gene_b"], "dgidb_druggable").fillna(False).astype(bool)
    edges["druggable_a"] = drug_a
    edges["druggable_b"] = drug_b
    edges["druggable_pair"] = drug_a & drug_b

    edges["sex_class_a"] = gene_col(edges["gene_a"], "sex_class")
    edges["sex_class_b"] = gene_col(edges["gene_b"], "sex_class")

    edges["conserved_mouse_a"] = gene_col(edges["gene_a"], "is_conserved")
    edges["conserved_mouse_b"] = gene_col(edges["gene_b"], "is_conserved")

    edges["ferroptosis_a"] = gene_col(edges["gene_a"], "ferroptosis_class")
    edges["ferroptosis_b"] = gene_col(edges["gene_b"], "ferroptosis_class")

    edges["zonation_a"] = gene_col(edges["gene_a"], "zonation_class")
    edges["zonation_b"] = gene_col(edges["gene_b"], "zonation_class")

    return edges


def propagate_communities(edges, comm):
    for suffix in ["a", "b"]:
        for col in ["community_macro", "community_meso", "community_micro"]:
            edges[f"{col}_{suffix}"] = np.nan
    if comm is None:
        return edges
    cmap = comm.set_index("symbol")
    for suffix, gcol in [("a", "gene_a"), ("b", "gene_b")]:
        for col in ["community_macro", "community_meso", "community_micro"]:
            if col in cmap.columns:
                edges[f"{col}_{suffix}"] = edges[gcol].map(cmap[col])
    return edges


SRC_LAYER_MAP = {
    "S":        "__src_S",
    "D-F2":     "__src_F2",
    "D-COLOC":  "__src_COLOC",
    "D-LR":     "__src_LR",
    "D-ceRNA":  "__src_ceRNA",
    # New layers added 2026-04-23:
    "D-REG":    "__src_REG",
    "D-STABLE": "__src_STABLE",
    "D-XS":     "__src_XS",
}


def compute_type_labels(edges):
    """Vectorized type-label + contested flag. Row-wise apply() at 15M rows
    is unusable (hours of Python iteration); build the comma-joined tag
    string with numpy char ops instead."""
    src_cols = dict(SRC_LAYER_MAP)
    for c in src_cols.values():
        if c not in edges.columns:
            edges[c] = False
        edges[c] = edges[c].fillna(False).astype(bool)

    n = len(edges)
    # Comma-join the set of layers whose source flag is True, purely with
    # numpy string ops. For each layer: prefix with "," if the flag is on,
    # else "". Then concatenate all per-layer string arrays and lstrip.
    joined = np.full(n, "", dtype=object)
    for layer, col in src_cols.items():
        mask = edges[col].to_numpy()
        addend = np.where(mask, "," + layer, "")
        joined = joined + addend   # numpy object string concat, vectorized
    # Strip leading comma on non-empty entries
    joined = np.where(joined.astype(str) != "",
                       np.char.lstrip(joined.astype(str), ","),
                       "")
    edges["type"] = joined

    edges["in_string_ge700"] = edges["__src_S"]
    d_cols = [v for k, v in src_cols.items() if k != "S"]
    d_any = edges[d_cols].any(axis=1)
    edges["is_contested"] = edges["__src_S"] & d_any
    return edges


def compute_provenance_hash(edges):
    """Vectorized hash: build one concatenated key column then hash in a
    list comprehension (fast enough on 15M rows, ~1M hashes/sec)."""
    version = SCRIPT_VERSION
    src_cols = list(SRC_LAYER_MAP.values())

    n = len(edges)
    # Compact source tag: one char per layer, "1" or "0". Concatenate via
    # a Python loop over columns (ufunc.reduce on np.char.add hits dtype
    # mismatches; a plain loop is ~1s for 15M rows x 8 layers).
    src_tag = np.full(n, "", dtype=object)
    for c in src_cols:
        bit = np.where(edges[c].to_numpy(), "1", "0")
        src_tag = src_tag + bit

    ga = edges["gene_a"].astype(str).to_numpy()
    gb = edges["gene_b"].astype(str).to_numpy()
    # Build the key with object-dtype string ops (numpy char ops stumble on
    # variable widths across 15M rows).
    version_tail = "|" + version
    keys = ga + "|" + gb + "|" + src_tag.astype(str) + version_tail

    edges["provenance_hash"] = [hashlib.md5(k.encode()).hexdigest() for k in keys]
    return edges


def annotation_completeness(edges):
    annotation_fields = [
        "string_score", "r_F2", "delta_max", "emergence_stage", "fdr_f2",
        "loco_replication_fraction", "pp4_min", "gwas_list", "locus_id",
        "score_diff_lr", "cell_type_pairs", "emergence_stage_lr",
        "n_shared_mirnas", "shared_mirnas",
        "celltype_driver_a", "celltype_driver_b",
        "coloc_linked", "druggable_a", "druggable_b", "druggable_pair",
        "sex_class_a", "sex_class_b",
        "conserved_mouse_a", "conserved_mouse_b",
        "ferroptosis_a", "ferroptosis_b",
        "zonation_a", "zonation_b",
        "community_macro_a", "community_macro_b",
        "community_meso_a", "community_meso_b",
        "community_micro_a", "community_micro_b",
    ]
    rows = []
    n = len(edges)
    trivial = {"unknown", "none", "na", "nan", "", "false", "0"}

    def is_non_trivial(val):
        if pd.isna(val):
            return False
        if isinstance(val, bool):
            return bool(val)
        s = str(val).strip().lower()
        if s in trivial:
            return False
        return True

    per_edge_counts = np.zeros(n, dtype=int)
    trivial_strs = np.array(list(trivial))
    for f in annotation_fields:
        if f not in edges.columns:
            rows.append({"field_name": f, "n_non_trivial_edges": 0, "pct_non_trivial": 0.0})
            continue
        col = edges[f]
        if col.dtype == bool:
            mask = col.fillna(False).astype(bool).to_numpy()
        else:
            notna = col.notna().to_numpy()
            # Stringify non-null values and strip+lower-case en-masse
            strs = col.astype(object).fillna("").astype(str).str.strip().str.lower().to_numpy()
            trivial_mask = np.isin(strs, trivial_strs)
            mask = notna & ~trivial_mask
        per_edge_counts = per_edge_counts + mask.astype(int)
        rows.append({"field_name": f, "n_non_trivial_edges": int(mask.sum()),
                     "pct_non_trivial": float(mask.mean() * 100.0)})
    summary = pd.DataFrame(rows)

    median_count = float(np.median(per_edge_counts)) if n else 0.0
    pct_ge3 = float((per_edge_counts >= 3).mean() * 100.0) if n else 0.0
    summary_row = pd.DataFrame([{
        "field_name": "__SUMMARY__",
        "n_non_trivial_edges": n,
        "pct_non_trivial": pct_ge3,
        "median_annotation_count_per_edge": median_count,
        "pct_edges_with_ge_3_annotations": pct_ge3,
    }])
    summary = pd.concat([summary, summary_row], ignore_index=True)
    return summary, pct_ge3, median_count


def main():
    t_start = time.time()
    print(f"=== 290_edge_annotation_atlas.py (v{SCRIPT_VERSION}) ===")
    NETWORK_DIR.mkdir(parents=True, exist_ok=True)

    string_df = load_string()
    f2_df = load_d_f2()
    coloc_df = load_d_coloc()
    lr_df = load_d_lr()
    cerna_df = load_d_cerna()
    # New Phase 1/3 layers:
    reg_df = load_d_regulon()
    stable_df = load_d_stable()
    xs_df = load_d_xs()
    # Annotation-only (merges onto D-F2 rows, no new edges):
    f2_deconv_df = load_d_f2_deconv_flag()

    base = None
    if string_df is not None:
        string_df["__src_S"] = True
        base = merge_edges(None, string_df, "S")
    if f2_df is not None:
        f2_df["__src_F2"] = True
        base = merge_edges(base, f2_df, "F2")
    if coloc_df is not None:
        coloc_df["__src_COLOC"] = True
        base = merge_edges(base, coloc_df, "COLOC")
    if lr_df is not None:
        lr_df["__src_LR"] = True
        base = merge_edges(base, lr_df, "LR")
    if cerna_df is not None:
        cerna_df["__src_ceRNA"] = True
        base = merge_edges(base, cerna_df, "ceRNA")
    if reg_df is not None:
        reg_df["__src_REG"] = True
        base = merge_edges(base, reg_df, "REG")
    if stable_df is not None:
        stable_df["__src_STABLE"] = True
        base = merge_edges(base, stable_df, "STABLE")
    if xs_df is not None:
        xs_df["__src_XS"] = True
        base = merge_edges(base, xs_df, "XS")

    if base is None or len(base) == 0:
        print("[FATAL] no edges loaded from any source")
        sys.exit(1)

    # Annotation-only merge for the F2 deconv flag (left-join).
    if f2_deconv_df is not None:
        before_cols = set(base.columns)
        base = pd.merge(base, f2_deconv_df, on=["gene_a", "gene_b"], how="left")
        added = sorted(set(base.columns) - before_cols)
        print(f"[merge] F2 deconv flag: added cols {added}")

    base = base.drop_duplicates(["gene_a", "gene_b"]).reset_index(drop=True)
    print(f"[merge] union edges: {len(base):,}")

    # Per-layer edge counts (post-merge). Emitted to layer_counts.csv so
    # 294 and the portal UI can apply the collapse-to-badge rule.
    layer_counts = {}
    for layer, col in SRC_LAYER_MAP.items():
        if col in base.columns:
            n = int(base[col].fillna(False).astype(bool).sum())
        else:
            n = 0
        layer_counts[layer] = n
    collapse_to_badge = [layer for layer, n in layer_counts.items()
                         if n > 0 and n < LAYER_COUNT_FLOOR and layer != "S"]
    print("[layers] edge counts by source layer:")
    for layer, n in layer_counts.items():
        marker = " -> BADGE" if layer in collapse_to_badge else ""
        print(f"    {layer:<10s} {n:>10,d}{marker}")
    layer_counts_df = pd.DataFrame([
        {"layer": layer, "n_edges": n,
         "collapse_to_badge": layer in collapse_to_badge}
        for layer, n in layer_counts.items()
    ])
    layer_counts_df.to_csv(NETWORK_DIR / "layer_counts.csv", index=False)
    print(f"[write] {(NETWORK_DIR / 'layer_counts.csv').name}")

    base = compute_type_labels(base)

    atlas = load_node_atlas()
    base = propagate_gene_annotations(base, atlas)

    comm = load_communities()
    base = propagate_communities(base, comm)

    base = compute_provenance_hash(base)

    internal_cols = [c for c in base.columns if c.startswith("__src_")]
    front = ["gene_a", "gene_b", "type", "in_string_ge700", "is_contested", "string_score"]
    rest = [c for c in base.columns if c not in front and c not in internal_cols and c != "provenance_hash"]
    final_cols = front + rest + ["provenance_hash"]
    final_cols = [c for c in final_cols if c in base.columns]
    out = base[final_cols].copy()

    t0 = time.time()
    out.to_parquet(OUT_PARQUET, index=False)
    print(f"[write] {OUT_PARQUET.name}: {len(out):,} rows in {time.time()-t0:.1f}s")
    t0 = time.time()
    out.to_csv(OUT_CSV, index=False)
    print(f"[write] {OUT_CSV.name}: {len(out):,} rows in {time.time()-t0:.1f}s")

    summary, pct_ge3, median_count = annotation_completeness(out)
    summary.to_csv(OUT_SUMMARY, index=False)
    print(f"[write] {OUT_SUMMARY.name}: {len(summary):,} rows")

    print(f"\n[summary] median annotations/edge: {median_count:.1f}")
    print(f"[summary] % edges with >= 3 annotations: {pct_ge3:.1f}%")
    verdict = "PASS" if pct_ge3 >= 80.0 else "FAIL"
    print(f"CRITERION 6: {verdict} (threshold: 80% edges with >=3 non-trivial annotations)")

    n_contested = int(out["is_contested"].sum())
    n_string = int(out["in_string_ge700"].sum())
    n_d_only = len(out) - n_string
    print(f"[breakdown] Type-S only: {n_string - n_contested:,} | "
          f"Contested (S+D): {n_contested:,} | D-only: {n_d_only:,}")
    print(f"[done] total runtime: {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    main()
