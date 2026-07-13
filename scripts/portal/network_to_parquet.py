#!/usr/bin/env python3
"""Collapse the 15.3M-edge network annotation atlas -> two web-scale parquets.

Replaces the retired 16,601 per-gene-graph JSONs consumed by the /network page
with:
  network_edges.parquet   (filtered subset, data contract Section 6.1)
  network_nodes.parquet    (data contract Section 6.2)

Env: spatial (rnaseq has a NumPy2/pandas ABI break).
HEAVY read: source parquet is 1.5 GB / 15,295,605 rows -> COMPUTE NODE ONLY.

--------------------------------------------------------------------------------
WEB EDGE FILTER  (records data-contract Section 8 item 7 decision)
--------------------------------------------------------------------------------
An edge is KEPT iff:

    ( in_string_ge700 == True
      OR coloc_linked == True
      OR is_contested == True
      OR abs(delta_max) >= TAU )
    AND gene_a in NODE_SET
    AND gene_b in NODE_SET

then a per-node top-K cap (K=30 strongest edges/node, edge kept if it is in
gene_a's OR gene_b's top-K) preserves force-graph connectivity.

NODE_SET (bounded, from network_nodes.csv + atlas):
    union of {
        is_deg == True                         (network_nodes.csv, legacy flag),
        coloc_susie_best_pp4 >= 0.5            (network_nodes.csv),
        dgidb_druggable == True                (network_nodes.csv),
        top TOP_LAYERS by layers_active        (multi_evidence_atlas.csv)
    }

TAU is auto-tuned so the written network_edges.parquet is < SIZE_TARGET_MB.
The chosen TAU, the exact predicate, and the final edge count are PRINTED and
recorded in the run log.

--------------------------------------------------------------------------------
NODE PARQUET  (data contract Section 6.2)
--------------------------------------------------------------------------------
network_nodes.csv predates the 2026-06-08 C2 canonical swap: it carries retired
`dream_logFC/dream_padj/dream_tstat`. Per the contract these columns are RENAMED
to `bulk_logFC/bulk_padj/bulk_tstat`. To also resolve the value-staleness trap
(contract Section 8 item 9 -- "ideally regenerate from the current atlas"), the
renamed columns and `is_deg` are REFRESHED from the canonical
multi_evidence_atlas.csv (bulk_* + effect-size-aware interval-null FDR gate
`bulk_treat_fdr < 0.05`), falling back to the CSV value only for genes absent
from the atlas. No `dream_*` column name survives.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# ---- tunables ----------------------------------------------------------------
K_PER_NODE = 30          # per-node top-K strongest-edge cap
TOP_LAYERS = 3000        # top-N genes by atlas layers_active added to NODE_SET
SIZE_TARGET_MB = 9.0     # write budget for network_edges.parquet (< ~10 MB)
COLOC_PP4_MIN = 0.5      # node-set inclusion threshold on coloc_susie_best_pp4
TAU_FLOOR = 0.25         # streaming floor on |delta_max| (memory bound; TAU >= this)
TAU_GRID = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]

# §6.1 output columns (also exactly the columns needed for filtering).
EDGE_COLS = [
    "gene_a", "gene_b", "type",
    "in_string_ge700", "is_contested", "string_score",
    "r_F01", "r_F2", "r_F34", "delta_max", "emergence_stage",
    "fdr_f2", "loco_replication_fraction",
    "pp4_min", "gwas_list",
    "score_diff_lr", "cell_type_pairs",
    "reg_tf", "reg_target", "regulon_activity_diff",
    "coloc_linked", "druggable_pair",
    "community_macro_a", "community_macro_b",
]

# §6.2 node output columns (post-rename).
NODE_COLS = [
    "symbol", "ensembl_id", "gene_biotype",
    "bulk_logFC", "bulk_padj", "bulk_tstat",
    "is_deg", "coloc_susie_best_pp4", "is_conserved_core",
    "sex_class", "zonation_class", "ferroptosis_class",
    "dgidb_druggable", "n_coloc_sources", "attribution_class", "mouse_ortholog",
]


def default_root() -> Path:
    return Path(os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    ))


def as_bool(s: pd.Series) -> pd.Series:
    """R-style 'TRUE'/'FALSE'/'' strings (or real bools) -> Python bool."""
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.upper().eq("TRUE")


# Raw `type` is a comma-joined SET of channel tokens (an edge can belong to
# several channels), e.g. "S,D-XS", "D-COLOC,D-STABLE". Tokens: S (STRING/PPI),
# D-STABLE / D-F2 (differential coexpr), D-COLOC, D-LR, D-XS (cross-species),
# D-REG (regulon), D-ceRNA. We surface ONE primary channel per edge, most
# specific wins over generic coexpression.
_CHANNEL_PRIORITY = ["coloc", "lr", "regulon", "cerna", "xspecies", "coexpr"]
_TOKEN_MAP = {
    "COLOC": "coloc",
    "LR": "lr",
    "REG": "regulon", "REGULON": "regulon",
    "CERNA": "cerna", "MIRNA": "cerna",
    "XS": "xspecies", "XSPECIES": "xspecies", "CONSERVED": "xspecies",
    "STABLE": "coexpr", "F2": "coexpr", "S": "coexpr",
    "COEXPR": "coexpr", "PPI": "coexpr",
}


def map_channel(t) -> str | None:
    """Reduce a composite raw `type` string to one canonical channel (priority-ranked)."""
    if t is None:
        return None
    u = str(t).strip().upper()
    if u in ("", "NAN", "NONE"):
        return None
    found: set[str] = set()
    for tok in u.split(","):
        tok = tok.strip()
        if tok.startswith("D-"):
            tok = tok[2:]
        ch = _TOKEN_MAP.get(tok)
        if ch is None:
            if "COLOC" in tok:
                ch = "coloc"
            elif "CERNA" in tok:
                ch = "cerna"
            elif tok.startswith("XS"):
                ch = "xspecies"
            elif tok.startswith("REG"):
                ch = "regulon"
            elif tok == "LR":
                ch = "lr"
            else:
                ch = "coexpr"
        found.add(ch)
    for ch in _CHANNEL_PRIORITY:
        if ch in found:
            return ch
    return "coexpr"


def build_node_set(nodes_raw: pd.DataFrame, atlas: pd.DataFrame) -> set[str]:
    ns: set[str] = set()
    ns |= set(nodes_raw.loc[as_bool(nodes_raw["is_deg"]), "human_symbol"])
    pp4 = pd.to_numeric(nodes_raw["coloc_susie_best_pp4"], errors="coerce")
    ns |= set(nodes_raw.loc[pp4 >= COLOC_PP4_MIN, "human_symbol"])
    ns |= set(nodes_raw.loc[as_bool(nodes_raw["dgidb_druggable"]), "human_symbol"])
    la = atlas[["human_symbol", "layers_active"]].copy()
    la["layers_active"] = pd.to_numeric(la["layers_active"], errors="coerce")
    top = la.dropna(subset=["layers_active"]).sort_values(
        "layers_active", ascending=False).head(TOP_LAYERS)["human_symbol"]
    ns |= set(top)
    ns.discard(np.nan)
    ns = {s for s in ns if isinstance(s, str) and s}
    # Referential integrity: the shipped node table IS network_nodes.csv, so a
    # NODE_SET member that is not a node (atlas-only top-layers gene) would leave
    # dangling edge endpoints. Intersect to the node universe.
    ns &= set(nodes_raw["human_symbol"].dropna())
    print(f"  NODE_SET size = {len(ns):,}  "
          f"(is_deg={int(as_bool(nodes_raw['is_deg']).sum())}, "
          f"coloc>={COLOC_PP4_MIN}:{int((pp4>=COLOC_PP4_MIN).sum())}, "
          f"druggable={int(as_bool(nodes_raw['dgidb_druggable']).sum())}, "
          f"top{TOP_LAYERS}_layers)")
    return ns


def stream_candidate_edges(src: Path, node_set: set[str]) -> pd.DataFrame:
    """One projected pass over the 15.3M edges.

    Keeps rows with both endpoints in NODE_SET that carry ANY web signal, i.e.
    (in_string_ge700 | coloc_linked | is_contested | |delta_max| >= TAU_FLOOR).
    Rows with no signal can never pass at any TAU >= TAU_FLOOR, so dropping them
    here is lossless for the web filter and bounds memory.
    """
    pf = pq.ParquetFile(src)
    n_groups = pf.num_row_groups
    parts: list[pd.DataFrame] = []
    total = 0
    for i in range(n_groups):
        df = pf.read_row_group(i, columns=EDGE_COLS).to_pandas()
        total += len(df)
        both = df["gene_a"].isin(node_set) & df["gene_b"].isin(node_set)
        df = df[both]
        if df.empty:
            continue
        dmax = df["delta_max"].abs()
        keep = (
            df["in_string_ge700"].fillna(False).astype(bool)
            | df["coloc_linked"].fillna(False).astype(bool)
            | df["is_contested"].fillna(False).astype(bool)
            | (dmax >= TAU_FLOOR)
        )
        df = df[keep]
        if not df.empty:
            parts.append(df)
        print(f"    row-group {i+1}/{n_groups}: read {len(df):,} kept "
              f"(cumulative source {total:,})", flush=True)
    cand = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=EDGE_COLS)
    print(f"  candidate edges (both-in-set & any-signal, |dmax|>={TAU_FLOOR}) = {len(cand):,}")
    return cand


def edge_strength(df: pd.DataFrame) -> pd.Series:
    dmax = df["delta_max"].abs().fillna(0.0)
    sstr = (pd.to_numeric(df["string_score"], errors="coerce") / 1000.0).fillna(0.0)
    pp4 = pd.to_numeric(df["pp4_min"], errors="coerce").fillna(0.0)
    coloc = df["coloc_linked"].fillna(False).astype(bool).astype(float)
    cont = df["is_contested"].fillna(False).astype(bool).astype(float)
    return np.maximum.reduce([dmax.values, sstr.values, pp4.values,
                              coloc.values, cont.values])


def apply_topk(df: pd.DataFrame, k: int) -> pd.DataFrame:
    """Keep an edge if it is in gene_a's OR gene_b's top-k by strength."""
    df = df.reset_index(drop=True)
    df["_eid"] = np.arange(len(df))
    a = df[["_eid", "gene_a", "_strength"]].rename(columns={"gene_a": "node"})
    b = df[["_eid", "gene_b", "_strength"]].rename(columns={"gene_b": "node"})
    long = pd.concat([a, b], ignore_index=True)
    long["_rk"] = long.groupby("node")["_strength"].rank(method="first", ascending=False)
    keep_eids = pd.unique(long.loc[long["_rk"] <= k, "_eid"])
    out = df[df["_eid"].isin(keep_eids)].drop(columns=["_eid", "_strength"])
    return out


def finalize_edges(cand: pd.DataFrame, tau: float, k: int) -> pd.DataFrame:
    dmax = cand["delta_max"].abs()
    keep = (
        cand["in_string_ge700"].fillna(False).astype(bool)
        | cand["coloc_linked"].fillna(False).astype(bool)
        | cand["is_contested"].fillna(False).astype(bool)
        | (dmax >= tau)
    )
    sub = cand[keep].copy()
    sub["_strength"] = edge_strength(sub)
    sub = apply_topk(sub, k)
    return sub.sort_values(["gene_a", "gene_b"]).reset_index(drop=True)


def write_edges(df: pd.DataFrame, path: Path) -> float:
    out = df.copy()
    out["type"] = out["type"].map(map_channel)
    # community IDs -> nullable Int64 (they are integer community labels)
    for c in ("community_macro_a", "community_macro_b"):
        v = pd.to_numeric(out[c], errors="coerce")
        out[c] = v.round().astype("Int64")
    # +/-Inf -> NaN for all float columns (contract Section 0 item 4)
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].replace([np.inf, -np.inf], np.nan)
    out = out[EDGE_COLS]
    out.to_parquet(path, engine="pyarrow", index=False)
    return path.stat().st_size / 1e6


def build_nodes(nodes_raw: pd.DataFrame, atlas: pd.DataFrame) -> pd.DataFrame:
    df = nodes_raw.rename(columns={
        "human_symbol": "symbol",
        "dream_logFC": "bulk_logFC",
        "dream_padj": "bulk_padj",
        "dream_tstat": "bulk_tstat",
    }).copy()

    # Refresh from the canonical atlas (resolves the value-staleness trap, not
    # just the column names). A gene absent from the atlas gets NaN DE rather
    # than a retired dream number relabeled as bulk_*.
    a = atlas.drop_duplicates("human_symbol").set_index("human_symbol")
    n_absent = int((~df["symbol"].isin(a.index)).sum())
    for c in ("bulk_logFC", "bulk_padj", "bulk_tstat"):
        df[c] = pd.to_numeric(a[c], errors="coerce").reindex(df["symbol"]).to_numpy()

    # canonical is_deg = effect-size-aware interval-null FDR gate (< 0.05); NaN -> False.
    treat = pd.to_numeric(a["bulk_treat_fdr"], errors="coerce")
    df["is_deg"] = (treat < 0.05).reindex(df["symbol"]).fillna(False).to_numpy().astype(bool)
    print(f"  [nodes] {n_absent:,} network genes absent from atlas -> NaN bulk_* (no stale dream fallback)")

    df["is_conserved_core"] = as_bool(df["is_conserved_core"]).astype(bool)
    df["dgidb_druggable"] = as_bool(df["dgidb_druggable"]).astype(bool)
    df["coloc_susie_best_pp4"] = pd.to_numeric(df["coloc_susie_best_pp4"], errors="coerce")
    df["n_coloc_sources"] = pd.to_numeric(df["n_coloc_sources"], errors="coerce").round().astype("Int64")
    for c in ("sex_class", "zonation_class", "ferroptosis_class",
              "attribution_class", "mouse_ortholog", "ensembl_id", "gene_biotype"):
        df[c] = df[c].where(df[c].notna(), None)
        df[c] = df[c].replace({"": None})

    for c in ("bulk_logFC", "bulk_padj", "bulk_tstat", "coloc_susie_best_pp4"):
        df[c] = df[c].replace([np.inf, -np.inf], np.nan)

    df = df[NODE_COLS].sort_values("symbol").reset_index(drop=True)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=None,
                    help="default: <root>/masld-atlas-v2/public/data")
    ap.add_argument("--tau", type=float, default=None,
                    help="override auto-tuned |delta_max| threshold")
    ap.add_argument("--k", type=int, default=K_PER_NODE)
    args = ap.parse_args()

    root = default_root()
    net = root / "RNA-seq/results/network"
    src = net / "edge_annotation_atlas.parquet"
    nodes_csv = net / "network_nodes.csv"
    atlas_csv = root / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
    out_dir = Path(args.output_dir) if args.output_dir else root / "masld-atlas-v2/public/data"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[network_to_parquet] root={root}")
    print(f"  edge src   = {src}")
    print(f"  nodes csv  = {nodes_csv}")
    print(f"  atlas csv  = {atlas_csv}")
    print(f"  out dir    = {out_dir}")

    nodes_raw = pd.read_csv(nodes_csv)
    atlas = pd.read_csv(
        atlas_csv,
        usecols=["human_symbol", "bulk_logFC", "bulk_padj", "bulk_tstat",
                 "bulk_treat_fdr", "layers_active"],
        low_memory=False,
    )

    # ---------------- nodes ----------------
    print("\n[nodes] building network_nodes.parquet")
    node_df = build_nodes(nodes_raw, atlas)
    node_path = out_dir / "network_nodes.parquet"
    node_df.to_parquet(node_path, engine="pyarrow", index=False)
    assert not any(c.startswith("dream_") for c in node_df.columns), "dream_* leaked into node parquet"
    print(f"  wrote {node_path.name}: {len(node_df):,} rows x {len(node_df.columns)} cols "
          f"({node_path.stat().st_size/1e6:.2f} MB); is_deg TRUE={int(node_df['is_deg'].sum())}")

    # ---------------- edges ----------------
    print("\n[edges] building NODE_SET + streaming candidate edges")
    node_set = build_node_set(nodes_raw, atlas)
    cand = stream_candidate_edges(src, node_set)

    print("\n[edges] tuning TAU to fit <"
          f" {SIZE_TARGET_MB} MB (K={args.k})")
    edge_path = out_dir / "network_edges.parquet"
    chosen = None
    tau_list = [args.tau] if args.tau is not None else TAU_GRID
    for tau in tau_list:
        fin = finalize_edges(cand, tau, args.k)
        size_mb = write_edges(fin, edge_path)
        print(f"    TAU={tau:>4}: edges={len(fin):>8,}  ->  {size_mb:6.2f} MB")
        if args.tau is not None or size_mb <= SIZE_TARGET_MB:
            chosen = (tau, len(fin), size_mb)
            break
        chosen = (tau, len(fin), size_mb)  # keep last as fallback
    tau, n_edges, size_mb = chosen

    # channel distribution + raw->mapped verification
    fin = finalize_edges(cand, tau, args.k)
    raw_counts = fin["type"].value_counts(dropna=False)
    mapping = {t: map_channel(t) for t in raw_counts.index}
    write_edges(fin, edge_path)

    print("\n==================== RESULT ====================")
    print(f"FILTER: (in_string_ge700 OR coloc_linked OR is_contested OR "
          f"|delta_max|>={tau}) AND both endpoints in NODE_SET; then per-node top-{args.k}")
    print(f"NODE_SET = union(is_deg, coloc_susie_best_pp4>={COLOC_PP4_MIN}, "
          f"dgidb_druggable, top{TOP_LAYERS} by layers_active) = {len(node_set):,} genes")
    print(f"TAU = {tau}")
    print(f"network_edges.parquet : {n_edges:,} rows x {len(EDGE_COLS)} cols  ->  {size_mb:.2f} MB")
    print(f"network_nodes.parquet : {len(node_df):,} rows x {len(NODE_COLS)} cols")
    print("\nedge channel (type) raw -> mapped counts:")
    for raw, cnt in raw_counts.items():
        print(f"   {str(raw):<16} -> {str(mapping[raw]):<10} {cnt:,}")
    print("===============================================")


if __name__ == "__main__":
    main()
