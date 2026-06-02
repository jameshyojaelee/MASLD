"""Network Explorer.

Force-directed multi-evidence gene network. All controls live on the page
itself (not the sidebar) because the graph is the primary workspace.
Deep-linkable as ``/Network?genes=THRB,RORA``.

Controls:
    * Gene multiselect (max 100) -- URL-synced via ``?genes=SYM1,SYM2``
    * Layer toggles for the six v4 edge types (S, D-F2, D-COLOC, D-LR,
      D-ceRNA, D-XS).
    * Three thresholds (STRING, LOCO replication, PP.H4) that dim rather
      than remove below-threshold edges -- preserves context.
    * F-stage filter (``all | F0-F1 | F2 | F3-F4 | F2_emerging |
      F2_dissolving``) applied to D-F2 edges via ``emergence_stage``.
    * Label top-K slider (0-100).
"""
from __future__ import annotations

import math
import time
from typing import Any, Iterable

import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from components.layout import inject_css, page_header, sidebar_header
from data.loaders import (
    load_gene_graph,
    load_gene_index,
    load_layer_metadata,
)
from theme.plotly_theme import LAYER_COLOR_MAP, PLOTLY_THEME, apply_theme

inject_css()
sidebar_header()

MAX_GENES = 100

V4_LAYERS: list[str] = [
    "ppi", "coexpr", "coexpr_stable", "regulon",
    "genetic", "lr", "cerna", "xspecies",
]

LAYER_TO_EDGE_TYPE: dict[str, str] = {
    "ppi":           "S",
    "coexpr":        "D-F2",
    "coexpr_stable": "D-STABLE",
    "regulon":       "D-REG",
    "genetic":       "D-COLOC",
    "lr":            "D-LR",
    "cerna":         "D-ceRNA",
    "xspecies":      "D-XS",
}
EDGE_TYPE_TO_LAYER: dict[str, str] = {v: k for k, v in LAYER_TO_EDGE_TYPE.items()}

LAYER_COLOR_FALLBACK = dict(LAYER_COLOR_MAP)

# User-facing layer names — keep internal keys stable so downstream
# filtering code doesn't change, but relabel everything visible.
LAYER_LABEL: dict[str, str] = {
    "ppi":           "Protein interaction (STRING)",
    "coexpr":        "Co-expression — F2-dynamic",
    "coexpr_stable": "Co-expression — stable across stages",
    "regulon":       "TF regulon (SCENIC+)",
    "genetic":       "Genetic (GWAS colocalization)",
    "lr":            "Ligand–receptor signalling",
    "cerna":         "miRNA-mediated (ceRNA)",
    "xspecies":      "Cross-species conserved",
}

LAYER_DESC_FALLBACK: dict[str, str] = {
    "ppi":           "Protein–protein interactions from STRING v12, high confidence (score ≥ 0.7).",
    "coexpr":        "Stage-dynamic co-expression: edges that change with fibrosis stage. Threshold with the |Δr| slider below.",
    "coexpr_stable": "Stage-invariant co-expression: pairs with strong correlation (|r_mean| ≥ 0.6) that barely move across F01 / F2 / F34. Baseline against which the F2-dynamic layer measures change.",
    "regulon":       "Directed TF → target edges from SCENIC+ disease and hepatocyte regulons.",
    "genetic":       "GWAS–eQTL colocalizations — gene pairs with shared causal variants (PP.H4 ≥ 0.3 shared locus or eQTL sentinel sharing).",
    "lr":            "Ligand–receptor signalling pairs (LIANA consensus across cell-type pairs).",
    "cerna":         "Competing endogenous RNA — pairs sharing microRNA regulators.",
    "xspecies":      "Edges conserved between human and mouse liver (WGCNA preservation + cell-type concordance).",
}

# Fibrosis-stage filter options. Keys here are what the user sees; values map
# to the `emergence_stage` labels on D-F2 edges. Order matters (rendered as a
# horizontal radio).
FSTAGE_MAP: dict[str, set[str] | None] = {
    "All stages": None,
    "F0–F1 (no/mild fibrosis)": {"F0_specific", "F2_dissolving"},
    "F2 (transition point)": {"F2_emerging", "F2_dissolving", "transient_F2"},
    "F3–F4 (advanced fibrosis)": {"F34_specific", "F2_emerging", "progressive_up", "progressive_down"},
    "F2 onset only": {"F2_emerging"},
    "F2 resolution only": {"F2_dissolving"},
}

NEIGHBOR_BUCKETS: list[str] = [
    "neighbors_string",
    "neighbors_d_f2_emerging",
    "neighbors_d_f2_dissolving",
    "neighbors_d_f2_transient",
    "neighbors_d_f2_progressive_up",
    "neighbors_d_f2_progressive_down",
    "neighbors_d_f2_f01_specific",
    "neighbors_d_f2_f34_specific",
    "neighbors_d_stable",
    "neighbors_d_regulon_out",
    "neighbors_d_regulon_in",
    "neighbors_d_regulon_other",
    "neighbors_d_coloc",
    "neighbors_d_lr",
    "top_d_cerna",
    "top_d_xs",
]

MAX_EDGES_PER_TYPE = 1000


def _parse_genes_param(raw: Any) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, list):
        raw = raw[0] if raw else ""
    if not isinstance(raw, str) or not raw.strip():
        return []
    seen: list[str] = []
    for tok in raw.split(","):
        s = tok.strip().upper()
        if s and s not in seen:
            seen.append(s)
        if len(seen) >= MAX_GENES:
            break
    return seen


_url_genes = _parse_genes_param(st.query_params.get("genes"))


# ---------------------------------------------------------------------------
# Header + on-page controls
# ---------------------------------------------------------------------------

page_header(
    "Network Explorer",
    "Force-directed multi-evidence gene network. Pick one gene to load its "
    "neighborhood, or add up to 100 to see each micronetwork with shared bridges highlighted.",
)

# --- Layer metadata (colours / counts / descriptions) needed by both the
# --- control row below and the legend further down the page.
layer_meta_raw = load_layer_metadata()
edge_type_counts = {row.get("type"): row.get("n_edges", 0) for row in layer_meta_raw}
layer_color: dict[str, str] = dict(LAYER_COLOR_FALLBACK)
layer_desc: dict[str, str] = dict(LAYER_DESC_FALLBACK)
for row in layer_meta_raw:
    et = row.get("type")
    if et in EDGE_TYPE_TO_LAYER:
        lk = EDGE_TYPE_TO_LAYER[et]
        if "color" in row:
            layer_color[lk] = row["color"]
        if "description" in row:
            layer_desc[lk] = row["description"]

try:
    gene_symbols = sorted(load_gene_index()["symbol"].dropna().unique().tolist())
except Exception:
    gene_symbols = []

# Pre-select RORA on a cold landing (no ?genes= URL param) so the user sees
# a populated graph immediately instead of an empty-state info banner.
DEFAULT_GENE = "RORA"
if _url_genes:
    default_genes = [g for g in _url_genes if not gene_symbols or g in gene_symbols][:MAX_GENES]
elif not gene_symbols or DEFAULT_GENE in gene_symbols:
    default_genes = [DEFAULT_GENE]
else:
    default_genes = []

# Primary action: gene picker (full width, prominent).
selected_genes: list[str] = st.multiselect(
    "Genes (max 100)",
    options=gene_symbols or default_genes,
    default=default_genes,
    max_selections=MAX_GENES,
    help="Pick up to 100 genes. Multi-gene mode shows each gene's full neighborhood with shared bridges highlighted. Layouts for 30+ genes may take a few seconds.",
)

# URL-sync the selection so sharing / refreshing preserves state.
if selected_genes != _url_genes:
    if selected_genes:
        st.query_params["genes"] = ",".join(selected_genes)
    else:
        try:
            del st.query_params["genes"]
        except KeyError:
            pass

# --- Filters: layer multiselect, F-stage radio, label slider -----------------
layer_labels = {
    lk: f"{LAYER_LABEL.get(lk, lk)}  ({edge_type_counts.get(LAYER_TO_EDGE_TYPE[lk], 0):,})"
    for lk in V4_LAYERS
}
ctrl_a, ctrl_b, ctrl_c = st.columns([3, 2, 1])
with ctrl_a:
    selected_layers = st.multiselect(
        "Edge layers",
        options=V4_LAYERS,
        default=V4_LAYERS,
        format_func=lambda k: layer_labels.get(k, k),
    )
with ctrl_b:
    fstage = st.radio(
        "Fibrosis stage",
        options=list(FSTAGE_MAP.keys()),
        index=0,
        horizontal=True,
        help=(
            "Restricts the co-expression (stage-dynamic) edges to a specific "
            "fibrosis-stage window. F0 = no fibrosis, F4 = cirrhosis; F2 is the "
            "clinical inflection point. Other edge types pass through unfiltered."
        ),
    )
with ctrl_c:
    label_top_k = st.slider("Label top-K nodes", 0, 100, 20, 5)

# --- Thresholds: 3 sliders side by side --------------------------------------
th_a, th_b, th_c = st.columns(3)
with th_a:
    th_string = st.slider(
        "Protein interaction score ≥", 0.40, 1.00, 0.70, 0.01,
        help=(
            "STRING v12 combined score. 0.7 is the standard 'high-confidence' "
            "cutoff. Below-threshold edges are dimmed, not removed."
        ),
    )
with th_b:
    th_loco = st.slider(
        "Co-expression replication ≥", 0.00, 1.00, 0.70, 0.01,
        help=(
            "Leave-one-cohort-out replication fraction: proportion of the 10 "
            "cohorts in which the co-expression edge is re-detected."
        ),
    )
with th_c:
    th_pp4 = st.slider(
        "GWAS colocalization ≥", 0.00, 1.00, 0.50, 0.01,
        help=(
            "PP.H4 = posterior probability that the gene and the GWAS trait "
            "share a causal variant. ≥0.5 = suggestive, ≥0.8 = strong."
        ),
    )

# Phase 2: continuous |Δr| slider for D-F2 edges + deconv-adjusted toggle.
th_d, th_e, _pad = st.columns([2, 2, 2])
with th_d:
    th_delta = st.slider(
        "F2 co-expression |Δr| ≥", 0.00, 1.00, 0.10, 0.01,
        help=(
            "Continuous threshold on |Δr| — the largest stage-to-stage "
            "change in Pearson correlation across F01 / F2 / F34. Replaces "
            "the old 0.3 / 0.5 hard cutpoints. Higher = only the largest "
            "shifts at F2."
        ),
    )
with th_e:
    deconv_only = st.checkbox(
        "Deconv-adjusted co-expression only",
        value=False,
        help=(
            "Hide D-F2 edges flagged as composition-driven: pairs whose "
            "|Δr| drops >50% once sample cell-type proportions are "
            "regressed out. When on, you see regulatory rewiring; when "
            "off, you see raw bulk correlation (which includes macrophage "
            "infiltration artefact)."
        ),
    )

st.markdown('<div class="masld-spacer-sm"></div>', unsafe_allow_html=True)

if not selected_genes:
    st.info("Search for a gene above to explore its network.")
    st.stop()


# ---------------------------------------------------------------------------
# Graph assembly + render
# ---------------------------------------------------------------------------

def _bucket_layer(bucket: str, type_hint: str | None = None) -> str:
    # Bucket names are layer-authoritative. `type_hint` used to be a
    # fallback but for multi-source edges (e.g. type="S,D-STABLE") it
    # mis-routed every row of the bucket to the first layer in the tag
    # list, leaving D-STABLE / D-REG neighbors stamped as PPI.
    if bucket == "neighbors_string":
        return "ppi"
    if bucket.startswith("neighbors_d_f2_"):
        return "coexpr"
    if bucket == "neighbors_d_stable":
        return "coexpr_stable"
    if bucket.startswith("neighbors_d_regulon"):
        return "regulon"
    if bucket == "neighbors_d_coloc":
        return "genetic"
    if bucket == "neighbors_d_lr":
        return "lr"
    if bucket == "top_d_cerna":
        return "cerna"
    if bucket == "top_d_xs":
        return "xspecies"
    # Unknown bucket -- fall back on any recognizable tag in the type hint.
    if type_hint:
        for tag in str(type_hint).split(","):
            if tag in EDGE_TYPE_TO_LAYER:
                return EDGE_TYPE_TO_LAYER[tag]
    return "ppi"


def _edge_score(layer: str, n: dict) -> float:
    if layer == "ppi":
        return float(n.get("string_score", 0.0) or 0.0)
    if layer == "coexpr":
        return float(n.get("delta_max", 0.0) or 0.0)
    if layer == "coexpr_stable":
        return float(n.get("stable_r_mean", 0.0) or 0.0)
    if layer == "regulon":
        padj = n.get("regulon_activity_padj")
        return float(1.0 - padj) if isinstance(padj, (int, float)) else 0.0
    if layer == "genetic":
        return float(n.get("pp4_min", 0.0) or 0.0)
    if layer == "lr":
        return float(n.get("score_diff_lr", 0.0) or 0.0)
    if layer == "xspecies":
        return float(n.get("xs_preservation_Z", 0.0) or 0.0)
    return float(
        n.get("posterior", n.get("string_score", n.get("pp4_min", 0.0))) or 0.0
    )


def _edge_visible(
    layer: str, n: dict, fstage_set: set[str] | None,
    th_string_v: float, th_loco_v: float, th_pp4_v: float,
    th_delta_v: float = 0.0,
    deconv_only: bool = False,
) -> tuple[bool, bool]:
    stage = n.get("emergence_stage")
    if layer == "coexpr" and fstage_set is not None:
        if not stage or stage not in fstage_set:
            return False, False
    if layer == "coexpr" and deconv_only:
        # "Deconv-adjusted only" checkbox: hide edges flagged as
        # composition-driven by 281d.
        if bool(n.get("f2_composition_driven")):
            return False, False
    dim = False
    if layer == "ppi":
        if float(n.get("string_score", 0.0) or 0.0) < th_string_v:
            dim = True
    elif layer == "coexpr":
        dm = float(n.get("delta_max", 0.0) or 0.0)
        if dm < th_delta_v:
            return False, False  # continuous slider is a HARD cut, not a dim
        rep = float(n.get("loco_replication_fraction", 0.0) or 0.0)
        # rep may not exist on the new continuous D-F2 table; treat NaN as
        # not-dimmed so the slider alone drives visibility.
        if rep > 0 and rep < th_loco_v:
            dim = True
    elif layer == "genetic":
        if float(n.get("pp4_min", 0.0) or 0.0) < th_pp4_v:
            dim = True
    return True, dim


def _collect_edges(
    graph_json: dict, center: str, selected: Iterable[str],
    fstage_set: set[str] | None,
    th_string_v: float, th_loco_v: float, th_pp4_v: float,
    th_delta_v: float = 0.0,
    deconv_only_v: bool = False,
) -> list[dict]:
    selected_set = set(selected)
    neighbors = graph_json.get("neighbors", {})
    out: list[dict] = []
    for bucket in NEIGHBOR_BUCKETS:
        rows = neighbors.get(bucket, []) or []
        layer = _bucket_layer(bucket, rows[0].get("type") if rows else None)
        if layer not in selected_set:
            continue
        rows_sorted = sorted(rows, key=lambda r: _edge_score(layer, r), reverse=True)[
            :MAX_EDGES_PER_TYPE
        ]
        for n in rows_sorted:
            keep, dim = _edge_visible(
                layer, n, fstage_set, th_string_v, th_loco_v, th_pp4_v,
                th_delta_v, deconv_only_v,
            )
            if not keep:
                continue
            partner = n.get("partner")
            if not partner:
                continue
            out.append({
                "source": center, "target": partner, "layer": layer,
                "score": _edge_score(layer, n), "dim": dim,
                "emergence_stage": n.get("emergence_stage"),
                "loco_replication_fraction": n.get("loco_replication_fraction"),
                "pp4_min": n.get("pp4_min"), "string_score": n.get("string_score"),
                "druggable_pair": n.get("druggable_pair"),
                "conserved_mouse_a": n.get("conserved_mouse_a"),
                "conserved_mouse_b": n.get("conserved_mouse_b"),
                "bucket": bucket,
            })
    return out


@st.cache_data(ttl=600, max_entries=32, show_spinner="Building graph...")
def build_graph(
    genes_key: tuple[str, ...], selected_layers_key: tuple[str, ...],
    fstage_key: str, th_string_key: float, th_loco_key: float, th_pp4_key: float,
    th_delta_key: float = 0.0, deconv_only_key: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    fstage_set = FSTAGE_MAP.get(fstage_key)
    selected_set = set(selected_layers_key)

    graphs: dict[str, dict] = {}
    missing: list[str] = []
    for sym in genes_key:
        g = load_gene_graph(sym)
        if g is None:
            missing.append(sym)
        else:
            graphs[sym] = g

    all_edges: list[dict] = []
    for sym, g in graphs.items():
        all_edges.extend(_collect_edges(
            g, sym, selected_set, fstage_set,
            th_string_key, th_loco_key, th_pp4_key,
            th_delta_key, deconv_only_key,
        ))

    seen: set[tuple] = set()
    dedup_edges: list[dict] = []
    for e in all_edges:
        k = (frozenset({e["source"], e["target"]}), e["layer"])
        if k in seen:
            continue
        seen.add(k)
        dedup_edges.append(e)

    queried_set = set(genes_key)
    neighbor_owners: dict[str, set[str]] = {}
    for e in all_edges:
        src, tgt = e["source"], e["target"]
        if src in queried_set and tgt not in queried_set:
            neighbor_owners.setdefault(tgt, set()).add(src)
        elif tgt in queried_set and src not in queried_set:
            neighbor_owners.setdefault(src, set()).add(tgt)

    # Multi-gene mode: render the UNION of each gene's micronetwork and
    # highlight shared neighbors (sharedCount >= 2) as bridges. Pairwise-
    # specific neighbors stay visible (faded) so each gene's own neighborhood
    # is preserved — this is the "union view" the user expects.
    multi_mode = len(queried_set) > 1

    node_ids: set[str] = set(queried_set)
    for e in dedup_edges:
        node_ids.add(e["source"])
        node_ids.add(e["target"])

    edges_df = pd.DataFrame(dedup_edges) if dedup_edges else pd.DataFrame(
        columns=["source", "target", "layer", "score", "dim"]
    )

    degree: dict[str, int] = {nid: 0 for nid in node_ids}
    for e in dedup_edges:
        degree[e["source"]] = degree.get(e["source"], 0) + 1
        degree[e["target"]] = degree.get(e["target"], 0) + 1

    primary_layer: dict[str, str] = {}
    max_score: dict[str, float] = {}
    for e in dedup_edges:
        for endpoint in (e["source"], e["target"]):
            if endpoint in queried_set:
                continue
            if e["score"] > max_score.get(endpoint, -1.0):
                max_score[endpoint] = e["score"]
                primary_layer[endpoint] = e["layer"]

    shared_count: dict[str, int] = {
        nid: len(neighbor_owners.get(nid, set())) for nid in node_ids
    }

    G = nx.Graph()
    G.add_nodes_from(node_ids)
    for e in dedup_edges:
        weight = max(0.1, e["score"]) if not e["dim"] else 0.05
        G.add_edge(e["source"], e["target"], weight=weight)

    n = max(1, G.number_of_nodes())
    k_val = 1.5 / math.sqrt(n)
    try:
        pos = nx.spring_layout(G, seed=42, k=k_val, iterations=30 if n > 400 else 50)
    except Exception:
        pos = nx.random_layout(G, seed=42)

    nodes_df = pd.DataFrame({"id": list(node_ids)})
    nodes_df["x"] = nodes_df["id"].map(lambda i: pos.get(i, (0.0, 0.0))[0])
    nodes_df["y"] = nodes_df["id"].map(lambda i: pos.get(i, (0.0, 0.0))[1])
    nodes_df["degree"] = nodes_df["id"].map(lambda i: degree.get(i, 0))
    nodes_df["shared_count"] = nodes_df["id"].map(lambda i: shared_count.get(i, 0))
    nodes_df["is_queried"] = nodes_df["id"].isin(queried_set)
    nodes_df["primary_layer"] = nodes_df.apply(
        lambda r: "queried" if r["is_queried"] else primary_layer.get(r["id"], "ppi"),
        axis=1,
    )
    return nodes_df, edges_df, {"missing": missing}


t0 = time.perf_counter()
nodes_df, edges_df, meta = build_graph(
    tuple(selected_genes), tuple(sorted(selected_layers)),
    fstage, round(th_string, 3), round(th_loco, 3), round(th_pp4, 3),
    round(th_delta, 3), bool(deconv_only),
)
missing = meta.get("missing", []) if isinstance(meta, dict) else []
build_ms = (time.perf_counter() - t0) * 1000.0

if missing:
    st.warning(
        f"No network graph found for: {', '.join(missing)}. "
        "These may not have exported neighborhoods yet."
    )
if edges_df.empty:
    st.warning(
        "No edges survive the current thresholds / layer selection. "
        "Try lowering a threshold or re-enabling layers."
    )


# ---- Plotly -----------------------------------------------------------------

fig = go.Figure()
if not edges_df.empty:
    for layer in selected_layers:
        sub = edges_df[edges_df["layer"] == layer]
        if sub.empty:
            continue
        for dim_flag in (True, False):
            rows = sub[sub["dim"] == dim_flag]
            if rows.empty:
                continue
            xs: list[float | None] = []
            ys: list[float | None] = []
            x_lookup = dict(zip(nodes_df["id"], nodes_df["x"]))
            y_lookup = dict(zip(nodes_df["id"], nodes_df["y"]))
            for _, e in rows.iterrows():
                if e["source"] not in x_lookup or e["target"] not in x_lookup:
                    continue
                xs.extend([x_lookup[e["source"]], x_lookup[e["target"]], None])
                ys.extend([y_lookup[e["source"]], y_lookup[e["target"]], None])
            if not xs:
                continue
            alpha = 0.12 if dim_flag else 0.45
            fig.add_trace(go.Scattergl(
                x=xs, y=ys, mode="lines",
                line=dict(color=layer_color.get(layer, "#888"), width=1.2),
                opacity=alpha, hoverinfo="skip",
                name=f"{layer}{' (dim)' if dim_flag else ''}",
                legendgroup=layer, showlegend=(not dim_flag),
            ))


def _node_color(r: pd.Series) -> str:
    if r["is_queried"]:
        return "#c0392b"
    return layer_color.get(r["primary_layer"], "#95a5a6")


def _node_size(r: pd.Series) -> float:
    # Lighter marker sizing to match the Next.js force-graph aesthetic
    # (Next.js uses ~2 + log2(degree)*0.9; we scale up modestly for Plotly
    # which interprets marker.size differently). Queried genes always pop.
    base = 4.0 + math.log2(max(1, int(r["degree"]))) * 1.1
    return base * (1.6 if r.get("is_queried") else 1.0)


if not nodes_df.empty:
    nodes_df = nodes_df.copy()
    nodes_df["color"] = nodes_df.apply(_node_color, axis=1)
    nodes_df["size"] = nodes_df.apply(_node_size, axis=1)
    # In multi-gene mode, shared bridges get a larger size bump + white ring
    # so they stand out against the pairwise-specific neighbors.
    multi_mode_view = len(selected_genes) > 1
    if multi_mode_view:
        bridge_mask = nodes_df["shared_count"] >= 2
        nodes_df.loc[bridge_mask, "size"] = nodes_df.loc[bridge_mask, "size"] * 1.5
    top_k = (
        nodes_df[~nodes_df["is_queried"]]
        .sort_values(["shared_count", "degree"], ascending=[False, False])
        .head(label_top_k)["id"].tolist()
    )
    # Always label queried genes AND shared bridges (sharedCount>=2).
    bridge_ids = set(nodes_df.loc[nodes_df["shared_count"] >= 2, "id"])
    labelled_ids = set(top_k) | set(nodes_df.loc[nodes_df["is_queried"], "id"]) | bridge_ids
    nodes_df["label"] = nodes_df["id"].where(nodes_df["id"].isin(labelled_ids), "")
    # Outline: 3px white ring on queried; 2px on shared bridges; 0.6px on others.
    outline_width = nodes_df.apply(
        lambda r: 3.0 if r["is_queried"] else (2.0 if r["shared_count"] >= 2 else 0.6),
        axis=1,
    )
    # Opacity fades pairwise-specific neighbors when multi-gene mode is active;
    # queried + bridges stay at full opacity so they pop out.
    if multi_mode_view:
        opacity_arr = nodes_df.apply(
            lambda r: 1.0 if (r["is_queried"] or r["shared_count"] >= 2) else 0.35,
            axis=1,
        ).values
    else:
        opacity_arr = 0.92

    fig.add_trace(go.Scattergl(
        x=nodes_df["x"], y=nodes_df["y"], mode="markers+text",
        marker=dict(
            size=nodes_df["size"], color=nodes_df["color"],
            line=dict(color="white", width=outline_width), opacity=opacity_arr,
        ),
        text=nodes_df["label"], textposition="top center",
        textfont=dict(size=10, color="#222"),
        customdata=nodes_df[["id", "degree", "shared_count", "primary_layer"]].values,
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>Degree: %{customdata[1]}<br>"
            "Shared across queries: %{customdata[2]}<br>"
            "Primary layer: %{customdata[3]}<extra></extra>"
        ),
        showlegend=False, name="nodes",
    ))

apply_theme(fig, margin=dict(l=10, r=10, t=10, b=10))
fig.update_layout(
    height=640,
    xaxis=dict(visible=False),
    yaxis=dict(visible=False, scaleanchor="x", scaleratio=1),
    legend=dict(
        orientation="h", yanchor="bottom", y=-0.06,
        xanchor="right", x=1,
        font=dict(size=11, color="#737782"),
        itemsizing="constant",
    ),
    hovermode="closest",
)

if len(selected_genes) > 1:
    shared_nodes = int((nodes_df["shared_count"] >= 2).sum())
    direct_edges = 0
    if not edges_df.empty:
        q = set(selected_genes)
        direct_edges = int(
            edges_df.apply(lambda r: (r["source"] in q) and (r["target"] in q), axis=1).sum()
        )
    st.markdown(
        f"**{len(selected_genes)}** queried  "
        f"&middot;  **{shared_nodes}** shared neighbors  "
        f"&middot;  **{direct_edges}** direct edges between queries"
    )

# Per-queried-gene breakdown — explains what is/isn't rendered for each gene.
# Helps users interpret cases where a gene looks "missing" because its edges
# are below a threshold or excluded by the F-stage filter.
if selected_genes and not edges_df.empty:
    per_gene_rows = []
    for sym in selected_genes:
        mask = (edges_df["source"] == sym) | (edges_df["target"] == sym)
        g_edges = edges_df[mask]
        visible = int((~g_edges["dim"]).sum()) if "dim" in g_edges.columns else len(g_edges)
        dimmed = int(g_edges["dim"].sum()) if "dim" in g_edges.columns else 0
        by_layer = g_edges["layer"].value_counts().to_dict() if len(g_edges) else {}
        per_gene_rows.append({
            "gene": sym,
            "visible edges": visible,
            "dimmed (below threshold)": dimmed,
            **{f"{lyr}": int(n) for lyr, n in by_layer.items()},
        })
    with st.expander("Per-gene edge breakdown", expanded=False):
        st.dataframe(
            pd.DataFrame(per_gene_rows).fillna(0),
            use_container_width=True, hide_index=True,
        )
        st.caption(
            "If a queried gene has 0 visible edges, it may be (a) below one of "
            "the three threshold sliders above, (b) excluded by the fibrosis-"
            "stage filter, or (c) only has legacy-layer edges (regulon / "
            "pathway / spatial — not shown in v4)."
        )

st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False})
st.caption(f"Rendered {len(nodes_df)} nodes / {len(edges_df)} edges in {build_ms:.0f} ms.")


# ---- Layer legend -----------------------------------------------------------

legend_cols = st.columns(3)
for i, layer in enumerate(V4_LAYERS):
    with legend_cols[i % 3]:
        active = layer in selected_layers
        sub = edges_df[edges_df["layer"] == layer] if not edges_df.empty else edges_df
        n_visible = len(sub)
        n_dim = int(sub["dim"].sum()) if not sub.empty else 0
        swatch = (
            f"<span class='masld-layer-swatch' "
            f"style='background:{layer_color.get(layer, '#888')};'></span>"
        )
        state = "on" if active else "off"
        st.markdown(
            f"{swatch}**{LAYER_LABEL.get(layer, layer)}** _{state}_ &middot; "
            f"{n_visible} edges ({n_dim} dim)",
            unsafe_allow_html=True,
        )
        st.caption(layer_desc.get(layer, ""))


# ---- Queried gene summaries -------------------------------------------------

st.markdown("### Queried gene summaries")
for sym in selected_genes:
    g = load_gene_graph(sym)
    if g is None:
        continue
    attrs = g.get("attributes", {})
    ec = g.get("edge_counts", {})
    trajectory = g.get("f_stage_trajectory", {})
    with st.expander(f"{sym} -- {g.get('biotype', '')}", expanded=False):
        c1, c2, c3, c4 = st.columns(4)
        lfc = attrs.get("dream_logFC")
        padj = attrs.get("dream_padj")
        c1.metric("Integrated logFC", f"{lfc:+.3f}" if isinstance(lfc, (int, float)) else "--")
        c2.metric("padj", f"{padj:.2e}" if isinstance(padj, (int, float)) and padj else "--")
        c3.metric("COLOC PP4 max", f"{attrs.get('coloc_susie_best_pp4', 0):.3f}")
        c4.metric("Edges (total)", sum(ec.values()) if ec else 0)
        st.markdown(
            f"Emergence stage: **{trajectory.get('emergence_stage', '--')}**  "
            f"&middot; sex_class: **{attrs.get('sex_class', '--')}**  "
            f"&middot; attribution: **{attrs.get('attribution_class', '--')}**"
        )
        st.markdown(f"[Open gene detail for {sym} \u2192](/Gene?symbol={sym})")


# ---- Edge provenance table --------------------------------------------------

st.markdown("### Visible edges (provenance)")
if edges_df.empty:
    st.info("No edges to list.")
else:
    display_cols = [
        "source", "target", "layer", "score", "emergence_stage",
        "loco_replication_fraction", "pp4_min", "string_score",
        "druggable_pair", "conserved_mouse_a", "conserved_mouse_b", "dim",
    ]
    present_cols = [c for c in display_cols if c in edges_df.columns]
    st.dataframe(
        edges_df[present_cols],
        hide_index=True, use_container_width=True, height=320,
    )
