"""Prioritize -- Convergence Matrix + Bayesian Ranking + Gene-Drug-TF Map.

Three internal tabs for prioritizing genes across the 8 modalities:
    * Convergence Matrix -- per-gene modality scores, heatmap + table.
    * Bayesian Ranking -- composite posterior ordering.
    * Gene-Drug-TF Map (formerly "Knowledge Graph") -- a small curated
      network linking top MASLD genes to their drug targets, regulatory
      TFs, and pathway context. Complements the Network page, which is
      the fully gene-centric force graph.
"""
from __future__ import annotations

import math

import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from components.layout import inject_css, page_header, sidebar_header
from data.loaders import (
    load_bayesian_ranking,
    load_convergence_matrix,
    load_knowledge_graph,
)
from theme.plotly_theme import PLOTLY_THEME, apply_theme

inject_css()
sidebar_header()


# User-facing modality order (2026-04-23).
# Data columns `m1..m8` in convergence_matrix.json are in the OLD positional
# scheme (human / mouse / genetic / essentiality / epigenomic / spatial /
# singlecell / proteomics) — the third tuple element remaps to the new
# display sequence: human → genetic → epigenomic → spatial → sc → proteomics
# → mouse → essentiality.
MODALITY_DISPLAY: list[tuple[str, str, str]] = [
    ("M1", "Human bulk DE",              "m1"),
    ("M2", "Genetic (GWAS+eQTL)",        "m3"),
    ("M3", "Epigenomic (ATAC/SCENIC+)",  "m5"),
    ("M4", "Spatial (Visium)",           "m6"),
    ("M5", "Single-cell (sc + LIANA)",   "m7"),
    ("M6", "Proteomics",                 "m8"),
    ("M7", "Mouse bulk DE",              "m2"),
    ("M8", "Essentiality (DepMap)",      "m4"),
]
# Legacy alias used by a couple of call sites + the "Sort by" dropdown.
MODALITY_HEADERS: list[tuple[str, str]] = [(m[0], m[1]) for m in MODALITY_DISPLAY]


page_header(
    "Prioritize",
    "Rank genes across the 8 modalities. Convergence Matrix shows per-modality "
    "scores; Bayesian Ranking gives a composite posterior; Gene-Drug-TF Map "
    "places top genes in their therapeutic + regulatory context.",
)

tabs = st.tabs(["Convergence Matrix", "Bayesian Ranking", "Gene-Drug-TF Map"])


# ---------------------------------------------------------------------------
# Tab 1 -- Convergence Matrix
# ---------------------------------------------------------------------------

with tabs[0]:
    st.caption("Primary modality-convergence view. Rows = genes, columns = 8 modalities.")
    try:
        df = load_convergence_matrix()
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    # Display order: remap data columns m1..m8 through MODALITY_DISPLAY so the
    # user-facing column sequence matches the new M1..M8 convention.
    display_cols = [m[2] for m in MODALITY_DISPLAY]      # e.g. ["m1","m3","m5", ...]
    display_labels = [f"{m[0]} {m[1]}" for m in MODALITY_DISPLAY]
    display_short = [m[0] for m in MODALITY_DISPLAY]

    c1, c2, c3, c4 = st.columns([2, 1, 1.5, 2])
    with c1:
        search = st.text_input("Filter by gene symbol", value="", placeholder="e.g. THRB", key="cm_search")
    with c2:
        max_count = int(df["count"].max()) if not df.empty else 8
        min_count = st.slider("Min modality count", 0, max(max_count, 8), 3)
    with c3:
        sort_options = ["count (desc)", "gene A-Z"] + [f"{m[0]} {m[1].split(' (')[0]}" for m in MODALITY_DISPLAY]
        sort_choice = st.selectbox("Sort by", sort_options)
    with c4:
        required_modalities = st.multiselect(
            "Require modality ≥ 0.5",
            options=display_cols,
            format_func=lambda k: f"{MODALITY_DISPLAY[display_cols.index(k)][0]}  "
                                  f"{MODALITY_DISPLAY[display_cols.index(k)][1]}",
        )

    filtered = df.copy()
    if search:
        filtered = filtered[
            filtered["gene"].astype(str).str.upper().str.contains(search.strip().upper())
        ]
    filtered = filtered[filtered["count"] >= min_count]
    for mcol in required_modalities:
        filtered = filtered[filtered[mcol] >= 0.5]

    if sort_choice == "count (desc)":
        filtered = filtered.sort_values(["count", "gene"], ascending=[False, True])
    elif sort_choice == "gene A-Z":
        filtered = filtered.sort_values("gene")
    else:
        # Match the prefix "M1 ", "M2 " … to the corresponding data column.
        for code, _long, dcol in MODALITY_DISPLAY:
            if sort_choice.startswith(code + " "):
                filtered = filtered.sort_values(dcol, ascending=False)
                break

    # Cap at 2000 genes (was 500 in earlier builds).
    shown = filtered.head(2000)
    max_mod = int(filtered["count"].max()) if not filtered.empty else 0
    st.markdown(
        f"Showing **{len(shown):,}** of **{len(df):,}** genes "
        f"(filtered: {len(filtered):,}). Highest convergence: **{max_mod}/8** modalities."
    )

    if not shown.empty:
        # Build z in the display order so columns line up with display_labels.
        z = shown[display_cols].values.tolist()
        y = shown["gene"].tolist()
        heatmap = go.Figure(data=go.Heatmap(
            z=z, x=display_labels, y=y,
            colorscale=[[0, "#fff7ed"], [1, "#c2410c"]],
            zmin=0, zmax=1, colorbar=dict(title="Score"),
            hovertemplate="<b>%{y}</b><br>%{x}: %{z:.3f}<extra></extra>",
        ))
        # Put column headers on TOP and give them room with a bigger top margin.
        apply_theme(heatmap, margin=dict(l=70, r=20, t=90, b=20))
        # Row height shrinks as the matrix grows so 2K genes stay tractable;
        # hover picks individual rows when the pixel height drops below ~6px.
        row_px = 14 if len(shown) < 300 else (8 if len(shown) < 800 else 4)
        heatmap.update_layout(
            height=max(400, min(row_px * len(shown) + 140, 6000)),
            xaxis=dict(
                side="top",
                tickangle=-35,
                tickfont=dict(size=11),
                ticklabelposition="outside top",
            ),
            yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
        )
        st.plotly_chart(heatmap, use_container_width=True)

    link_df = shown.copy()
    link_df["Gene page"] = link_df["gene"].apply(lambda s: f"/Gene?symbol={s}")
    # Build a rename map so the table shows M1…M8 column headers in the
    # requested order. display_cols is already in M1→M8 order above.
    table_cols = ["gene", "count", *display_cols, "is_deg", "is_coloc", "is_druggable", "Gene page"]
    col_config: dict = {
        "Gene page": st.column_config.LinkColumn("Gene page", display_text="open"),
    }
    for code, long_label, dcol in MODALITY_DISPLAY:
        col_config[dcol] = st.column_config.NumberColumn(
            code, format="%.3f", help=long_label,
        )
    st.dataframe(
        link_df[table_cols],
        column_config=col_config,
        use_container_width=True, hide_index=True,
    )


# ---------------------------------------------------------------------------
# Tab 2 -- Bayesian Ranking (swapped to 2nd; KG moves to 3rd)
# ---------------------------------------------------------------------------

with tabs[1]:
    st.caption(
        "Genes ordered by **posterior probability of convergent multi-modal "
        "biological signal** across 8 modalities (Script 46d). "
        "Tier 1 = genetically validated (PP4 > 0.7 + ≥2 active modalities); "
        "Tier 2 = ≥5 active modalities concordant; "
        "Tier 3 = 3-4 active modalities. "
        "Concordance state classifies each gene's evidence pattern "
        "(Concordant-up/down, Protective-LOF, Conflicted)."
    )

    try:
        br = load_bayesian_ranking()
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    highlights = ["THRB", "HKDC1", "GAS6", "CHI3L1", "ADH4", "TM6SF2"]
    flags = br[br["symbol"].isin(highlights)].copy() if "symbol" in br.columns else pd.DataFrame()
    if not flags.empty:
        bits = [f"**{row['symbol']}** (rank {int(row['rank'])}, {row.get('tier','')})"
                for _, row in flags.iterrows()]
        st.info("Highlighted targets — " + "  &middot;  ".join(bits))

    # Tier filter
    if "tier" in br.columns:
        tier_options = sorted(br["tier"].dropna().unique().tolist())
        tier_pick = st.multiselect(
            "Filter by tier",
            options=tier_options,
            default=[t for t in tier_options if t.startswith(("1_", "2_", "3_"))],
            key="br_tier",
        )
        if tier_pick:
            br = br[br["tier"].isin(tier_pick)]

    q = st.text_input("Filter by gene symbol", value="", placeholder="e.g. THRB", key="br_search")
    view = br.copy()
    if q and "symbol" in view.columns:
        view = view[view["symbol"].astype(str).str.upper().str.contains(q.strip().upper())]

    if "symbol" in view.columns:
        view["Gene page"] = view["symbol"].apply(lambda s: f"/Gene?symbol={s}")

    preferred = [
        "rank", "symbol", "posterior_prob", "tier", "concordance_state",
        "dominant_stage_S1", "n_modalities_active", "druggability_tier",
        "coloc_best_pp4", "coloc_best_gwas", "cross_species_concordant",
        "is_deg", "is_coloc", "is_druggable", "Gene page",
    ]
    display_cols = [c for c in preferred if c in view.columns]
    for c in view.columns:
        if c not in display_cols:
            display_cols.append(c)

    st.dataframe(
        view[display_cols],
        column_config={
            "rank": st.column_config.NumberColumn("rank", format="%d"),
            "posterior_prob": st.column_config.NumberColumn(
                "P(convergent)", format="%.4f",
                help="Empirically-calibrated posterior probability of convergent multi-modal signal.",
            ),
            "tier": st.column_config.TextColumn(
                "tier",
                help=(
                    "1_Genetic_validated = strong COLOC or TWAS + ≥2 active modalities | "
                    "2_Convergent = ≥5 active modalities, concordant direction | "
                    "3_Suggestive = 3-4 active modalities | "
                    "4_Weak = background"
                ),
            ),
            "concordance_state": st.column_config.TextColumn(
                "concordance state",
                help=(
                    "Concordant-up/down = ≥70% direction agreement | "
                    "Protective-LOF = genetic + expression-down (inhibitor-target signature) | "
                    "Conflicted = mixed directions | "
                    "Insufficient-evidence = <2 active signed modalities"
                ),
            ),
            "dominant_stage_S1": st.column_config.TextColumn(
                "dominant stage",
                help="Which S1 sub-contrast drives the signal: early_disease / nafl_to_nash / f2_switch / advanced_fibrosis / sex_dimorphic.",
            ),
            "n_modalities_active": st.column_config.NumberColumn(
                "n mod", format="%d",
                help="Number of modalities (of 8) with log-BF > log(3).",
            ),
            "druggability_tier": st.column_config.TextColumn("druggability"),
            "coloc_best_pp4": st.column_config.NumberColumn(
                "best PP4", format="%.3f",
                help="Highest SuSiE PP.H4 across the 28-GWAS portfolio.",
            ),
            "coloc_best_gwas": st.column_config.TextColumn("best GWAS"),
            "cross_species_concordant": st.column_config.CheckboxColumn(
                "mouse↔human concordant",
                help="Mouse meta and human dream agree in direction at padj < 0.05.",
            ),
            "Gene page": st.column_config.LinkColumn("Gene page", display_text="open"),
        },
        hide_index=True, use_container_width=True, height=600,
    )


# ---------------------------------------------------------------------------
# Tab 3 -- Gene-Drug-TF Map (renamed from Knowledge Graph)
# ---------------------------------------------------------------------------

with tabs[2]:
    st.caption(
        "Curated network linking top MASLD genes to their drug targets, regulatory "
        "TFs, and pathway context. For the gene-centric force graph use the dedicated "
        "**Network** page."
    )
    kg = load_knowledge_graph()
    if kg is None:
        st.warning("knowledge_graph.json not found in DATA_DIR.")
    else:
        kg_nodes = kg.get("nodes", [])
        kg_edges = kg.get("edges", kg.get("links", []))

        type_mix = pd.Series([n.get("type", "gene") for n in kg_nodes]).value_counts()
        gene_nodes_raw = [n for n in kg_nodes if n.get("type", "gene") == "gene"]
        n_genes_loaded = len(gene_nodes_raw)

        # Explain what's in the default graph ------------------------------
        with st.expander("How are genes selected for this map?", expanded=False):
            st.markdown(
                f"""
                The map is pre-built by
                `masld-atlas-v2/scripts/generate_phase4_data.py` from the
                multi-evidence atlas. Default population:

                1. **Top 200 genes by modality convergence** (`layers_active`
                   in the atlas — how many of the 8 independent modalities
                   show signal for that gene).
                2. **~69 extra genes** added only because they are a drug
                   target, SCENIC+ TF target, or pathway member — so every
                   rendered edge has both endpoints in the graph.

                You're looking at **{n_genes_loaded:,} genes · {int(type_mix.get('drug', 0))} drugs
                · {int(type_mix.get('tf', 0))} TFs · {int(type_mix.get('pathway', 0))} pathways**
                loaded from `knowledge_graph.json`. The thresholds below
                filter *within* this set — they let you tighten (show only
                high-convergence / DEG / COLOC-supported genes) or relax (show
                everything) without rebuilding. To pull in genes that aren't
                in the loaded 269, the generator's `head(200)` cap would need
                to be raised (offline rebuild + data-tarball re-push).
                """
            )

        # Threshold row ----------------------------------------------------
        tcol1, tcol2, tcol3 = st.columns([1.3, 1, 1])
        with tcol1:
            min_modalities = st.slider(
                "Min modality convergence (layers_active ≥)",
                min_value=0, max_value=8, value=0, step=1,
                help=(
                    "Restrict gene nodes to those active in at least N of the 8 "
                    "independent modalities. 0 = show every gene in the loaded "
                    "set; 4+ reveals the heavy-convergence hubs."
                ),
            )
        with tcol2:
            require_deg = st.checkbox(
                "DEG only",
                value=False,
                help="Restrict to genes that are differentially expressed in the "
                     "human bulk RNA-seq mega-analysis (padj < 0.05).",
            )
        with tcol3:
            require_coloc = st.checkbox(
                "Genetic evidence only",
                value=False,
                help="Restrict to genes that colocalize with a GWAS liver-trait "
                     "signal (any COLOC PP.H4 ≥ 0.5).",
            )

        # Display caps / entity filter
        mcol1, mcol2 = st.columns([1, 2])
        with mcol1:
            max_nodes = st.slider(
                "Max nodes to display",
                min_value=50, max_value=1000, value=900, step=50,
            )
        with mcol2:
            node_type_filter = st.multiselect(
                "Entity types",
                options=sorted(type_mix.index.tolist()),
                default=sorted(type_mix.index.tolist()),
                help="Uncheck to hide drugs, TFs, or pathways.",
            )

        # Apply the biology filters ONLY to gene nodes — drugs/TFs/pathways
        # don't have layers_active / is_deg / is_coloc to filter on.
        def _gene_passes(n: dict) -> bool:
            if min_modalities > 0 and int(n.get("evidence_score", 0) or 0) < min_modalities:
                return False
            if require_deg and not bool(n.get("is_deg", False)):
                return False
            if require_coloc and not bool(n.get("is_coloc", False)):
                return False
            return True

        typ_set = set(node_type_filter)
        filtered_nodes = []
        for n in kg_nodes:
            t = n.get("type", "gene")
            if t not in typ_set:
                continue
            if t == "gene" and not _gene_passes(n):
                continue
            filtered_nodes.append(n)

        # Sort genes by convergence so the top-N cap keeps the best-connected.
        filtered_nodes.sort(
            key=lambda n: (n.get("type") == "gene", int(n.get("evidence_score", 0) or 0)),
            reverse=True,
        )
        keep_ids = {n["id"] for n in filtered_nodes[:max_nodes]}

        # Second pass: drop non-gene nodes that lost all edges after filtering.
        edge_subset = [
            e for e in kg_edges
            if e.get("source") in keep_ids and e.get("target") in keep_ids
        ]
        endpoints_in_edges = set()
        for e in edge_subset:
            endpoints_in_edges.add(e.get("source"))
            endpoints_in_edges.add(e.get("target"))
        node_subset = [
            n for n in filtered_nodes
            if n["id"] in keep_ids and (
                n.get("type") == "gene" or n["id"] in endpoints_in_edges
            )
        ]
        keep_ids = {n["id"] for n in node_subset}
        edge_subset = [
            e for e in edge_subset
            if e.get("source") in keep_ids and e.get("target") in keep_ids
        ]

        # Active-filter summary line
        n_gene_kept = sum(1 for n in node_subset if n.get("type") == "gene")
        filter_bits = []
        if min_modalities > 0:
            filter_bits.append(f"≥{min_modalities} modalities")
        if require_deg:
            filter_bits.append("DEG-only")
        if require_coloc:
            filter_bits.append("genetic-only")
        filter_str = " · ".join(filter_bits) if filter_bits else "no biology filters"
        st.caption(
            f"**{n_gene_kept}/{n_genes_loaded}** genes visible  ·  "
            f"{len(node_subset)} total nodes  ·  {len(edge_subset)} edges  ·  "
            f"filters: {filter_str}"
        )

        G = nx.Graph()
        G.add_nodes_from(keep_ids)
        for e in edge_subset:
            G.add_edge(e["source"], e["target"])
        k_val = 1.5 / math.sqrt(max(1, G.number_of_nodes()))
        pos = nx.spring_layout(G, seed=42, k=k_val, iterations=40)

        ex: list[float | None] = []
        ey: list[float | None] = []
        for e in edge_subset:
            if e["source"] not in pos or e["target"] not in pos:
                continue
            ex.extend([pos[e["source"]][0], pos[e["target"]][0], None])
            ey.extend([pos[e["source"]][1], pos[e["target"]][1], None])

        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=ex, y=ey, mode="lines",
            line=dict(color="#b0b7c1", width=0.6),
            opacity=0.35, hoverinfo="skip", showlegend=False,
        ))

        type_colors = {
            "gene": "#3b5ad2",
            "drug": "#d67f3c",
            "pathway": "#23a369",
            "tf": "#9b59b6",
        }
        xs, ys, cs, ts, ls = [], [], [], [], []
        for n in node_subset:
            xs.append(pos.get(n["id"], (0.0, 0.0))[0])
            ys.append(pos.get(n["id"], (0.0, 0.0))[1])
            cs.append(type_colors.get(n.get("type", "gene"), "#737782"))
            ts.append(n.get("label", n.get("id", "")))
            ls.append(n.get("type", "gene"))
        fig.add_trace(go.Scattergl(
            x=xs, y=ys, mode="markers",
            marker=dict(size=7, color=cs, line=dict(color="white", width=0.5)),
            hovertemplate="<b>%{text}</b><br>type: %{customdata}<extra></extra>",
            text=ts, customdata=ls, showlegend=False,
        ))
        apply_theme(fig, margin=dict(l=10, r=10, t=10, b=10))
        fig.update_layout(
            height=600,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False, scaleanchor="x", scaleratio=1),
        )
        st.plotly_chart(fig, use_container_width=True)
        st.caption(
            f"Showing {len(node_subset)} of {len(kg_nodes)} nodes; "
            f"{len(edge_subset)} of {len(kg_edges)} edges."
        )


