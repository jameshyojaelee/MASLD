"""Atlas view -- Overview + Drugs + Stratifications (3 internal tabs).

Absorbs the previous Atlas + Explorer + Drugs + Resolution pages.

- **Overview**: dataset inventory, stat cards, volcano plot, filterable
  34K-gene DEG table (toggle "DEGs only | All genes").
- **Drugs**: clinical MASLD drug table + evidence-tier funnel.
- **Stratifications**: 4 stacked-bar summaries for sex / attribution /
  ferroptosis / zonation classes.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from components.layout import eyebrow, inject_css, page_header, sidebar_header, stat_card
from components.modality_radar import make_radar_wall_svg
from data.loaders import (
    load_atlas_parquet,
    load_atlas_summary,
    load_clinical_drugs,
    load_drugs_data,
    load_gene_index,
)
from theme.plotly_theme import DEG_DOWN, DEG_NONSIG, DEG_UP, PLOTLY_THEME, apply_theme

inject_css()
sidebar_header()

page_header(
    "MASLD Atlas",
    "Human RNA-seq mega-analysis across 10 patient cohorts (1,444 QC-passing "
    "samples). Volcano, 34K-gene evidence table, clinical drug evidence, and "
    "four stratification axes. Other modalities (scRNA, ATAC, proteomics, "
    "spatial, mouse, GWAS) live in the Prioritize / Validation / Network tabs.",
)


# ---------------------------------------------------------------------------
# Top-of-page stat cards (shared across tabs)
# ---------------------------------------------------------------------------

try:
    summary = load_atlas_summary()
except FileNotFoundError as exc:
    st.error(str(exc))
    st.stop()

top_stats = st.columns(4)
with top_stats[0]:
    stat_card("Human RNA-seq cohorts", summary.get("total_cohorts", 10))
with top_stats[1]:
    stat_card("Human RNA-seq samples", summary.get("total_samples", 1444))
with top_stats[2]:
    stat_card("Integrated DEGs", summary.get("total_degs", 5484))
with top_stats[3]:
    stat_card("Human-Mouse Concordant", summary.get("conserved_count", "--"))


eyebrow("Cross-modal coverage")
xmcols = st.columns(6)
with xmcols[0]:
    stat_card("Bulk RNA-seq", "19", "10 human · 9 mouse")
with xmcols[1]:
    stat_card("scRNA-seq", "7", "liver single-cell")
with xmcols[2]:
    stat_card("Spatial", "3", "Visium tissue arrays")
with xmcols[3]:
    stat_card("ATAC-seq", "2", "bulk + scATAC")
with xmcols[4]:
    stat_card("Proteomics", "3", "plasma + tissue DIA-MS")
with xmcols[5]:
    stat_card("GWAS", "13", "+ 3 liver eQTL panels")
st.caption(
    "This Atlas page focuses on the human bulk RNA-seq layer. For per-modality "
    "detail see **Prioritize** (modality convergence), **Network** (edge layers), "
    "or **Validation** (cross-species + proteomics)."
)


tabs = st.tabs(["Overview", "Drugs", "Stratifications"])


# ---------------------------------------------------------------------------
# Tab 1 — Overview
# ---------------------------------------------------------------------------

# Per-cohort QC-pass counts + NAS-harmonised diagnosis breakdown.
# Numbers derived from:
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv
#   RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv
# These are the same QC-pass denominators the figure 1b treemap uses
# (scripts/figures/fig1_atlas_overview_v2.R). Total = 1,444 QC-pass samples.
# Update by re-running the grouping query in that script whenever the QC
# report changes.
COHORTS = [
    {"dataset": "GSE213621",   "reference": "Chen 2023",       "condition": "Fibrosis staging",       "samples": 359, "control": 67, "nafl": 0,  "borderline": 0,   "nash": 0,   "fibrosis_only": 292, "unknown": 0},
    {"dataset": "GSE135251",   "reference": "Govaere 2020",    "condition": "MASLD staging (NAS)",    "samples": 215, "control": 10, "nafl": 31, "borderline": 64,  "nash": 110, "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "PRJNA512027", "reference": "Gerhard 2018",    "condition": "NASH/Healthy",           "samples": 185, "control": 34, "nafl": 49, "borderline": 0,   "nash": 102, "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "GSE193066",   "reference": "Hoshida 2022",    "condition": "MASLD staging",          "samples": 160, "control": 0,  "nafl": 14, "borderline": 101, "nash": 45,  "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "GSE162694",   "reference": "Bril 2021",       "condition": "Steatosis/NASH",         "samples": 142, "control": 31, "nafl": 21, "borderline": 24,  "nash": 40,  "fibrosis_only": 0,   "unknown": 26},
    {"dataset": "GSE167523",   "reference": "Pantano 2021",    "condition": "NAFLD progression",      "samples": 96,  "control": 0,  "nafl": 51, "borderline": 0,   "nash": 45,  "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "GSE174478",   "reference": "Kawamura 2021",   "condition": "NAFL/NASH",              "samples": 93,  "control": 0,  "nafl": 12, "borderline": 26,  "nash": 55,  "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "GSE130970",   "reference": "Hoang 2019",      "condition": "NASH vs Healthy",        "samples": 76,  "control": 4,  "nafl": 13, "borderline": 33,  "nash": 26,  "fibrosis_only": 0,   "unknown": 0},
    {"dataset": "GSE240729",   "reference": "Verschuren 2024", "condition": "Fibrosis",               "samples": 63,  "control": 0,  "nafl": 0,  "borderline": 0,   "nash": 0,   "fibrosis_only": 63,  "unknown": 0},
    {"dataset": "GSE126848",   "reference": "Suppli 2019",     "condition": "NAFL/NASH/Fibrosis",     "samples": 55,  "control": 25, "nafl": 15, "borderline": 0,   "nash": 15,  "fibrosis_only": 0,   "unknown": 0},
]


ATLAS_COLS = (
    "human_symbol",
    "ensembl_id",
    "gene_biotype",
    "dream_logFC",
    "dream_padj",
    "dream_tstat",
    "is_conserved",
    "dgidb_druggable",
)


def _load_full_gene_index_for_table() -> pd.DataFrame:
    """Prefer the pre-built gene_index (34K genes w/ evidence scores).

    Returns a DataFrame with ``symbol``, ``biotype``, ``dream_logfc``,
    ``dream_padj`` and evidence scores s1..s8. Callers filter / sort further.
    """
    try:
        idx = load_gene_index()
    except FileNotFoundError:
        return pd.DataFrame()
    return idx


with tabs[0]:
    # Load the atlas first so the two-column section below can use it.
    try:
        atlas = load_atlas_parquet(columns=ATLAS_COLS)
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    atlas = atlas.dropna(subset=["dream_logFC", "dream_padj"]).copy()
    atlas["neglog10_padj"] = -np.log10(atlas["dream_padj"].clip(lower=1e-300))
    atlas["neglog10_padj"] = atlas["neglog10_padj"].clip(upper=300)

    sig_mask = (atlas["dream_padj"] < 0.05) & (atlas["dream_logFC"].abs() > 0.2)
    atlas["is_sig"] = sig_mask
    atlas["is_deg"] = (atlas["dream_padj"] < 0.05) & (atlas["dream_logFC"].abs() > 0.3)

    nonsig = atlas.loc[~sig_mask]
    sig_up = atlas.loc[sig_mask & (atlas["dream_logFC"] > 0)]
    sig_down = atlas.loc[sig_mask & (atlas["dream_logFC"] < 0)]

    # Side-by-side: left column = dataset inventory, right column = compact volcano.
    # Ratio [1.25, 1] gives the table a bit more breathing room since it has
    # more columns than the plot has useful pixels.
    inv_col, volc_col = st.columns([1.25, 1], gap="large")

    with inv_col:
        st.markdown("### Human RNA-seq dataset inventory")
        st.caption(
            "10 cohorts · 1,444 QC-pass samples · NAS-harmonised diagnosis "
            "(source: fig1b of the manuscript)."
        )
        cohort_df = pd.DataFrame(COHORTS)
        st.dataframe(
            cohort_df,
            hide_index=True,
            use_container_width=True,
            height=430,
            column_config={
                "dataset":       st.column_config.TextColumn("Dataset", width="small"),
                "reference":     st.column_config.TextColumn("Reference"),
                "condition":     st.column_config.TextColumn("Condition"),
                "samples":       st.column_config.NumberColumn("QC-pass n", format="%d", width="small"),
                "control":       st.column_config.NumberColumn("Ctrl",  format="%d", help="Non-MASLD controls"),
                "nafl":          st.column_config.NumberColumn("NAFL",  format="%d", help="NAFL / MASL"),
                "borderline":    st.column_config.NumberColumn("Border.", format="%d", help="Borderline NASH (NAS 3–4, Kleiner 2005)"),
                "nash":          st.column_config.NumberColumn("NASH",  format="%d", help="NASH / MASH (NAS ≥ 5)"),
                "fibrosis_only": st.column_config.NumberColumn("Fib-only", format="%d", help="Fibrosis-staging cohorts with no NAS label"),
                "unknown":       st.column_config.NumberColumn("Unk.",  format="%d", help="QC-pass samples with missing diagnosis"),
            },
        )

    with volc_col:
        st.markdown("### Integrated volcano")
        st.caption(
            "Human RNA-seq mega-analysis (dream). "
            "Significant: padj < 0.05 AND |logFC| > 0.2. y-axis clamped to 300."
        )
        vfig = go.Figure()
        vfig.add_trace(go.Scattergl(
            x=nonsig["dream_logFC"], y=nonsig["neglog10_padj"], mode="markers", name="n.s.",
            marker=dict(color=DEG_NONSIG, size=4, opacity=0.35, line=dict(width=0)),
            text=nonsig["human_symbol"],
            hovertemplate="<b>%{text}</b><br>logFC=%{x:.2f}<br>-log10 padj=%{y:.1f}<extra></extra>",
        ))
        vfig.add_trace(go.Scattergl(
            x=sig_up["dream_logFC"], y=sig_up["neglog10_padj"], mode="markers", name="Up",
            marker=dict(color=DEG_UP, size=5, opacity=0.85, line=dict(width=0)),
            text=sig_up["human_symbol"],
            hovertemplate="<b>%{text}</b><br>logFC=%{x:.2f}<br>-log10 padj=%{y:.1f}<extra></extra>",
        ))
        vfig.add_trace(go.Scattergl(
            x=sig_down["dream_logFC"], y=sig_down["neglog10_padj"], mode="markers", name="Down",
            marker=dict(color=DEG_DOWN, size=5, opacity=0.85, line=dict(width=0)),
            text=sig_down["human_symbol"],
            hovertemplate="<b>%{text}</b><br>logFC=%{x:.2f}<br>-log10 padj=%{y:.1f}<extra></extra>",
        ))
        vfig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="#9aa3ad", line_width=1, opacity=0.55)
        vfig.add_vline(x=0.2, line_dash="dash", line_color="#9aa3ad", line_width=1, opacity=0.55)
        vfig.add_vline(x=-0.2, line_dash="dash", line_color="#9aa3ad", line_width=1, opacity=0.55)
        apply_theme(vfig, margin=dict(l=40, r=10, t=10, b=35))
        vfig.update_layout(
            height=410,
            xaxis_title="Integrated logFC",
            yaxis_title="-log10(padj)",
            legend=dict(
                orientation="h", yanchor="bottom", y=1.01, xanchor="right", x=1,
                font=dict(size=11),
            ),
        )
        st.plotly_chart(vfig, use_container_width=True, config={"displayModeBar": False})

    # Load the full gene index once — shared by the radar wall + gene table.
    gene_idx = _load_full_gene_index_for_table()

    # -- Radar wall: small-multiples evidence fingerprints -------------------
    # Lets users scan the *shape* of each gene's modality profile side-by-side
    # — a view that neither the convergence matrix (one cell per modality) nor
    # the gene table (numeric columns) gives you. Biological use: spot genes
    # with unusual combinations ("spatial + genetic only, no expression"),
    # triage a top-N list by modality, compare candidates side-by-side.
    st.markdown("### Evidence radar wall")
    st.caption(
        "Each tile = one gene's 8-modality polygon (M1–M8). Shape, not colour, "
        "carries the signal — a balanced octagon means all modalities converge; "
        "a spike means a single modality dominates. Click any tile to open the "
        "gene card."
    )

    if not gene_idx.empty:
        rw_c1, rw_c2, rw_c3, rw_c4 = st.columns([1.3, 1.1, 1.0, 1.0])
        # Sort modes surface different "discovery" angles. "Balanced" rewards
        # genes high across many modalities; the per-M codes let you drill into
        # "top-M3 (epigenomic)" or "top-M4 (spatial)" specifically.
        sort_options = {
            "Most modalities (balanced convergence)": "_balanced",
            "|Integrated logFC| (M1 effect size)":    "_abs_lfc",
            "M1 Human bulk DE":                        "s1_human",
            "M2 Genetic (GWAS + eQTL)":               "s3_genetic",
            "M3 Epigenomic (ATAC / SCENIC+)":          "s5_epigenomic",
            "M4 Spatial (Visium)":                     "s6_spatial",
            "M5 Single-cell":                          "s7_singlecell",
            "M6 Proteomics":                           "s8_proteomics",
            "M7 Mouse bulk DE":                        "s2_mouse",
            "M8 Essentiality (DepMap)":                "s4_essential",
        }
        sort_label = rw_c1.selectbox(
            "Sort by", list(sort_options.keys()), index=0,
            help="Controls which genes surface first in the wall. Change this "
                 "to triage by a specific modality.",
        )
        min_layers = rw_c2.slider(
            "Min modalities active", 0, 8, 3,
            help="Each gene's `layers_active` is the count of modalities with "
                 "non-zero evidence. Raise this to hide single-modality hits.",
        )
        deg_only_wall = rw_c3.checkbox(
            "DEGs only", value=False, key="radarwall_deg",
            help="Restrict to human bulk RNA-seq DEGs (padj < 0.05, |logFC| > 0.3).",
        )
        protein_only = rw_c4.checkbox(
            "Protein-coding", value=True,
            help="Hide lncRNAs / pseudogenes. Uncheck to include all biotypes.",
        )

        rw_view = gene_idx.copy()
        # Recompute "layers active" from the 8 canonical modality columns so
        # that the displayed N/8 always matches what the polygon shows.
        # (The upstream `layers_active` field in gene_index.json counts more
        # than 8 sub-modalities — e.g. FAH = 9 — which is confusing here.)
        _ev_cols_for_count = [c for c in ["s1_human", "s2_mouse", "s3_genetic",
                                           "s4_essential", "s5_epigenomic",
                                           "s6_spatial", "s7_singlecell",
                                           "s8_proteomics"] if c in rw_view.columns]
        if _ev_cols_for_count:
            rw_view["active_n8"] = (rw_view[_ev_cols_for_count] > 0).sum(axis=1)
        else:
            rw_view["active_n8"] = 0
        if min_layers > 0:
            rw_view = rw_view[rw_view["active_n8"] >= min_layers]
        if deg_only_wall and "is_deg" in rw_view.columns:
            rw_view = rw_view[rw_view["is_deg"]]
        if protein_only and "biotype" in rw_view.columns:
            rw_view = rw_view[rw_view["biotype"] == "protein_coding"]

        sort_key = sort_options[sort_label]
        if sort_key == "_balanced":
            # Secondary sort by summed evidence so ties break meaningfully.
            rw_view["_ev_sum"] = rw_view[_ev_cols_for_count].sum(axis=1) if _ev_cols_for_count else 0
            rw_view = rw_view.sort_values(
                ["active_n8", "_ev_sum"], ascending=[False, False]
            )
        elif sort_key == "_abs_lfc" and "dream_logfc" in rw_view.columns:
            rw_view = rw_view.assign(_abs=rw_view["dream_logfc"].abs())
            rw_view = rw_view.sort_values("_abs", ascending=False)
        elif sort_key in rw_view.columns:
            rw_view = rw_view.sort_values(sort_key, ascending=False)

        N_TILES = 24
        top = rw_view.head(N_TILES)
        st.caption(
            f"Showing top **{len(top)}** of {len(rw_view):,} genes matching "
            f"filters (sorted by *{sort_label}*)."
        )

        ev_cols = ["s1_human", "s2_mouse", "s3_genetic", "s4_essential",
                   "s5_epigenomic", "s6_spatial", "s7_singlecell", "s8_proteomics"]

        tiles: list[str] = []
        for _, row in top.iterrows():
            sym = str(row.get("symbol", "")) or ""
            if not sym:
                continue
            evidence = {c: float(row.get(c, 0.0) or 0.0) for c in ev_cols if c in top.columns}
            svg = make_radar_wall_svg(evidence, size=92)
            lfc = row.get("dream_logfc", None)
            lfc_html = ""
            if lfc is not None and pd.notna(lfc):
                cls = "up" if float(lfc) > 0 else "down"
                lfc_html = f'<span class="{cls}">{float(lfc):+.2f}</span>'
            layers = row.get("active_n8", None)
            layers_html = f" · {int(layers)}/8" if layers is not None and pd.notna(layers) else ""
            tiles.append(
                f'<a class="masld-radar-tile" href="/Gene?symbol={sym}" target="_self" '
                f'title="{sym} — open gene card">'
                f'{svg}'
                f'<div class="masld-radar-tile__sym">{sym}</div>'
                f'<div class="masld-radar-tile__meta">{lfc_html}{layers_html}</div>'
                f'</a>'
            )

        if tiles:
            st.markdown(
                f'<div class="masld-radar-wall">{"".join(tiles)}</div>',
                unsafe_allow_html=True,
            )
            # Compact reference legend so a first-time visitor can decode
            # the 8 axis codes without navigating away.
            legend_items = [
                ("M1", "Human bulk DE"),
                ("M2", "Genetic (GWAS+eQTL)"),
                ("M3", "Epigenomic (ATAC/SCENIC+)"),
                ("M4", "Spatial (Visium)"),
                ("M5", "Single-cell"),
                ("M6", "Proteomics"),
                ("M7", "Mouse bulk DE"),
                ("M8", "Essentiality (DepMap)"),
            ]
            legend_html = "".join(f"<div><b>{k}</b>{v}</div>" for k, v in legend_items)
            st.markdown(
                f'<div class="masld-radar-legend">{legend_html}</div>',
                unsafe_allow_html=True,
            )
        else:
            st.info("No genes match the current filters — relax the threshold or biotype.")

    # -- Gene table (single filterable table, replaces old Explorer) -------
    st.markdown("### Gene table")
    st.caption(
        "Filterable across all 34K genes (toggle 'DEGs only' to restrict to the 5,484 integrated "
        "DEGs). **Sorted by |integrated logFC|** — for composite-score ordering see "
        "*Translate → Bayesian Ranking*."
    )

    if gene_idx.empty:
        st.warning("gene_index.json not found; falling back to atlas.parquet.")
        base = atlas.rename(columns={
            "human_symbol": "symbol",
            "gene_biotype": "biotype",
            "dream_logFC": "dream_logfc",
            "dream_padj": "dream_padj",
        })
    else:
        base = gene_idx

    # Filter row
    fc1, fc2, fc3, fc4 = st.columns([2, 1.3, 1.3, 2])
    query = fc1.text_input("Search symbol", value="", placeholder="e.g. THRB, PNPLA3")
    degs_only = fc2.checkbox(
        "DEGs only",
        value=False,
        help="Differentially expressed in human bulk RNA-seq (padj < 0.05, |logFC| > 0.3).",
    )
    coloc_only = fc3.checkbox(
        "Genetic evidence",
        value=False,
        help="Gene colocalizes with a GWAS liver-trait signal (COLOC PP.H4 ≥ 0.5).",
    )
    biotype_opts = sorted(base["biotype"].dropna().unique().tolist()) if "biotype" in base.columns else []
    biotypes = fc4.multiselect("Biotype", options=biotype_opts, default=[])

    view = base
    if query:
        q = query.strip().upper()
        view = view[view["symbol"].astype(str).str.upper().str.contains(q, na=False)]
    if degs_only and "is_deg" in view.columns:
        view = view[view["is_deg"]]
    if coloc_only and "s3_genetic" in view.columns:
        view = view[view["s3_genetic"].astype(float) > 0.5]
    if biotypes and "biotype" in view.columns:
        view = view[view["biotype"].isin(biotypes)]

    view = view.copy()
    if "dream_logfc" in view.columns:
        view["abs_lfc"] = view["dream_logfc"].abs()
    elif "dream_logFC" in view.columns:
        view["abs_lfc"] = view["dream_logFC"].abs()
    else:
        view["abs_lfc"] = 0.0
    view = view.sort_values("abs_lfc", ascending=False).drop(columns=["abs_lfc"])

    st.caption(f"{len(view):,} genes shown (sorted by |logFC|).")

    view["Gene page"] = view["symbol"].apply(lambda s: f"/Gene?symbol={s}")

    # Modality column ordering for the gene table — must match the M1..M8
    # labels used in the Convergence Matrix. Each tuple:
    #   (data column in gene_index.json,
    #    compact header with modality suffix [shown in column name],
    #    full label [shown in header tooltip]).
    ATLAS_MODALITIES: list[tuple[str, str, str]] = [
        ("s1_human",      "M1 Human",    "M1 · Human bulk DE"),
        ("s3_genetic",    "M2 Genetic",  "M2 · Genetic (GWAS + eQTL colocalization)"),
        ("s5_epigenomic", "M3 Epigen",   "M3 · Epigenomic (ATAC / SCENIC+)"),
        ("s6_spatial",    "M4 Spatial",  "M4 · Spatial (Visium)"),
        ("s7_singlecell", "M5 scRNA",    "M5 · Single-cell (pseudobulk + LIANA)"),
        ("s8_proteomics", "M6 Protein",  "M6 · Proteomics (tissue + plasma)"),
        ("s2_mouse",      "M7 Mouse",    "M7 · Mouse bulk DE"),
        ("s4_essential",  "M8 Essent",   "M8 · Essentiality (DepMap)"),
    ]
    modality_data_cols = [d for d, _, _ in ATLAS_MODALITIES if d in view.columns]

    display_cols = [c for c in [
        "symbol", "biotype", "dream_logfc", "dream_padj", "layers_active",
        *modality_data_cols,
        "is_deg", "Gene page",
    ] if c in view.columns]
    col_config = {
        "symbol": st.column_config.TextColumn("Symbol", width="small"),
        "biotype": st.column_config.TextColumn("Biotype", width="small"),
        "dream_logfc": st.column_config.NumberColumn("Integrated logFC", format="%+.3f"),
        "dream_padj": st.column_config.NumberColumn("padj", format="%.2e"),
        "layers_active": st.column_config.NumberColumn("Layers", format="%d"),
        "is_deg": st.column_config.CheckboxColumn("DEG"),
        "Gene page": st.column_config.LinkColumn("Open", display_text="View"),
    }
    for data_col, code, tip in ATLAS_MODALITIES:
        if data_col in view.columns:
            col_config[data_col] = st.column_config.NumberColumn(
                code, format="%.2f", help=tip,
            )

    PAGE_SIZE = 100
    n_pages = max(1, (len(view) + PAGE_SIZE - 1) // PAGE_SIZE)
    page = st.number_input("Page", min_value=1, max_value=n_pages, value=1, step=1)
    start = (page - 1) * PAGE_SIZE
    end = start + PAGE_SIZE
    st.dataframe(
        view.iloc[start:end][display_cols],
        hide_index=True,
        use_container_width=True,
        column_config=col_config,
    )


# ---------------------------------------------------------------------------
# Tab 2 — Drugs
# ---------------------------------------------------------------------------

with tabs[1]:
    pipeline = load_drugs_data()
    if pipeline and pipeline.get("funnel"):
        st.markdown("### Drug repurposing funnel")
        funnel_df = pd.DataFrame(pipeline["funnel"])
        fcols = st.columns(len(funnel_df))
        for col, row in zip(fcols, funnel_df.itertuples(index=False)):
            with col:
                stat_card(row.stage, f"{int(row.count):,}", row.description)

    df = load_clinical_drugs()
    if df.empty:
        st.warning("clinical_drug_validation_table.csv not found.")
    else:
        st.markdown("### Clinical MASLD drugs -- atlas evidence")

        if "atlas_support" in df.columns:
            tiers = sorted(df["atlas_support"].dropna().unique().tolist())
            selected_tiers = st.multiselect("Evidence tier", options=tiers, default=tiers)
            df = df[df["atlas_support"].isin(selected_tiers)]

        primary_cols = [
            "drug", "target_gene", "stage", "moa", "atlas_support",
            "dream_logFC", "dream_padj", "best_liver_enzyme_pp4",
            "is_conserved", "is_deg", "sex_class",
            "dgidb_druggable", "opentargets_drug", "lincs_reversal",
        ]
        present = [c for c in primary_cols if c in df.columns]
        view_df = df[present].copy()
        if "target_gene" in view_df.columns:
            view_df["Gene page"] = view_df["target_gene"].apply(
                lambda s: f"/Gene?symbol={s}" if isinstance(s, str) and s else ""
            )

        st.dataframe(
            view_df,
            hide_index=True,
            use_container_width=True,
            column_config={
                "dream_logFC": st.column_config.NumberColumn("logFC", format="%+.3f"),
                "dream_padj": st.column_config.NumberColumn("padj", format="%.2e"),
                "best_liver_enzyme_pp4": st.column_config.NumberColumn(
                    "liver enzyme PP4", format="%.3f"
                ),
                "Gene page": st.column_config.LinkColumn("Gene page", display_text="open"),
            },
        )
        st.caption(
            f"{len(df)} drugs. Source: "
            "`RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv`."
        )


# ---------------------------------------------------------------------------
# Tab 3 — Stratifications
# ---------------------------------------------------------------------------

with tabs[2]:
    st.caption(
        "Four orthogonal resolution axes from the atlas: sex-dimorphic class, "
        "deconvolution attribution (hepatocyte-intrinsic vs composition-driven), "
        "ferroptosis program, and pericentral/periportal zonation."
    )

    strat_cols = (
        "human_symbol", "dream_logFC", "dream_padj",
        "sex_class", "attribution_class", "ferroptosis_class", "zonation_class",
    )
    try:
        sdf = load_atlas_parquet(columns=strat_cols)
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    if "dream_padj" in sdf.columns and "dream_logFC" in sdf.columns:
        sdf["is_deg"] = (sdf["dream_padj"] < 0.05) & (sdf["dream_logFC"].abs() > 0.3)
    else:
        sdf["is_deg"] = False

    def _stacked_bar(col: str, title: str) -> None:
        if col not in sdf.columns:
            st.info(f"{col} not available.")
            return
        all_ct = sdf[col].fillna("None").value_counts().rename_axis("class").reset_index(name="count")
        all_ct["scope"] = "All genes"
        deg_ct = (
            sdf.loc[sdf["is_deg"], col]
            .fillna("None").value_counts()
            .rename_axis("class").reset_index(name="count")
        )
        deg_ct["scope"] = "DEGs only"
        combined = pd.concat([all_ct, deg_ct], ignore_index=True)

        # Build stacked bar via graph_objects (plotly.express unavailable in some envs)
        bfig = go.Figure()
        for cls in combined["class"].unique():
            sub = combined[combined["class"] == cls]
            bfig.add_trace(go.Bar(
                x=sub["scope"], y=sub["count"], name=str(cls),
                hovertemplate=f"<b>{cls}</b><br>%{{x}}: %{{y:,}}<extra></extra>",
            ))
        apply_theme(bfig, margin=dict(l=10, r=10, t=40, b=10))
        bfig.update_layout(
            barmode="stack",
            title=dict(text=title, x=0.02, font=dict(size=13)),
            height=340,
            legend=dict(orientation="h", y=-0.18, yanchor="top"),
            yaxis_title="Genes",
        )
        st.plotly_chart(bfig, use_container_width=True)

    r1 = st.columns(2)
    with r1[0]:
        _stacked_bar("sex_class", "Sex-dimorphic class")
    with r1[1]:
        _stacked_bar("attribution_class", "Deconvolution attribution")

    r2 = st.columns(2)
    with r2[0]:
        _stacked_bar("ferroptosis_class", "Ferroptosis program")
    with r2[1]:
        _stacked_bar("zonation_class", "Pericentral / periportal zonation")

    st.caption(
        "Class definitions: see `RNA-seq/Human/.../results/integration/` for sex_class "
        "(Concordant / Male_biased / Female_biased / Divergent) and attribution_class "
        "(Hepatocyte_intrinsic / Composition_driven / Not_significant)."
    )
