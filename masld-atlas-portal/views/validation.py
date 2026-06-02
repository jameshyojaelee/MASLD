"""Validation -- Cross-Species + Proteomics.

Two internal tabs:

- **Cross-Species**: searchable human × mouse logFC scatter with category
  filters + a unified filterable / searchable gene table across all 8
  concordance categories (not just Conserved).
- **Proteomics**: searchable mRNA-protein concordance with stratum bar chart,
  direction filter, and a single combined protein table.

The friendly-label layer (``CATEGORY_LABELS``) keeps UI copy in plain
English -- the upstream cross_species.json uses engineering labels like
``Conserved`` which confuse first-time viewers.
"""
from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from components.layout import inject_css, page_header, sidebar_header, stat_card
from data.loaders import (
    load_atlas_parquet,
    load_cross_species_data,
    load_proteomics_summary,
)
from theme.plotly_theme import PLOTLY_THEME, apply_theme

inject_css()
sidebar_header()

page_header(
    "Validation",
    "Cross-species concordance (human vs mouse DE) and tissue + plasma "
    "proteomics that validate the transcriptomic atlas.",
)


# ---------------------------------------------------------------------------
# Category label dictionary.
# Upstream `cross_species.json` uses engineering tags (`Conserved`,
# `Species_Discordant` ...) — we surface plain-English equivalents so
# first-time viewers understand the category without reading methods.
# Colour palette is consistent across scatter + table + bar plots.
# ---------------------------------------------------------------------------
CATEGORY_LABELS = {
    "Conserved":        "Human-Mouse Concordant",
    "Moderate_Concordance":  "Moderate Concordance",
    "Species_Discordant":    "Species Discordant",
    "Human_Enriched":        "Human-Only",
    "Mouse_Specific":        "Mouse-Only",
    "Diet_Selective":        "Diet-Selective (mouse)",
    "Not_Significant":       "Not Significant",
    "Unclassified":          "Unclassified",
}
CATEGORY_COLORS = {
    "Human-Mouse Concordant": "#2c9e50",   # green — the hero class
    "Moderate Concordance":   "#5ea9d8",
    "Species Discordant":     "#c85a5a",
    "Human-Only":             "#d48f3b",
    "Mouse-Only":             "#7a6bbf",
    "Diet-Selective (mouse)": "#b36fa2",
    "Not Significant":        "#a2a8b2",
    "Unclassified":           "#8a8f98",
}
CATEGORY_ORDER = list(CATEGORY_LABELS.values())


tabs = st.tabs(["Cross-Species", "Proteomics"])


# ===========================================================================
# Tab 1 — Cross-Species
# ===========================================================================
with tabs[0]:
    st.caption(
        "Each gene's human bulk DE (dream mega-analysis) vs mouse meta-analysis "
        "logFC is classified into 8 categories based on direction + significance "
        "in both species. **Concordance-tested** means the gene has a measurable "
        "logFC in both human and mouse arms; only these enter the scatter."
    )

    raw = load_cross_species_data()
    if raw is not None:
        df = pd.DataFrame(raw.get("genes", []))
    else:
        st.info("cross_species.json absent — synthesising view from atlas.parquet.")
        try:
            raw_df = load_atlas_parquet(columns=(
                "human_symbol", "dream_logFC", "dream_padj",
                "mouse_meta_logFC", "mouse_meta_padj",
                "is_conserved", "n_diets_sig",
            ))
        except FileNotFoundError as exc:
            st.error(str(exc))
            st.stop()
        df = raw_df.rename(columns={
            "human_symbol": "symbol",
            "dream_logFC": "human_logfc",
            "mouse_meta_logFC": "mouse_logfc",
        })
        df["category"] = df["is_conserved"].map(
            {True: "Conserved", False: "Unclassified"}
        ).fillna("Unclassified")

    if df.empty:
        st.info("No cross-species data available.")
        st.stop()

    # Apply friendly label to a new column — keep the raw `category` for audit.
    df["category_label"] = df["category"].map(CATEGORY_LABELS).fillna(df["category"])

    # ---------------------------------------------------------------------
    # Expand coverage: cross_species.json only carries ~3,800 curated genes
    # (mostly protein-coding). Merge in every atlas row that has measurable
    # human + mouse logFC so users can query lncRNAs / MT-* / LINC / MIR etc.
    # On-the-fly categories (simple significance rules) are assigned to any
    # gene not already classified upstream.
    # ---------------------------------------------------------------------
    try:
        atlas_full = load_atlas_parquet(columns=(
            "human_symbol", "gene_biotype", "dream_logFC", "dream_padj",
            "mouse_meta_logFC", "mouse_meta_padj", "n_diets_sig",
        )).rename(columns={
            "human_symbol": "symbol",
            "gene_biotype": "biotype",
            "dream_logFC": "human_logfc",
            "dream_padj": "human_padj",
            "mouse_meta_logFC": "mouse_logfc",
            "mouse_meta_padj": "mouse_padj",
        })
    except FileNotFoundError:
        atlas_full = pd.DataFrame()

    # Index curated labels for quick lookup.
    curated = df.set_index("symbol")[["category_label"]].to_dict()["category_label"] if not df.empty else {}
    # Translatability + n_concordant are only in cross_species.json — preserve.
    curated_meta = df.set_index("symbol")[[c for c in ("translatability_score", "n_concordant") if c in df.columns]].to_dict("index") if not df.empty else {}

    def _classify(row):
        sym = row["symbol"]
        if sym in curated and pd.notna(curated.get(sym)):
            return curated[sym]
        h_lfc, m_lfc = row.get("human_logfc"), row.get("mouse_logfc")
        h_p, m_p = row.get("human_padj"), row.get("mouse_padj")
        if pd.isna(h_lfc) and pd.isna(m_lfc):
            return "Not tested"
        if pd.isna(h_lfc):
            return "Mouse-Only"
        if pd.isna(m_lfc):
            return "Human-Only"
        h_sig = pd.notna(h_p) and float(h_p) < 0.05
        m_sig = pd.notna(m_p) and float(m_p) < 0.05
        same_dir = (float(h_lfc) * float(m_lfc)) > 0
        if h_sig and m_sig and same_dir:
            return "Human-Mouse Concordant"
        if h_sig and m_sig and not same_dir:
            return "Species Discordant"
        if h_sig and not m_sig:
            return "Human-Only"
        if m_sig and not h_sig:
            return "Mouse-Only"
        return "Not Significant"

    if not atlas_full.empty:
        atlas_full["category_label"] = atlas_full.apply(_classify, axis=1)
        # Attach curated metadata where available
        atlas_full["translatability_score"] = atlas_full["symbol"].map(
            lambda s: (curated_meta.get(s) or {}).get("translatability_score", None)
        )
        atlas_full["n_concordant"] = atlas_full["symbol"].map(
            lambda s: (curated_meta.get(s) or {}).get("n_concordant", None)
        )
        merged = atlas_full
    else:
        merged = df.assign(biotype="protein_coding", human_padj=None, mouse_padj=None)

    # Restrict to genes measurable in BOTH species — "concordance test"
    # only makes sense when both arms have a logFC. Keeping only these
    # rows also stops Human-Only / Mouse-Only from being inflated by
    # genes where the other species just wasn't measured (no ortholog,
    # too-low mouse expression, etc.) which previously gave a misleading
    # 26K Human-Only count.
    merged = merged.dropna(subset=["human_logfc", "mouse_logfc"]).copy()

    has_both = merged
    n_tested = len(has_both)
    n_conserved = int((merged["category_label"] == "Human-Mouse Concordant").sum())
    n_discordant = int((merged["category_label"] == "Species Discordant").sum())
    cons_frac = (n_conserved / n_tested * 100) if n_tested else 0.0

    # Stat cards
    c1, c2, c3, c4 = st.columns(4)
    with c1:
        stat_card(
            "Concordance-tested",
            f"{n_tested:,}",
            "genes with measurable logFC in both human and mouse",
        )
    with c2:
        stat_card(
            "Human-Mouse Concordant",
            f"{n_conserved:,}",
            f"{cons_frac:.1f}% of tested — same direction + sig in both",
        )
    with c3:
        stat_card("Species Discordant", f"{n_discordant:,}",
                  "opposite direction, both sig — flagged for followup")
    with c4:
        # Strongest concordant up by (h+m) sum among Human-Mouse Concordant genes
        cc = has_both[has_both["category_label"] == "Human-Mouse Concordant"]
        if len(cc):
            top_row = cc.loc[(cc["human_logfc"] + cc["mouse_logfc"]).idxmax()]
            stat_card(
                "Strongest concordant ↑",
                str(top_row.get("symbol", "--")),
                f"H {top_row['human_logfc']:+.2f} · M {top_row['mouse_logfc']:+.2f}",
            )
        else:
            stat_card("Strongest concordant ↑", "--")

    # Category count chart — at-a-glance breakdown replacing the old scatter.
    cat_counts = (
        merged["category_label"]
        .value_counts()
        .reindex(CATEGORY_ORDER + ["Not tested"])
        .dropna()
        .astype(int)
    )
    cfig = go.Figure()
    cfig.add_trace(go.Bar(
        x=cat_counts.values, y=cat_counts.index, orientation="h",
        marker=dict(color=[CATEGORY_COLORS.get(c, "#888") for c in cat_counts.index]),
        text=[f"{v:,}" for v in cat_counts.values], textposition="outside",
        hovertemplate="<b>%{y}</b><br>%{x:,} genes<extra></extra>",
    ))
    apply_theme(cfig, margin=dict(l=180, r=50, t=10, b=40))
    cfig.update_layout(
        height=320, xaxis_title="Genes", yaxis=dict(autorange="reversed"),
        showlegend=False,
    )
    st.plotly_chart(cfig, use_container_width=True)
    st.caption(
        "**Note on Human-Only vs Mouse-Only counts.** The human arm "
        "(1,444-sample mega-analysis) has ~3× the statistical power of "
        "the mouse arm (~463 samples across 5 diet models), so many "
        "genes pass padj<0.05 in human but not mouse at the same effect size."
    )

    # -- Searchable concordance table ---------------------------------------
    # A small mouse SVG sits next to the section header — it's a single
    # inline <svg> (no external assets) with a subtle CSS `wiggle`
    # animation defined in theme/custom.css. Pure decoration, signals the
    # "mouse side" of the concordance test visually.
    st.markdown(
        """
        <div class="masld-section-flex">
          <h3 class="masld-section-heading" style="margin:0;">Human-Mouse Concordance Test</h3>
          <svg class="masld-mouse-wiggle" xmlns="http://www.w3.org/2000/svg"
               width="42" height="32" viewBox="0 0 60 44" aria-label="mouse">
            <!-- body -->
            <ellipse cx="30" cy="28" rx="18" ry="12" fill="#c9ccd2"
                     stroke="#6a7280" stroke-width="1.2"/>
            <!-- ears -->
            <circle cx="18" cy="16" r="6" fill="#c9ccd2"
                    stroke="#6a7280" stroke-width="1.2"/>
            <circle cx="42" cy="16" r="6" fill="#c9ccd2"
                    stroke="#6a7280" stroke-width="1.2"/>
            <circle cx="18" cy="16" r="3" fill="#e8a5b7"/>
            <circle cx="42" cy="16" r="3" fill="#e8a5b7"/>
            <!-- eyes -->
            <circle cx="24" cy="26" r="1.4" fill="#1a1a1a"/>
            <circle cx="36" cy="26" r="1.4" fill="#1a1a1a"/>
            <!-- nose -->
            <circle cx="30" cy="32" r="1.6" fill="#e8a5b7"/>
            <!-- whiskers -->
            <path d="M22 33 L12 32 M22 35 L12 36 M38 33 L48 32 M38 35 L48 36"
                  stroke="#6a7280" stroke-width="0.9" fill="none"/>
            <!-- tail -->
            <path d="M48 30 Q58 28 57 38" stroke="#6a7280" stroke-width="1.4"
                  fill="none" stroke-linecap="round"/>
          </svg>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.caption(
        "Search **any gene** measurable in both human and mouse — including "
        "lncRNAs, miRNAs, and MT-* — to see its concordance status. Categories "
        "come from the curated cross-species analysis where available; "
        "otherwise classified on the fly from atlas logFC + padj."
    )

    present_cats = [c for c in CATEGORY_ORDER if c in has_both["category_label"].unique()]

    tc1, tc2, tc3 = st.columns([2, 2.4, 1.2])
    t_query = tc1.text_input(
        "Search symbol", value="",
        placeholder="e.g. THRB, MALAT1, MT-ND1", key="xspec_table_q",
    )
    t_cats = tc2.multiselect(
        "Categories", options=present_cats, default=present_cats,
        key="xspec_table_cats",
    )
    sort_choice = tc3.selectbox(
        "Sort by",
        ["Translatability score",
         "|Human logFC|",
         "|Mouse logFC|",
         "|Human + Mouse logFC|"],
        index=0,
    )

    tview = has_both[has_both["category_label"].isin(t_cats)].copy()
    if t_query.strip():
        q = t_query.strip().upper()
        tview = tview[tview["symbol"].astype(str).str.upper().str.contains(q, na=False)]

    if sort_choice == "Translatability score" and "translatability_score" in tview.columns:
        tview = tview.sort_values("translatability_score", ascending=False)
    elif sort_choice == "|Human logFC|":
        tview = tview.assign(_k=tview["human_logfc"].abs()).sort_values("_k", ascending=False).drop(columns=["_k"])
    elif sort_choice == "|Mouse logFC|":
        tview = tview.assign(_k=tview["mouse_logfc"].abs()).sort_values("_k", ascending=False).drop(columns=["_k"])
    else:
        tview = tview.assign(
            _k=(tview["human_logfc"] + tview["mouse_logfc"]).abs()
        ).sort_values("_k", ascending=False).drop(columns=["_k"])

    if not tview.empty:
        tview["Gene page"] = tview["symbol"].apply(lambda s: f"/Gene?symbol={s}")
        display_cols = [c for c in [
            "symbol", "biotype", "category_label", "human_logfc", "mouse_logfc",
            "translatability_score", "n_concordant", "n_diets_sig", "Gene page",
        ] if c in tview.columns]
        st.dataframe(
            tview.head(1000)[display_cols],
            hide_index=True, use_container_width=True, height=460,
            column_config={
                "symbol": st.column_config.TextColumn("Symbol", width="small"),
                "biotype": st.column_config.TextColumn("Biotype", width="small"),
                "category_label": st.column_config.TextColumn("Category"),
                "human_logfc": st.column_config.NumberColumn("Human logFC", format="%+.3f"),
                "mouse_logfc": st.column_config.NumberColumn("Mouse logFC", format="%+.3f"),
                "translatability_score": st.column_config.NumberColumn(
                    "Translatability", format="%.2f",
                    help="Composite score rewarding concordant direction and "
                         "significance across both species; higher = more translatable.",
                ),
                "n_concordant": st.column_config.NumberColumn("# concordant", format="%d"),
                "n_diets_sig": st.column_config.NumberColumn(
                    "Mouse diets sig.", format="%d",
                    help="Number of mouse diet cohorts in which the gene is DE.",
                ),
                "Gene page": st.column_config.LinkColumn("Gene page", display_text="open"),
            },
        )
        st.caption(f"Top 1,000 of {len(tview):,} matching genes shown.")
    else:
        st.info("No genes match the current filters.")


# ===========================================================================
# Tab 2 — Proteomics
# ===========================================================================
with tabs[1]:
    try:
        summary = load_proteomics_summary()
    except FileNotFoundError as exc:
        st.error(str(exc))
        st.stop()

    st.caption(
        "Tissue + plasma proteomics cross-validate transcriptomic signals "
        "at the protein level. A matching direction + significance at the "
        "protein level substantially raises confidence in a target."
    )

    # Stat cards — renamed for clarity
    mc1, mc2, mc3, mc4, mc5 = st.columns(5)
    with mc1:
        stat_card(
            "Proteins tested",
            f"{summary.get('n_proteins_tested', 0):,}",
            "measurable in at least one proteomics dataset",
        )
    with mc2:
        stat_card(
            "Significant proteins",
            f"{summary.get('n_proteins_significant', 0):,}",
            "padj < 0.05 in tissue DIA-MS",
        )
    with mc3:
        stat_card(
            "mRNA-protein ρ (overall)",
            f"{summary.get('overall_rho', 0):.3f}",
            "Spearman correlation, all proteins",
        )
    with mc4:
        stat_card(
            "ρ in Human-Mouse Concordant",
            f"{summary.get('conserved_rho', 0):.3f}",
            "subset with strongest cross-species signal",
        )
    with mc5:
        stat_card(
            "Validated drug targets",
            f"{summary.get('n_validated_targets', 0):,}",
            "clinical MASLD targets confirmed at protein level",
        )

    # -- Unified, searchable protein table ----------------------------------
    st.markdown("### Protein-level atlas")
    st.caption(
        "Unified top-regulated proteins from tissue DIA-MS + plasma Olink. "
        "Filter by direction or search a specific gene."
    )

    up_rows = summary.get("top_upregulated", []) or []
    dn_rows = summary.get("top_downregulated", []) or []
    for r in up_rows:
        r["direction"] = "Up"
    for r in dn_rows:
        r["direction"] = "Down"
    all_rows = up_rows + dn_rows

    if not all_rows:
        st.info("No protein table rows available in proteomics_summary.json.")
    else:
        pdf = pd.DataFrame(all_rows)

        pc1, pc2, pc3 = st.columns([2, 1.2, 1.2])
        p_query = pc1.text_input(
            "Search protein / gene", value="",
            placeholder="e.g. COL1A1, THRB", key="prot_q",
        )
        p_dir = pc2.multiselect(
            "Direction", options=["Up", "Down"], default=["Up", "Down"],
        )
        p_sort = pc3.selectbox(
            "Sort by", ["|logFC|", "tstat", "padj (ascending)"], index=0,
        )

        pview = pdf[pdf["direction"].isin(p_dir)].copy()
        if p_query.strip():
            q = p_query.strip().upper()
            if "gene" in pview.columns:
                pview = pview[pview["gene"].astype(str).str.upper().str.contains(q, na=False)]

        if p_sort == "|logFC|" and "logfc" in pview.columns:
            pview = pview.assign(_k=pview["logfc"].abs()).sort_values("_k", ascending=False).drop(columns=["_k"])
        elif p_sort == "tstat" and "tstat" in pview.columns:
            pview = pview.assign(_k=pview["tstat"].abs()).sort_values("_k", ascending=False).drop(columns=["_k"])
        elif p_sort == "padj (ascending)" and "padj" in pview.columns:
            pview = pview.sort_values("padj", ascending=True)

        if "gene" in pview.columns:
            pview["Gene page"] = pview["gene"].apply(lambda s: f"/Gene?symbol={s}")

        display_cols = [c for c in [
            "gene", "direction", "logfc", "tstat", "padj", "dataset", "Gene page",
        ] if c in pview.columns]
        st.dataframe(
            pview[display_cols], hide_index=True, use_container_width=True,
            height=460,
            column_config={
                "gene":   st.column_config.TextColumn("Gene", width="small"),
                "direction": st.column_config.TextColumn("Direction", width="small"),
                "logfc":  st.column_config.NumberColumn("logFC", format="%+.3f"),
                "tstat":  st.column_config.NumberColumn("tstat", format="%.2f"),
                "padj":   st.column_config.NumberColumn("padj", format="%.2e"),
                "Gene page": st.column_config.LinkColumn("Gene page", display_text="open"),
            },
        )
        st.caption(f"{len(pview):,} proteins match filters.")

    st.divider()
    st.markdown(
        "**Provenance.** Tissue proteomics: PXD052937 (DIA-MS, 72 samples) + "
        "GSE276114_fibrosis. Plasma proteomics: Olink Explore 3072."
    )
