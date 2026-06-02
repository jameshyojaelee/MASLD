"""Gene detail view -- per-gene evidence card.

Reorganised into four ``section_card`` blocks to mirror the Next.js
``src/app/gene/[symbol]/gene-evidence-card.tsx`` layout:

1. **Overview** — modality radar + Integrated-DEG badge + biotype + ext links.
2. **Expression** — 2x2 trajectory grid + per-cohort forest plot w/ contrast.
3. **Causal + Therapeutic** — COLOC table + druggability + drugs + proteomics.
4. **Cell-type + Sex** — pseudobulk DE + sex/NMF subtype.

The page is registered via ``st.Page(url_path="Gene")`` so the URL stays
``/Gene?symbol=THRB``. The sidebar entry is hidden in CSS.
"""
from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from components.layout import (
    badge,
    inject_css,
    page_header,
    section_card,
    sidebar_header,
)
from components.modality_radar import make_radar
from data.loaders import load_gene_graph, load_gene_index, load_gene_profile
from theme.plotly_theme import DEG_DOWN, DEG_UP, PLOTLY_THEME, apply_theme

inject_css()
sidebar_header()


EXCLUDED_COHORTS = {"PRJNA512027"}


# ---------------------------------------------------------------------------
# Resolve symbol from URL or picker
# ---------------------------------------------------------------------------

symbol = st.query_params.get("symbol", None)
if isinstance(symbol, list):
    symbol = symbol[0] if symbol else None

if not symbol:
    try:
        symbols = sorted(load_gene_index()["symbol"].dropna().unique().tolist())
    except Exception:
        symbols = []
    symbol = st.selectbox(
        "Pick a gene",
        options=symbols or ["THRB"],
        index=(symbols.index("THRB") if "THRB" in symbols else 0) if symbols else 0,
    )
    if symbol:
        st.query_params["symbol"] = symbol

if not symbol:
    st.stop()

profile = load_gene_profile(symbol)
if profile is None:
    from data import paths
    expected = paths.gene_profile_path(symbol)
    genes_dir = expected.parent
    if not genes_dir.exists():
        # `genes/` subdir missing entirely — data hasn't been fetched yet, or
        # the gene-profile payload is missing from the data repo.
        st.info(
            f"**Gene data still initializing.** The per-gene profiles "
            f"(`genes/{symbol}.json`) aren't available yet at `{genes_dir}`. "
            "If this is a fresh deploy, wait for the data bootstrap to finish "
            "(first boot can take a few minutes) and reload. If the Space has "
            "been running for a while, the data repo may be missing the "
            "`genes/` tree — check the dataset repo."
        )
    else:
        # Directory exists but this specific symbol isn't in it — most likely
        # the user typed something that doesn't map to a profile we emitted.
        st.warning(
            f"No profile found for `{symbol}`. The Atlas covers 33,943 genes "
            "(including lncRNAs and pseudogenes); confirm the symbol is a "
            "HGNC-style human gene symbol."
        )
    st.stop()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _evidence_from_profile(p: dict) -> dict[str, float]:
    e: dict[str, float] = {k: 0.0 for k in (
        "s1_human", "s2_genetic", "s3_essential", "s4_epigenomic",
        "s5_spatial", "s6_singlecell", "s7_mouse", "s8_proteomics",
    )}
    expr_ = p.get("expression") or {}
    if expr_.get("dream_padj") is not None:
        padj_ = float(expr_["dream_padj"])
        abs_lfc = abs(float(expr_.get("dream_logfc", 0) or 0))
        if padj_ <= 0.05 and abs_lfc >= 0.2:
            e["s1_human"] = min(1.0, 0.5 + abs_lfc)
        elif padj_ <= 0.05:
            e["s1_human"] = 0.3
        else:
            e["s1_human"] = min(0.15, abs_lfc * 0.5)
    causal_ = p.get("causal") or {}
    if causal_:
        pp4_ = float(causal_.get("coloc_pp4_max", 0) or 0)
        e["s2_genetic"] = min(1.0, pp4_)
        if causal_.get("twas_pval") is not None and float(causal_["twas_pval"]) < 0.05:
            e["s2_genetic"] = min(1.0, e["s2_genetic"] + 0.2)
    ct_ = (p.get("celltype") or {}).get("pseudobulk_de") or []
    n_sig = sum(1 for r in ct_ if float(r.get("padj", 1)) < 0.05)
    if n_sig:
        e["s6_singlecell"] = min(1.0, n_sig * 0.25)
    prot_ = p.get("proteomics") or {}
    if prot_.get("logfc") is not None:
        lfc_ = float(prot_["logfc"])
        padj_ = float(prot_.get("padj", 1) or 1)
        if padj_ < 0.05:
            e["s8_proteomics"] = min(1.0, abs(lfc_) / 2)
        else:
            e["s8_proteomics"] = min(0.15, abs(lfc_) * 0.3)
    return e


def _stage_bar(rows: list[dict], label_key: str, title: str) -> None:
    if not rows:
        st.info(f"{title}: no data.")
        return
    dfp = pd.DataFrame(rows)
    dfp["label"] = dfp[label_key]
    dfp["is_sig"] = dfp["padj"] < 0.05
    dfp["color"] = dfp.apply(
        lambda r: (DEG_UP if r["is_sig"] else "#f3b4b4")
        if r["logfc"] >= 0
        else (DEG_DOWN if r["is_sig"] else "#b0c9f1"),
        axis=1,
    )
    fig = go.Figure(go.Bar(
        x=dfp["label"], y=dfp["logfc"],
        marker_color=dfp["color"],
        hovertemplate="<b>%{x}</b><br>logFC: %{y:+.3f}<extra></extra>",
    ))
    apply_theme(fig, margin=dict(l=30, r=10, t=30, b=30))
    fig.update_layout(
        height=260,
        title=dict(text=title, x=0.02, font=dict(size=13)),
        yaxis_title="logFC",
    )
    fig.add_hline(y=0, line_color="#9aa3ad", line_width=1)
    st.plotly_chart(fig, use_container_width=True)


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

expr = profile.get("expression") or {}
lfc_top = expr.get("dream_logfc")
padj_top = expr.get("dream_padj")
is_deg = (
    isinstance(padj_top, (int, float)) and padj_top < 0.05
    and isinstance(lfc_top, (int, float)) and abs(lfc_top) > 0.3
)

page_header(
    symbol,
    f"Ensembl: {profile.get('ensembl_id', '')} &middot; "
    f"biotype: {profile.get('biotype', '--').replace('_', ' ')}",
)


# ---------------------------------------------------------------------------
# Card 1 -- Overview
# ---------------------------------------------------------------------------

with section_card("Overview"):
    radar_col, meta_col = st.columns([1, 2])
    with radar_col:
        # Use gene_index evidence scores so the radar here matches the
        # Atlas radar wall + Home featured-gene fingerprints exactly.
        # (Previously this recomputed from gene_profile and diverged.)
        _idx_row = None
        try:
            _idx = load_gene_index()
            _hit = _idx[_idx["symbol"] == symbol]
            if len(_hit):
                _idx_row = _hit.iloc[0]
        except Exception:
            _idx_row = None

        if _idx_row is not None:
            evidence = {
                k: float(_idx_row.get(k, 0.0) or 0.0)
                for k in ("s1_human", "s2_mouse", "s3_genetic", "s4_essential",
                          "s5_epigenomic", "s6_spatial", "s7_singlecell",
                          "s8_proteomics")
            }
        else:
            evidence = _evidence_from_profile(profile)

        st.plotly_chart(
            make_radar(evidence, size="large"),
            use_container_width=False,
            config={"displayModeBar": False},
        )
    with meta_col:
        badge_html = badge(
            "Integrated DEG" if is_deg else "Not a DEG (atlas threshold)",
            "success" if is_deg else "neutral",
        )
        biotype_html = badge(profile.get("biotype", "--").replace("_", " "), "primary")
        st.markdown(
            f"<div class='masld-badge-row'>{badge_html} {biotype_html}</div>",
            unsafe_allow_html=True,
        )
        mrow = st.columns(3)
        tstat_top = expr.get("dream_tstat")
        mrow[0].metric("Integrated logFC", f"{lfc_top:+.3f}" if isinstance(lfc_top, (int, float)) else "--")
        mrow[1].metric(
            "padj",
            f"{padj_top:.2e}" if isinstance(padj_top, (int, float)) and padj_top > 0
            else ("0" if padj_top == 0 else "--"),
        )
        mrow[2].metric("tstat", f"{tstat_top:.2f}" if isinstance(tstat_top, (int, float)) else "--")

        ext = profile.get("external") or {}
        if ext:
            st.markdown(
                "<div class='masld-badge-row' style='margin-top:10px;'>"
                + "".join(
                    f"<a class='masld-extlink' href='{v}' target='_blank'>"
                    f"{k.capitalize()} &rarr;</a>"
                    for k, v in ext.items() if v
                )
                + "</div>",
                unsafe_allow_html=True,
            )

        # Phase-3 collapsed-layer badges: sparse layers (D-REG, D-ceRNA...)
        # that the network view hides appear here as node-level badges
        # when the gene participates in them. Gene graph JSONs carry a
        # `badges` dict (set in 294_export_portal_v2.py).
        graph = load_gene_graph(symbol) or {}
        badge_map = (graph.get("badges") or {}) if isinstance(graph, dict) else {}
        if badge_map:
            _labels = {
                "D-REG":    "TF regulon participant",
                "D-ceRNA":  "ceRNA participant",
                "D-XS":     "Cross-species conserved",
                "D-LR":     "Ligand / receptor",
                "D-COLOC":  "GWAS-colocalized",
                "D-STABLE": "Stable co-expression hub",
            }
            parts = []
            for k, n in sorted(badge_map.items()):
                label = _labels.get(k, k) + (f" ({int(n)})" if isinstance(n, (int, float)) and n else "")
                parts.append(badge(label, "info"))
            st.markdown(
                "<div class='masld-badge-row' style='margin-top:10px;'>"
                + "".join(parts) + "</div>",
                unsafe_allow_html=True,
            )


# ---------------------------------------------------------------------------
# Card 2 -- Expression
# ---------------------------------------------------------------------------

with section_card("Expression evidence"):
    fib_variant = st.radio("Fibrosis reference", options=["vs F0", "vs Healthy"], horizontal=True, key="fib_variant")
    nas_variant = st.radio("NAS reference", options=["vs NAS0", "vs Healthy"], horizontal=True, key="nas_variant")

    fib_rows = (
        expr.get("stage_trajectory_vs_control")
        if fib_variant == "vs F0"
        else expr.get("stage_trajectory_vs_healthy")
    ) or []
    fib_rest = expr.get("stage_trajectory_vs_rest") or []
    nas_rows = (
        expr.get("nas_trajectory_vs_nas0")
        if nas_variant == "vs NAS0"
        else expr.get("nas_trajectory_vs_healthy")
    ) or []
    nas_rest = expr.get("nas_trajectory") or []

    r1 = st.columns(2)
    with r1[0]:
        _stage_bar(fib_rows, "stage", f"Fibrosis stage ({fib_variant})")
    with r1[1]:
        _stage_bar(fib_rest, "stage", "Fibrosis one-vs-rest")
    r2 = st.columns(2)
    with r2[0]:
        _stage_bar(nas_rows, "nas", f"NAS ({nas_variant})")
    with r2[1]:
        _stage_bar(nas_rest, "nas", "NAS one-vs-rest")

    CONTRAST_LABELS = {
        "disease_vs_control": "Disease vs Healthy control",
        "MASH_vs_control": "MASH vs Healthy control",
        "MASL_vs_control": "MASL vs Healthy control",
        "MASH_vs_MASL": "MASH vs MASL",
    }
    per_cohort = expr.get("per_cohort") or {}
    available_contrasts = [k for k in CONTRAST_LABELS if per_cohort.get(k)]
    if not available_contrasts and expr.get("per_cohort_dream"):
        per_cohort = {"disease_vs_control": expr["per_cohort_dream"]}
        available_contrasts = ["disease_vs_control"]

    if available_contrasts:
        st.markdown("#### Per-cohort forest plot")
        contrast = st.selectbox(
            "Contrast", options=available_contrasts,
            format_func=lambda k: CONTRAST_LABELS[k],
        )
        rows = per_cohort.get(contrast, [])
        rows = [r for r in rows if r.get("dataset") not in EXCLUDED_COHORTS]
        if rows:
            dff = pd.DataFrame(rows).sort_values("logfc")
            dff["is_sig"] = dff["pval"] < 0.05 if "pval" in dff.columns else False
            dff["color"] = dff.apply(
                lambda r: (DEG_UP if r["is_sig"] else "#f3b4b4")
                if r["logfc"] >= 0
                else (DEG_DOWN if r["is_sig"] else "#b0c9f1"),
                axis=1,
            )
            fig = go.Figure(go.Bar(
                y=dff["dataset"], x=dff["logfc"], orientation="h",
                marker_color=dff["color"],
                hovertemplate="<b>%{y}</b><br>logFC %{x:+.3f}<extra></extra>",
            ))
            apply_theme(fig, margin=dict(l=10, r=10, t=10, b=30))
            fig.update_layout(
                height=max(200, 22 * len(dff) + 60),
                xaxis_title="logFC",
            )
            fig.add_vline(x=0, line_color="#9aa3ad", line_width=1)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No cohort data for this contrast.")


# ---------------------------------------------------------------------------
# Card 3 -- Causal + Therapeutic
# ---------------------------------------------------------------------------

with section_card("Causal + Therapeutic"):
    causal = profile.get("causal") or {}
    therapeutic = profile.get("therapeutic") or {}
    prot = profile.get("proteomics") or {}

    # 2-column split: left = COLOC (genetic), right = druggability + drugs
    left, right = st.columns(2)

    with left:
        st.markdown("#### Genetic causal")
        if not causal:
            st.info("No genetic causal data.")
        else:
            g1, g2 = st.columns(2)
            g1.metric("COLOC PP4 max", f"{causal.get('coloc_pp4_max', 0):.3f}")
            g2.metric("N GWAS PP4>=0.5", causal.get("coloc_n_gwas_05", 0))
            st.caption(f"Best GWAS: **{causal.get('coloc_best_gwas', '--')}**")

            by_gwas = causal.get("coloc_by_gwas") or {}
            if by_gwas:
                cdf = (
                    pd.DataFrame({"gwas": list(by_gwas.keys()), "pp4": list(by_gwas.values())})
                    .sort_values("pp4", ascending=False)
                )
                st.dataframe(
                    cdf, hide_index=True, use_container_width=True,
                    column_config={"pp4": st.column_config.NumberColumn("PP.H4", format="%.4f")},
                )

            twas_p = causal.get("twas_pval")
            if twas_p is not None:
                z_val = causal.get("twas_z")
                if isinstance(z_val, (int, float)):
                    st.markdown(f"**TWAS** &middot; z = {z_val:.2f}  &middot; p = {twas_p:.2e}")
                else:
                    st.markdown(f"**TWAS** &middot; p = {twas_p}")

    with right:
        st.markdown("#### Therapeutic")
        if not therapeutic:
            st.info("No therapeutic data.")
        else:
            druggable = therapeutic.get("dgidb_druggable")
            st.markdown(
                f"Druggable (DGIdb): **{'Yes' if druggable else 'No' if druggable is not None else '--'}**"
            )
            drugs = therapeutic.get("drugs") or []
            if drugs:
                st.markdown("**Approved/investigational drugs:**")
                st.dataframe(pd.DataFrame(drugs), hide_index=True, use_container_width=True)
            lincs = therapeutic.get("lincs_compounds") or []
            if lincs:
                st.markdown("**LINCS reversal compounds:**")
                ldf = pd.DataFrame(lincs)
                st.dataframe(
                    ldf, hide_index=True, use_container_width=True,
                    column_config={"score": st.column_config.NumberColumn("score", format="%.3f")},
                )
            if therapeutic.get("progression_class") or therapeutic.get("progression_drug"):
                st.markdown(
                    f"Progression class: **{therapeutic.get('progression_class', '--')}**  "
                    f"&middot;  drug: **{therapeutic.get('progression_drug', '--')}**"
                )

    st.markdown("#### Proteomics")
    if not prot:
        st.info("No proteomics data.")
    else:
        p1, p2, p3, p4 = st.columns(4)
        p1.metric("logFC", f"{prot['logfc']:+.3f}" if isinstance(prot.get("logfc"), (int, float)) else "--")
        p2.metric("padj", f"{prot['padj']:.2e}" if isinstance(prot.get("padj"), (int, float)) else "--")
        p3.metric("tstat", f"{prot['tstat']:.2f}" if isinstance(prot.get("tstat"), (int, float)) else "--")
        p4.metric("dataset", prot.get("dataset", "--"))


# ---------------------------------------------------------------------------
# Card 4 -- Cell-type + Sex
# ---------------------------------------------------------------------------

with section_card("Cell-type + Sex"):
    # 2-column split: left = cell-type pseudobulk, right = sex/NMF
    ct_col, sex_col = st.columns(2)

    with ct_col:
        ct = profile.get("celltype") or {}
        pb = ct.get("pseudobulk_de") or []
        st.markdown("#### Cell-type pseudobulk DE")
        if pb:
            pdf = pd.DataFrame(pb).sort_values("logfc")
            st.dataframe(
                pdf, hide_index=True, use_container_width=True,
                column_config={
                    "logfc": st.column_config.NumberColumn("logFC", format="%+.3f"),
                    "padj": st.column_config.NumberColumn("padj", format="%.2e"),
                },
            )
        else:
            st.info("No pseudobulk DE data.")

    with sex_col:
        st.markdown("#### Sex & NMF subtype")
        sx = profile.get("sex_subtype") or {}
        if sx:
            s1, s2 = st.columns(2)
            s1.metric("sex_class", sx.get("sex_class", "--"))
            s2.metric("NMF subtype", sx.get("nmf_subtype", "--"))
            s3, s4 = st.columns(2)
            s3.metric(
                "logFC F",
                f"{sx['logfc_female']:+.3f}"
                if isinstance(sx.get("logfc_female"), (int, float)) else "--",
            )
            s4.metric(
                "logFC M",
                f"{sx['logfc_male']:+.3f}"
                if isinstance(sx.get("logfc_male"), (int, float)) else "--",
            )
        else:
            st.info("No sex/NMF data.")


# Raw JSON for debugging
with st.expander("Raw profile JSON", expanded=False):
    st.json(profile)
