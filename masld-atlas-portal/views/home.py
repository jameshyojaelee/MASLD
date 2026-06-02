"""MASLD Multi-Evidence Atlas -- Home view.

Mirrors the Next.js landing page: compact hero, 4 quick-link cards, and a
3-column featured-gene gallery. Downloads / modality reference are collapsed
into the bottom expander so the landing area stays lean.
"""
from __future__ import annotations

import streamlit as st

from components.layout import (
    badge,
    eyebrow,
    inject_css,
    sidebar_header,
)
from components.modality_radar import make_fingerprint_svg
from data.loaders import load_featured_genes

inject_css()
sidebar_header()


# ---------------------------------------------------------------------------
# Hero (centered, Next.js pattern)
# ---------------------------------------------------------------------------

st.markdown(
    """
    <div class="masld-hero">
      <h1 class="masld-hero-title">MASLD Atlas</h1>
      <p class="masld-hero-subtitle">
        Multi-modal atlas for metabolic dysfunction-associated steatotic
        liver disease.
      </p>
    </div>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# Stat cards -- 4 headline stats (the most important ones)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Quick links -- 4 primary routes
# ---------------------------------------------------------------------------

eyebrow("Explore")

QUICK_LINKS = [
    (
        "Atlas",
        "Gene-by-gene multi-modal evidence inventory across the 34K-gene atlas.",
        "Atlas",
    ),
    (
        "Network",
        "Force-directed gene network linking genes by protein, co-expression, genetic, and ligand-receptor edges.",
        "Network",
    ),
    (
        "Prioritize",
        "Convergence matrix, Bayesian ranking, and a gene–drug–TF target map.",
        "Prioritize",
    ),
    (
        "Validation",
        "Cross-species concordance and proteomics validation of atlas targets.",
        "Validation",
    ),
]

qcols = st.columns(4)
for col, (title, desc, path) in zip(qcols, QUICK_LINKS):
    with col:
        st.markdown(
            f"""
            <a class="masld-quicklink" href="/{path}" target="_self">
              <div class="masld-quicklink-title">{title} <span>&rarr;</span></div>
              <p class="masld-quicklink-desc">{desc}</p>
            </a>
            """,
            unsafe_allow_html=True,
        )


# ---------------------------------------------------------------------------
# Featured genes -- 3-column responsive grid
# ---------------------------------------------------------------------------

st.markdown(
    '<h2 class="masld-section-heading">Featured genes</h2>'
    '<p class="masld-section-heading-sub">'
    "Curated targets spanning approved drug candidates, GWAS loci, and progression drivers.</p>",
    unsafe_allow_html=True,
)

featured = load_featured_genes()

# Curated subset — keep only six hero genes so the grid stays compact.
# Order controls display order on the page.
_FEATURED_ALLOW = ["THRB", "PNPLA3", "RORA", "HSD17B13", "PPARA", "DGAT2"]
_by_sym = {g.get("symbol"): g for g in featured}
featured = [_by_sym[s] for s in _FEATURED_ALLOW if s in _by_sym]

GRID_COLS = 3
for row_start in range(0, len(featured), GRID_COLS):
    cols = st.columns(GRID_COLS, gap="small")
    for offset, gene in enumerate(featured[row_start : row_start + GRID_COLS]):
        with cols[offset]:
            symbol = gene.get("symbol", "?")
            category = gene.get("category", "")
            tagline = gene.get("tagline", "")
            lfc = gene.get("dream_logfc")
            padj = gene.get("dream_padj")
            evidence = gene.get("evidence", {})

            lfc_txt = f"{lfc:+.3f}" if isinstance(lfc, (int, float)) else "--"
            padj_txt = (
                f"{padj:.1e}"
                if isinstance(padj, (int, float)) and padj > 0
                else ("0" if padj == 0 else "--")
            )
            fp_svg = make_fingerprint_svg(evidence, size=60)
            st.markdown(
                f"""
                <a class="masld-gene-card" href="/Gene?symbol={symbol}" target="_self">
                  <div class="masld-gene-card__row">
                    <div class="masld-gene-card__body">
                      <div class="masld-gene-card__header">
                        <span class="masld-gene-card__symbol">{symbol}</span>
                        {badge(category, "neutral")}
                      </div>
                      <p class="masld-gene-card__tagline">{tagline}</p>
                      <div class="masld-gene-card__metric">
                        logFC <b>{lfc_txt}</b> &nbsp;|&nbsp; padj <b>{padj_txt}</b>
                      </div>
                    </div>
                    <div class="masld-gene-card__fp">{fp_svg}</div>
                  </div>
                </a>
                """,
                unsafe_allow_html=True,
            )


# ---------------------------------------------------------------------------
# More (downloads + modality reference), collapsed by default
# ---------------------------------------------------------------------------

st.markdown('<div class="masld-spacer-lg"></div>', unsafe_allow_html=True)

with st.expander("Downloads & modality reference", expanded=False):
    st.markdown("##### Data resources")
    st.markdown(
        "- **Multi-Evidence Atlas** (Parquet): 33,943 genes × 197 columns — `atlas.parquet`\n"
        "- **Gene index** (JSON): compact search index — `gene_index.json`\n"
        "- **Per-gene profiles** (JSON, 34K files): full evidence cards — `genes/{SYMBOL}.json`\n"
        "- **Convergence matrix** (JSON): modality-level scores — `convergence_matrix.json`\n"
        "- **Bayesian ranking** (JSON): composite posterior per gene — `bayesian_ranking.json`\n"
        "- **Proteomics summary** (JSON): mRNA-protein concordance — `proteomics_summary.json`\n"
    )
    st.markdown("##### Modalities")
    st.markdown(
        "- **S1 Human Bulk** — Integrated mega-analysis, 10 cohorts, 5,484 DEGs (LOO mean recovery 88.1%).\n"
        "- **S2 Mouse Bulk** — metafor random-effects across 5 diet models (463 samples).\n"
        "- **S3 Genetic Causal** — TWAS + SuSiE/ABF COLOC across 24 EUR GWAS + FinnGen + BBJ.\n"
        "- **S4 Essentiality** — DepMap Chronos scores.\n"
        "- **S5 Epigenomic** — scATAC + SCENIC+ regulons + chromVAR + motif disruption (380 GWAS variants).\n"
        "- **S6 Spatial** — Visium GSE192741, 150 SVGs, pericentral/periportal zonation.\n"
        "- **S7 Single-cell** — pseudobulk DE + LIANA across 5 cell types.\n"
        "- **S8 Proteomics** — tissue DIA-MS + plasma Olink."
    )
    st.caption(
        "Data served from masld-atlas-v2/public/data (local dev) or an LFS-tracked "
        "data subtree on Hugging Face Spaces. Set ``MASLD_DATA_DIR`` to override."
    )
