"""Plotly radar chart for the 8-modality evidence profile.

The home page uses :func:`make_radar` at a compact 250x250 size; detail pages
can pass ``size='large'`` for a 400x400 render. ``normalize_evidence`` maps the
per-modality keys used by ``gene_index.json`` (``s1_human``, ``s2_mouse`` …)
into the display order defined by :data:`MODALITY_DISPLAY`.

Modality ordering
-----------------
The portal uses a single M1–M8 display order everywhere (radar, convergence
matrix, atlas gene table). Order is bio-evidence first, then cross-species
/ essentiality derived layers last, per user feedback (2026-04-23):

    M1  Human bulk DE
    M2  Genetic (GWAS+eQTL)
    M3  Epigenomic (ATAC / SCENIC+)
    M4  Spatial (Visium)
    M5  Single-cell (pseudobulk + LIANA)
    M6  Proteomics
    M7  Mouse bulk DE              (cross-species — second last)
    M8  Essentiality (DepMap)      (non-MASLD-specific — last)
"""
from __future__ import annotations

from typing import Mapping

import plotly.graph_objects as go


# Single source of truth for modality display order + labels across the portal.
# Each entry: (M-code, long label, gene_index JSON key, legacy s*-style key).
# gene_index.json uses s1_human, s2_mouse, s3_genetic ... per old schema;
# both the new Mn and the legacy sn keys are preserved for backward compat.
MODALITY_DISPLAY: tuple[tuple[str, str, str, str], ...] = (
    ("M1", "Human bulk DE",                  "s1_human",       "s1_human"),
    ("M2", "Genetic (GWAS+eQTL)",            "s3_genetic",     "s2_genetic"),
    ("M3", "Epigenomic (ATAC/SCENIC+)",      "s5_epigenomic",  "s4_epigenomic"),
    ("M4", "Spatial (Visium)",               "s6_spatial",     "s5_spatial"),
    ("M5", "Single-cell (sc + LIANA)",       "s7_singlecell",  "s6_singlecell"),
    ("M6", "Proteomics",                     "s8_proteomics",  "s8_proteomics"),
    ("M7", "Mouse bulk DE",                  "s2_mouse",       "s7_mouse"),
    ("M8", "Essentiality (DepMap)",          "s4_essential",   "s3_essential"),
)

# Radar uses the legacy 4th-element keys (for the existing normalize_evidence
# pipeline) — order here controls the angular layout.
LEGACY_ORDER: tuple[str, ...] = tuple(m[3] for m in MODALITY_DISPLAY)
LEGACY_LABELS: dict[str, str] = {
    m[3]: f"{m[0]} {m[1].split(' (')[0]}"  # "M1 Human bulk DE"
    for m in MODALITY_DISPLAY
}


def normalize_evidence(evidence: Mapping[str, float] | None) -> dict[str, float]:
    """Map gene_index's M1-M7 aligned keys -> legacy 8-modality scheme."""
    raw = dict(evidence or {})
    return {
        "s1_human": float(raw.get("s1_human", 0.0) or 0.0),
        "s2_genetic": float(raw.get("s3_genetic", raw.get("s2_genetic", 0.0)) or 0.0),
        "s3_essential": float(raw.get("s4_essential", raw.get("s3_essential", 0.0)) or 0.0),
        "s4_epigenomic": float(raw.get("s5_epigenomic", raw.get("s4_epigenomic", 0.0)) or 0.0),
        "s5_spatial": float(raw.get("s6_spatial", raw.get("s5_spatial", 0.0)) or 0.0),
        "s6_singlecell": float(raw.get("s7_singlecell", raw.get("s6_singlecell", 0.0)) or 0.0),
        "s7_mouse": float(raw.get("s2_mouse", raw.get("s7_mouse", 0.0)) or 0.0),
        "s8_proteomics": float(raw.get("s8_proteomics", 0.0) or 0.0),
    }


def make_radar(
    evidence: Mapping[str, float] | None,
    *,
    color: str = "#0d62b8",
    size: str = "small",
    name: str | None = None,
) -> go.Figure:
    """Return a Plotly Scatterpolar radar chart.

    Evidence values are expected to be in ``[0, 1]``. Accepts either the
    legacy (s1..s8) keys or the M1-M7 aligned keys -- ``normalize_evidence``
    handles both.
    """
    normed = normalize_evidence(evidence)
    theta = [LEGACY_LABELS[k] for k in LEGACY_ORDER]
    r = [normed[k] for k in LEGACY_ORDER]

    # Close the polygon
    theta_closed = list(theta) + [theta[0]]
    r_closed = list(r) + [r[0]]

    fig = go.Figure()
    fig.add_trace(
        go.Scatterpolar(
            r=r_closed,
            theta=theta_closed,
            fill="toself",
            fillcolor=_hex_to_rgba(color, 0.18),
            line=dict(color=color, width=1.5),
            name=name or "Evidence",
            hovertemplate="<b>%{theta}</b>: %{r:.2f}<extra></extra>",
        )
    )

    # Labels like "M1 Human bulk DE" are ~16 chars; the 3/9-o'clock ticks
    # need generous horizontal margin so Plotly doesn't clip them. Small
    # size stays compact (used in tight grids); large expands horizontally
    # and shrinks the polar domain so the polygon leaves room for text.
    dims_map = {"small": (260, 250), "large": (520, 400)}  # (width, height)
    width, height = dims_map.get(size, (260, 250))
    tick_font = 9 if size == "small" else 11
    # Horizontal margin + polar domain shrink prevents left/right clipping.
    lr_margin = 50 if size == "small" else 90
    polar_x = [0.18, 0.82] if size == "large" else [0.12, 0.88]

    fig.update_layout(
        polar=dict(
            domain=dict(x=polar_x, y=[0.05, 0.95]),
            radialaxis=dict(
                visible=True,
                range=[0, 1],
                tickfont=dict(size=tick_font - 1),
                tickvals=[0.25, 0.5, 0.75, 1.0],
                showline=False,
                gridcolor="rgba(0,0,0,0.08)",
            ),
            angularaxis=dict(
                tickfont=dict(size=tick_font),
                gridcolor="rgba(0,0,0,0.08)",
                rotation=90,          # M1 at 12 o'clock (top)
                direction="clockwise",  # M2 → M3 → ... clockwise
            ),
            bgcolor="rgba(0,0,0,0)",
        ),
        showlegend=False,
        margin=dict(l=lr_margin, r=lr_margin, t=20, b=20),
        width=width,
        height=height,
        paper_bgcolor="rgba(0,0,0,0)",
    )
    # Inherit the shared Plotly typography (Inter, muted axes) so the radar
    # matches every other chart in the portal. Local polar/margin settings
    # above take precedence over anything apply_theme would set.
    try:
        from theme.plotly_theme import apply_theme  # local import avoids cycle
        apply_theme(fig, margin=dict(l=lr_margin, r=lr_margin, t=20, b=20))
    except Exception:
        pass
    return fig


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# ---------------------------------------------------------------------------
# Lightweight SVG "evidence fingerprint" — no Plotly, renders instantly.
# Used on the Home featured-gene grid where we had 9 simultaneous Plotly
# radars; replacing them cuts first-paint by several seconds on mobile and
# a noticeable fraction on desktop. Mirrors the Next.js EvidenceFingerprint
# component: 8 pie slices, opacity proportional to the score for each of
# the 8 modalities.
# ---------------------------------------------------------------------------

import math as _math


def make_fingerprint_svg(
    evidence: Mapping[str, float] | None,
    *,
    size: int = 56,
    color: str = "#0d62b8",
) -> str:
    """Return inline SVG for an 8-slice evidence donut.

    Each slice's fill-opacity encodes the modality score (0..1). Slice
    ordering matches ``LEGACY_ORDER`` (starts at 12 o'clock and proceeds
    clockwise). Returns a self-contained ``<svg>`` string that can be
    embedded directly via ``st.markdown(..., unsafe_allow_html=True)``.
    """
    normed = normalize_evidence(evidence)
    values = [max(0.0, min(1.0, normed[k])) for k in LEGACY_ORDER]
    n = len(values)
    cx, cy = size / 2, size / 2
    r_outer = size / 2 - 1
    r_inner = r_outer * 0.55  # donut hole for a crisper, ring-shaped look

    slices: list[str] = []
    # Start at -90° so slice 0 is at 12 o'clock, then clockwise.
    for i, v in enumerate(values):
        a0 = -_math.pi / 2 + (i / n) * 2 * _math.pi
        a1 = -_math.pi / 2 + ((i + 1) / n) * 2 * _math.pi
        x0_out, y0_out = cx + r_outer * _math.cos(a0), cy + r_outer * _math.sin(a0)
        x1_out, y1_out = cx + r_outer * _math.cos(a1), cy + r_outer * _math.sin(a1)
        x0_in, y0_in = cx + r_inner * _math.cos(a0), cy + r_inner * _math.sin(a0)
        x1_in, y1_in = cx + r_inner * _math.cos(a1), cy + r_inner * _math.sin(a1)
        # Large-arc-flag is always 0 here because each slice spans 45°.
        d = (
            f"M {x0_out:.2f},{y0_out:.2f} "
            f"A {r_outer:.2f},{r_outer:.2f} 0 0 1 {x1_out:.2f},{y1_out:.2f} "
            f"L {x1_in:.2f},{y1_in:.2f} "
            f"A {r_inner:.2f},{r_inner:.2f} 0 0 0 {x0_in:.2f},{y0_in:.2f} Z"
        )
        op = 0.12 + 0.88 * v  # floor of 0.12 so zero-evidence slices still register
        slices.append(
            f'<path d="{d}" fill="{color}" fill-opacity="{op:.3f}" '
            'stroke="#ffffff" stroke-width="0.6" />'
        )

    # Count "active" modalities (score > 0) and put the count in the hole.
    active = sum(1 for v in values if v > 0)
    label = f'<text x="{cx}" y="{cy + 3.5}" text-anchor="middle" ' \
            f'font-size="{size * 0.26:.1f}" font-weight="600" ' \
            f'fill="var(--masld-fg, #1a1a1a)" font-family="ui-monospace, monospace">' \
            f'{active}/{n}</text>'

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
        f'viewBox="0 0 {size} {size}" role="img" '
        f'aria-label="Evidence fingerprint: {active}/{n} modalities active">'
        + "".join(slices)
        + label
        + "</svg>"
    )


# ---------------------------------------------------------------------------
# Small-multiples radar (octagon polygon). Used by the Atlas "Radar wall"
# to let users scan hundreds of evidence fingerprints at a glance and
# spot unusual modality shapes that a heatmap or table can't convey.
# Pure SVG so a wall of 60 radars renders instantly vs. 60 Plotly figures.
# ---------------------------------------------------------------------------


def make_radar_wall_svg(
    evidence: Mapping[str, float] | None,
    *,
    size: int = 92,
    color: str = "#3b5ad2",
    show_ticks: bool = True,
) -> str:
    """Return inline SVG for an 8-axis polygon radar (small-multiples mode).

    Unlike ``make_fingerprint_svg`` (donut slices), this draws an actual
    polygon: axes at 8 angles, filled shape inside, grid rings at 0.25 /
    0.5 / 0.75 / 1.0. Optimised for compact (~80-100px) rendering in a
    grid where readers compare polygon *shapes* across genes.

    Axis order matches :data:`LEGACY_ORDER` (M1 Human at top, clockwise).
    """
    normed = normalize_evidence(evidence)
    values = [max(0.0, min(1.0, float(normed[k]))) for k in LEGACY_ORDER]
    n = len(values)
    cx, cy = size / 2, size / 2
    r_max = size / 2 - size * 0.14  # leave room for M-code ticks

    def _pt(angle: float, radius: float) -> tuple[float, float]:
        return cx + radius * _math.cos(angle), cy + radius * _math.sin(angle)

    angles = [-_math.pi / 2 + (i / n) * 2 * _math.pi for i in range(n)]

    # Background grid: concentric octagons at 0.25, 0.5, 0.75, 1.0
    grid_paths: list[str] = []
    for frac in (0.25, 0.5, 0.75, 1.0):
        pts = [_pt(a, r_max * frac) for a in angles]
        d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
        grid_paths.append(
            f'<path d="{d}" fill="none" stroke="var(--masld-border,#e3e6ec)" '
            f'stroke-width="0.5" opacity="{0.9 if frac == 1.0 else 0.5:.2f}" />'
        )

    # Axis spokes
    spokes: list[str] = []
    for a in angles:
        xe, ye = _pt(a, r_max)
        spokes.append(
            f'<line x1="{cx:.2f}" y1="{cy:.2f}" x2="{xe:.2f}" y2="{ye:.2f}" '
            f'stroke="var(--masld-border,#e3e6ec)" stroke-width="0.4" opacity="0.6" />'
        )

    # Filled polygon for this gene's evidence
    poly_pts = [_pt(a, r_max * v) for a, v in zip(angles, values)]
    poly_d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in poly_pts) + " Z"
    poly = (
        f'<path d="{poly_d}" fill="{color}" fill-opacity="0.28" '
        f'stroke="{color}" stroke-width="1.2" stroke-linejoin="round" />'
    )

    # Value dots on each axis (only non-zero, so sparse evidence reads clean)
    dots: list[str] = []
    for (x, y), v in zip(poly_pts, values):
        if v > 0:
            dots.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.6" fill="{color}" />'
            )

    # M1..M8 tick labels outside the polygon
    ticks: list[str] = []
    if show_ticks:
        for i, a in enumerate(angles):
            tx, ty = _pt(a, r_max * 1.16)
            ticks.append(
                f'<text x="{tx:.2f}" y="{ty + 2.4:.2f}" text-anchor="middle" '
                f'font-size="{size * 0.11:.1f}" font-weight="600" '
                f'fill="var(--masld-muted-fg,#6a7280)" '
                f'font-family="ui-monospace, monospace">M{i + 1}</text>'
            )

    active = sum(1 for v in values if v > 0)

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" '
        f'viewBox="0 0 {size} {size}" role="img" '
        f'aria-label="Evidence radar: {active}/{n} modalities active">'
        + "".join(grid_paths)
        + "".join(spokes)
        + poly
        + "".join(dots)
        + "".join(ticks)
        + "</svg>"
    )
