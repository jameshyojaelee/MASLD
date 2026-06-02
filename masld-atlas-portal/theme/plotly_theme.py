"""Plotly theme matching the Next.js portal.

Import ``PLOTLY_THEME`` and apply via ``fig.update_layout(**PLOTLY_THEME)``,
or use :func:`apply_theme(fig)` for convenience.

Color palettes mirror ``src/components/network/*`` LAYER_COLOR_MAP and the
evidence-source palette used on the atlas and convergence pages.
"""
from __future__ import annotations

from typing import Any

# ---- Categorical palettes ---------------------------------------------------

# Network layer colors (mirrors layer_metadata.json / sidebar config).
LAYER_COLOR_MAP: dict[str, str] = {
    "ppi":           "#34495e",
    "coexpr":        "#e74c3c",
    "coexpr_stable": "#1abc9c",
    "regulon":       "#d35400",
    "genetic":       "#9b59b6",
    "lr":            "#f39c12",
    "cerna":         "#e91e63",
    "xspecies":      "#2ecc71",
}

# Evidence-source (modality) palette — oklch(0.65 0.18 ...) approximations.
MODALITY_COLORS: list[str] = [
    "#23a369",  # s1 human bulk
    "#7b5cff",  # s2 genetic / mouse
    "#d67f3c",  # s3 essential
    "#c65dba",  # s4 epigenomic
    "#c8942a",  # s5 spatial
    "#2e93b0",  # s6 single-cell
    "#d24b7a",  # s7 mouse
    "#9ba33e",  # s8 proteomics
]

# Sequential (up/down DE) colors
DEG_UP = "#d13b3b"
DEG_DOWN = "#2566c2"
DEG_NONSIG = "#9aa3ad"

# ---- Theme dict -------------------------------------------------------------

_AXIS = dict(
    showgrid=True,
    gridcolor="#eceff3",
    gridwidth=1,
    zerolinecolor="#d9dde4",
    linecolor="#d9dde4",
    tickfont=dict(size=11, color="#4a4f5a"),
    title=dict(font=dict(size=12, color="#4a4f5a")),
)

PLOTLY_THEME: dict[str, Any] = {
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
    "font": dict(
        family=(
            "Inter, 'Segoe UI', -apple-system, BlinkMacSystemFont, "
            "'Helvetica Neue', Arial, sans-serif"
        ),
        size=12,
        color="#1a1a1a",
    ),
    "margin": dict(l=10, r=10, t=10, b=10),
    "colorway": MODALITY_COLORS,
    "xaxis": _AXIS,
    "yaxis": _AXIS,
    "legend": dict(
        bgcolor="rgba(0,0,0,0)",
        font=dict(size=11, color="#4a4f5a"),
        bordercolor="rgba(0,0,0,0)",
    ),
    "hoverlabel": dict(
        bgcolor="#ffffff",
        bordercolor="#d9dde4",
        font=dict(size=12, family="Inter, sans-serif", color="#1a1a1a"),
    ),
}


def apply_theme(fig, *, margin: dict | None = None) -> None:
    """In-place theme application for a Plotly figure.

    Callers sometimes want to override the default tight margins (e.g. when
    plotting long vertical bar charts) — pass ``margin=dict(...)`` to merge.
    """
    layout_update = dict(PLOTLY_THEME)
    if margin:
        layout_update = {**layout_update, "margin": margin}
    fig.update_layout(**layout_update)
