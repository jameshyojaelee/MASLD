"""Shared layout helpers for the MASLD Streamlit portal.

All primitives target visual parity with the Next.js portal at
``masld-atlas-v2``. The CSS they depend on lives in ``theme/custom.css``
(injected by :func:`inject_css`).
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable

import streamlit as st


# ---------------------------------------------------------------------------
# CSS injection
# ---------------------------------------------------------------------------

_CSS_PATH = Path(__file__).resolve().parent.parent / "theme" / "custom.css"


def inject_css() -> None:
    """Inject ``theme/custom.css`` into the current page.

    Safe to call multiple times within a single Streamlit run; repeated
    ``st.markdown`` writes are idempotent at the DOM level.
    """
    if _CSS_PATH.exists():
        st.markdown(
            f"<style>{_CSS_PATH.read_text()}</style>",
            unsafe_allow_html=True,
        )


# ---------------------------------------------------------------------------
# Headers / brand
# ---------------------------------------------------------------------------

def page_header(title: str, subtitle: str | None = None) -> None:
    """Render a page-level title + optional subtitle (Next.js H1 + muted p)."""
    st.markdown(
        f"""
        <div style="margin:0 0 16px 0;">
            <h1 style="margin:0;">{title}</h1>
            {f'<p style="margin:4px 0 0 0;color:var(--masld-muted-fg);font-size:14px;line-height:1.5;">{subtitle}</p>' if subtitle else ''}
        </div>
        """,
        unsafe_allow_html=True,
    )


def sidebar_header() -> None:
    """No-op. Brand block removed per user request — sidebar starts directly
    with the nav. Kept as a stub so existing callers don't break."""
    return


def eyebrow(text: str) -> None:
    """Render a small uppercase section label (the 'ATLAS AT A GLANCE' style)."""
    st.markdown(
        f'<div class="masld-eyebrow">{text}</div>',
        unsafe_allow_html=True,
    )


def spacer(size: str = "md") -> None:
    """Emit a vertical spacer div. ``size`` is one of ``sm`` (8px),
    ``md`` (16px), or ``lg`` (32px). Keeps ad-hoc ``<div style='margin-top...'>``
    out of view code so rhythm stays consistent across pages."""
    size = size if size in ("sm", "md", "lg") else "md"
    st.markdown(
        f'<div class="masld-spacer-{size}"></div>',
        unsafe_allow_html=True,
    )


# ---------------------------------------------------------------------------
# Cards
# ---------------------------------------------------------------------------

def stat_card(label: str, value: str | int | float, hint: str | None = None) -> None:
    """Compact stat card (matches Next.js ``StatCard``).

    - ~140px min-width, ~80px height
    - label: 11px uppercase tracked muted
    - value: 22px bold tabular-nums
    - hint: 11px muted (optional)
    """
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        try:
            shown = f"{int(value):,}"
        except (TypeError, ValueError):
            shown = str(value)
    else:
        shown = str(value)

    st.markdown(
        f"""
        <div class="masld-stat">
          <div class="masld-stat-value">{shown}</div>
          <div class="masld-stat-label">{label}</div>
          {f'<div class="masld-stat-hint">{hint}</div>' if hint else ''}
        </div>
        """,
        unsafe_allow_html=True,
    )


def badge(text: str, variant: str = "neutral") -> str:
    """Return HTML for a pill-style badge. Variants: primary|success|neutral."""
    variant = variant if variant in ("primary", "success", "neutral") else "neutral"
    return f'<span class="masld-badge {variant}">{text}</span>'


@contextmanager
def section_card(title: str | None = None):
    """Context manager that wraps content in a bordered card.

    Uses ``st.container(border=True)`` (Streamlit >=1.35). Falls back silently
    to a plain container if the flag is unsupported. Title is rendered as a
    separator-divided heading inside the card (matches Next.js ``<Card>`` +
    ``<CardHeader>``).
    """
    try:
        container = st.container(border=True)
    except TypeError:  # pragma: no cover — only on very old Streamlit
        container = st.container()
    with container:
        if title:
            st.markdown(
                f'<div class="masld-section-title">{title}</div>',
                unsafe_allow_html=True,
            )
        yield


def render_section_card(title: str, body: Callable[[], None]) -> None:
    """Functional wrapper: call ``body()`` inside a bordered card with a title."""
    with section_card(title):
        body()


# ---------------------------------------------------------------------------
# Plotly helper
# ---------------------------------------------------------------------------

# Default Plotly config: hide modebar entirely (matches the Next.js aesthetic
# where charts are clean, non-interactive-looking). Individual callers can
# override with ``config={...}`` if they need the modebar for a specific chart.
_DEFAULT_PLOTLY_CONFIG: dict[str, Any] = {
    "displayModeBar": False,
    "displaylogo": False,
    "responsive": True,
}


def plot_figure(
    fig,
    *,
    use_container_width: bool = True,
    config: dict | None = None,
    **kwargs,
) -> None:
    """Wrapper around ``st.plotly_chart`` that defaults to hidden modebar.

    Use this in place of ``st.plotly_chart`` on pages that want the Next.js
    chart aesthetic. Additional kwargs are forwarded verbatim.
    """
    merged_config = dict(_DEFAULT_PLOTLY_CONFIG)
    if config:
        merged_config.update(config)
    st.plotly_chart(
        fig,
        use_container_width=use_container_width,
        config=merged_config,
        **kwargs,
    )
