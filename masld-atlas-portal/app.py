"""MASLD Multi-Evidence Atlas -- Streamlit entrypoint.

Declares the 5 visible nav entries (Home, Atlas, Network, Translate,
Validation) plus a hidden ``Gene`` page that stays deep-linkable as
``/Gene?symbol=THRB``.

The sidebar uses the **dict form** of ``st.navigation`` (Streamlit 1.36+)
to group entries under section labels -- this mirrors the Next.js
sidebar's DISCOVER / TRANSLATE / TOOLS structure.

Hiding strategy: ``st.navigation`` has no per-page ``hidden`` flag. The
Gene page is registered normally, then its sidebar entry is hidden via a
single CSS rule in ``theme/custom.css`` that targets its stable
``href$="/Gene"`` attribute. Deep linking still works because
``st.Page(url_path="Gene")`` registers the route in the Streamlit router.

When ``st.navigation`` is used, Streamlit's automatic ``pages/`` discovery
is disabled — we keep page code under ``views/`` so there is no risk of
accidental duplicate routing.
"""
from __future__ import annotations

import streamlit as st

from data.bootstrap import ensure_data

ensure_data()

st.set_page_config(
    page_title="MASLD Multi-Evidence Atlas",
    page_icon="\U0001f9ec",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ---------------------------------------------------------------------------
# Page registry
# ---------------------------------------------------------------------------

HOME_PAGE = st.Page(
    "views/home.py",
    title="Home",
    icon=":material/home:",
    default=True,
)

ATLAS_PAGE = st.Page(
    "views/atlas.py",
    title="Atlas",
    icon=":material/grid_view:",
    url_path="Atlas",
)

NETWORK_PAGE = st.Page(
    "views/network.py",
    title="Network",
    icon=":material/hub:",
    url_path="Network",
)

TRANSLATE_PAGE = st.Page(
    "views/translate.py",
    title="Prioritize",
    icon=":material/sync_alt:",
    url_path="Prioritize",
)

VALIDATION_PAGE = st.Page(
    "views/validation.py",
    title="Validation",
    icon=":material/verified:",
    url_path="Validation",
)

# Registered so it stays deep-linkable as /Gene?symbol=THRB. The sidebar
# entry is hidden by CSS (see theme/custom.css — `[href$="/Gene"]`).
GENE_PAGE = st.Page(
    "views/gene.py",
    title="Gene",
    icon=":material/science:",
    url_path="Gene",
)


# ---------------------------------------------------------------------------
# Sidebar quick gene lookup (rendered BEFORE st.navigation so it sits above
# the nav). Submitting jumps to /Gene?symbol=<SYM>. Press "/" anywhere to
# focus the input — the keyboard handler is injected once and debounces
# itself via a window-level flag.
# ---------------------------------------------------------------------------

import streamlit.components.v1 as components  # noqa: E402

with st.sidebar:
    with st.form("sidebar_gene_search", clear_on_submit=True, border=False):
        gene_q = st.text_input(
            "Jump to gene",
            placeholder="gene symbol (e.g. THRB)  —  press  /",
            label_visibility="collapsed",
            key="sidebar_gene_search_input",
        )
        submitted = st.form_submit_button(
            "Open gene card  →",
            use_container_width=True,
        )

    if submitted and gene_q:
        sym = gene_q.strip().upper()
        if sym:
            # JS navigation: st.switch_page doesn't accept query params, and a
            # meta-refresh inside Streamlit's iframe is unreliable. Navigating
            # window.top bubbles out of the iframe to the actual browser URL.
            components.html(
                f'<script>window.top.location.href="/Gene?symbol={sym}";</script>',
                height=0,
            )
            st.stop()

    # Keyboard shortcut: press "/" anywhere (except while typing in an input
    # already) to focus the gene search box. Idempotent — re-attaches cleanly
    # on every Streamlit rerun by replacing the previous handler.
    components.html(
        """
        <script>
        (function() {
          const doc = window.parent.document;
          if (window.__masldSearchHandler) {
            doc.removeEventListener('keydown', window.__masldSearchHandler);
          }
          const handler = function(e) {
            if (e.key !== '/' || e.ctrlKey || e.metaKey || e.altKey) return;
            const t = e.target;
            if (t && ['INPUT', 'TEXTAREA', 'SELECT'].includes(t.tagName)) return;
            if (t && t.isContentEditable) return;
            const inp = doc.querySelector('input[aria-label="Jump to gene"]');
            if (inp) { e.preventDefault(); inp.focus(); inp.select(); }
          };
          window.__masldSearchHandler = handler;
          doc.addEventListener('keydown', handler);
        })();
        </script>
        """,
        height=0,
    )



# ---------------------------------------------------------------------------
# Navigation — grouped by section (matches Next.js sidebar.tsx NAV_SECTIONS).
# The dict-form of st.navigation renders each key as a section heading.
# ---------------------------------------------------------------------------

pg = st.navigation(
    [
        HOME_PAGE,
        ATLAS_PAGE,
        TRANSLATE_PAGE,
        VALIDATION_PAGE,
        NETWORK_PAGE,
        GENE_PAGE,  # hidden from sidebar via CSS (see theme/custom.css)
    ],
    position="sidebar",
    expanded=True,
)

# ---------------------------------------------------------------------------
# Floating dark/light mode toggle — renders on every page.
# Injects a fixed-position button into the parent document's <body> so it
# sits in the top-right corner regardless of which view the user is on.
# Preference persists in localStorage and survives navigation / refresh.
# ---------------------------------------------------------------------------

components.html(
    """
    <script>
    (function() {
      const doc = window.parent.document;
      const root = doc.documentElement;

      // Remove any stale button from a previous rerun so we don't stack.
      const prev = doc.getElementById('masld-theme-float');
      if (prev) prev.remove();

      // Default is ALWAYS light unless the user has explicitly toggled to dark.
      // We deliberately do not follow prefers-color-scheme — dark mode is
      // opt-in only.
      const STORED = localStorage.getItem('masld-theme');  // 'light' | 'dark' | null
      let current = (STORED === 'dark') ? 'dark' : 'light';
      root.setAttribute('data-masld-theme', current);

      const btn = doc.createElement('button');
      btn.id = 'masld-theme-float';
      btn.setAttribute('aria-label', 'Toggle dark mode');
      btn.style.cssText = [
        'position:fixed', 'top:14px', 'right:14px', 'z-index:9999',
        'height:38px', 'padding:0 14px', 'border-radius:999px',
        'border:1px solid var(--masld-border, rgba(128,128,128,0.25))',
        'background:var(--masld-card, #ffffff)',
        'color:var(--masld-fg, #1a1a1a)',
        'cursor:pointer', 'font-size:13px', 'font-weight:500',
        'display:inline-flex', 'align-items:center', 'gap:8px',
        'box-shadow:0 2px 6px rgba(0,0,0,0.08)',
        'transition:transform 0.15s ease, background 0.15s ease',
        'font-family:inherit', 'line-height:1',
      ].join(';');
      function render() {
        // In dark mode: clicking switches to light → label "Light Mode" + sun.
        // In light mode: clicking switches to dark  → label "Dark Mode"  + moon.
        const icon = current === 'dark' ? '☀' : '☾';
        const label = current === 'dark' ? 'Light Mode' : 'Dark Mode';
        btn.innerHTML = '<span style="font-size:15px;line-height:1;">' + icon +
                        '</span><span>' + label + '</span>';
        btn.title = 'Switch to ' + (current === 'dark' ? 'light' : 'dark') + ' mode';
      }
      render();
      btn.addEventListener('mouseenter', function() {
        btn.style.transform = 'scale(1.08)';
      });
      btn.addEventListener('mouseleave', function() {
        btn.style.transform = 'scale(1)';
      });
      btn.addEventListener('click', function() {
        current = (current === 'dark') ? 'light' : 'dark';
        localStorage.setItem('masld-theme', current);
        root.setAttribute('data-masld-theme', current);
        render();
      });
      doc.body.appendChild(btn);
    })();
    </script>
    """,
    height=0,
)

pg.run()
