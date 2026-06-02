"""Gene compare panel -- Phase 2 stub.

Will render up to 5 genes' modality radars on the same polar axes for direct
comparison. For Phase 1 this is a placeholder that explains the intent.
"""
from __future__ import annotations

import streamlit as st


def render_compare_placeholder() -> None:
    st.info(
        "Gene comparison panel arrives in Phase 2. It will mirror the Next.js "
        "`gene-compare-panel.tsx`: up to five genes plotted on a shared 8-modality "
        "radar with LFC/padj side-by-side."
    )
