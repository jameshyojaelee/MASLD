#!/usr/bin/env python3
"""
Canonical Okabe-Ito modality palette for MASLD atlas Figure 1 (locked 2026-05-29).

SINGLE SOURCE OF TRUTH for panel colors. All Fig-1 Python panels import from here so the
figure (and Figs 2-5 downstream) share one modality->color key. Okabe-Ito = colorblind-safe.

Usage (works regardless of cwd):
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from fig1_palette import MODALITY_COLORS, SOURCE_COLORS, SPECIES, CONTROL_GRAY
"""

# --- 7 data modalities (A1 schematic + B1 sunburst share these keys) ---
MODALITY_COLORS = {
    "bulk":    "#0072B2",  # bulk RNA-seq          - Okabe-Ito blue
    "scrna":   "#56B4E9",  # single-cell / snRNA   - sky blue
    "spatial": "#CC79A7",  # spatial               - reddish purple
    "atac":    "#E69F00",  # ATAC / epigenomic     - orange
    "proteo":  "#D55E00",  # proteomics            - vermillion
    "gwas":    "#C2185B",  # genetics (GWAS/eQTL)  - deep crimson (matches S2 causal in figS_convergence/figS_network)
    "pharma":  "#999999",  # pharmacological profiling - neutral gray
}

# extra data types used only by the evidence-source panels (not A1/B1 modalities)
MOUSE_BULK   = "#009E73"   # mouse bulk RNA-seq (cross-species)  - bluish green
ESSENTIALITY = "#F0E442"   # DepMap essentiality                 - yellow

# secondary axis: species (desaturated tints, deliberately NOT modality hues)
SPECIES = {"human": "#6BAED6", "mouse": "#74C476"}

CONTROL_GRAY = "#9E9E9E"   # control/healthy / "absent" everywhere

# --- evidence sources (S1-S8) inherit the color of their underlying modality ---
# Channels per RNA-seq/results/multi_evidence/sources_active_definition.md (S6 dropped).
SOURCE_COLORS = {
    "S1": MODALITY_COLORS["bulk"],     # human bulk RNA-seq
    "S8": MOUSE_BULK,                  # mouse bulk RNA-seq
    "S2": MODALITY_COLORS["gwas"],     # genetic causal (TWAS+COLOC)
    "S4": MODALITY_COLORS["atac"],     # epigenomic (ATAC+SCENIC)
    "S5": MODALITY_COLORS["spatial"],  # spatial
    "S7": MODALITY_COLORS["proteo"],   # proteomics
    "S3": ESSENTIALITY,                # DepMap essentiality
}
