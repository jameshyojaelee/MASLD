#!/usr/bin/env python3
"""
Figure 1b (alt) — Multi-Modal Evidence Convergence Wheel, Modular Sunburst Variation
Matches the exact style of the screenshot with proper layout, right-hand legend,
radiating outer bars, and layered structure.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import numpy as np
import matplotlib.colors as mc
import os

# ─── Sanjana Lab Publication Standards ───────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":        10,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
    "figure.dpi":       150,
})

TEXT_DARK  = "#2B2B2B" # warm near-black (Liang-style ink, not pure jet)
LABEL_GRAY = "#7A7A7A"
BG_COLOR   = "white"
LABEL_SIZE = 9  # unified font size for all text except the center "MASLD Atlas"

# Liang et al. 2025 cell-line palette extended to 7 modalities. Each hue
# sits in its own family with ~50% chroma — distinct at panel scale but
# softer than primary jewel tones.
# Canonical Okabe-Ito modality palette (single source of truth)
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS

# Species ring: matched to the modality saturation so it sits in the same
# visual register without dominating.
BINARY_COLORS = {
    "Human": "#6FA2C8",  # cornflower
    "Mouse": "#8FBE6F",  # spring sage
}

# Disease-vs-control composition per modality, RESOLVED BY SPECIES — drives the
# proportional INNER ring (human block | mouse block, each split disease/control).
# Format {"human": (disease, control), "mouse": (disease, control)}; omit a species
# with no data. Counts sourced 2026-06-18 from harmonized metadata (proportions render;
# the human|mouse split angle is taken from the species `context` so it aligns with the
# species ring, and disease/control within each block uses these species-specific n):
#   bulk    = human 1124/160 (9 cohorts, unified_metadata.csv) + mouse 179/114 (unified_mouse_metadata.csv)
#   scrna   = human 196/64 (donor condition_binary; 9 NA excl). Atlas scRNA is HUMAN-ONLY:
#             all 7 datasets are human liver studies (269 donors, donor_metadata.tsv). The
#             prior (283,114) species split was spurious — corrected to (269,0) 2026-06-18.
#   spatial = human 22/5 (Govaere 12 + Vu 10 disease vs GSE192741 5 healthy/steatotic; dataset-level)
#   atac    = human 13/5 (donor_metadata_curated) + mouse 9/3 (config SCD_Control = control)
#   proteo  = human 111/19 (plasma 65/7 is_masld + liver 46/12 PXD051911 group)
#   pharma / gwas = None → not a case/control design (compounds / population genetics)
SPECIES_DISEASE = {
    "bulk":    {"human": (1124, 160), "mouse": (179, 114)},
    "scrna":   {"human": (196, 64)},
    "spatial": {"human": (22, 5)},
    "atac":    {"human": (13, 5), "mouse": (9, 3)},
    "proteo":  {"human": (111, 19)},
    "pharma":  None,
    "gwas":    None,
}
CONTROL_GRAY  = "#9E9E9E"  # canonical control gray (FIGURE_GUIDELINES)
DISEASE_COLOR = "#C0524E"  # single semantic disease red — deliberately NOT a modality/blue hue
NA_GRAY       = "#E2E2E2"  # neutral fill for non-case/control modalities
SPECIES_LIGHTEN = 0.45     # mouse arcs = this much lighter than human (lightness = species)

# GWAS inner ring = ANCESTRY of the 23 GWAS (by study count). The 3 eQTL
# (GTEx/Broadaway/sc-eQTL, EUR reference panels) are EXCLUDED — not ancestry-stratified
# GWAS. Violet gradient (dark EUR → pale SAS). Source: gwas_registry.tsv.
ANCESTRY = {
    "gwas": [("EUR", 14, "#3A2259"), ("EAS", 3, "#7B4FB0"),
             ("AFR", 3, "#C08AC9"), ("SAS", 3, "#EBD4E8")],
}

# (key, title, l1, l2, (n_human, n_mouse), weight, datasets)
# Verified from unified metadata & processed data (updated 2026-06-03)
# scRNA: human donors + mouse donors across 7 source datasets
# GWAS: 23 GWAS (14 EUR + 3 EAS + 3 AFR + 3 SAS Pan-UKBB) + 3 eQTL (GTEx v8, Broadaway, sc-eQTL)
#   Registry: gwas_registry.tsv; sex-stratified Pan-UKBB F/M arms + cirrhosis/HCC GWAS excluded (2026 portfolio refactor)
# bulk: 9 human QC-pass (1,259) + 8 mouse QC-pass (463) = 1,722 total; 17 datasets
# pharma: LINCS L1000, network proximity, DGIdb/OT, ClinicalTrials MASH pipeline, multi-layer (5 sources)
# spatial: Govaere 2026 (GeoMx+CosMx, integrated 2026-05-21), GSE192741 (Visium), Vu_et_al_2025 (3 active)
SECTORS = [
    ("bulk",    "RNA-seq",                    "17\nRNA-seq",      "1,722\nsamples",  (1259, 463),  18.0, 17),
    ("scrna",   "scRNA-seq",                  "7\nscRNA-seq",     "2.29M+\ncells",   (269, 0),      8.0, 7),
    ("spatial", "Spatial Omics",              "3\nSpatial",       "27\nsamples",     (27, 0),        4.5, 3),
    ("atac",    "ATAC-seq",                   "2\nATAC-seq",      "30\nsamples",     (18, 12),       4.5, 2),
    ("proteo",  "Proteomics",                 "2\nProteomics",    "130\nsamples",    (130, 0),       4.5, 2),
    ("pharma",  "Therapeutics",               "5\nTherapeutics",  "1,107\ncmpds",    (1, 0),         4.5, 5),
    ("gwas",    "Genomics",                   "23 GWAS\n3 eQTL",  "1.5M+\nsubj.",   (1, 0),        16.0, 26),
]

# Per-dataset sample counts for proportional outer-bar heights
_DATASET_SIZES_RAW = {
    "bulk": [  # 9 Human (QC-pass from sample_qc_report.csv) + 8 Mouse; STAR -s2 canonical 2026-05-28
        358, 215, 160, 142, 97, 93, 76, 64, 55,       # Human QC-pass: Chen(GSE213621), Suppli(GSE135251), Govaere(GSE193066), Bril(GSE162694), Kawamura(GSE167523), Kozumi(GSE174478), Hoshida(GSE130970), Verschuren(GSE240729), Hoang(GSE126848)
        213, 151, 29, 22, 12, 11, 10, 6,               # Mouse (unified_mouse_metadata.csv): CDAHFD/FPC(GSE162876), LIDPAD(GSE159911), HFD(GSE224069), HFD-long(GSE274914), InHouse_MCD, MCD-SE(GSE156918), GAN(GSE225616), MCD-PE(GSE205974)
    ],
    "scrna": [  # 7 source datasets (unique donors; SRR technical runs excluded)
        38, 117, 67, 21, 20, 2, 4,  # Liver_Atlas(38 donors), GSE244832, GSE202379, GSE185477, GSE136103, GSE189600(2 human donors), GSE174748
    ],
    "spatial": [12, 10, 5],  # Govaere2026(GeoMx 8pt+CosMx 4pt), Vu_et_al_2025(10 arrays), GSE192741(5 Visium sections)
    "atac": [18, 12],                   # Human_Multiome (18 donors), Mouse_Bulk (12 samples)
    "proteo": [72, 58],                 # PXD052937 plasma DIA-MS, PXD051911 liver DIA-MS (GSE276114 removed 2026-06: GEO confirms bulk RNA-seq, not SomaScan proteomics)
    "pharma": [1107, 1173, 58, 23, 20],  # LINCS(1107 compounds), network proximity(1173 screened), DGIdb/OT, ClinicalTrials MASH, multi-layer
    "gwas": [  # 23 GWAS + 3 eQTL (subjects in thousands); from gwas_registry.tsv (2026 portfolio refactor)
        # EUR (14): N_tot in thousands
        778.6,   # 2021_34841290 NAFLD EUR (Ghodsian meta)
        438.9,   # FinnGen NAFLD
        435.0,   # FinnGen NASH
        343.9,   # UKBB_ALT
        343.9,   # UKBB_AST
        343.9,   # UKBB_GGT
        397.0,   # 2023 UKBB NAFLD (Sveinbjornsson cohort)
        359.2,   # deCODE NAFLD
        44.9,    # 2022 PDFF EUR
        36.1,    # 2021 PDFF EUR (34128465)
        32.9,    # 2021 PDFF EUR (34957434)
        32.7,    # Intermountain NAFLD EUR
        9.5,     # 2020 NAFLD EUR
        8.4,     # 2019 NAFLD EUR
        # EAS (3)
        160.0,   # BBJ_ALT
        160.0,   # BBJ_AST
        160.0,   # BBJ_GGT
        # AFR Pan-UKBB (3)
        6.6, 6.6, 6.6,
        # SAS/CSA Pan-UKBB (3)
        8.9, 8.9, 8.9,
        # eQTL (3, in thousands of donors)
        1.183, 0.312, 0.208,
    ],
}

# Shuffle within each sector (seeded) so bars aren't monotonically sorted
_rng_shuffle = np.random.RandomState(99)
DATASET_SIZES = {}
for _k, _v in _DATASET_SIZES_RAW.items():
    _arr = list(_v)
    _rng_shuffle.shuffle(_arr)
    DATASET_SIZES[_k] = _arr

BAR_LEN_MIN = 0.08   # shortest bar (smallest dataset)
BAR_LEN_MAX = 0.45   # tallest bar  (largest dataset)

# Per-sector multiplier so stored values map to real-world units for labels
# GWAS values are stored in thousands → multiply by 1000 for display
_UNIT_MULT = {"gwas": 1000}

# Per-sector units for axis labels
_SECTOR_UNITS = {
    "bulk":    "samples",
    "scrna":   "donors",
    "spatial": "sections",
    "atac":    "samples",
    "proteo":  "samples",
    "pharma":  "cmpds",
    "gwas":    "subjects",
}

def _fmt_axis_val(key, raw_val):
    """Format a value for a radial axis tick, respecting unit multipliers."""
    real = raw_val * _UNIT_MULT.get(key, 1)
    if real >= 1_000_000:
        return f"{real/1e6:.1f}M"
    elif real >= 10_000:
        return f"{real/1e3:.0f}K"
    elif real >= 1_000:
        return f"{real/1e3:.1f}K"
    elif real < 1:
        return f"{int(real*1000)}"
    else:
        return f"{int(real)}"

def get_text_color(bg_color):
    rgb = mc.to_rgb(bg_color)
    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return "white" if lum < 0.65 else TEXT_DARK

def blend_white(color, factor):
    rgb = np.array(mc.to_rgb(color))
    white = np.array([1.0, 1.0, 1.0])
    blended = rgb * (1 - factor) + white * factor
    return mc.to_hex(blended)

def main():
    fig = plt.figure(figsize=(10.5, 7.5))
    fig.patch.set_facecolor(BG_COLOR)

    # Left axes: Main sunburst
    ax = fig.add_axes([0.02, 0.05, 0.65, 0.90])
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.65, 1.65)
    ax.set_aspect("equal")
    ax.axis("off")

    # Right axes: Half-sunburst legend (vertically centred; no modality legend above)
    ax_leg = fig.add_axes([0.72, 0.30, 0.25, 0.50])
    ax_leg.set_xlim(-0.3, 1.3)
    ax_leg.set_ylim(-1.35, 1.5)
    ax_leg.set_aspect("equal")
    ax_leg.axis("off")

    # =========================================================================
    # MAIN SUNBURST 
    # =========================================================================

    # Central hub
    hub_r = 0.28
    hub = plt.Circle((0, 0), hub_r, facecolor="white",
                      edgecolor=TEXT_DARK, linewidth=1.5, zorder=10)
    ax.add_patch(hub)
    
    # Decorator ring for hub
    hub_outer = plt.Circle((0, 0), hub_r + 0.015, facecolor="none",
                           edgecolor=LABEL_GRAY, linewidth=0.8, zorder=9, alpha=0.5)
    ax.add_patch(hub_outer)

    ax.text(0, 0, "MASLD\nAtlas", fontsize=13, fontweight="bold",
            color=TEXT_DARK, ha="center", va="center", zorder=11)

    # Rings (Hierarchical thicknesses)
    # Layer 0: Disease vs control composition (inner, proportional arcs)
    # Layer 1: Datasets — modality name + n datasets
    # Layer 2: Binary species flag (outermost slim ring)
    # (The old inner per-modality "Counts" totals ring was removed 2026-06-18 — it
    #  duplicated the outer bars; replaced here by disease/control composition.)
    r_starts = [0.31, 0.62, 0.95]
    r_ends   = [0.62, 0.95, 1.05]

    MOD_SHADE = 0.15   # modality (name) ring — softened so it doesn't punch
    
    total_weight = sum(s[5] for s in SECTORS)
    gap_deg = 3.5  # Distinct sector gaps
    total_arc = 360 - len(SECTORS) * gap_deg
    
    current_angle = 185.0
    rng = np.random.RandomState(42)

    
    for i, (key, title, l1, l2, context, weight, n_data) in enumerate(SECTORS):
        span = (weight / total_weight) * total_arc
        end_angle = current_angle - span
        
        base_color = MODALITY_COLORS[key]

        # 1. Inner ring — disease vs control composition, split by species.
        #    Each sector is a human block then a mouse block (angles from `context`
        #    so they line up with the species ring); within each block, disease
        #    (red) vs control (gray). Mouse arcs are lightened so species reads off
        #    lightness without adding another hue.
        r_in0, r_out0 = r_starts[0], r_ends[0]
        sd = SPECIES_DISEASE.get(key)
        if sd is not None:
            sp_total = sum(context) if sum(context) > 0 else 1
            all_d = sum(v[0] for v in sd.values())
            all_c = sum(v[1] for v in sd.values())
            overall_dfrac = all_d / (all_d + all_c)
            seg_start = current_angle
            for sp_idx, sp_name in enumerate(("human", "mouse")):
                sp_count = context[sp_idx]
                if sp_count <= 0:
                    continue
                sp_span = (sp_count / sp_total) * span
                if sp_name in sd:
                    d, c = sd[sp_name]
                    dfrac = d / (d + c)
                else:
                    dfrac = overall_dfrac  # mouse block w/o per-species data → overall
                lighten = 0.0 if sp_name == "human" else SPECIES_LIGHTEN
                for cond_frac, base in [(dfrac, DISEASE_COLOR), (1 - dfrac, CONTROL_GRAY)]:
                    if cond_frac > 0:
                        seg_end = seg_start - cond_frac * sp_span
                        ax.add_patch(Wedge((0, 0), r_out0, seg_end, seg_start, width=r_out0 - r_in0,
                                           facecolor=blend_white(base, lighten),
                                           edgecolor=BG_COLOR, linewidth=2.0, zorder=5))
                        seg_start = seg_end
        elif ANCESTRY.get(key) is not None:
            # GWAS/eQTL — ancestry composition (violet gradient, dark EUR → pale SAS)
            anc = ANCESTRY[key]
            tot = sum(n for _, n, _ in anc)
            seg_start = current_angle
            for _lab, n, col in anc:
                if n > 0:
                    seg_end = seg_start - (n / tot) * span
                    ax.add_patch(Wedge((0, 0), r_out0, seg_end, seg_start, width=r_out0 - r_in0,
                                       facecolor=col, edgecolor=BG_COLOR, linewidth=2.0, zorder=5))
                    seg_start = seg_end
        else:
            # pharma — sources not commensurable as case/control or ancestry → neutral fill
            ax.add_patch(Wedge((0, 0), r_out0, end_angle, current_angle, width=r_out0 - r_in0,
                               facecolor=NA_GRAY, edgecolor=BG_COLOR, linewidth=2.0, zorder=5))

        # 2. Modality ring — name + n datasets
        r_in1, r_out1 = r_starts[1], r_ends[1]
        color = blend_white(base_color, MOD_SHADE)
        ax.add_patch(Wedge((0, 0), r_out1, end_angle, current_angle, width=r_out1 - r_in1,
                           facecolor=color, edgecolor=BG_COLOR, linewidth=2.0, zorder=5))
        # Text placement — adaptive for narrow wedges
        r_mid = (r_in1 + r_out1) / 2.0
        mid_angle = (current_angle + end_angle) / 2.0
        rot = mid_angle - 90
        if rot < -90:
            rot += 180
        elif rot > 90:
            rot -= 180
        x = r_mid * np.cos(np.radians(mid_angle))
        y = r_mid * np.sin(np.radians(mid_angle))
        t_color = get_text_color(color)
        base_size = LABEL_SIZE
        if span < 18:
            f_size = base_size * max(0.55, span / 22.0)
        elif span < 25:
            f_size = base_size * max(0.7, span / 25.0)
        else:
            f_size = base_size
        ax.text(x, y, l1, rotation=rot, color=t_color,
                fontsize=f_size, fontweight="medium",
                ha="center", va="center", zorder=6)

        # 3. Outer slim layer: Binary Species Flag (Proportionally Scaled)
        total_samples = sum(context)
        if total_samples > 0:
            h_frac = context[0] / total_samples
            m_frac = context[1] / total_samples
            
            # Start from outer bound (current_angle)
            seg_start = current_angle
            
            for f, col_name in zip([h_frac, m_frac], ["Human", "Mouse"]):
                if f > 0:
                    seg_span = f * span
                    seg_end = seg_start - seg_span
                    
                    sp_color = BINARY_COLORS[col_name]
                    w_species = Wedge((0, 0), r_ends[2], seg_end, seg_start, width=r_ends[2] - r_starts[2],
                                      facecolor=sp_color, edgecolor=BG_COLOR, linewidth=1.5, zorder=5)
                    ax.add_patch(w_species)
                    
                    seg_start = seg_end
                          
        # 2b. Outer bars (Layer 3) — one bar per dataset, height ∝ dataset size (linear % of max)
        margin = 1.0
        usable_span = span - 2 * margin
        r_base = r_ends[-1] + 0.02
        if n_data > 0:
            sizes = DATASET_SIZES.get(key, [1] * n_data)
            s_max = max(sizes)
            bar_span = usable_span / n_data
            for j in range(n_data):
                b_start = current_angle - margin - j * bar_span
                b_end = b_start - bar_span

                b_gap = max(0.6, bar_span * 0.18)
                # Linear scaling: bar height = fraction of sector max
                frac = sizes[j] / s_max if s_max > 0 else 0.5
                length = BAR_LEN_MIN + frac * (BAR_LEN_MAX - BAR_LEN_MIN)

                spike_color = "#B0B0B0"  # uniform muted gray for all bars

                w_bar = Wedge((0, 0), r_base + length, b_end + b_gap, b_start - b_gap,
                              width=length, facecolor=spike_color, edgecolor=BG_COLOR,
                              linewidth=0.8, alpha=0.55, zorder=4)
                ax.add_patch(w_bar)

            # Per-sector radial axis: line at the leading edge of the sector
            # with a single tick at the per-modality max value, so bar heights
            # are interpretable in the modality's own units (samples / donors /
            # subjects / compounds).
            axis_ang = current_angle - margin
            axis_rad = np.radians(axis_ang)
            r_axis0  = r_base
            r_axis1  = r_base + BAR_LEN_MAX
            # Axis stem
            ax.plot([r_axis0 * np.cos(axis_rad), r_axis1 * np.cos(axis_rad)],
                    [r_axis0 * np.sin(axis_rad), r_axis1 * np.sin(axis_rad)],
                    color="#555555", linewidth=0.7, zorder=6)
            # Outward tick at max, perpendicular to the radial axis
            tick_len = 0.025
            tang_x = -np.sin(axis_rad) * tick_len
            tang_y =  np.cos(axis_rad) * tick_len
            tx0 = r_axis1 * np.cos(axis_rad)
            ty0 = r_axis1 * np.sin(axis_rad)
            ax.plot([tx0 - tang_x, tx0 + tang_x],
                    [ty0 - tang_y, ty0 + tang_y],
                    color="#555555", linewidth=0.7, zorder=6)
            # Max-value tick label, sitting just beyond the axis line
            s_max = max(sizes)
            unit  = _SECTOR_UNITS.get(key, "")
            max_label = f"{_fmt_axis_val(key, s_max)} {unit}".strip()
            r_label = r_axis1 + 0.06
            lab_x = r_label * np.cos(axis_rad) - tang_x * 0.6
            lab_y = r_label * np.sin(axis_rad) - tang_y * 0.6
            # Tangential rotation so label reads along the axis
            rot = (axis_ang - 90) % 360
            if rot > 90 and rot < 270:
                rot -= 180
            ax.text(lab_x, lab_y, max_label, rotation=rot,
                    fontsize=LABEL_SIZE, color="#444444",
                    ha="center", va="center", zorder=6)

        current_angle = end_angle - gap_deg

    # =========================================================================
    # HALF-SUNBURST LEGEND
    # =========================================================================
    
    cx, cy = 0.0, 0.0
    
    leg_colors = [MODALITY_COLORS["scrna"], MODALITY_COLORS["bulk"]]
    leg_angles = [(45, 90), (0, 45)]
    # Inner ring — Disease vs Control (disease red + control gray)
    r_in, r_out = r_starts[0], r_ends[0]
    dc_demo = [DISEASE_COLOR, CONTROL_GRAY]
    for ang_idx, (th1, th2) in enumerate(leg_angles):
        ax_leg.add_patch(Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                               facecolor=dc_demo[ang_idx], edgecolor=BG_COLOR, lw=2.0))
    r_mid = (r_in + r_out) / 2.0
    ax_leg.text(-0.06, r_mid, "Condition", ha="right", va="center",
                fontsize=LABEL_SIZE, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)

    # Modality (Datasets) ring
    r_in, r_out = r_starts[1], r_ends[1]
    for ang_idx, (th1, th2) in enumerate(leg_angles):
        c = blend_white(leg_colors[ang_idx], MOD_SHADE)
        ax_leg.add_patch(Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                               facecolor=c, edgecolor=BG_COLOR, lw=2.0))
    r_mid = (r_in + r_out) / 2.0
    ax_leg.text(-0.06, r_mid, "Modality", ha="right", va="center",
                fontsize=LABEL_SIZE, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)

    # Species ring
    r_in, r_out = r_starts[2], r_ends[2]
    sp_legend_colors = [BINARY_COLORS["Human"], BINARY_COLORS["Mouse"]]
    for ang_idx, (th1, th2) in enumerate(leg_angles):
        ax_leg.add_patch(Wedge((cx, cy), r_out, th1, th2, width=r_out - r_in,
                               facecolor=sp_legend_colors[ang_idx], edgecolor=BG_COLOR, lw=2.0))
    r_mid = (r_in + r_out) / 2.0
    ax_leg.text(-0.06, r_mid, "Species", ha="right", va="center",
                fontsize=LABEL_SIZE, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_mid, r_mid], color="black", lw=1.0)

    # Outer bars legend
    ax_leg.text(-0.06, r_ends[-1] + 0.2, "Datasets", ha="right", va="center",
                fontsize=LABEL_SIZE, color=TEXT_DARK, fontweight="medium")
    ax_leg.plot([-0.04, -0.01], [r_ends[-1] + 0.2, r_ends[-1] + 0.2], color="black", lw=1.0)

    # Dummy outer bars (varying heights to illustrate proportional sizing)
    dummy_lengths = [0.42, 0.25, 0.10, 0.38, 0.15, 0.30]
    bar_idx = 0
    for ang, c_idx in zip([30, 60], [1, 0]):
        for offset in [-6, 0, 6]:
            length = dummy_lengths[bar_idx]
            bar_idx += 1
            w = Wedge((cx, cy), r_ends[-1] + length + 0.02, ang+offset-1.8, ang+offset+1.8,
                      width=length, facecolor="#B0B0B0", edgecolor=BG_COLOR, lw=0.8, alpha=0.55)
            ax_leg.add_patch(w)

    # Radial axis line in legend (at right edge of dummy bars)
    leg_axis_ang = 70
    leg_axis_rad = np.radians(leg_axis_ang)
    leg_r0 = r_ends[-1] + 0.02
    leg_r1 = r_ends[-1] + 0.02 + BAR_LEN_MAX
    ax_leg.plot([leg_r0 * np.cos(leg_axis_rad) + cx, leg_r1 * np.cos(leg_axis_rad) + cx],
                [leg_r0 * np.sin(leg_axis_rad) + cy, leg_r1 * np.sin(leg_axis_rad) + cy],
                color="#888888", linewidth=0.5)
    # "max" label at top of legend axis
    ax_leg.text((leg_r1 + 0.04) * np.cos(leg_axis_rad) + cx,
                (leg_r1 + 0.04) * np.sin(leg_axis_rad) + cy,
                "max", fontsize=LABEL_SIZE, color="#888888", ha="center", va="center")

    # Bottom swatch legends — two columns: Species (left), Disease/Control (right)
    sp_y = -0.5
    ax_leg.text(0, sp_y + 0.15, "Species", fontsize=LABEL_SIZE, fontweight="bold", color=TEXT_DARK)
    for key, color in BINARY_COLORS.items():
        ax_leg.add_patch(plt.Rectangle((0, sp_y - 0.05), 0.1, 0.1, facecolor=color, edgecolor=BG_COLOR))
        ax_leg.text(0.15, sp_y, key, fontsize=LABEL_SIZE, color=TEXT_DARK, va="center")
        sp_y -= 0.18

    dc_x, dc_y = 0.62, -0.5
    ax_leg.text(dc_x, dc_y + 0.15, "Condition", fontsize=LABEL_SIZE, fontweight="bold", color=TEXT_DARK)
    for lab, col in [("Disease", DISEASE_COLOR), ("Control", CONTROL_GRAY)]:
        ax_leg.add_patch(plt.Rectangle((dc_x, dc_y - 0.05), 0.1, 0.1, facecolor=col, edgecolor=BG_COLOR))
        ax_leg.text(dc_x + 0.15, dc_y, lab, fontsize=LABEL_SIZE, color=TEXT_DARK, va="center")
        dc_y -= 0.18

    # Ancestry legend (GWAS inner ring) — 2×2 grid under the Species column
    ax_leg.text(0, -0.95, "Ancestry", fontsize=LABEL_SIZE, fontweight="bold", color=TEXT_DARK)
    anc_pos = [(0.0, -1.10), (0.30, -1.10), (0.0, -1.26), (0.30, -1.26)]
    for (label, _n, col), (ax_, ay_) in zip(ANCESTRY["gwas"], anc_pos):
        ax_leg.add_patch(plt.Rectangle((ax_, ay_ - 0.05), 0.09, 0.1, facecolor=col, edgecolor=BG_COLOR))
        ax_leg.text(ax_ + 0.12, ay_, label, fontsize=LABEL_SIZE - 1, color=TEXT_DARK, va="center")

    # ── Save ─────────────────────────────────────────────────────────────────
    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)))),
        "figures", "main", "fig1_atlas_overview", "panels",
    )
    os.makedirs(out_dir, exist_ok=True)

    pdf_path = os.path.join(out_dir, "fig1b.pdf")

    fig.savefig(pdf_path, bbox_inches="tight", dpi=600, facecolor=BG_COLOR)
    plt.close(fig)
    print(f"Saved: {pdf_path}")

if __name__ == "__main__":
    main()
