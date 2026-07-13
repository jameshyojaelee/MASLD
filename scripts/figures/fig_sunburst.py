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
    "bulk":    {"human": (1124, 160)},
    "scrna":   {"human": (196, 64)},
    "spatial": {"human": (22, 5)},
    "atac":    {"human": (13, 5)},
    "proteo":  {"human": (111, 19)},
    "pharma":  None,
    "gwas":    None,
}
CONTROL_GRAY  = "#9E9E9E"  # canonical control gray (FIGURE_GUIDELINES)
DISEASE_COLOR = "#C0524E"  # single semantic disease red — deliberately NOT a modality/blue hue
NA_GRAY       = "#E2E2E2"  # neutral fill for non-case/control modalities
SPECIES_LIGHTEN = 0.45     # mouse arcs = this much lighter than human (lightness = species)

# GWAS inner ring = ANCESTRY of the 50 GWAS strata (by dataset count). The 3 eQTL
# (GTEx/Broadaway/sc-eQTL, EUR reference panels) are EXCLUDED — not ancestry-stratified
# GWAS. Violet gradient (dark EUR → pale SAS). Source: gwas_registry.tsv (MVP 50-GWAS
# portfolio, 2026-07-05: EUR 21 + AFR 10 + EAS 9 + AMR 7 + SAS 3 = 50; AMR is the new
# MVP ancestry). Counts verified on disk against gwas_registry.tsv.
ANCESTRY = {
    "gwas": [("EUR", 21, "#3A2259"), ("AFR", 10, "#5F3A8F"),
             ("EAS", 9, "#8B5FC0"), ("AMR", 7, "#B592D4"),
             ("SAS", 3, "#E4CDEE")],
}

# (key, title, l1, l2, (n_human, n_mouse), weight, n_bars)
# Verified from unified metadata & processed data (GWAS updated 2026-07-08 to MVP 50-GWAS portfolio)
# scRNA: human donors + mouse donors across 7 source datasets
# GWAS: 50 GWAS strata (EUR 21 + AFR 10 + EAS 9 + AMR 7 + SAS 3) + 3 eQTL (GTEx v8, Broadaway, sc-eQTL)
#   Registry: GWAS/finemapping/config/gwas_registry.tsv (MVP R4 added 2026-07-05; 23 legacy + 27 MVP strata
#   that colocalized). 50 datasets = 48 distinct studies (Sveinbjornsson split into deCODE/Intermountain/UKBB).
#   Fig 2 uses the same 50-GWAS portfolio (35 are Tier-1/2 liver-specific). eQTL EXCLUDED from ancestry ring.
#   n_bars=13 (NOT 50): the outer per-dataset bars are AGGREGATED BY TRAIT (10 traits: NAFLD/NASH/
#   PDFF/ALT/AST/GGT + MVP-added Cirrhosis/ChronLiver/Albumin/Platelet) + 3 eQTL bars, 2026-07-08 —
#   53 individual-study bars degraded into an illegible sub-pixel fringe at this angular weight; trait
#   aggregation keeps the bars legible while the l1 text label still states the true 50-study count.
# bulk: 9 human QC-pass (1,259); 9 datasets (mouse removed 2026-06-25)
# pharma: LINCS L1000, network proximity, DGIdb/OT, ClinicalTrials MASH pipeline, multi-layer (5 sources)
# spatial: Govaere 2026 (GeoMx+CosMx, integrated 2026-05-21), GSE192741 (Visium), Vu_et_al_2025 (3 active)
SECTORS = [
    ("bulk",    "RNA-seq",                    "9\nRNA-seq",       "1,259\nsamples",  (1259, 0),    18.0, 9),
    ("scrna",   "scRNA-seq",                  "7\nscRNA-seq",     "269\ndonors",     (269, 0),      8.0, 7),
    ("spatial", "Spatial Omics",              "3\nSpatial",       "27\nsamples",     (27, 0),        4.5, 3),
    ("atac",    "ATAC-seq",                   "1\nATAC-seq",      "18\nsamples",     (18, 0),        4.5, 1),
    ("proteo",  "Proteomics",                 "2\nProteomics",    "130\nsamples",    (130, 0),       4.5, 2),
    ("pharma",  "Clinical",                   "5\nClinical",      "1,107\ncmpds",    (1, 0),         4.5, 5),
    ("gwas",    "Genomics",                   "50 GWAS\n3 eQTL",  "1.5M+\nsubj.",   (1, 0),        16.0, 13),
]

# Per-dataset sample counts for proportional outer-bar heights
_DATASET_SIZES_RAW = {
    "bulk": [  # 9 Human QC-pass from sample_qc_report.csv; STAR -s2 canonical 2026-05-28 (mouse removed 2026-06-25)
        358, 215, 160, 142, 97, 93, 76, 64, 55,       # Human QC-pass: Chen(GSE213621), Suppli(GSE135251), Govaere(GSE193066), Bril(GSE162694), Kawamura(GSE167523), Kozumi(GSE174478), Hoshida(GSE130970), Verschuren(GSE240729), Hoang(GSE126848)
    ],
    "scrna": [  # 7 source datasets (unique donors; SRR technical runs excluded)
        38, 117, 67, 21, 20, 2, 4,  # Liver_Atlas(38 donors), GSE244832, GSE202379, GSE185477, GSE136103, GSE189600(2 human donors), GSE174748
    ],
    "spatial": [12, 10, 5],  # Govaere2026(GeoMx 8pt+CosMx 4pt), Vu_et_al_2025(10 arrays), GSE192741(5 Visium sections)
    "atac": [18],                       # Human_Multiome (18 donors); Mouse_Bulk removed 2026-06-25
    "proteo": [72, 58],                 # PXD052937 plasma DIA-MS, PXD051911 liver DIA-MS (GSE276114 removed 2026-06: GEO confirms bulk RNA-seq, not SomaScan proteomics)
    "pharma": [1107, 1173, 58, 23, 20],  # LINCS(1107 compounds), network proximity(1173 screened), DGIdb/OT, ClinicalTrials MASH, multi-layer
    "gwas": [  # AGGREGATED BY TRAIT (N_tot summed across all 5 ancestries per trait, in thousands;
               # from gwas_registry.tsv, MVP 50-GWAS portfolio 2026-07-05). 10 traits + 3 eQTL = 13
               # bars — the prior 53 individual-study bars were sub-pixel/illegible at this angular
               # weight. Verified: trait sums total 8,253.719K, exactly matching the registry's full
               # 50-row N_tot sum (binary 4,310.369K + quantitative 3,943.350K), so no study double-
               # counted or dropped. Direct Tier-1/2 liver traits first, then MVP-added Tier-3/4 traits.
        2637.014,  # NAFLD (11 studies, all 5 ancestries; incl. Ghodsian 778.6K meta)
        435.000,   # NASH (1 study: FinnGen)
        113.841,   # PDFF (3 studies, EUR only)
        1103.344,  # ALT (8 studies, all 5 ancestries)
        1081.068,  # AST (8 studies, all 5 ancestries)
        519.362,   # GGT (4 studies: UKBB/BBJ/PanUKBB-AFR/PanUKBB-CSA, no MVP GGT stratum)
        626.711,   # Cirrhosis (3 MVP studies: EUR/AFR/AMR — no EAS stratum)
        611.644,   # ChronLiver (4 MVP studies, all non-SAS ancestries)
        542.276,   # Albumin (4 MVP studies, all non-SAS ancestries)
        583.459,   # Platelet (4 MVP studies, all non-SAS ancestries)
        # eQTL (3, in thousands of donors) — unchanged, not trait-aggregated
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
    # (Species ring removed 2026-06-25 — human-only atlas, so the flag was constant.
    #  The old inner per-modality "Counts" totals ring was removed 2026-06-18 — it
    #  duplicated the outer bars; replaced by disease/control composition.)
    r_starts = [0.31, 0.62]
    r_ends   = [0.62, 0.95]

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

        # (Species ring removed 2026-06-25 — human-only atlas.)
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

            # Per-sector TOTAL label — modality-level n (from l2), CENTERED on the
            # sector mid-angle just beyond the tallest bar, matching the modality-name
            # orientation so it reads as that wedge's total. Shows the modality total
            # (not the max single dataset); bars stay proportional to the per-modality
            # max as a relative-size texture. GWAS/Clinical totals are descriptive
            # (GWAS subjects overlap biobanks; compound sources not disjoint) — NOT sums.
            label_ang = (current_angle + end_angle) / 2.0
            label_rad = np.radians(label_ang)
            r_label = r_base + BAR_LEN_MAX + 0.11
            lab_x = r_label * np.cos(label_rad)
            lab_y = r_label * np.sin(label_rad)
            rot = label_ang - 90
            if rot < -90:
                rot += 180
            elif rot > 90:
                rot -= 180
            ax.text(lab_x, lab_y, l2, rotation=rot,
                    fontsize=LABEL_SIZE, color=TEXT_DARK, fontweight="medium",
                    ha="center", va="center", zorder=6, linespacing=0.9)

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

    # Species ring removed 2026-06-25 (human-only atlas).

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

    # (Bar axis + "max" annotation removed 2026-06-25 — per-modality totals are
    #  labeled directly on each sector; bars just show relative dataset sizes.)

    # Bottom swatch legends — two columns: Species (left), Disease/Control (right)
    # Species swatch removed 2026-06-25 (human-only atlas).
    dc_x, dc_y = 0.0, -0.5
    ax_leg.text(dc_x, dc_y + 0.15, "Condition", fontsize=LABEL_SIZE, fontweight="bold", color=TEXT_DARK)
    for lab, col in [("Disease", DISEASE_COLOR), ("Control", CONTROL_GRAY)]:
        ax_leg.add_patch(plt.Rectangle((dc_x, dc_y - 0.05), 0.1, 0.1, facecolor=col, edgecolor=BG_COLOR))
        ax_leg.text(dc_x + 0.15, dc_y, lab, fontsize=LABEL_SIZE, color=TEXT_DARK, va="center")
        dc_y -= 0.18

    # Ancestry legend (GWAS inner ring) — 5 ancestries (3 top / 2 bottom) under the Condition column
    ax_leg.text(0, -0.95, "Ancestry", fontsize=LABEL_SIZE, fontweight="bold", color=TEXT_DARK)
    anc_pos = [(0.0, -1.10), (0.30, -1.10), (0.60, -1.10), (0.0, -1.26), (0.30, -1.26)]
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
