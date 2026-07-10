#!/usr/bin/env python3
"""
Fig 1 Panel C3 — Canonical-driver evidence portfolio (PROTOTYPE -> figures/misc/)
KEY MESSAGE: The atlas recovers the field's established MASLD genes, each supported by a
             DISTINCT combination of independent sources -- none by all -- so multi-modal
             integration is required to capture them (face validity + the necessity argument).

HONEST FRAMING (data-driven, 2026-05-29): the convergence *ranking* does NOT place canonical
drivers at the apex (THRB rank ~3091, PNPLA3 ~5409; score saturates), so this is NOT a
"top-of-ranking" panel. It is a curated recovery/portfolio panel. Per-source support uses the
correct metric for each channel:
  - Genetic (S2): COLOC-validated  (tier startswith '1_Genetic_validated') -- the S2-INTACT
                  log-BF is conservative for COLOC-only hits (THRB/TM6SF2 read ~0 there).
  - All other channels: calibrated log Bayes factor > Jeffreys log(3) = 1.099.

Data: RNA-seq/results/multi_evidence/convergence_evidence.csv
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import csv, math, os

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.grid": False, "legend.frameon": False,
})

THR = math.log(3)
ABSENT = "#D9D9D9"

# canonical Okabe-Ito palette (single source of truth)
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import SOURCE_COLORS

# 6 sources shown (essentiality dropped — empty for these drivers); colors from locked palette
SOURCES = [
    ("S1", "log_BF_S1",     "Human\nbulk",      SOURCE_COLORS["S1"]),
    ("S8", "log_BF_S8",     "Mouse\nbulk",      SOURCE_COLORS["S8"]),
    ("S2", "GENETIC",       "Genetic\n(COLOC)", SOURCE_COLORS["S2"]),
    ("S4", "log_BF_S4",     "Epigen-\nomic",    SOURCE_COLORS["S4"]),
    ("S5", "log_BF_S5",     "Spatial",          SOURCE_COLORS["S5"]),
    ("S7", "log_BF_S7",     "Proteo-\nmics",    SOURCE_COLORS["S7"]),
]

DRIVERS = ["THRB", "PNPLA3", "TM6SF2", "GCKR", "MBOAT7", "HSD17B13", "MTTP",
           "HNF4A", "HKDC1", "NR1H4", "PPARA", "PPARG", "TREM2", "GPAM", "DGAT2", "FABP4"]


def _root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load():
    f = os.path.join(_root(), "RNA-seq/results/multi_evidence/convergence_evidence.csv")
    rows = {r["human_symbol"]: r for r in csv.DictReader(open(f))}
    out = {}
    for g in DRIVERS:
        r = rows.get(g)
        if r is None:
            continue
        present = {}
        genetic = r.get("tier", "").startswith("1_Genetic_validated")
        for code, col, _, _ in SOURCES:
            if col == "GENETIC":
                present[code] = genetic
            else:
                try:
                    present[code] = float(r.get(col, "nan")) > THR
                except ValueError:
                    present[code] = False
        out[g] = present
    return out


def main():
    data = load()
    genes = sorted(data.keys(), key=lambda g: -sum(data[g].values()))  # portfolio breadth desc
    n_g, n_s = len(genes), len(SOURCES)

    fig, ax = plt.subplots(figsize=(4.0, 4.9))
    fig.patch.set_facecolor("white")

    BOX = 0.66
    for gi, g in enumerate(genes):
        y = n_g - 1 - gi
        if gi % 2 == 0:
            ax.axhspan(y - 0.5, y + 0.5, color="#F7F7F7", zorder=0)
        for si, (code, _, _, color) in enumerate(SOURCES):
            if data[g][code]:
                ax.add_patch(mpatches.FancyBboxPatch((si - BOX / 2, y - BOX / 2), BOX, BOX,
                             boxstyle="round,pad=0.02", facecolor=color, edgecolor="white",
                             lw=0.8, alpha=0.92, zorder=2))
            else:
                ax.add_patch(mpatches.Circle((si, y), 0.07, facecolor="none",
                             edgecolor=ABSENT, lw=1.0, zorder=2))
        ax.text(-0.75, y, g, ha="right", va="center", fontsize=6,
                color="#2D3436", fontstyle="italic")

    # source column headers, colored
    for si, (code, _, lab, color) in enumerate(SOURCES):
        ax.text(si, n_g - 0.32, lab, ha="center", va="bottom", fontsize=6,
                color=color, linespacing=0.9)

    ax.set_xlim(-2.0, n_s - 0.45)
    ax.set_ylim(-0.7, n_g + 0.6)
    ax.axis("off")

    print("[caption] Canonical MASLD drivers recovered across complementary sources")

    fig.savefig(os.path.join(_root(), "figures/misc/fig1_C3_driver_portfolio.pdf"),
                bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print("Saved figures/misc/fig1_C3_driver_portfolio.pdf  (%d drivers, %d sources)" % (n_g, n_s))


if __name__ == "__main__":
    main()
