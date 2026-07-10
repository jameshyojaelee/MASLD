#!/usr/bin/env python
"""Fig 4f-i — MASLD trait families dissociate their spatial genetic risk onto
different cell types (depth-corrected, replicated across two Visium cohorts).

gsMap per-spot GWAS-heritability −log10P is strongly depth-confounded, so this
uses the DEPTH-CORRECTED partial correlation (control total_counts) between
per-spot risk and each cell type's cell2location proportion (from 15m). Headline:
  NAFLD-diagnosis risk -> Hepatocytes (anti-Fibroblast)
  PDFF (liver-fat)     risk -> Fibroblasts (anti-Hepatocyte)
The hepatocyte<->fibroblast crossover replicates in GSE192741 and Vu (15 sections).
Liver-enzyme traits (weak hepatocyte lean, no dissociation) are the reference.

Input:  Analysis/Spatial/results/gsmap/gsmap_celltype_localization.csv  (15m)
Output: figures/main/fig4_validation/panels/fig4f_gsmap_celltype_dissociation.pdf
Env: rnaseq or spatial (matplotlib/pandas).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 6,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.4, "xtick.major.width": 0.4, "ytick.major.width": 0.4,
})

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CSV = os.path.join(BASE, "Analysis/Spatial/results/gsmap/gsmap_celltype_localization.csv")
OUT = os.path.join(BASE, "figures/main/fig4_validation/panels/fig4f_gsmap_celltype_dissociation.pdf")

COH = {"gse192741": ("GSE192741", "#1565C0"), "vu": ("Vu et al. 2025", "#C9265E")}
# rows top->bottom, grouped by cell type
GROUPS = [("Hepatocytes", ["NAFLD", "PDFF"]),
          ("Fibroblasts", ["NAFLD", "PDFF"])]
TRAIT_LAB = {"NAFLD": "NAFLD", "PDFF": "liver fat (PDFF)"}


def main():
    d = pd.read_csv(CSV)
    d["se"] = d["sd"] / np.sqrt(d["n_sections"])

    # build row layout
    rows = []          # (y, cell_type, trait)
    y = 0.0
    yticks, ylabs, group_spans = [], [], []
    for ct, traits in GROUPS:
        g0 = y
        for tr in traits:
            rows.append((y, ct, tr)); yticks.append(y); ylabs.append(TRAIT_LAB[tr])
            y -= 1.0
        group_spans.append((ct, g0, y + 1.0))
        y -= 0.7        # gap between cell-type groups

    fig, ax = plt.subplots(figsize=(3.5, 2.2))
    ax.axvline(0, color="#9E9E9E", lw=0.4, ls="--", zorder=1)
    dodge = 0.17
    for (yy, ct, tr) in rows:
        for i, (coh, (lab, col)) in enumerate(COH.items()):
            r = d[(d.cohort == coh) & (d.cell_type == ct) & (d.trait_family == tr)]
            if r.empty:
                continue
            m, se = r["mean_partial_r"].iloc[0], r["se"].iloc[0]
            yp = yy + (dodge if i == 0 else -dodge)
            ax.plot([m - se, m + se], [yp, yp], color=col, lw=0.9, zorder=2)
            ax.plot(m, yp, "o", color=col, ms=3.2, zorder=3)

    ax.set_yticks(yticks); ax.set_yticklabels(ylabs, fontsize=6)
    ax.set_ylim(min(yy for yy, _, _ in rows) - 0.5, max(yy for yy, _, _ in rows) + 0.5)
    # cell-type group labels at left (rotated), + faint separators
    xlo = d["mean_partial_r"].min() - 0.06
    for ct, y0, y1 in group_spans:
        ax.text(-0.34, (y0 + y1) / 2, ct, transform=ax.get_yaxis_transform(),
                rotation=90, va="center", ha="center", fontsize=6.2, color="black")
    ax.set_xlabel("depth-corrected partial r\n(gsMap risk vs cell-type abundance)", fontsize=6)
    ax.tick_params(labelsize=6, width=0.4, length=2)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(min(-0.28, xlo), 0.30)

    ax.legend(handles=[Line2D([0], [0], marker="o", color=c, lw=0, markersize=3.2, label=l)
                       for l, c in COH.values()],
              loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=5,
              frameon=False, handletextpad=0.3, columnspacing=1.0)

    fig.subplots_adjust(left=0.32, right=0.97, top=0.86, bottom=0.24)
    fig.savefig(OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig4f-i] saved: {OUT}")

    # report + caption
    print("\n[fig4f-i] values (mean partial r | cross-cohort concordant):")
    for ct, traits in GROUPS:
        for tr in traits:
            g = d[(d.cohort == "gse192741") & (d.cell_type == ct) & (d.trait_family == tr)]
            v = d[(d.cohort == "vu") & (d.cell_type == ct) & (d.trait_family == tr)]
            print(f"  {ct:12} {tr:6} GSE {g.mean_partial_r.iloc[0]:+.2f}  Vu {v.mean_partial_r.iloc[0]:+.2f}  "
                  f"{'concordant' if bool(g.cross_cohort_concordant.iloc[0]) else 'DISCORDANT'}")

    print("\nCAPTION (Fig 4f-i): MASLD trait families dissociate their spatial genetic risk onto "
          "different cell types. Depth-corrected partial Spearman correlation (control per-spot "
          "total counts) between gsMap per-spot GWAS-heritability −log10P and cell2location cell-type "
          "proportion, per section, mean ± SE across sections, in two independent Visium cohorts "
          "(GSE192741 n=5; Vu et al. 2025 n=10). Most liver-trait GWAS load onto a hepatocyte↔"
          "fibroblast axis; NAFLD-diagnosis (and enzyme) risk associates with hepatocytes and against "
          "fibroblasts, whereas liver-fat (PDFF) risk uniquely REVERSES sign — associating with "
          "fibroblasts and against hepatocytes. Both patterns exceed a 1000× spot-permutation null "
          "(Stouffer cross-cohort P≈4.5×10⁻⁶ each) and STRENGTHEN under double-partialling that also "
          "controls the two traits' shared polygenic −log10P (PDFF→Fib 0.15→0.19; NAFLD→Hep 0.19→0.23), "
          "so the crossover is not an artifact of their shared steatosis loci. Consistent with the gene "
          "level: fibroblast/ECM genes (COL1A1/COL3A1/ACTA2/PDGFRA) are the most reproducibly flagged "
          "PDFF spatial-risk genes, hepatocyte genes (HNF4A/APOB/TTR) the NAFLD ones. Raw −log10P is "
          "depth-confounded (r≈0.7 with total counts) so only the depth-corrected association is shown; "
          "effect sizes are modest (r 0.11–0.23). CAVEATS: the pattern is directional not PDFF-exclusive "
          "(deCODE-NAFLD also reverses); the four PDFF GWAS are overlapping UKBB releases (one signal, "
          "two tissue cohorts); and the cell2location reference signature is shared across cohorts.")


if __name__ == "__main__":
    main()
