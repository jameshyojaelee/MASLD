#!/usr/bin/env python
"""Supplementary exploratory gsMap cell-type dissociation panel.

MASLD trait groups partition their spatial genetic risk across cell types
(grouped bars; two Visium cohorts pooled).

Bar = does a trait's genetic risk concentrate in spots where a cell type lives?
i.e. the spatial partial correlation between gsMap per-spot GWAS-heritability
−log10P and each cell type's cell2location abundance, corrected for (a) sequencing
depth and (b) compositionality (CLR). + = risk concentrates in that cell type,
− = avoids it. Estimate pools all 15 sections (GSE192741 5 + Vu 10); the ± is the
pooled SE; cross-cohort sign-concordance is annotated.

Story: liver-enzyme & NAFLD/NASH (injury/diagnosis, BLUES) risk → Hepatocytes &
Macrophages, away from Fibroblasts & B cells; liver-fat PDFF (MAGENTA) REVERSES on
the hepatocyte↔fibroblast axis.

Input:  Analysis/Spatial/results/gsmap/gsmap_celltype_localization.csv  (15m, CLR)
Output: figures/main/fig5_molecular_context/panels/figS4_gsmap_celltype_dissociation_exploratory.pdf
Env: rnaseq or spatial (matplotlib/pandas).
Demoted from main Figure 4 on 2026-08-06 because these descriptive partial
correlations do not validate the frozen Figure 3 programs.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 6,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.4, "xtick.major.width": 0.4, "ytick.major.width": 0.4,
})

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CSV = os.path.join(BASE, "Analysis/Spatial/results/gsmap/gsmap_celltype_localization.csv")
OUT = os.path.join(
    BASE,
    "figures/main/fig5_molecular_context/panels/figS4_gsmap_celltype_dissociation_exploratory.pdf",
)

CELLS = ["Hepatocytes", "Fibroblasts", "Macrophages"]           # cell2location atlas labels
# biologically-precise display names (verified by reference markers):
#   Fibroblasts = hepatic stellate cells (RGS5/PDGFRB/ACTA2/COL1A1/PDGFRA)
#   Macrophages = Kupffer cells (MARCO/VSIG4/TIMD4/CD5L/FOLR2); monocyte-derived are separate
DISPLAY = {"Hepatocytes": "Hepatocytes", "Fibroblasts": "Stellate cells",
           "Macrophages": "Kupffer cells"}
# trait groups (left->right within a cell) + colors: injury/diagnosis = blues, PDFF = magenta
TRAITS = [("Liver enzyme", "enzyme", "#90CAF9"),
          ("NAFLD/NASH", "NAFLD/NASH", "#1565C0"),
          ("PDFF", "liver fat", "#C9265E")]


def pooled(g, v):
    """Pool two per-cohort summaries (mean, pop-sd, n) into mean + SE over all sections."""
    n1, m1, s1 = g["n_sections"], g["mean_partial_r"], g["sd"]
    n2, m2, s2 = v["n_sections"], v["mean_partial_r"], v["sd"]
    N = n1 + n2
    mp = (n1 * m1 + n2 * m2) / N
    ss = n1 * s1**2 + n2 * s2**2 + n1 * (m1 - mp)**2 + n2 * (m2 - mp)**2   # 15m sd is ddof=0
    se = np.sqrt((ss / N) / N)
    return mp, se


def main():
    d = pd.read_csv(CSV)

    fig, ax = plt.subplots(figsize=(2.55, 1.75))
    ax.axhline(0, color="black", lw=0.4, zorder=3)
    nT = len(TRAITS)
    bw = 0.8 / nT
    xc = np.arange(len(CELLS))
    for j, (tr, lab, col) in enumerate(TRAITS):
        xs, ys, es, hatch = [], [], [], []
        for i, ct in enumerate(CELLS):
            g = d[(d.cohort == "gse192741") & (d.cell_type == ct) & (d.trait_family == tr)]
            v = d[(d.cohort == "vu") & (d.cell_type == ct) & (d.trait_family == tr)]
            if g.empty or v.empty:
                continue
            m, se = pooled(g.iloc[0], v.iloc[0])
            concord = bool(g["cross_cohort_concordant"].iloc[0])
            x = xc[i] + (j - (nT - 1) / 2) * bw
            xs.append(x); ys.append(m); es.append(se)
            hatch.append(None if concord else "////")
        bars = ax.bar(xs, ys, bw * 0.92, color=col, edgecolor="white", linewidth=0.3,
                      zorder=2, label=lab)
        for b, h in zip(bars, hatch):
            if h:                                    # non-concordant across cohorts
                b.set_hatch(h); b.set_alpha(0.45); b.set_edgecolor(col)
        ax.errorbar(xs, ys, yerr=es, fmt="none", ecolor="#424242", elinewidth=0.5,
                    capsize=1.2, zorder=4)

    ax.set_xticks(xc)
    ax.set_xticklabels([DISPLAY.get(c, c) for c in CELLS], fontsize=5.5)
    ax.set_ylabel("genetic-risk enrichment\n(spatial partial r)", fontsize=5.5)
    ax.tick_params(labelsize=5.2, width=0.35, length=1.8, pad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_ylim(-0.26, 0.26)
    ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])

    ax.legend(handles=[Patch(facecolor=c, label=l) for _, l, c in TRAITS],
              loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=3, fontsize=4.6,
              frameon=False, handlelength=0.8, handleheight=0.8, columnspacing=0.7,
              handletextpad=0.25, borderaxespad=0.0)

    fig.subplots_adjust(left=0.17, right=0.99, top=0.90, bottom=0.11)
    fig.savefig(OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig4f-i] saved: {OUT}")

    print("\n[fig4f-i] pooled estimate (mean ± SE over 15 sections | cross-cohort concordant):")
    for ct in CELLS:
        for tr, lab, _ in TRAITS:
            g = d[(d.cohort == "gse192741") & (d.cell_type == ct) & (d.trait_family == tr)].iloc[0]
            v = d[(d.cohort == "vu") & (d.cell_type == ct) & (d.trait_family == tr)].iloc[0]
            m, se = pooled(g, v)
            print(f"  {ct:12} {tr:13} {m:+.3f} ± {se:.3f}  "
                  f"{'concord' if bool(g['cross_cohort_concordant']) else 'n.c.'}")

    print("\nCAPTION (Fig 4f-i): MASLD trait groups partition their spatial genetic risk across liver "
          "cell types. Each bar is the spatial partial correlation between gsMap per-spot GWAS-"
          "heritability −log10P and a cell type's cell2location abundance — corrected for sequencing "
          "depth and for compositionality (centered-log-ratio) — pooled across two Visium cohorts "
          "(GSE192741 n=5 + Vu et al. 2025 n=10 sections; error bars pooled SE). Positive = a trait's "
          "genetic risk concentrates in spots rich in that cell type, negative = avoids them. Cell-type "
          "annotation (cell2location reference; top cluster-specific genes): Hepatocytes ALB, APOB, TTR, "
          "CYP2E1, HNF4A, CPS1, ASGR1; Stellate cells (atlas ‘Fibroblasts’) RGS5, PDGFRB, ACTA2, COL1A1, "
          "COL3A1, DCN, PDGFRA; Kupffer cells (atlas ‘Macrophages’) MARCO, VSIG4, TIMD4, CD5L, FOLR2, "
          "SIGLEC1, C1QA/B/C — monocyte-derived macrophages form a separate cluster; excluded sinusoidal "
          "endothelium/LSEC (atlas ‘Endothelial cells’) STAB2, CLEC4G, FCN2, OIT3. Liver-enzyme and "
          "NAFLD/NASH-diagnosis risk (blues) concentrate in hepatocytes and Kupffer cells and are "
          "depleted from stellate cells; liver-fat (PDFF, magenta) reverses the hepatocyte↔stellate "
          "axis. Associations shown are directionally concordant across both cohorts (hatched/faded = "
          "not concordant); the PDFF→stellate reversal additionally survives a 1000× spot-permutation "
          "null (Stouffer P≈4.5×10⁻⁶) and control for the traits' shared polygenic signal. Sinusoidal "
          "endothelial (LSEC), cholangiocyte and the sparse immune cell types were tested but excluded "
          "as not cross-cohort concordant. Effect sizes are modest (|r| 0.05–0.23); claims rest on the "
          "compositional treatment, cross-cohort replication and sign-consistency. The PDFF GWAS are "
          "overlapping UKBB releases (one signal, two tissue cohorts) and the cell2location reference "
          "is shared across cohorts.")


if __name__ == "__main__":
    main()
