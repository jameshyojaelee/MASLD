#!/usr/bin/env python
"""Fig S4l — spatial GWAS-risk is convergence-specific (non-circularity control).

Companion to the cell-type dissociation panel. gsMap spatial-risk-gene enrichment
of the Fig4 prioritized set, DECOMPOSED into evidence arms (15h ->
arm_decomposition_spatial_risk.csv). The prioritized set = transcriptomic arm ∪
genetic arm; gsMap is itself genetics-based, so the circularity control is: only
the CONVERGENT (transcriptomic ∩ genetic) genes are enriched, while
transcriptomic-only is a flat null and genetic-only is DEPLETED. So the signal is
neither a genetics tautology (genetic-alone depletes) nor genetics-independent
(transcriptomic-alone null) — it is transcriptomic × genetic convergence.

Summarized over the liver-enzyme traits (ALT/AST/GGT), median OR + across-trait
range, per Visium cohort.

Input:  Analysis/Spatial/results/gsmap/arm_decomposition_spatial_risk.csv
Output: figures/main/fig4_validation/panels/figS4l.pdf
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
CSV = os.path.join(BASE, "Analysis/Spatial/results/gsmap/arm_decomposition_spatial_risk.csv")
OUT = os.path.join(BASE, "figures/main/fig4_validation/panels/figS4l.pdf")

ENZYME = ["ukbb_alt", "ukbb_ast", "ukbb_ggt"]
ARMS = ["Convergent", "Transcriptomic-only", "Genetic-only"]   # top -> bottom
ARM_LAB = {"Convergent": "Convergent\n(T×G)", "Transcriptomic-only": "Transcriptomic-only",
           "Genetic-only": "Genetic-only"}
COH = {"gse192741": ("GSE192741", "#1565C0"), "vu": ("Vu et al. 2025", "#C9265E")}


def main():
    d = pd.read_csv(CSV)
    d = d[d["trait"].isin(ENZYME) & d["group"].isin(ARMS)]

    fig, ax = plt.subplots(figsize=(3.4, 1.9))
    ax.axvline(1, color="#9E9E9E", lw=0.4, ls="--", zorder=1)
    ymap = {arm: -i for i, arm in enumerate(ARMS)}
    dodge = 0.16
    for arm in ARMS:
        for i, (coh, (lab, col)) in enumerate(COH.items()):
            s = d[(d.group == arm) & (d.cohort == coh)]["fisher_or"]
            if s.empty:
                continue
            med, lo, hi = np.median(s), s.min(), s.max()
            yp = ymap[arm] + (dodge if i == 0 else -dodge)
            ax.plot([lo, hi], [yp, yp], color=col, lw=0.9, zorder=2)
            ax.plot(med, yp, "o", color=col, ms=3.4, zorder=3)

    ax.set_yticks(list(ymap.values()))
    ax.set_yticklabels([ARM_LAB[a] for a in ARMS], fontsize=6)
    ax.set_ylim(min(ymap.values()) - 0.5, max(ymap.values()) + 0.5)
    ax.set_xscale("log")
    ax.set_xticks([0.4, 0.6, 1.0, 1.6])
    ax.set_xticklabels(["0.4", "0.6", "1.0", "1.6"])
    ax.set_xlim(0.38, 2.0)
    ax.set_xlabel("spatial GWAS-risk enrichment OR\n(prioritized-gene arm, ALT/AST/GGT)", fontsize=6)
    ax.tick_params(labelsize=6, width=0.4, length=2)
    ax.tick_params(axis="y", length=0)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)

    ax.legend(handles=[Line2D([0], [0], marker="o", color=c, lw=0, markersize=3.2, label=l)
                       for l, c in COH.values()],
              loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=5,
              frameon=False, handletextpad=0.3, columnspacing=1.0)

    fig.subplots_adjust(left=0.30, right=0.97, top=0.85, bottom=0.28)
    fig.savefig(OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig4f-ii] saved: {OUT}")

    print("\n[fig4f-ii] median OR (range) over ALT/AST/GGT by arm x cohort:")
    for arm in ARMS:
        for coh, (lab, _) in COH.items():
            s = d[(d.group == arm) & (d.cohort == coh)]["fisher_or"]
            print(f"  {arm:20} {lab:14} {np.median(s):.2f} ({s.min():.2f}-{s.max():.2f})")

    print("\nCAPTION (Fig S4l): spatial GWAS-risk is convergence-specific. Enrichment (Fisher OR) "
          "of the Fig4 prioritized gene set among gsMap spatial GWAS-risk genes, decomposed into "
          "evidence arms and summarized over liver-enzyme GWAS (ALT/AST/GGT; median OR, across-trait "
          "range) in two Visium cohorts. Only Convergent genes (transcriptomic ∩ genetic, n=1,244) "
          "are enriched (both cohorts); Transcriptomic-only (n=6,844) is a flat null and Genetic-only "
          "(n=1,794) is significantly DEPLETED (OR≈0.5). Because gsMap is itself genetics-based, this "
          "rules out both a genetics tautology (genetic-alone depletes) and a genetics-independent "
          "artifact (transcriptomic-alone null): spatial genetic risk lands specifically on the "
          "transcriptomic×genetic convergent targets. Dashed line, OR=1.")


if __name__ == "__main__":
    main()
