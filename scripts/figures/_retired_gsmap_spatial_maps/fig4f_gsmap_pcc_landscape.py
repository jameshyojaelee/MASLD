#!/usr/bin/env python
"""Fig 4f (main) — gsMap spatially narrows GWAS risk genes across all MASLD traits.

Comprehensive across-trait view built on gsMap's native gene-narrowing statistic
PCC (across-spots correlation of a gene's spatial specificity GSS with the trait's
spatial heritability -log10P; gsMap/diagnosis.py). Gene x phenotype heatmap:
  rows  = GWAS/COLOC risk genes gsMap spatially narrows to (top by mean PCC),
  cols  = the 6 MASLD phenotypes (18 gsMap GWAS clumped across studies),
  color = spatial concordance (PCC) in the GSE192741 Visium cohort (5 sections).
Genes bright across many columns = pan-MASLD spatially-supported risk genes.
Prioritized (multi-evidence) genes are marked.

Why GSE192741 only: the Vu cohort's PCC magnitudes are systematically lower
(different Visium platform), so a pooled color scale would be dominated by cohort,
not gene; the cross-cohort replication lives in Fig S4 / the phenotype forest.

Input (15l): Analysis/Spatial/results/gsmap/gsmap_pcc_by_trait.csv
Output: figures/main/fig4_validation/panels/fig4f_gsmap_risk_in_tissue.pdf
Env: rnaseq or spatial (pandas/numpy/matplotlib).
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
GS = os.path.join(BASE, "Analysis/Spatial/results/gsmap")
FIG_OUT = os.path.join(BASE, "figures/main/fig4_validation/panels/fig4f_gsmap_risk_in_tissue.pdf")

COHORT = "gse192741"
N_GENES = int(os.environ.get("FIG4F_NGENES", "28"))
CMAP = "magma"
COL_PRIOR = "#C9265E"

PHENO = {
    "ukbb_alt": "ALT", "mvp_alt": "ALT",
    "ukbb_ast": "AST", "mvp_ast": "AST",
    "ukbb_ggt": "GGT",
    "finngen_nafld": "NAFLD", "ghodsian_nafld": "NAFLD", "mvp_nafld": "NAFLD",
    "decode_nafld": "NAFLD", "ukbb2023_nafld": "NAFLD", "intermtn_nafld": "NAFLD",
    "anstee2020_nafld": "NAFLD", "nafld_2019": "NAFLD",
    "finngen_nash": "NASH",
    "pdff": "PDFF", "pdff_2021a": "PDFF", "pdff_2021b": "PDFF", "pdff_2022": "PDFF",
}
PHENO_ORDER = ["ALT", "AST", "GGT", "NAFLD", "NASH", "PDFF"]


def main():
    d = pd.read_csv(os.path.join(GS, "gsmap_pcc_by_trait.csv"))
    d = d[d.cohort == COHORT].copy()
    d["phenotype"] = d["trait"].map(PHENO)

    # per (gene, phenotype): mean PCC over the phenotype's traits (studies)
    gp = (d.groupby(["gene", "phenotype"])
            .agg(pcc=("pcc_mean", "mean"),
                 is_prioritized=("is_prioritized", "max"),
                 is_genetic=("is_genetic", "max"))
            .reset_index())

    mat = gp.pivot(index="gene", columns="phenotype", values="pcc").reindex(columns=PHENO_ORDER)
    flags = gp.groupby("gene")[["is_prioritized", "is_genetic"]].max()

    # rows = GWAS/COLOC (genetic) genes, ranked by mean PCC across phenotypes
    genetic_genes = flags.index[flags["is_genetic"]]
    mat = mat.loc[mat.index.isin(genetic_genes)]
    mat["mean_pcc"] = mat[PHENO_ORDER].mean(axis=1)
    mat = mat.sort_values("mean_pcc", ascending=False).head(N_GENES)
    genes = mat.index.tolist()
    M = mat[PHENO_ORDER].values                      # top at row 0

    prio = flags.loc[genes, "is_prioritized"].values

    # ── heatmap ──
    nrow, ncol = M.shape
    fig, ax = plt.subplots(figsize=(2.9, 0.135 * nrow + 0.75))
    vmax = np.nanpercentile(M, 98)
    im = ax.imshow(M, aspect="auto", cmap=CMAP, vmin=0, vmax=vmax)

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(PHENO_ORDER, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(genes, fontstyle="italic", fontsize=5.6)
    ax.tick_params(length=1.5, width=0.4)
    # family separators (enzymes | diagnosis | imaging)
    for x in (2.5, 4.5):
        ax.axvline(x, color="white", linewidth=1.1)

    # mark prioritized (multi-evidence) genes with a dot left of the row
    for i, p in enumerate(prio):
        if p:
            ax.plot(-0.85, i, marker="o", ms=2.2, color=COL_PRIOR,
                    clip_on=False, zorder=5)
    ax.set_xlim(-1.3, ncol - 0.5)

    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("gsMap spatial concordance (PCC)", size=6)
    cb.ax.tick_params(labelsize=5, width=0.4, length=2)
    cb.outline.set_linewidth(0.3)

    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=COL_PRIOR,
                              markersize=3, label="also prioritized")],
              loc="lower left", bbox_to_anchor=(0, 1.01), fontsize=5,
              frameon=False, handletextpad=0.3)

    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)

    fig.savefig(FIG_OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig4f] saved: {FIG_OUT}  ({nrow} genes x {ncol} phenotypes)")
    print(f"[fig4f] top narrowed GWAS-risk genes: {', '.join(genes[:12])} ...")
    print(f"[fig4f] prioritized among shown: {int(prio.sum())}/{nrow}")

    print("\nCAPTION (Fig 4f): gsMap spatially narrows GWAS risk genes across MASLD "
          "traits. Heatmap of gsMap spatial concordance (PCC — the across-spots correlation "
          "between a gene's spatial specificity and a trait's spatial heritability −log10P) in "
          f"the GSE192741 Visium cohort (5 sections), for the {nrow} GWAS/COLOC risk genes with "
          "the highest mean PCC (rows) across the 6 MASLD phenotypes (columns; 18 GWAS clumped "
          "across studies). Bright = the gene's tissue expression tracks where that trait's genetic "
          "risk concentrates. Genes bright across many phenotypes are pan-MASLD spatially-supported "
          "risk genes; magenta dots mark genes also in the multi-evidence prioritized set. Single "
          "spatial annotation (condition) precludes per-niche claims; cross-cohort replication in "
          "Fig S4.")


if __name__ == "__main__":
    main()
