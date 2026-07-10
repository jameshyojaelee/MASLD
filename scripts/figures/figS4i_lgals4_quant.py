#!/usr/bin/env python
"""
figS4i — LGALS4 cross-assay concordance (companion to figS4h spatial images).

Single panel: LGALS4 log2FC (MASLD vs control) across four independent assays, all human,
all using the CANONICAL POOLED (cohort-adjusted) analyses — NOT per-cohort breakdowns:

  Bulk RNA-seq (pooled)   +0.97  padj 7.6e-12   (846 samples / 5 control-bearing cohorts, UP both sexes)
  Proteomics (DIA-MS)     +1.64  padj 0.050
  scRNA hepatocyte        +0.81  n.s.
  Spatial (GeoMx)         +0.53  n.s.

Bars shaded by significance (dark = padj<0.05, light = n.s.); * marks padj<0.05.
Gene name "LGALS4" centered at top.

Deliberately NOT shown:
  - Individual RNA-seq cohorts — the canonical analysis is the POOLED cohort-adjusted result;
    per-cohort bars misrepresent it and the per-study contrasts are not a uniform D-vs-C.
  - PRJNA512027 (Gerhard) — excluded from the paper for the L0/S0 library-prep batch confound.
  - Visium violin — its spatial signal is in figS4h's images (would duplicate the GeoMx bar).
  - COLOC — atlas 0.92 is an OFF-TRAIT platelet-GWAS colocalization (MVP_Platelet_AFR);
    LGALS4 does not colocalize with any MASLD/liver GWAS (PP.H4 approx 0).

Convention: 6 pt text, black, no bold, gene symbol italic, caption to stdout. up/sig = #E65100.

Output: figures/main/fig4_validation/panels/figS4i.pdf
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE  = os.environ.get("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ATLAS = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
OUT   = os.path.join(BASE, "figures/main/fig4_validation/panels/figS4i.pdf")
GENE  = "LGALS4"
DARK, LIGHT = "#E65100", "#FFB74D"   # significant / n.s.

plt.rcParams.update({
    "font.size": 6, "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5, "savefig.dpi": 300,
})

a = pd.read_csv(ATLAS, low_memory=False)
r = a[a.human_symbol == GENE].iloc[0]
def num(c): return float(r[c])
assays = [  # (label, log2FC, padj) — all pooled / canonical
    ("Bulk\nRNA-seq",    num("bulk_logFC"),                              num("bulk_padj")),
    ("Proteo-\nmics",    num("best_protein_logFC"),                      num("best_protein_padj")),
    ("scRNA\nhep.",      num("sc_hepatocyte_logFC"),                     num("sc_hepatocyte_padj")),
    ("Spatial\n(GeoMx)", num("spatial_govaere2026_geomx_sh_vs_pt_logfc"),num("spatial_govaere2026_geomx_sh_vs_pt_padj")),
]
labs = [x[0] for x in assays]; lfcs = [x[1] for x in assays]; padjs = [x[2] for x in assays]

fig, ax = plt.subplots(figsize=(2.5, 2.0))
xs = np.arange(len(assays))
ax.bar(xs, lfcs, width=0.66, color=[DARK if p < 0.05 else LIGHT for p in padjs])
ax.axhline(0, color="gray", lw=0.4)
for x, lf, p in zip(xs, lfcs, padjs):
    if p < 0.05:
        ax.text(x, lf + max(lfcs) * 0.03, "*", ha="center", va="bottom", fontsize=7)
ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=5.2)
ax.set_ylabel("log$_2$FC (MASLD vs control)", fontsize=6)
ax.set_ylim(0, max(lfcs) * 1.22)
ax.tick_params(axis="y", labelsize=6, length=2, width=0.4); ax.tick_params(axis="x", length=0)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.set_title(GENE, fontstyle="italic", fontsize=8, color="black", pad=3)  # gene label, top-center

fig.subplots_adjust(left=0.20, right=0.97, top=0.86, bottom=0.20)
fig.savefig(OUT, bbox_inches="tight")
print(f"saved {OUT}")
print(f"[caption] {GENE} cross-assay concordance (human, pooled/canonical). LGALS4 log2FC "
      f"(MASLD vs control): bulk RNA-seq {assays[0][1]:+.2f} (padj {assays[0][2]:.1e}, pooled 846 "
      f"samples/5 cohorts, UP both sexes), proteomics {assays[1][1]:+.2f} (padj {assays[1][2]:.3f}), "
      f"scRNA hepatocyte {assays[2][1]:+.2f} (n.s.), GeoMx spatial {assays[3][1]:+.2f} (n.s.). "
      f"Dark = padj<0.05 (*). Spatial Visium images in figS4h. No per-cohort bars (canonical = pooled), "
      f"no COLOC (atlas 0.92 is off-trait platelet GWAS).")
