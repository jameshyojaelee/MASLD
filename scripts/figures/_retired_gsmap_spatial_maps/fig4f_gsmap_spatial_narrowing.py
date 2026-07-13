#!/usr/bin/env python
"""Fig 4f — gsMap narrows MASLD GWAS risk genes using spatial context.

Three-part composite that shows gsMap's actual workflow, not a set-level OR:
  (a) Spatial heritability map: a steatotic GSE192741 Visium section (H&E) with
      each spot colored by gsMap -log10P of GWAS (NAFLD) heritability enrichment
      -> genetic risk is spatially STRUCTURED in the lobule.
  (b) gsMap-narrowed risk genes: among GWAS/COLOC candidate genes, the top ones
      by PCC (across-spots correlation of a gene's spatial specificity GSS with
      the trait's -log10P) -> spatial context narrows to specific causal genes.
  (c) GSS colocalization: the top narrowed gene's spatial specificity on the same
      section overlaps the heritability hotspot (gsMap's signature diagnosis).

Data (all on disk from 15e + 15l; no recompute of gsMap):
  report/{trait}/gsMap_plot/{sample}_{trait}_gsMap_plot.csv   (per-spot logp)
  Analysis/Spatial/results/gsmap/gsmap_pcc_by_trait.csv        (15l: PCC + flags)
  latent_to_gene/{sample}_gene_marker_score.feather           (GSS gene x spot)
  results/spaceranger/GSE192741/{JBO*}/outs/                   (H&E + coords)

Honesty note: the run used --annotation condition (single-valued), so no per-niche
claims; the spatial -log10P map + PCC ranking are annotation-independent and valid.

Output: figures/main/fig4_validation/panels/fig4f_gsmap_risk_in_tissue.pdf
Env: spatial (scanpy). Run on a compute node.
"""
import json
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import gridspec
import scanpy as sc

sc.settings.verbosity = 0
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 6,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.4, "xtick.major.width": 0.4, "ytick.major.width": 0.4,
})

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SR = os.path.join(BASE, "Analysis/Spatial/results/spaceranger")
GS = os.path.join(BASE, "Analysis/Spatial/results/gsmap")
FIG_OUT = os.path.join(BASE, "figures/main/fig4_validation/panels/figS4_gsmap_spatial_narrowing.pdf")

# ── configuration (chosen after focality + PCC verification) ──
SAMPLE = os.environ.get("FIG4F_SAMPLE", "gse192741_JBO019")   # steatotic, join 0.98
TRAIT = os.environ.get("FIG4F_TRAIT", "mvp_nafld")
TRAIT_LABEL = "NAFLD"                    # short, biomarker-agnostic
COHORT = "gse192741"
N_GENES = 10                             # narrowed genes shown
CMAP = "magma"
COL_PRIOR = "#C9265E"                    # prioritized+COLOC marker (Liang magenta)


# ── helpers (self-contained; adapted from fig4_spatial_candidates.py) ──
def strip_bc(bcs):
    out = []
    for b in bcs:
        m = re.match(r"^([ACGTN]{16}-\d+)", str(b))
        out.append(m.group(1) if m else str(b))
    return out


def load_he(gsmap_sample):
    """H&E + fullres coords for a GSE192741 sample. gsMap sample name carries the
    'gse192741_' prefix; spaceranger dirs do not."""
    sr_sample = gsmap_sample.replace("gse192741_", "")
    outs = os.path.join(SR, "GSE192741", sr_sample, "outs")
    ad = sc.read_10x_h5(os.path.join(outs, "filtered_feature_bc_matrix.h5"))
    ad.var_names_make_unique()
    spdir = os.path.join(outs, "spatial")
    pos_path = next(os.path.join(spdir, fn) for fn in
                    ("tissue_positions_list.csv", "tissue_positions.csv")
                    if os.path.exists(os.path.join(spdir, fn)))
    hdr = 0 if pos_path.endswith("tissue_positions.csv") else None
    pos = pd.read_csv(pos_path, header=hdr)
    if hdr is None:
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]
    pos = pos.set_index("barcode").reindex(ad.obs_names)
    sf = json.load(open(os.path.join(spdir, "scalefactors_json.json")))
    img = plt.imread(os.path.join(spdir, "tissue_hires_image.png"))
    ad.obs["x"] = pos["pxl_col_in_fullres"].values
    ad.obs["y"] = pos["pxl_row_in_fullres"].values
    return ad, img, sf["tissue_hires_scalef"]


def per_spot_series(sample, trait, col):
    """Return a barcode-indexed series (col from gsMap_plot.csv)."""
    f = os.path.join(GS, COHORT, sample, "report", trait, "gsMap_plot",
                     f"{sample}_{trait}_gsMap_plot.csv")
    d = pd.read_csv(f, index_col=0)
    s = d[col].copy()
    s.index = strip_bc(s.index)
    return s


def gss_row(sample, gene):
    """GSS values (barcode-indexed) for one gene from the marker-score feather."""
    fe = os.path.join(GS, COHORT, sample, "latent_to_gene",
                      f"{sample}_gene_marker_score.feather")
    mk = pd.read_feather(fe)
    row = mk.loc[mk["HUMAN_GENE_SYM"] == gene]
    if row.empty:
        return None
    vals = row.drop(columns=["HUMAN_GENE_SYM"]).iloc[0]
    vals.index = strip_bc(vals.index)
    return vals.astype(float)


def draw_spatial(ax, ad, img, scalef, values, cbar_label, panel_label="",
                 label_italic=False, vlo=2, vhi=98):
    bc = strip_bc(ad.obs_names)
    v = pd.Series(values).reindex(bc).values.astype(float)
    ax.imshow(img, origin="upper")
    px, py = ad.obs["x"].values * scalef, ad.obs["y"].values * scalef
    finite = v[np.isfinite(v)]
    vmin, vmax = np.nanpercentile(finite, [vlo, vhi])
    rng = np.nanmax(px) - np.nanmin(px)
    ssz = max((2.2 * 72 / rng * (rng / np.sqrt(len(px))) * 0.62) ** 2 * np.pi, 0.5)
    sca = ax.scatter(px, py, c=v, s=ssz, cmap=CMAP, vmin=vmin, vmax=vmax,
                     linewidths=0, rasterized=True)
    # crop tightly to the tissue (drop the grey fiducial frame)
    pad = 0.04 * rng
    ax.set_xlim(np.nanmin(px) - pad, np.nanmax(px) + pad)
    ax.set_ylim(np.nanmax(py) + pad, np.nanmin(py) - pad)   # inverted (origin upper)
    ax.set_aspect("equal"); ax.axis("off")
    cb = ax.figure.colorbar(sca, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label(cbar_label, size=6)
    cb.ax.tick_params(labelsize=5, width=0.4, length=2)
    cb.outline.set_linewidth(0.3)
    if panel_label:
        ax.set_title(panel_label, fontsize=6.5, pad=2, color="black",
                     fontstyle=("italic" if label_italic else "normal"))
    return np.isfinite(v).mean()


def main():
    # ── data: narrowed genes (top genetic/COLOC by PCC) ──
    pcc = pd.read_csv(os.path.join(GS, "gsmap_pcc_by_trait.csv"))
    sub = pcc[(pcc.trait == TRAIT) & (pcc.cohort == COHORT)].copy()
    genetic = sub[sub.is_genetic].sort_values("pcc_mean", ascending=False).head(N_GENES)
    genetic = genetic.iloc[::-1]                       # ascending for horizontal plot (top at top)
    top_gene = sub[sub.is_genetic].sort_values("pcc_mean", ascending=False).iloc[0]["gene"]

    # discrimination stat for caption
    from scipy.stats import mannwhitneyu
    gen_all = sub.loc[sub.is_genetic, "pcc_mean"]; bg = sub.loc[~sub.is_genetic, "pcc_mean"]
    mwu_p = mannwhitneyu(gen_all, bg, alternative="greater").pvalue

    med_all = float(sub["pcc_mean"].median())          # typical-gene reference

    # ── figure scaffold (2 panels: map + narrowed genes) ──
    fig = plt.figure(figsize=(4.9, 2.55))
    gsp = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.95], wspace=0.55,
                            left=0.02, right=0.94, top=0.9, bottom=0.16)
    ax_map = fig.add_subplot(gsp[0])
    ax_gene = fig.add_subplot(gsp[1])

    # (a) heritability map
    ad, img, sf = load_he(SAMPLE)
    logp = per_spot_series(SAMPLE, TRAIT, "logp")
    cov_a = draw_spatial(ax_map, ad, img, sf, logp,
                         cbar_label="−log$_{10}$P", panel_label=TRAIT_LABEL)

    # (b) narrowed genes — horizontal dots by PCC, referenced to the typical gene
    y = np.arange(len(genetic))
    ax_gene.axvline(med_all, color="#9E9E9E", linewidth=0.5, linestyle="--", zorder=1)
    ax_gene.hlines(y, med_all, genetic["pcc_mean"], color="#E0AEC0", linewidth=0.6, zorder=2)
    ax_gene.scatter(genetic["pcc_mean"], y, s=11, color=COL_PRIOR, zorder=3, linewidths=0)
    ax_gene.set_yticks(y)
    ax_gene.set_yticklabels(genetic["gene"], fontstyle="italic", fontsize=6)
    ax_gene.set_ylim(-0.6, len(genetic) - 0.4)
    ax_gene.set_xlim(0, max(0.62, genetic["pcc_mean"].max() * 1.08))
    ax_gene.set_xlabel("gsMap spatial concordance (PCC)", fontsize=6)
    ax_gene.tick_params(labelsize=6, width=0.4, length=2)
    ax_gene.annotate("median\ngene", xy=(med_all, len(genetic) - 0.5),
                     xytext=(med_all + 0.02, len(genetic) - 1.4), fontsize=5,
                     color="#616161", ha="left", va="top")
    for sp in ("top", "right"):
        ax_gene.spines[sp].set_visible(False)

    fig.savefig(FIG_OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[figS4] saved: {FIG_OUT}")
    print(f"[figS4] sample={SAMPLE} trait={TRAIT} map_cov={cov_a:.2f}")
    print(f"[figS4] narrowed genetic/COLOC genes (by PCC): "
          f"{', '.join(genetic['gene'].iloc[::-1])}")

    print("\nCAPTION (Fig S4): gsMap narrows MASLD GWAS risk genes using spatial "
          f"context (worked example, one GWAS). LEFT — a steatotic GSE192741 Visium section "
          f"({SAMPLE.replace('gse192741_','')}) with each spot colored by gsMap −log10P of "
          f"{TRAIT_LABEL} ({TRAIT}) GWAS heritability enrichment; genetic risk is spatially "
          "structured across the tissue. RIGHT — among GWAS/COLOC candidate genes, the top "
          f"{N_GENES} by PCC, the across-spots correlation between a gene's spatial specificity "
          "(GSS) and the trait's spatial heritability, i.e. the genes whose tissue expression "
          "tracks where risk concentrates (dashed line, genome-wide median gene). Genetic-risk "
          f"genes have significantly higher PCC than background (Mann-Whitney P={mwu_p:.1e}). The "
          "run used a single spatial annotation (condition), so no per-niche localization is "
          "claimed. The comprehensive across-trait view is main Fig 4f.")


if __name__ == "__main__":
    main()
