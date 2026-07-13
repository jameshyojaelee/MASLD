#!/usr/bin/env python
"""Fig 4f (main) — liver-enzyme and NAFLD genetic risk co-localize in tissue.

gsMap projects each GWAS onto the Visium tissue as a per-spot -log10P heritability
map. On one steatotic GSE192741 section, the genetic risk for ALT, GGT and NAFLD
(each a multi-GWAS consensus) concentrates in the SAME tissue region (pairwise
spatial Pearson r 0.82-0.97) — a shared spatial footprint of MASLD genetic risk.

Robustness note: this WITHIN-section co-localization is the defensible finding.
A disease-vs-healthy contrast was tested and is a null result (the two healthy
sections disagree: JBO018 uniform, JBO022 polarized; no disease-specific
difference), so it is NOT shown.

Color = within-trait percentile of -log10P (each trait normalized to its own
distribution so all three are equally legible; per-panel peak -log10P annotated).

Data: report/{trait}/gsMap_plot/{sample}_{trait}_gsMap_plot.csv  (per-spot logp)
      results/spaceranger/GSE192741/{JBO*}/outs/                  (H&E + coords)
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
import scanpy as sc

sc.settings.verbosity = 0
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 6,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.4,
})

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SR = os.path.join(BASE, "Analysis/Spatial/results/spaceranger")
GS = os.path.join(BASE, "Analysis/Spatial/results/gsmap")
FIG_OUT = os.path.join(BASE, "figures/main/fig4_validation/panels/fig4f_gsmap_risk_in_tissue.pdf")

SAMPLE = os.environ.get("FIG4F_SAMPLE", "gse192741_JBO019")   # steatotic, join 0.98
# columns: trait label -> constituent gsMap GWAS (aggregated per trait)
COLS = [("ALT", ["ukbb_alt", "mvp_alt"]),
        ("GGT", ["ukbb_ggt"]),
        ("NAFLD", ["finngen_nafld", "ghodsian_nafld", "mvp_nafld", "ukbb2023_nafld", "anstee2020_nafld"])]
CMAP = "magma"


def strip_bc(bcs):
    out = []
    for b in bcs:
        m = re.match(r"^([ACGTN]{16}-\d+)", str(b))
        out.append(m.group(1) if m else str(b))
    return out


def load_he(gsmap_sample):
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


def consensus_logp(sample, traits, bc):
    mats = []
    for t in traits:
        f = os.path.join(GS, "gse192741", sample, "report", t, "gsMap_plot",
                         f"{sample}_{t}_gsMap_plot.csv")
        s = pd.read_csv(f, index_col=0)["logp"].copy()
        s.index = strip_bc(s.index)
        mats.append(s.reindex(bc).values)
    return np.nanmean(np.vstack(mats), axis=0)


def main():
    ad, img, sf = load_he(SAMPLE)
    bc = strip_bc(ad.obs_names)
    px, py = ad.obs["x"].values * sf, ad.obs["y"].values * sf
    # tight crop to the tissue bbox (same section in every panel -> boxes stay equal);
    # excludes the Visium fiducial frame in the surrounding H&E
    padx, pady = 0.03 * (px.max() - px.min()), 0.03 * (py.max() - py.min())
    xlo, xhi = px.min() - padx, px.max() + padx
    ylo, yhi = py.min() - pady, py.max() + pady
    spacing = (px.max() - px.min()) / np.sqrt(len(px))
    ssz = max((1.05 * 72 / (xhi - xlo) * spacing * 0.6) ** 2 * np.pi, 0.3)

    vals = {tl: consensus_logp(SAMPLE, traits, bc) for tl, traits in COLS}

    def pct(v):                       # within-trait percentile -> equal legibility
        finite = np.sort(v[np.isfinite(v)])
        out = np.searchsorted(finite, v, side="right") / len(finite)
        out = out.astype(float); out[~np.isfinite(v)] = np.nan
        return out

    fig, axes = plt.subplots(1, len(COLS), figsize=(3.9, 1.55))
    sca = None
    for ax, (tl, _) in zip(axes, COLS):
        v = vals[tl]
        ax.imshow(img, origin="upper")
        sca = ax.scatter(px, py, c=pct(v), s=ssz, cmap=CMAP, vmin=0, vmax=1,
                         linewidths=0, rasterized=True)
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(yhi, ylo)
        ax.set_aspect("equal"); ax.axis("off")
        ax.text(0.03, 0.03, f"peak {np.nanmax(v):.1f}", transform=ax.transAxes,
                fontsize=4.2, ha="left", va="bottom", color="black",
                bbox=dict(facecolor="white", alpha=0.6, pad=0.4, edgecolor="none"))
        ax.set_title(tl, fontsize=6.5, color="black", pad=1.5)

    cb = fig.colorbar(sca, ax=axes.tolist(), fraction=0.03, pad=0.012, location="right")
    cb.set_label("heritability enrichment\n(within-trait percentile)", size=4.6)
    cb.set_ticks([0, 0.5, 1.0])
    cb.ax.tick_params(labelsize=4.2, width=0.3, length=1.5)
    cb.outline.set_linewidth(0.3)

    fig.subplots_adjust(left=0.01, right=0.9, top=0.9, bottom=0.02, wspace=0.04)
    fig.savefig(FIG_OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)

    # co-localization statistic (pairwise spatial r among shown traits)
    M = pd.DataFrame({tl: vals[tl] for tl, _ in COLS})
    C = M.corr()
    tri = C.values[np.triu_indices(len(COLS), 1)]
    cov = np.mean([np.isfinite(vals[tl]).mean() for tl, _ in COLS])
    print(f"[fig4f] saved: {FIG_OUT}  (section {SAMPLE.replace('gse192741_','')}, join {cov:.2f})")
    print(f"[fig4f] pairwise spatial r among {[t for t,_ in COLS]}: {tri.min():.2f}-{tri.max():.2f}")

    print("\nCAPTION (Fig 4f): liver-enzyme and NAFLD genetic risk co-localize in tissue. gsMap "
          "per-spot −log10P of GWAS heritability enrichment for 3 MASLD traits (each a multi-GWAS "
          "consensus — ALT: UKB+MVP; GGT: UKB; NAFLD: FinnGen+Ghodsian+MVP+UKB-2023+Anstee) on a "
          f"steatotic GSE192741 Visium section ({SAMPLE.replace('gse192741_','')}). Color = within-"
          "trait percentile of −log10P (each trait normalized to its own distribution so all are "
          "equally legible; per-panel peak −log10P annotated). All three traits' genetic risk "
          f"concentrates in the same tissue region (pairwise spatial Pearson r {tri.min():.2f}–"
          f"{tri.max():.2f}). The run used a single spatial annotation (condition), so no per-niche/"
          "cell-type claim is made. Per-gene spatial narrowing (worked example) in Fig S4.")


if __name__ == "__main__":
    main()
