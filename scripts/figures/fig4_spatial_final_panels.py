#!/usr/bin/env python3
"""
Fig 4 spatial — finalize the chosen candidate panels to publication spec,
written into the spatial panel dir:  figures/main/fig4_validation/panels/spatial/

  1. fig4_spatial_cyp3a4_he      — CYP3A4 over H&E, Healthy (JBO018) | Steatotic
                                    (JBO019), shared color scale (drug-target loss)
  2. fig4_spatial_cosmx_celltype — Govaere CosMx single-cell cell-type architecture
                                    (Leuven_1)
  3. fig4_spatial_zonationscore  — periportal→pericentral zonation backbone over H&E
                                    (JBO018)
  4. fig4_spatial_serpine1_fibro — SERPINE1 expression | Fibroblast abundance,
                                    Steatotic (JBO019), clean scatter (fibroblast-
                                    ligand reorganization)

Reuses loaders/constants from fig4_spatial_candidates.py so the data path is
identical to the gallery. Env: spatial. Run on a compute node (sbatch).

House rules: PDF only, one panel per PDF, all text black, short/italic gene
titles, control/healthy gray #9E9E9E, no subtitles/lollipops/3D.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import TwoSlopeNorm
import scanpy as sc

sc.settings.verbosity = 0
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.5, "savefig.dpi": 300,
})

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig4_spatial_candidates as G   # loaders + constants (module-level is safe)

BASE = G.BASE
SPA  = G.SPA
OUT  = os.path.join(BASE, "figures/main/fig4_validation/panels/spatial")
os.makedirs(OUT, exist_ok=True)

CMAP_EXPR  = G.CMAP_EXPR     # magma
CMAP_ABUND = G.CMAP_ABUND    # rocket/magma
CMAP_ZON   = G.CMAP_ZON      # RdBu_r
C2L_PREFIX = G.C2L_PREFIX
CAT_PAL    = G.CAT_PAL
ZONE_COLORS = G.ZONE_COLORS  # PP1 blue ... PC1 orange
ZONE_ORDER  = G.ZONE_ORDER
load_visium_he = G.load_visium_he

# Single-column publication width
W1 = 88 / 25.4        # 3.46 in
W2 = 120 / 25.4       # 4.72 in (two-up)


def _gene_vec(ad, g):
    x = ad[:, g].X
    return np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()


def _style_axis(ax):
    ax.set_aspect("equal"); ax.axis("off")


def _tissue_bbox(px, py, pad=0.04):
    """Tight bounding box around the spots (hires-px), with fractional padding."""
    x0, x1 = np.nanmin(px), np.nanmax(px)
    y0, y1 = np.nanmin(py), np.nanmax(py)
    dx, dy = (x1 - x0) * pad, (y1 - y0) * pad
    return x0 - dx, x1 + dx, y0 - dy, y1 + dy


def _apply_crop(ax, bbox):
    x0, x1, y0, y1 = bbox
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)          # origin='upper' → small y at top


def _match_bbox_aspect(bbox_a, bbox_b):
    """Pad the shorter dimension of each bbox (about its own center) so both
    share the wider of the two aspect ratios. Two independently-cropped tissue
    sections otherwise have different width/height ratios, so `set_aspect
    ('equal')` letterboxes them by different amounts and the two boxes render
    at visibly different sizes side by side — this equalizes them without
    cropping either section."""
    def aspect(b):
        x0, x1, y0, y1 = b
        return (x1 - x0) / (y1 - y0)

    target = max(aspect(bbox_a), aspect(bbox_b))

    def pad(b):
        x0, x1, y0, y1 = b
        w, h = x1 - x0, y1 - y0
        if w / h < target - 1e-9:
            new_w = h * target
            cx = (x0 + x1) / 2
            return (cx - new_w / 2, cx + new_w / 2, y0, y1)
        return b

    return pad(bbox_a), pad(bbox_b)


def _he_inset(ax, V, bbox, loc=(0.015, 0.015, 0.30, 0.30)):
    """Bare cropped-H&E thumbnail (no spots) as a bordered corner inset, so the
    histology is visible without spot occlusion."""
    ins = ax.inset_axes(loc)
    ins.imshow(V["img"], origin="upper")
    _apply_crop(ins, bbox)
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.6); sp.set_edgecolor("black")
    return ins


def _he_xy(V):
    sf = V["scalef"]
    return (V["adata"].obs["x_fullres"].values * sf,
            V["adata"].obs["y_fullres"].values * sf)


def _spot_size(x, s_scale):
    rng = (np.nanmax(x) - np.nanmin(x)) or 1.0
    radius = (rng / np.sqrt(len(x))) * (2.4 * 72 / rng) * s_scale
    return max(np.pi * radius ** 2, 0.3)


def _scatter_he(ax, V, c, *, cmap, vmin=None, vmax=None, norm=None,
                s_scale=0.50, alpha=0.88, crop=True, inset=True, bbox=None):
    ax.imshow(V["img"], origin="upper")
    x, y = _he_xy(V)
    kw = dict(s=_spot_size(x, s_scale), linewidths=0, alpha=alpha, rasterized=True)
    if norm is not None:
        sca = ax.scatter(x, y, c=c, cmap=cmap, norm=norm, **kw)
    else:
        sca = ax.scatter(x, y, c=c, cmap=cmap, vmin=vmin, vmax=vmax, **kw)
    bbox = bbox if bbox is not None else _tissue_bbox(x, y)
    _style_axis(ax)
    if crop:
        _apply_crop(ax, bbox)
    if inset:
        _he_inset(ax, V, bbox)
    return sca


def _scatter_he_cat(ax, V, labels, color_map, *, s_scale=0.50, alpha=0.92,
                    crop=True, inset=True, na="#E8E8E8"):
    ax.imshow(V["img"], origin="upper")
    x, y = _he_xy(V)
    cols = pd.Series(labels).astype(object).map(color_map).fillna(na).values
    ax.scatter(x, y, c=cols, s=_spot_size(x, s_scale), linewidths=0,
               alpha=alpha, rasterized=True)
    bbox = _tissue_bbox(x, y)
    _style_axis(ax)
    if crop:
        _apply_crop(ax, bbox)
    if inset:
        _he_inset(ax, V, bbox)


# ── feature table (zonation_score + c2l fibroblast) from the zonation object ──
print("[load] zonation feature object", flush=True)
_zad = sc.read_h5ad(os.path.join(SPA, "zonation/spatial_with_zonation.h5ad"))
_zobs = _zad.obs.copy(); _zobs["bc"] = _zobs.index.astype(str)
del _zad


def feat_for(sample, ad):
    # barcode-normalized join (strips anndata make_unique suffixes); see
    # G.sample_feature — without it the per-sample join matches almost nothing.
    return G.sample_feature(_zobs, "sample_id", sample, ad.obs_names.astype(str))


# ══════════════════════════════════════════════════════════════════════════════
# Panel 1 — CYP3A4 over H&E, Healthy | Steatotic, shared scale
# ══════════════════════════════════════════════════════════════════════════════
def panel_cyp3a4():
    print("[panel] CYP3A4 H&E pair", flush=True)
    Vh = load_visium_he("JBO018"); Vs = load_visium_he("JBO019")
    g = "CYP3A4"
    gh, gs = _gene_vec(Vh["adata"], g), _gene_vec(Vs["adata"], g)
    vmax = np.nanpercentile(np.concatenate([gh, gs]), 99)
    bbox_h = _tissue_bbox(*_he_xy(Vh))
    bbox_s = _tissue_bbox(*_he_xy(Vs))
    bbox_h, bbox_s = _match_bbox_aspect(bbox_h, bbox_s)
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.3))
    sca = _scatter_he(axes[0], Vh, gh, cmap=CMAP_EXPR, vmin=0, vmax=vmax, bbox=bbox_h)
    _scatter_he(axes[1], Vs, gs, cmap=CMAP_EXPR, vmin=0, vmax=vmax, bbox=bbox_s)
    cb = fig.colorbar(sca, ax=axes, fraction=0.030, pad=0.02)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
    cb.set_label(f"{g}  log1p", fontsize=6, fontstyle="italic")
    print(f"[caption] {g} expression over H&E, left=Healthy (JBO018), right=Steatotic (JBO019), shared color scale", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_cyp3a4_he.pdf"), bbox_inches="tight")
    # PROMOTED to a MAIN panel (2026-07-07); relettered e (2026-07-08, H&E leads the
    # Fig4 lineup's cyp3a4 pair, immediately followed by cyp3a4_zonation/f).
    fig.savefig(os.path.join(BASE, "figures/main/fig4_validation/panels", "fig4e_cyp3a4_he.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_cyp3a4_he.pdf + panels/fig4e_cyp3a4_he.pdf (main)", flush=True)


def panel_cyp3a4_single():
    """RECOMMENDED illustrative anchor: CYP3A4 on a single Steatotic section
    (JBO019). Frames the honest claim — CYP3A4 keeps its zonated pericentral
    expression in human MASLD liver — without an unsupported disease-difference
    claim (its spatial autocorrelation is unchanged: ΔMoran's I ~0)."""
    print("[panel] CYP3A4 single steatotic", flush=True)
    V = load_visium_he("JBO019")
    g = "CYP3A4"
    gv = _gene_vec(V["adata"], g)
    vmax = np.nanpercentile(gv[np.isfinite(gv)], 99)
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    sca = _scatter_he(ax, V, gv, cmap=CMAP_EXPR, vmin=0, vmax=vmax)
    cb = fig.colorbar(sca, ax=ax, fraction=0.040, pad=0.02)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
    cb.set_label(f"{g}  log1p", fontsize=6, fontstyle="italic")
    print(f"[caption] {g} expression over H&E, Steatotic section (JBO019)", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_cyp3a4_steatotic_he.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_cyp3a4_steatotic_he.pdf", flush=True)


def panel_visium_gene(sample, gene):
    """Generic single-section Visium gene map over H&E (crop + inset). Used for
    the multi-modal-convergent candidates (LYZ, KRT7) that also sit on the CosMx
    panel — so the same gene is validatable on BOTH spatial platforms."""
    print(f"[panel] Visium {gene} ({sample})", flush=True)
    V = load_visium_he(sample)
    if gene not in V["adata"].var_names:
        print(f"  [skip] {gene} not in {sample}", flush=True); return
    gv = _gene_vec(V["adata"], gene)
    vmax = np.nanpercentile(gv[np.isfinite(gv)], 99)
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    sca = _scatter_he(ax, V, gv, cmap=CMAP_EXPR, vmin=0, vmax=vmax)
    cb = fig.colorbar(sca, ax=ax, fraction=0.040, pad=0.02)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
    cb.set_label(f"{gene}  log1p", fontsize=6, fontstyle="italic")
    print(f"[caption] {gene} expression over H&E, {sample}", flush=True)
    fig.savefig(os.path.join(OUT, f"fig4_spatial_visium_{gene.lower()}_{sample}_he.pdf"),
                bbox_inches="tight")
    plt.close(fig)
    print(f"  saved fig4_spatial_visium_{gene.lower()}_{sample}_he.pdf", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# Panel 3 — zonation-score backbone over H&E (JBO018)
# ══════════════════════════════════════════════════════════════════════════════
def _zon_norm(z):
    """Diverging norm centered on the section's OWN median (the score is shifted
    positive because the pericentral marker set is higher-magnitude; centering on
    0 would paint ~76% of spots red). Honest about within-section relative
    zonation, not cross-section absolute level."""
    fin = z[np.isfinite(z)]
    center = float(np.median(fin))
    vmin, vmax = np.quantile(fin, [0.02, 0.98])
    return TwoSlopeNorm(vcenter=center, vmin=min(vmin, center - 1e-3), vmax=max(vmax, center + 1e-3))


def panel_zonation():
    """RECOMMENDED backbone: single Healthy (JBO018) zonation map, palette
    centered on the section median so periportal (blue) and pericentral (red) are
    balanced. Serves as the anatomical coordinate frame for the CYP3A4 panel."""
    print("[panel] zonation backbone (single Healthy)", flush=True)
    V = load_visium_he("JBO018")
    z = pd.to_numeric(feat_for("JBO018", V["adata"])["zonation_score"], errors="coerce").values
    sca = None
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    sca = _scatter_he(ax, V, z, cmap=CMAP_ZON, norm=_zon_norm(z))
    cb = fig.colorbar(sca, ax=ax, fraction=0.040, pad=0.02)
    cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
    cb.set_label("periportal to pericentral", fontsize=6)
    print("[caption] Zonation score (periportal to pericentral) over H&E, Healthy (JBO018)", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_zonationscore_he.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_zonationscore_he.pdf", flush=True)


def panel_zonation_pair():
    """HONEST paired option: Healthy | Steatotic, EACH centered on its own median
    (independent norms + colorbars). Compares spatial ORGANIZATION only — a module
    score's absolute baseline is not comparable across sections, so no
    cross-section magnitude claim is made."""
    print("[panel] zonation per-section pair", flush=True)
    Vh = load_visium_he("JBO018"); Vs = load_visium_he("JBO019")
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.4))
    for ax, V, s, lab in ((axes[0], Vh, "JBO018", "Healthy"),
                          (axes[1], Vs, "JBO019", "Steatotic")):
        z = pd.to_numeric(feat_for(s, V["adata"])["zonation_score"], errors="coerce").values
        sca = _scatter_he(ax, V, z, cmap=CMAP_ZON, norm=_zon_norm(z))
        cb = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
        cb.set_label("periportal to pericentral", fontsize=6)
    print("[caption] Zonation score over H&E, left=Healthy (JBO018), right=Steatotic (JBO019), each centered on its own section median", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_zonation_pair_he.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_zonation_pair_he.pdf", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# Panel 4 — SERPINE1 | Fibroblast abundance (Steatotic JBO019, clean)
# ══════════════════════════════════════════════════════════════════════════════
def panel_serpine1_fibro():
    print("[panel] SERPINE1 | Fibroblast", flush=True)
    V = load_visium_he("JBO019")
    ad = V["adata"]
    x = ad.obs["x_fullres"].values; y = -ad.obs["y_fullres"].values
    serp = _gene_vec(ad, "SERPINE1")
    fib = pd.to_numeric(feat_for("JBO019", ad)[C2L_PREFIX + "Fibroblasts"], errors="coerce").values
    rng = (np.nanmax(x) - np.nanmin(x)) or 1.0
    radius = (rng / np.sqrt(len(x))) * (2.0 * 72 / rng) * 0.55
    s = max(np.pi * radius ** 2, 0.3)
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.3))
    for ax, vals, cmap, title, ital in (
        (axes[0], serp, CMAP_EXPR, "SERPINE1", True),
        (axes[1], fib, CMAP_ABUND, "Fibroblast", False)):
        v = np.asarray(vals, float)
        finite = v[np.isfinite(v)]
        vmax = np.nanpercentile(finite, 99) if finite.size else 1
        sca = ax.scatter(x, y, c=v, cmap=cmap, vmin=np.nanpercentile(finite, 1) if finite.size else 0,
                         vmax=vmax, s=s, linewidths=0, rasterized=True)
        _style_axis(ax)
        cb = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.02)
        cb.ax.tick_params(labelsize=6, width=0.4, length=2); cb.outline.set_linewidth(0.3)
    print("[caption] SERPINE1 expression (left) and fibroblast abundance (right), Steatotic (JBO019)", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_serpine1_fibroblast.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_serpine1_fibroblast.pdf", flush=True)


# ══════════════════════════════════════════════════════════════════════════════
# Panel 2 — CosMx single-cell cell-type architecture (Leuven_1)
# ══════════════════════════════════════════════════════════════════════════════
def panel_cosmx():
    print("[panel] CosMx cell-type", flush=True)
    ad = sc.read_h5ad(os.path.join(BASE, "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad"))
    if "gene_symbol" in ad.var:
        ad.var_names = ad.var["gene_symbol"].astype(str).values
        ad.var_names_make_unique()
    slide = "Leuven_1"
    m = (ad.obs["sample_id"] == slide).values
    so = ad.obs[m]
    x = so["centerX_global_px"].values.astype(float)
    y = so["centerY_global_px"].values.astype(float)
    ct = so["cell_type"].astype(str).values
    uniq = sorted(pd.unique(ct))
    cc = {c: CAT_PAL[i % len(CAT_PAL)] for i, c in enumerate(uniq)}
    cols = pd.Series(ct).map(cc).values
    fig, ax = plt.subplots(figsize=(W2, 3.0))
    ax.scatter(x, -y, c=cols, s=0.6, linewidths=0, rasterized=True)
    _style_axis(ax)
    handles = [Patch(facecolor=cc[c], edgecolor="none", label=c) for c in uniq]
    ax.legend(handles=handles, fontsize=6, loc="center left", bbox_to_anchor=(1.0, 0.5),
              frameon=False, handlelength=0.9, handleheight=0.9, labelspacing=0.3,
              borderaxespad=0.1)
    print("[caption] CosMx single-cell cell-type architecture, Leuven_1", flush=True)
    fig.savefig(os.path.join(OUT, "fig4_spatial_cosmx_celltype.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("  saved fig4_spatial_cosmx_celltype.pdf", flush=True)

    # ── two-panel: cell-type | <gene> (a TARGET on an independent single-cell
    #    platform). COL1A1 = clean mesenchymal/septal pattern (fibrosis in situ);
    #    IL32 kept for the record though it is hepatocyte-broad (the macrophage
    #    niche is a sub-µm proximity effect, NOT a visible expression colocaln). ──
    if "lognorm" in ad.layers:
        ad.X = ad.layers["lognorm"]
    for gene in ["COL1A1", "LYZ", "KRT7", "IL32"]:
        if gene not in ad.var_names:
            continue
        gv = ad[m, gene].X
        gv = np.asarray(gv.todense()).ravel() if hasattr(gv, "todense") else np.asarray(gv).ravel()
        fig2, ax2 = plt.subplots(1, 2, figsize=(W2 * 1.18, 3.0))
        ax2[0].scatter(x, -y, c=cols, s=0.5, linewidths=0, rasterized=True)
        _style_axis(ax2[0])
        ax2[0].legend(handles=handles, fontsize=6, loc="upper center",
                      bbox_to_anchor=(0.5, -0.01), ncol=4, frameon=False,
                      handlelength=0.8, handleheight=0.8, columnspacing=0.8, labelspacing=0.2)
        vmax = np.nanpercentile(gv[np.isfinite(gv)], 99) if np.isfinite(gv).any() else 1.0
        sca2 = ax2[1].scatter(x, -y, c=gv, cmap=CMAP_EXPR, vmin=0, vmax=vmax,
                              s=0.5, linewidths=0, rasterized=True)
        _style_axis(ax2[1])
        cb2 = fig2.colorbar(sca2, ax=ax2[1], fraction=0.045, pad=0.02)
        cb2.ax.tick_params(labelsize=6, width=0.4, length=2); cb2.outline.set_linewidth(0.3)
        print(f"[caption] Cell type (left) and {gene} expression (right), CosMx Leuven_1", flush=True)
        fig2.savefig(os.path.join(OUT, f"fig4_spatial_cosmx_celltype_{gene.lower()}.pdf"),
                     bbox_inches="tight")
        plt.close(fig2)
        print(f"  saved fig4_spatial_cosmx_celltype_{gene.lower()}.pdf", flush=True)
    del ad


def main():
    panel_cyp3a4()
    panel_cyp3a4_single()
    panel_visium_gene("JBO019", "LYZ")     # convergent: bulk + COLOC 0.92 + proteomics
    panel_visium_gene("JBO019", "KRT7")    # convergent: bulk + COLOC + proteomics + mouse (ductular)
    panel_zonation_pair()                  # main zonation figure (per-section-centered pair)
    panel_serpine1_fibro()
    panel_cosmx()                          # also emits CosMx cell-type | {COL1A1,LYZ,KRT7,IL32}
    print("DONE.", flush=True)


if __name__ == "__main__":
    main()
