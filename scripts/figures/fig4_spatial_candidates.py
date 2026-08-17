#!/usr/bin/env python3
"""
Fig 4 spatial candidate gallery — render EVERY on-tissue / single-cell map option.

Exploration deliverable (NOT wired into the live figure). Produces a large set of
candidate spatial panels so the user can flip through and pick the strongest one(s)
for main Fig 4. Nothing in the live Fig 4 is modified; all output lands under
  figures/main/fig5_molecular_context/panels/spatial/_candidates/

Cohorts (user-selected this session): GSE192741 human Visium + Vu 2025 human Visium
+ Govaere CosMx single-cell. Visium maps are rendered BOTH over the registered H&E
(`*_he.pdf`) and as a clean white-background spot scatter (`*_clean.pdf`).

House rules honoured: PDF only (no PNG), one panel per PDF, all text black, short
gene-name titles OK / no subtitles/headers/in-plot sentences, control/healthy gray
#9E9E9E, gene names italic, no lollipop/3D. Descriptive captions go to
  spatial/_candidates/README.md  (not onto the plots).

Environment: spatial (scanpy). Run via run_fig4_spatial_candidates.sbatch on a
compute node — NOT the login node.
"""
import os, sys, json, traceback, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.patches import Patch
import scanpy as sc

sc.settings.verbosity = 0
# Editable vector text in Illustrator; never embed bitmap fonts.
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.5, "savefig.dpi": 300, "figure.dpi": 120,
})

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SPA  = os.path.join(BASE, "Analysis/Spatial/results")
SR   = os.path.join(SPA, "spaceranger")
GSM  = os.path.join(BASE, "Analysis/Spatial/data/gsmap_input")
OUT  = os.path.join(BASE, "figures/main/fig5_molecular_context/panels/spatial/_candidates")
DIR_GSE = os.path.join(OUT, "visium_gse192741")
DIR_VU  = os.path.join(OUT, "visium_vu")
DIR_COS = os.path.join(OUT, "cosmx")
DIR_QNT = os.path.join(OUT, "quantitative")
for d in (OUT, DIR_GSE, DIR_VU, DIR_COS, DIR_QNT):
    os.makedirs(d, exist_ok=True)

sys.path.insert(0, os.path.join(BASE, "Analysis/Spatial/scripts"))
try:
    from spatial_stats import ensure_lognorm
except Exception:
    ensure_lognorm = None

# ── Colors / semantics (mirror scripts/figures/publication_theme.R) ───────────
C_HEALTHY = "#9E9E9E"   # control gray (invariant)
C_DISEASE = "#C9265E"   # disease magenta
ZONE_COLORS = {"PP1": "#1565C0", "PP2": "#5E92C8", "Mid": "#9E9E9E",
               "PC2": "#E07B39", "PC1": "#E65100"}
ZONE_ORDER  = ["PP1", "PP2", "Mid", "PC2", "PC1"]
# Categorical cell-type palette (color-blind-safe, control-ish grays avoided for marks)
CAT_PAL = ["#C9265E", "#1565C0", "#2E7D32", "#E65100", "#6A1B9A", "#00838F",
           "#AD1457", "#558B2F", "#4527A0", "#EF6C00", "#00695C", "#9E9E9E",
           "#5D4037", "#283593", "#827717", "#BF360C"]
CMAP_EXPR = "magma"          # sequential expression
CMAP_ABUND = "rocket" if "rocket" in plt.colormaps() else "magma"
CMAP_ZON  = "RdBu_r"         # diverging zonation score (blue PP ↔ red PC)
CMAP_PT   = "viridis"        # pseudotime

# ── Hero / marker gene panels (presence-checked at render time) ───────────────
HERO_GENES = ["CYP3A4", "SERPINE1", "FADS2", "HKDC1", "RORA", "THRB"]
ZON_MARKERS = ["CYP2E1", "GLUL", "ASS1", "HAL", "SDS"]      # PC / PP anatomy
FIB_GENES  = ["COL1A1", "COL1A2", "ACTA2"]
VIS_GENES  = HERO_GENES + ZON_MARKERS + FIB_GENES

# c2l abundance columns of interest -> short label
C2L_PREFIX = "c2l_q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
C2L_TYPES = {  # full obs name suffix : short label
    "Hepatocytes": "Hepatocyte", "Fibroblasts": "Fibroblast",
    "Macrophages": "Macrophage", "Endothelial cells": "Endothelial",
    "Cholangiocytes": "Cholangiocyte", "T cells": "T cell", "B cells": "B cell",
}

GSE_COND = {"JBO018": "Healthy", "JBO022": "Healthy",
            "JBO014": "Steatotic", "JBO015": "Steatotic", "JBO019": "Steatotic"}

LOG = []          # caption / index lines for README
CONTACT = {}      # category -> list of (title, render_fn(ax)) for contact sheets


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def strip_bc(bcs):
    """Recover the base 10x Visium barcode (16nt + '-1') by stripping the
    anndata `obs_names_make_unique` suffix that the merged feature objects added
    (e.g. 'AAAC...-1-2' -> 'AAAC...-1'). Spaceranger barcodes pass through
    unchanged. Without this, per-sample joins onto spaceranger barcodes match
    almost nothing (the same Visium barcode recurs across samples)."""
    out = []
    for b in bcs:
        m = re.match(r"^([ACGTN]{16}-\d+)", str(b))
        out.append(m.group(1) if m else str(b))
    return out


def sample_feature(frame, sample_col, sample, barcodes):
    """Subset `frame` to one sample, normalize its barcodes, dedupe, and reindex
    onto `barcodes` (spaceranger order)."""
    sub = frame[frame[sample_col] == sample].copy()
    sub.index = strip_bc(sub.index.astype(str))
    sub = sub[~sub.index.duplicated()]
    return sub.reindex(pd.Index(barcodes).astype(str))


# ══════════════════════════════════════════════════════════════════════════════
#  Visium H&E loader (manual: robust across spaceranger versions)
# ══════════════════════════════════════════════════════════════════════════════
def load_visium_he(sample):
    """Return dict with counts AnnData (var=symbols, normalized log1p in .X),
    hires image, hires scalefactor, and a barcode->(x,y fullres px) frame.
    Coordinates are FULLRES pixels; multiply by scalef to overlay on hires image."""
    outs = os.path.join(SR, "GSE192741", sample, "outs")
    h5   = os.path.join(outs, "filtered_feature_bc_matrix.h5")
    spdir = os.path.join(outs, "spatial")
    ad = sc.read_10x_h5(h5)
    ad.var_names_make_unique()
    # positions
    pos_path = None
    for fn in ("tissue_positions_list.csv", "tissue_positions.csv"):
        p = os.path.join(spdir, fn)
        if os.path.exists(p):
            pos_path = p; break
    hdr = 0 if pos_path.endswith("tissue_positions.csv") else None
    pos = pd.read_csv(pos_path, header=hdr)
    if hdr is None:
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]
    pos = pos.set_index("barcode")
    # scalefactors + image
    with open(os.path.join(spdir, "scalefactors_json.json")) as fh:
        sf = json.load(fh)
    img = plt.imread(os.path.join(spdir, "tissue_hires_image.png"))
    scalef = sf["tissue_hires_scalef"]
    # align positions to matrix barcodes
    pos = pos.reindex(ad.obs_names)
    ad.obs["x_fullres"] = pos["pxl_col_in_fullres"].values  # col -> x
    ad.obs["y_fullres"] = pos["pxl_row_in_fullres"].values  # row -> y
    # normalize for expression maps (keep raw counts in a layer)
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return {"adata": ad, "img": img, "scalef": scalef}


# ══════════════════════════════════════════════════════════════════════════════
#  Core drawing primitive — draws ONE spot/cell map onto a provided Axes
# ══════════════════════════════════════════════════════════════════════════════
def draw_map(ax, x, y, values, *, categorical=False, cmap=CMAP_EXPR,
             cat_colors=None, img=None, scalef=None, title="", spot_scale=0.55,
             vmin=None, vmax=None, cbar=True, point_size=None, flip_y=True,
             cell_mode=False, na_color="#E8E8E8"):
    """If img+scalef given -> overlay on H&E (image coords). Else clean scatter.
    `values` may be numeric (cmap) or category labels (cat_colors dict)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    if img is not None and scalef is not None:
        ax.imshow(img, origin="upper")            # imshow fixes orientation
        px, py = x * scalef, y * scalef           # fullres -> hires px
    else:
        px, py = x, (-y if flip_y else y)         # clean scatter: flip y upright
    # auto spot size (cells smaller than Visium spots)
    if point_size is None:
        rng = (np.nanmax(px) - np.nanmin(px)) or 1.0
        n = max(len(px), 1)
        if img is not None:
            w_in = 2.4
        else:
            w_in = 2.4
        spacing = rng / np.sqrt(n)
        scale = (w_in * 72) / rng
        radius = spacing * scale * (0.30 if cell_mode else spot_scale)
        point_size = max(np.pi * radius ** 2, 0.2)

    if categorical:
        labs = pd.Series(values).astype(object)
        colmap = cat_colors or {}
        cols = labs.map(colmap).fillna(na_color).values
        ax.scatter(px, py, c=cols, s=point_size, linewidths=0, rasterized=True)
        handles = [Patch(facecolor=colmap[k], edgecolor="none", label=str(k))
                   for k in colmap if (labs == k).any()]
        if handles:
            ax.legend(handles=handles, fontsize=6, loc="center left",
                      bbox_to_anchor=(1.0, 0.5), frameon=False, handlelength=0.9,
                      handleheight=0.9, labelspacing=0.25, borderaxespad=0.1)
    else:
        v = np.asarray(values, float)
        finite = v[np.isfinite(v)]
        if vmin is None and finite.size:
            vmin = np.nanpercentile(finite, 1)
        if vmax is None and finite.size:
            vmax = np.nanpercentile(finite, 99)
        sca = ax.scatter(px, py, c=v, s=point_size, cmap=cmap, vmin=vmin,
                         vmax=vmax, linewidths=0, rasterized=True)
        if cbar:
            cb = ax.figure.colorbar(sca, ax=ax, fraction=0.040, pad=0.02)
            cb.ax.tick_params(labelsize=6, width=0.4, length=2)
            cb.outline.set_linewidth(0.3)
    ax.set_aspect("equal")
    ax.axis("off")
    if title:
        # gene symbols italic; descriptive words upright
        ax.set_title(title, fontsize=6, fontstyle=("italic" if title in VIS_GENES
                     or title.split()[0] in VIS_GENES else "normal"), pad=2,
                     color="black")


def emit(category, outdir, fname, render_fn, *, figsize=(2.8, 2.6)):
    """Render `render_fn(ax)` to an individual PDF AND queue it for the contact
    sheet. render_fn must accept a single Axes."""
    try:
        fig, ax = plt.subplots(figsize=figsize)
        render_fn(ax)
        fig.savefig(os.path.join(outdir, fname), bbox_inches="tight")
        plt.close(fig)
        CONTACT.setdefault(category, []).append((fname.replace(".pdf", ""), render_fn))
        log(f"  [ok] {category}/{fname}")
    except Exception as e:
        plt.close("all")
        log(f"  [SKIP] {category}/{fname}: {e}")


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION A+B — GSE192741 human Visium
# ══════════════════════════════════════════════════════════════════════════════
def section_gse192741():
    log("\n=== SECTION A/B: GSE192741 human Visium ===")
    # Feature object (zonation + deconv + dominant cell type), loaded once
    zad = sc.read_h5ad(os.path.join(SPA, "zonation/spatial_with_zonation.h5ad"))
    zobs = zad.obs.copy()        # barcode normalization handled in sample_feature()
    # optional domain + trajectory objects
    dom = traj = None
    try:
        d = sc.read_h5ad(os.path.join(SPA, "domains/spatial_with_domains.h5ad"))
        dom = d.obs[["sample_id"]].copy()
        dcol = next((c for c in d.obs.columns if "domain" in c.lower()), None)
        if dcol: dom["spatial_domain"] = d.obs[dcol].astype(str).values
        del d
    except Exception as e:
        log(f"  domains object unavailable: {e}")
    try:
        t = sc.read_h5ad(os.path.join(SPA, "trajectory/spatial_with_trajectory.h5ad"))
        ptcol = next((c for c in t.obs.columns if "pseudotime" in c.lower() or c == "dpt_pseudotime"), None)
        if ptcol:
            traj = t.obs[["sample_id"]].copy()
            traj["pt"] = t.obs[ptcol].values
        del t
    except Exception as e:
        log(f"  trajectory object unavailable: {e}")

    samples = [s for s in GSE_COND if os.path.isdir(os.path.join(SR, "GSE192741", s, "outs"))]
    montage_store = {g: {} for g in HERO_GENES}        # gene -> {sample: (vis dict)}
    for s in samples:
        cond = GSE_COND[s]
        log(f"\n-- {s} ({cond}) --")
        try:
            V = load_visium_he(s)
        except Exception as e:
            log(f"  [SKIP sample] {s}: {e}")
            continue
        ad = V["adata"]; img = V["img"]; scalef = V["scalef"]
        bc = ad.obs_names.astype(str)
        x, y = ad.obs["x_fullres"].values, ad.obs["y_fullres"].values
        ok = np.isfinite(x) & np.isfinite(y)
        # per-barcode feature frame from zonation object (barcode-normalized join)
        feat = sample_feature(zobs, "sample_id", s, bc)

        # --- 1. H&E context (no overlay) ---
        def _he(ax, img=img):
            ax.imshow(img, origin="upper"); ax.set_aspect("equal"); ax.axis("off")
            ax.set_title(f"{s} {cond} H&E", fontsize=6, color="black", pad=2)
        emit("gse", DIR_GSE, f"{s}_{cond}_HE.pdf", _he)

        # --- 2. zonation bin (categorical) + score (continuous) ---
        if "zonation_bin" in feat:
            zb = feat["zonation_bin"].astype(object).values
            for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                emit("gse", DIR_GSE, f"{s}_{cond}_zonationbin_{style}.pdf",
                     lambda ax, zb=zb, im=im, sca=sca: draw_map(
                         ax, x, y, zb, categorical=True, cat_colors=ZONE_COLORS,
                         img=im, scalef=sca, title="Zonation"))
        if "zonation_score" in feat:
            zs = pd.to_numeric(feat["zonation_score"], errors="coerce").values
            for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                emit("gse", DIR_GSE, f"{s}_{cond}_zonationscore_{style}.pdf",
                     lambda ax, zs=zs, im=im, sca=sca: draw_map(
                         ax, x, y, zs, cmap=CMAP_ZON, img=im, scalef=sca,
                         title="Zonation score"))

        # --- 3. dominant cell type ---
        if "cell_type_dominant" in feat:
            dctypes = feat["cell_type_dominant"].astype(object)
            uniq = [u for u in pd.unique(dctypes.dropna())]
            cc = {u: CAT_PAL[i % len(CAT_PAL)] for i, u in enumerate(sorted(map(str, uniq)))}
            for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                emit("gse", DIR_GSE, f"{s}_{cond}_dominantCT_{style}.pdf",
                     lambda ax, v=dctypes.astype(str).values, im=im, sca=sca, cc=cc: draw_map(
                         ax, x, y, v, categorical=True, cat_colors=cc, img=im,
                         scalef=sca, title="Dominant cell type"))

        # --- 4. cell2location abundance maps ---
        for suff, lab in C2L_TYPES.items():
            col = C2L_PREFIX + suff
            if col in feat:
                vals = pd.to_numeric(feat[col], errors="coerce").values
                for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                    emit("gse", DIR_GSE, f"{s}_{cond}_c2l_{lab.replace(' ','')}_{style}.pdf",
                         lambda ax, vals=vals, im=im, sca=sca, lab=lab: draw_map(
                             ax, x, y, vals, cmap=CMAP_ABUND, img=im, scalef=sca,
                             title=f"{lab} abundance"))

        # --- 5. spatial domain ---
        if dom is not None and "spatial_domain" in dom.columns:
            dsub = sample_feature(dom, "sample_id", s, bc)
            if True:
                dd = dsub["spatial_domain"].astype(object)
                if dd.notna().any():
                    uq = sorted(map(str, pd.unique(dd.dropna())))
                    cc = {u: CAT_PAL[i % len(CAT_PAL)] for i, u in enumerate(uq)}
                    for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                        emit("gse", DIR_GSE, f"{s}_{cond}_domain_{style}.pdf",
                             lambda ax, v=dd.astype(str).values, im=im, sca=sca, cc=cc: draw_map(
                                 ax, x, y, v, categorical=True, cat_colors=cc,
                                 img=im, scalef=sca, title="Spatial domain"))

        # --- 6. pseudotime (steatotic only) ---
        if traj is not None and cond == "Steatotic":
            tsub = sample_feature(traj, "sample_id", s, bc)
            if "pt" in tsub and tsub["pt"].notna().any():
                pv = pd.to_numeric(tsub["pt"], errors="coerce").values
                for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                    emit("gse", DIR_GSE, f"{s}_{cond}_pseudotime_{style}.pdf",
                         lambda ax, pv=pv, im=im, sca=sca: draw_map(
                             ax, x, y, pv, cmap=CMAP_PT, img=im, scalef=sca,
                             title="Disease pseudotime"))

        # --- 7. hero / marker / fibrosis gene expression ---
        for g in VIS_GENES:
            if g in ad.var_names:
                gv = np.asarray(ad[:, g].X.todense()).ravel() if hasattr(ad[:, g].X, "todense") \
                     else np.asarray(ad[:, g].X).ravel()
                for style, im, sca in (("he", img, scalef), ("clean", None, None)):
                    emit("gse", DIR_GSE, f"{s}_{cond}_gene_{g}_{style}.pdf",
                         lambda ax, gv=gv, im=im, sca=sca, g=g: draw_map(
                             ax, x, y, gv, cmap=CMAP_EXPR, img=im, scalef=sca, title=g))
                if g in HERO_GENES:
                    montage_store[g][s] = (x.copy(), y.copy(), gv.copy())

    # --- B. all-sample montage per hero gene (clean small-multiples) ---
    log("\n-- hero-gene montages --")
    for g, store in montage_store.items():
        if len(store) < 2:
            continue
        # build dedicated montage figure (small-multiples across samples)
        try:
            order = [s for s in GSE_COND if s in store]
            n = len(order)
            fig, axes = plt.subplots(1, n, figsize=(1.55 * n, 1.8))
            if n == 1: axes = [axes]
            vmax = np.nanpercentile(np.concatenate([store[s][2] for s in order]), 99)
            for ax, s in zip(axes, order):
                xx, yy, vv = store[s]
                ax.scatter(xx, -yy, c=vv, s=2.0, cmap=CMAP_EXPR, vmin=0, vmax=vmax,
                           linewidths=0, rasterized=True)
                ax.set_aspect("equal"); ax.axis("off")
                ax.set_title(f"{s}\n{GSE_COND[s]}", fontsize=6, color="black", pad=1)
            fig.suptitle(g, fontsize=6, fontstyle="italic", y=1.02)
            fig.savefig(os.path.join(DIR_GSE, f"MONTAGE_{g}_clean.pdf"), bbox_inches="tight")
            plt.close(fig)
            log(f"  [ok] gse/MONTAGE_{g}_clean.pdf")
        except Exception as e:
            plt.close("all"); log(f"  [SKIP] montage {g}: {e}")

    # --- E. quantitative companions from the zonation object ---
    try:
        quantitative_companions(zobs)
    except Exception as e:
        log(f"  [SKIP] quantitative companions: {e}")
    del zad


def quantitative_companions(zobs):
    log("\n-- quantitative companions --")
    # CYP3A4-style zonation gradient: mean abundance / score by zonation bin
    if "zonation_bin" in zobs and "condition" in zobs:
        df = zobs[zobs["zonation_bin"].notna()].copy()
        df["zonation_bin"] = pd.Categorical(df["zonation_bin"], ZONE_ORDER, ordered=True)
        # deconv composition stacked bar by condition
        c2l_cols = [C2L_PREFIX + s for s in C2L_TYPES if (C2L_PREFIX + s) in df]
        if c2l_cols:
            comp = df.groupby("condition")[c2l_cols].mean()
            comp = comp.div(comp.sum(axis=1), axis=0)
            comp.columns = [c.replace(C2L_PREFIX, "") for c in comp.columns]
            fig, ax = plt.subplots(figsize=(2.6, 2.4))
            bottom = np.zeros(len(comp))
            for i, ct in enumerate(comp.columns):
                ax.bar(comp.index, comp[ct].values, bottom=bottom,
                       color=CAT_PAL[i % len(CAT_PAL)], width=0.7, label=ct)
                bottom += comp[ct].values
            ax.set_ylabel("Mean cell-type fraction", fontsize=6)
            ax.tick_params(labelsize=6)
            ax.legend(fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5),
                      frameon=False, handlelength=0.9)
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            fig.savefig(os.path.join(DIR_QNT, "deconv_composition_by_condition.pdf"),
                        bbox_inches="tight"); plt.close(fig)
            log("  [ok] quantitative/deconv_composition_by_condition.pdf")


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION C — Vu 2025 human Visium (clean scatter only; pattern replication)
# ══════════════════════════════════════════════════════════════════════════════
def section_vu():
    log("\n=== SECTION C: Vu 2025 human Visium (clean scatter; pattern only) ===")
    files = sorted(f for f in os.listdir(GSM) if f.startswith("vu_") and f.endswith(".h5ad"))
    for fn in files:
        sid = fn.replace("vu_", "").replace(".h5ad", "")
        try:
            ad = sc.read_h5ad(os.path.join(GSM, fn))
            if "spatial" not in ad.obsm:
                log(f"  [SKIP] {sid}: no spatial coords"); continue
            xy = ad.obsm["spatial"]
            x, y = xy[:, 0], xy[:, 1]
            # normalize for expression
            if ad.X.max() > 50:  # looks like counts
                sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
            for g in HERO_GENES + ["COL1A1"]:
                if g in ad.var_names:
                    gv = np.asarray(ad[:, g].X.todense()).ravel() if hasattr(ad[:, g].X, "todense") \
                         else np.asarray(ad[:, g].X).ravel()
                    emit("vu", DIR_VU, f"{sid}_gene_{g}_clean.pdf",
                         lambda ax, gv=gv, x=x, y=y, g=g: draw_map(
                             ax, x, y, gv, cmap=CMAP_EXPR, title=g))
            del ad
        except Exception as e:
            log(f"  [SKIP] {sid}: {e}")


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION D — Govaere CosMx single-cell
# ══════════════════════════════════════════════════════════════════════════════
def section_cosmx():
    log("\n=== SECTION D: Govaere CosMx single-cell ===")
    fn = os.path.join(BASE, "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad")
    ad = sc.read_h5ad(fn)
    # gene symbols live in var['gene_symbol']
    if "gene_symbol" in ad.var:
        ad.var["_orig"] = ad.var_names
        ad.var_names = ad.var["gene_symbol"].astype(str).values
        ad.var_names_make_unique()
    # use lognorm layer for expression
    if "lognorm" in ad.layers:
        ad.X = ad.layers["lognorm"]
    obs = ad.obs
    xcol = "centerX_global_px" if "centerX_global_px" in obs else None
    ycol = "centerY_global_px" if "centerY_global_px" in obs else None
    if xcol is None and "spatial" in ad.obsm:
        coords_all = ad.obsm["spatial"]
    cell_col = "cell_type" if "cell_type" in obs else None
    slide_col = "sample_id" if "sample_id" in obs else None
    dis_col = "sample_disease" if "sample_disease" in obs else None
    slides = list(pd.unique(obs[slide_col])) if slide_col else [None]
    ctypes = sorted(map(str, pd.unique(obs[cell_col].dropna()))) if cell_col else []
    ct_colors = {c: CAT_PAL[i % len(CAT_PAL)] for i, c in enumerate(ctypes)}
    cos_genes = [g for g in ["IL32", "CD74", "SPP1", "TREM2", "COL1A1", "CD68"] if g in ad.var_names]
    prot_cols = [c for c in ["CD68_intensity", "PanCK_intensity", "CD45_intensity"] if c in obs]

    for sl in slides:
        m = (obs[slide_col] == sl).values if slide_col else np.ones(ad.n_obs, bool)
        sub = ad[m]
        so = sub.obs
        if xcol:
            x, y = so[xcol].values.astype(float), so[ycol].values.astype(float)
        else:
            x, y = sub.obsm["spatial"][:, 0], sub.obsm["spatial"][:, 1]
        dis = str(so[dis_col].iloc[0]) if dis_col and len(so) else ""
        tag = f"{sl}_{dis}".replace(" ", "").replace("/", "-")
        log(f"\n-- CosMx slide {sl} ({dis}) n={m.sum():,} --")

        # cell-type map
        if cell_col:
            cv = so[cell_col].astype(str).values
            emit("cosmx", DIR_COS, f"{tag}_celltype.pdf",
                 lambda ax, cv=cv, x=x, y=y: draw_map(
                     ax, x, y, cv, categorical=True, cat_colors=ct_colors,
                     cell_mode=True, title="Cell type"), figsize=(3.4, 2.8))
        # macrophage density (KC + Mac subsets)
        if cell_col:
            mac = np.isin(so[cell_col].astype(str).values,
                          ["KC", "MetMac", "TransMac", "Macrophage", "preMac", "Monocyte"])
            emit("cosmx", DIR_COS, f"{tag}_macrophage_mask.pdf",
                 lambda ax, mac=mac, x=x, y=y: draw_map(
                     ax, x, y, np.where(mac, "Macrophage", "other"),
                     categorical=True,
                     cat_colors={"Macrophage": C_DISEASE, "other": "#DDDDDD"},
                     cell_mode=True, title="Macrophages"))
        # gene expression per-cell
        for g in cos_genes:
            gv = np.asarray(sub[:, g].X.todense()).ravel() if hasattr(sub[:, g].X, "todense") \
                 else np.asarray(sub[:, g].X).ravel()
            emit("cosmx", DIR_COS, f"{tag}_gene_{g}.pdf",
                 lambda ax, gv=gv, x=x, y=y, g=g: draw_map(
                     ax, x, y, gv, cmap=CMAP_EXPR, cell_mode=True, title=g))
        # protein intensity
        for pc in prot_cols:
            pv = pd.to_numeric(so[pc], errors="coerce").values
            emit("cosmx", DIR_COS, f"{tag}_prot_{pc.replace('_intensity','')}.pdf",
                 lambda ax, pv=pv, x=x, y=y, pc=pc: draw_map(
                     ax, x, y, pv, cmap="cividis", cell_mode=True,
                     title=pc.replace("_intensity", " protein")))
        del sub
    del ad


# ══════════════════════════════════════════════════════════════════════════════
#  Contact sheets + README index
# ══════════════════════════════════════════════════════════════════════════════
def build_contact_sheets():
    log("\n=== contact sheets ===")
    for cat, items in CONTACT.items():
        if not items:
            continue
        path = os.path.join(OUT, f"INDEX_{cat}_contact_sheet.pdf")
        per_page = 12; ncol = 4
        with PdfPages(path) as pdf:
            for start in range(0, len(items), per_page):
                chunk = items[start:start + per_page]
                nrow = int(np.ceil(len(chunk) / ncol))
                fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.7 * nrow))
                axes = np.array(axes).reshape(-1)
                for ax in axes: ax.axis("off")
                for ax, (title, fn) in zip(axes, chunk):
                    try:
                        fn(ax)
                    except Exception as e:
                        ax.text(0.5, 0.5, f"[err]\n{title}", fontsize=6, ha="center")
                    ax.set_title(title, fontsize=6, color="black")
                log(f"  [caption] {cat} contact sheet — candidates {start+1}-{start+len(chunk)}")
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
        log(f"  [ok] INDEX_{cat}_contact_sheet.pdf ({len(items)} panels)")


def write_readme():
    pdfs = []
    for root, _, files in os.walk(OUT):
        for f in sorted(files):
            if f.endswith(".pdf"):
                pdfs.append(os.path.relpath(os.path.join(root, f), OUT))
    lines = [
        "# Fig 4 spatial candidate gallery", "",
        "Exhaustive sweep of on-tissue / single-cell spatial maps for main Fig 4.",
        "Exploration only — nothing here is wired into the live figure. Pick the",
        "strongest panel(s) and promote them out of `spatial/_candidates/`.", "",
        "## Cohorts",
        "- **visium_gse192741/** — 5 human Visium (JBO018/JBO022 Healthy; "
        "JBO014/JBO015/JBO019 Steatotic). Maps rendered over registered H&E "
        "(`*_he.pdf`) and as clean spot scatter (`*_clean.pdf`).",
        "- **visium_vu/** — Vu 2025 human Visium, clean scatter only. "
        "PATTERN REPLICATION ONLY — different cohort/batch, never a disease-direction claim.",
        "- **cosmx/** — Govaere CosMx single-cell (Leuven_1-4). Per-cell maps: cell type, "
        "IL32/CD74/SPP1/TREM2/COL1A1 expression, CD68/PanCK/CD45 protein.",
        "- **quantitative/** — companion summaries (deconv composition).", "",
        f"## Panels rendered: {len(pdfs)}", "",
        "Browse the `INDEX_*_contact_sheet.pdf` files first; each thumbnail maps to a",
        "same-named individual PDF in the category subdir.", "",
        "## File list", "",
    ]
    lines += [f"- `{p}`" for p in pdfs]
    with open(os.path.join(OUT, "README.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.join(OUT, "_render_log.txt"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    log(f"\nREADME + log written. Total PDFs: {len(pdfs)}")


def main():
    for fn in (section_gse192741, section_vu, section_cosmx):
        try:
            fn()
        except Exception as e:
            log(f"\n!!! {fn.__name__} failed: {e}\n{traceback.format_exc()}")
    build_contact_sheets()
    write_readme()
    log("\nDONE.")


if __name__ == "__main__":
    main()
