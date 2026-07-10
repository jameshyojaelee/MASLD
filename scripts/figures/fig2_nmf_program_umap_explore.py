#!/usr/bin/env python3
"""Explore bulk NMF k=6 + cNMF k=16 programs on the all-cell scRNA UMAP.

Outputs (figures/misc/nmf_umap_exploration/):
  - bulk_nmf_k6_umap_facets.pdf      (2x3 small-multiples, viridis intensity)
  - bulk_nmf_k6_argmax_umap.pdf      (single panel, categorical argmax)
  - cnmf_k16_umap_facets.pdf         (4x4 small-multiples, viridis intensity)
  - cnmf_k16_argmax_umap.pdf         (single panel, categorical argmax)
  - subsampled_cells.csv.gz          (subsample IDs + per-cell scores)

Inputs:
  - atlas_cnmf_global.h5ad         (895K cells, raw counts, no UMAP)
  - scored_atlas.h5ad              (1.23M cells, UMAP, no X)  -- aligned by barcode
  - global.usages.k_16.dt_0_03.consensus.txt
  - nmf_k6_top50_genes.tsv         (from export_bulk_nmf_top_genes.R)
"""

from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import colors as mcolors
from matplotlib.colors import ListedColormap

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))

ATLAS_X = BASE / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_global.h5ad"
ATLAS_UMAP = BASE / "Analysis/SingleCell/results_gpu_v2/disease_signatures/scored_atlas.h5ad"
CNMF_USAGES = BASE / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs/global/global.usages.k_16.dt_0_03.consensus.txt"
BULK_TOPGENES = BASE / "figures/misc/nmf_umap_exploration/nmf_k6_top50_genes.tsv"
BULK_ATLAS = BASE / "RNA-seq/results/subtypes/nmf_ksweep_program_atlas.csv"

OUTDIR = BASE / "figures/misc/nmf_umap_exploration"
OUTDIR.mkdir(parents=True, exist_ok=True)

SUBSAMPLE_N = 200_000
RNG_SEED = 42

# Labels + colors mirror fig2d (scripts/figures/fig2_progression_sex.R lines 101-109).
BULK_NMF_LABELS = {
    "P1": "P1 Pro-inflammatory",
    "P2": "P2 Innate-immune",
    "P3": "P3 Parenchymal",
    "P4": "P4 lncRNA",
    "P5": "P5 Hepatic-metabolic",
    "P6": "P6 Stellate-myofibroblast",
}
BULK_NMF_COLORS = {
    "P1": "#C2185B", "P2": "#7B1FA2", "P3": "#1565C0",
    "P4": "#00897B", "P5": "#558B2F", "P6": "#E65100",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def msg(*a):
    print(*a, flush=True)


def load_umap(path: Path) -> pd.DataFrame:
    """Load UMAP coords + cell_type from scored_atlas.h5ad keyed by barcode."""
    msg(f"Loading UMAP: {path}")
    with h5py.File(path, "r") as f:
        idx = f["obs/_index"][:]
        idx = np.array([b.decode() if isinstance(b, bytes) else b for b in idx])
        umap = f["obsm/X_umap"][:]
        # cell_type may be categorical
        ct_obj = f["obs/cell_type"]
        if "codes" in ct_obj:
            cats = ct_obj["categories"][:]
            cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
            ct = cats[ct_obj["codes"][:]]
        else:
            ct_raw = ct_obj[:]
            ct = np.array([c.decode() if isinstance(c, bytes) else c for c in ct_raw])
    df = pd.DataFrame(
        {"umap_1": umap[:, 0].astype(np.float32),
         "umap_2": umap[:, 1].astype(np.float32),
         "cell_type": ct},
        index=pd.Index(idx, name="barcode"))
    msg(f"  UMAP rows: {len(df):,}")
    return df


def load_cnmf_usages(path: Path) -> pd.DataFrame:
    """cNMF cell x program usage matrix; row-normalize so each cell sums to 1."""
    msg(f"Loading cNMF usages: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    # programs are columns 1..16
    df.columns = [f"cnmf_K{c}" for c in df.columns]
    msg(f"  usages: {df.shape[0]:,} cells x {df.shape[1]} programs")
    rs = df.sum(axis=1).replace(0, np.nan)
    df = df.div(rs, axis=0).fillna(0).astype(np.float32)
    return df


def score_bulk_programs(adata: ad.AnnData, topgenes: pd.DataFrame) -> pd.DataFrame:
    """Run sc.tl.score_genes for each P1..P6 using top-50 marker genes.

    Atlas var_names use gene symbols (with some ENSG fallbacks). Bulk gene IDs
    are ENSEMBL; topgenes table includes a `symbol` column we resolve to.
    """
    var_names = set(adata.var_names)
    msg(f"  atlas vars: {len(var_names):,}")

    scored_cols = []
    for prog in sorted(topgenes["program"].unique()):
        symbols = topgenes.loc[topgenes["program"] == prog, "symbol"].tolist()
        ensembls = topgenes.loc[topgenes["program"] == prog, "gene_id"].tolist()
        # Prefer symbols, fall back to ENSG
        gene_list = [s for s in symbols if s in var_names]
        if len(gene_list) < 10:
            gene_list += [e for e in ensembls if e in var_names and e not in gene_list]
        gene_list = list(dict.fromkeys(gene_list))  # dedup, preserve order
        msg(f"  {prog}: {len(gene_list)}/{len(symbols)} genes resolved in atlas")

        if len(gene_list) < 5:
            msg(f"  WARN: {prog} has only {len(gene_list)} resolvable genes, skipping")
            continue

        col = f"bulk_{prog}"
        sc.tl.score_genes(
            adata, gene_list=gene_list, score_name=col,
            random_state=RNG_SEED, n_bins=25, ctrl_size=50,
            use_raw=False)
        scored_cols.append(col)

    out = adata.obs[scored_cols].copy()
    return out


def stratified_subsample(df: pd.DataFrame, group_col: str, n: int, rng) -> pd.DataFrame:
    """Stratified subsample by group_col preserving relative proportions."""
    if len(df) <= n:
        return df.copy()
    # proportional allocation
    sizes = df.groupby(group_col, observed=True).size()
    alloc = (sizes / sizes.sum() * n).round().astype(int)
    # Make sure no group gets 0 if it has cells
    alloc = alloc.clip(lower=1)
    parts = []
    for grp, k in alloc.items():
        sub = df[df[group_col] == grp]
        k_use = min(k, len(sub))
        parts.append(sub.sample(n=k_use, random_state=rng))
    out = pd.concat(parts, axis=0).sample(frac=1.0, random_state=rng)
    msg(f"  subsampled {len(out):,} / {len(df):,} cells (target {n:,})")
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _scatter_panel(ax, df, color_vals, *, vmin=None, vmax=None, cmap="viridis",
                   title="", point_size=0.3, rasterized=True):
    sc_handle = ax.scatter(df["umap_1"], df["umap_2"], c=color_vals,
                           s=point_size, vmin=vmin, vmax=vmax, cmap=cmap,
                           linewidths=0, rasterized=rasterized)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.4)
    ax.set_title(title, fontsize=6, pad=2)
    return sc_handle


def render_facets_continuous(df, score_cols, labels, ncols, outpath, *,
                             panel_w=2.0, panel_h=2.0, title=None):
    """Small-multiples UMAP, one panel per program, viridis intensity."""
    nrows = int(np.ceil(len(score_cols) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(panel_w*ncols, panel_h*nrows),
                             squeeze=False)
    for k, col in enumerate(score_cols):
        r, c = divmod(k, ncols)
        ax = axes[r][c]
        vals = df[col].to_numpy()
        # Per-panel min-max so each program is its own gradient
        vmin = float(np.percentile(vals, 1))
        vmax = float(np.percentile(vals, 99))
        if vmax <= vmin:
            vmax = vmin + 1e-6
        # Plot low scores first so high scores are on top
        order = np.argsort(vals)
        sub = df.iloc[order]
        sub_vals = vals[order]
        sc_h = _scatter_panel(ax, sub, sub_vals, vmin=vmin, vmax=vmax,
                              title=labels.get(col, col))
        cb = fig.colorbar(sc_h, ax=ax, fraction=0.04, pad=0.02, shrink=0.8)
        cb.ax.tick_params(labelsize=6, length=2)
        cb.outline.set_linewidth(0.3)
    # blank any leftover axes
    for k in range(len(score_cols), nrows*ncols):
        r, c = divmod(k, ncols)
        axes[r][c].axis("off")
    if title:
        msg(f"[caption] {title}")
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    msg(f"  wrote {outpath}")


def render_argmax(df, score_cols, labels, palette, outpath, *, title=None):
    """Single UMAP, cells colored by their argmax program (categorical)."""
    arr = df[score_cols].to_numpy()
    am = np.argmax(arr, axis=1)
    am_labels = np.array([score_cols[i] for i in am])
    # Render
    fig, ax = plt.subplots(figsize=(6, 5))
    # plot per-program so the legend is built from real handles
    rng = np.random.default_rng(RNG_SEED)
    perm = rng.permutation(len(df))
    df_perm = df.iloc[perm]
    am_perm = am_labels[perm]
    for col in score_cols:
        m = am_perm == col
        if not m.any():
            continue
        ax.scatter(df_perm.loc[m, "umap_1"], df_perm.loc[m, "umap_2"],
                   s=0.3, c=palette[col], label=labels.get(col, col),
                   linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.4)
    leg = ax.legend(markerscale=10, fontsize=6, frameon=False,
                    loc="center left", bbox_to_anchor=(1.02, 0.5))
    if title:
        msg(f"[caption] {title}")
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    msg(f"  wrote {outpath}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def render_only(sub_csv: Path):
    """Re-render PDFs from a previously-saved subsampled_cells.csv.gz.

    Skips atlas loading + scoring (saves ~10 min) when only label / color
    tweaks are needed.
    """
    msg(f"==> Render-only: loading {sub_csv}")
    sub = pd.read_csv(sub_csv)
    if "barcode" in sub.columns:
        sub = sub.set_index("barcode")
    bulk_score_cols = [c for c in sub.columns if c.startswith("bulk_P")]
    cnmf_cols = [c for c in sub.columns if c.startswith("cnmf_K")]
    msg(f"  rows: {len(sub):,}  bulk cols: {bulk_score_cols}  cnmf cols: {len(cnmf_cols)}")

    bulk_label_for_col = {f"bulk_{p}": BULK_NMF_LABELS[p]
                          for p in sorted(BULK_NMF_LABELS)}
    render_facets_continuous(
        sub, bulk_score_cols, bulk_label_for_col, ncols=3,
        outpath=OUTDIR / "bulk_nmf_k6_umap_facets.pdf",
        panel_w=2.4, panel_h=2.2,
        title="Bulk NMF k=6 program scores (sc.tl.score_genes top-50 markers)")

    bulk_palette = {f"bulk_{p}": BULK_NMF_COLORS[p]
                    for p in sorted(BULK_NMF_COLORS)}
    render_argmax(
        sub, bulk_score_cols, bulk_label_for_col, bulk_palette,
        outpath=OUTDIR / "bulk_nmf_k6_argmax_umap.pdf",
        title="Bulk NMF k=6 dominant program per cell")

    cnmf_label_for_col = {c: c.replace("cnmf_", "") for c in cnmf_cols}
    render_facets_continuous(
        sub, cnmf_cols, cnmf_label_for_col, ncols=4,
        outpath=OUTDIR / "cnmf_k16_umap_facets.pdf",
        panel_w=1.85, panel_h=1.9,
        title="cNMF k=16 program usages")

    cnmf_palette_cmap = plt.get_cmap("husl") if "husl" in plt.colormaps() \
                        else plt.get_cmap("tab20")
    cnmf_palette = {c: cnmf_palette_cmap(i / max(len(cnmf_cols)-1, 1))
                    for i, c in enumerate(cnmf_cols)}
    render_argmax(
        sub, cnmf_cols, cnmf_label_for_col, cnmf_palette,
        outpath=OUTDIR / "cnmf_k16_argmax_umap.pdf",
        title="cNMF k=16 dominant program per cell")


def main():
    rng = np.random.default_rng(RNG_SEED)

    # 1. Load top-50 bulk genes + canonical labels
    msg("==> Loading bulk top-50 marker genes")
    topgenes = pd.read_csv(BULK_TOPGENES, sep="\t")
    msg(f"  {len(topgenes)} rows; programs: {sorted(topgenes['program'].unique())}")
    bulk_label_map = dict(BULK_NMF_LABELS)  # P1..P6 -> "P1 Pro-inflammatory" etc.
    msg(f"  bulk labels: {bulk_label_map}")

    # 2. Load atlas with X (raw counts)
    msg("==> Loading atlas with X")
    adata = sc.read_h5ad(ATLAS_X)
    msg(f"  atlas: {adata.shape}")
    msg(f"  X type: {type(adata.X).__name__}, dtype: {adata.X.dtype}")
    msg(f"  X max: {adata.X[:200,:].max():.1f} (raw counts expected)")

    # Normalize + log1p
    msg("==> Normalizing and log-transforming X")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    msg(f"  X max after lognorm: {adata.X[:200,:].max():.2f}")

    # Drop the counts layer (already in X) to free memory
    if "counts" in adata.layers:
        del adata.layers["counts"]

    # 3. Score bulk programs on cells
    msg("==> Scoring bulk programs")
    bulk_scores = score_bulk_programs(adata, topgenes)
    bulk_score_cols = list(bulk_scores.columns)
    msg(f"  bulk score cols: {bulk_score_cols}")

    # 4. Pull metadata, drop X to free RAM
    obs = adata.obs[["sample", "dataset", "cell_type"]].copy()
    obs = obs.join(bulk_scores)
    obs.index = obs.index.rename("barcode")
    del adata
    import gc; gc.collect()

    # 5. Load UMAP and align
    msg("==> Loading UMAP and aligning")
    umap_df = load_umap(ATLAS_UMAP)
    common = obs.index.intersection(umap_df.index)
    msg(f"  intersection (atlas X & UMAP): {len(common):,} cells")
    obs = obs.loc[common].join(umap_df.loc[common, ["umap_1", "umap_2"]])
    msg(f"  obs after UMAP join: {obs.shape}")

    # 6. Load cNMF usages and align
    msg("==> Loading cNMF usages and aligning")
    cnmf = load_cnmf_usages(CNMF_USAGES)
    cnmf_cols = list(cnmf.columns)
    common2 = obs.index.intersection(cnmf.index)
    msg(f"  intersection (above & cNMF): {len(common2):,} cells")
    obs = obs.loc[common2].join(cnmf.loc[common2, cnmf_cols])

    # 7. Subsample stratified by cell_type
    msg("==> Stratified subsampling")
    sub = stratified_subsample(obs, "cell_type", SUBSAMPLE_N, rng)

    # 8. Sanity-check log: argmax-program x cell-type contingency
    msg("==> Sanity-check argmax x cell_type contingencies")
    bulk_argmax = sub[bulk_score_cols].idxmax(axis=1)
    print("BULK argmax x cell_type:")
    print(pd.crosstab(bulk_argmax, sub["cell_type"]).to_string())
    cnmf_argmax = sub[cnmf_cols].idxmax(axis=1)
    print("\ncNMF argmax x cell_type:")
    print(pd.crosstab(cnmf_argmax, sub["cell_type"]).to_string())

    # 9. Save subsample reproducibility CSV
    out_csv = OUTDIR / "subsampled_cells.csv.gz"
    msg(f"==> Writing {out_csv}")
    sub.reset_index().to_csv(out_csv, index=False, compression="gzip")

    # 10. Render PDFs
    msg("==> Rendering PDFs")
    # Bulk facets (2x3) -- use score cols ordered by program number
    bulk_label_for_col = {f"bulk_{p}": bulk_label_map.get(p, p)
                          for p in sorted(topgenes["program"].unique())}
    render_facets_continuous(
        sub, bulk_score_cols, bulk_label_for_col, ncols=3,
        outpath=OUTDIR / "bulk_nmf_k6_umap_facets.pdf",
        panel_w=2.4, panel_h=2.2,
        title="Bulk NMF k=6 program scores (sc.tl.score_genes top-50 markers)")

    # Bulk argmax -- use fig2d palette
    bulk_palette = {f"bulk_{p}": BULK_NMF_COLORS[p]
                    for p in sorted(BULK_NMF_COLORS)}
    render_argmax(
        sub, bulk_score_cols, bulk_label_for_col, bulk_palette,
        outpath=OUTDIR / "bulk_nmf_k6_argmax_umap.pdf",
        title="Bulk NMF k=6 dominant program per cell")

    # cNMF facets (4x4)
    cnmf_label_for_col = {c: c.replace("cnmf_", "") for c in cnmf_cols}
    render_facets_continuous(
        sub, cnmf_cols, cnmf_label_for_col, ncols=4,
        outpath=OUTDIR / "cnmf_k16_umap_facets.pdf",
        panel_w=1.85, panel_h=1.9,
        title="cNMF k=16 program usages")

    # cNMF argmax (16-color husl palette)
    cnmf_palette_cmap = plt.get_cmap("husl") if "husl" in plt.colormaps() \
                        else plt.get_cmap("tab20")
    cnmf_palette = {c: cnmf_palette_cmap(i / max(len(cnmf_cols)-1, 1))
                    for i, c in enumerate(cnmf_cols)}
    render_argmax(
        sub, cnmf_cols, cnmf_label_for_col, cnmf_palette,
        outpath=OUTDIR / "cnmf_k16_argmax_umap.pdf",
        title="cNMF k=16 dominant program per cell")

    msg("==> Done")


if __name__ == "__main__":
    sub_csv = OUTDIR / "subsampled_cells.csv.gz"
    if "--render-only" in sys.argv and sub_csv.exists():
        render_only(sub_csv)
    else:
        main()
