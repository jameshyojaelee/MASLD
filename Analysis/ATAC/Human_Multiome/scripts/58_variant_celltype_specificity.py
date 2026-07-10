#!/usr/bin/env python3
"""58_variant_celltype_specificity.py — lineage-SPECIFICITY of the fig4d in-peak
GWAS credible-set variants.

MOTIVATION (2026-07-09)
-----------------------
Fig 4d counts a prioritized gene in a lineage's bar if its credible-set variant
falls inside a MACS3 peak CALLED in that lineage. That is confounded two ways:
  (1) broadly-open promoter-proximal elements are accessible in many lineages
      (a gene lands in many bars regardless of lineage-specific action), and
  (2) per-cell-type peak-set size tracks cell abundance (Hepatocytes have the
      most cells -> most peaks called -> most overlaps), so bar height is partly
      a coverage/power artifact, not biology.

This script re-scores each in-peak variant by DIFFERENTIAL/UNIQUE accessibility:
per lineage, the fraction of that lineage's cells with an insertion at the
variant's locus, normalized by the lineage's own genome-wide background
fraction (removes the depth/abundance offset). A variant is "lineage-specific"
if one lineage's normalized accessibility exceeds the others by >= FOLD. This
converts "where is this locus open at all" into "where is it PREFERENTIALLY
open", so shared promoters (e.g. GCKR) drop out and only lineage-enriched
regulatory elements remain.

This is CELL-TYPE marker accessibility (one-vs-rest, per-cell) — well-powered
(thousands of cells/type) and disease-agnostic; it does NOT reintroduce the
donor-level / disease-vs-control pseudoreplication issues.

Substrate: the label-transferred h5ad 500bp tile matrix (var_names = tiles), the
same matrix 07c uses. Each variant's hg38 position maps to its containing tile.

Env: snapatac2. Run via run_variant_specificity.sbatch (not the login node).
Output: results/gwas_atac/variant_celltype_specificity.csv
"""

import argparse, importlib.util, logging, os, sys, time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S", handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
_spec = importlib.util.spec_from_file_location(
    "da07b", os.path.join(SCRIPT_DIR, "07b_corrected_da_hepatocytes.py"))
da07b = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(da07b)
TILE_SIZE = da07b.TILE_SIZE  # 500

PROJECT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

# obs cell_type -> fig4d lineage (matches fig4d GRP groupings)
CT2LIN = {
    "Hepatocyte": "Hepatocytes",
    "Stellate_Cell": "Fibroblasts",
    "Macrophage": "Macrophages", "Kupffer_Cell": "Macrophages",
    "Cholangiocyte": "Cholangiocytes",
    "Endothelial": "Endothelial", "LSEC": "Endothelial",
    "Plasma_Cell": "Lymphoid", "NK_T_Cell": "Lymphoid", "B_Cell": "Lymphoid",
}
LINEAGES = ["Hepatocytes", "Fibroblasts", "Macrophages", "Cholangiocytes",
            "Endothelial", "Lymphoid"]


def pos_to_tile(chrom, pos, chrom_to_range):
    """hg38 (chrom,pos) -> tile column index (tiles are contiguous 500bp from 0)."""
    if chrom not in chrom_to_range:
        return None
    c0, c1 = chrom_to_range[chrom]
    t = c0 + pos // TILE_SIZE
    return int(t) if c0 <= t <= c1 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=os.path.join(
        PROJECT, "Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_label_transferred.h5ad"))
    ap.add_argument("--variants", default=os.path.join(
        PROJECT, "GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv"))
    ap.add_argument("--out", default=os.path.join(
        PROJECT, "GWAS/finemapping/results/gwas_atac/variant_celltype_specificity.csv"))
    ap.add_argument("--n-bg-tiles", type=int, default=20000,
                    help="random tiles for the per-lineage background fraction")
    args = ap.parse_args()
    t0 = time.time()

    # unique in-peak variants (hg38) + linked gene
    v = pd.read_csv(args.variants)
    v = v.dropna(subset=["chr_hg38", "pos_hg38"])
    v = v.drop_duplicates("variant_id").copy()
    v["pos_hg38"] = v["pos_hg38"].astype(int)
    log.info("in-peak variants: %d", len(v))

    log.info("loading h5ad (backed)...")
    adata = ad.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs
    ct = obs["cell_type"].astype(str).values
    lin = np.array([CT2LIN.get(x, "other") for x in ct])
    tile_names = list(adata.var_names)
    chrom_to_range = da07b.parse_chrom_sizes_from_tiles(tile_names)
    n_tiles = len(tile_names)
    log.info("cells=%d tiles=%d; per-lineage cells: %s",
             adata.n_obs, n_tiles,
             {L: int((lin == L).sum()) for L in LINEAGES})

    # map variants -> tile columns
    v["tile"] = [pos_to_tile(c, p, chrom_to_range)
                 for c, p in zip(v["chr_hg38"], v["pos_hg38"])]
    v = v[v["tile"].notna()].copy(); v["tile"] = v["tile"].astype(int)
    tiles = sorted(v["tile"].unique())
    log.info("variants mapped to %d unique tiles", len(tiles))

    # background tiles (evenly spaced) + variant tiles -> single column set read once
    bg_tiles = np.unique(np.linspace(0, n_tiles - 1, args.n_bg_tiles).astype(int))
    needed = np.array(sorted(set(tiles) | set(bg_tiles.tolist())))
    col_pos = {int(t): i for i, t in enumerate(needed)}   # tile -> index within `needed`
    log.info("single chunked pass over %d cells, keeping %d columns (%d variant + %d bg)...",
             adata.n_obs, len(needed), len(tiles), len(bg_tiles))

    # per-lineage binarized column-sums over `needed`, accumulated across cell chunks
    lin_codes = {L: i for i, L in enumerate(LINEAGES)}
    sums = np.zeros((len(LINEAGES), len(needed)), dtype=np.float64)
    counts = np.zeros(len(LINEAGES), dtype=np.int64)
    for L in LINEAGES:
        counts[lin_codes[L]] = int((lin == L).sum())
    chunk = 5000
    for s in range(0, adata.n_obs, chunk):
        e = min(s + chunk, adata.n_obs)
        Xc = adata.X[s:e, :]
        if not sparse.issparse(Xc):
            Xc = sparse.csr_matrix(Xc)
        Xb = (Xc[:, needed] > 0).astype(np.float32)   # chunk_cells x |needed|
        lc = lin[s:e]
        for L in LINEAGES:
            m = lc == L
            if m.any():
                sums[lin_codes[L]] += np.asarray(Xb[m, :].sum(axis=0)).ravel()
        if (s // chunk) % 4 == 0:
            log.info("  rows %d/%d", e, adata.n_obs)
    adata.file.close()

    frac = sums / np.clip(counts[:, None], 1, None)    # (n_lineage x |needed|) fraction accessible
    bg = np.clip(frac.mean(axis=1), 1e-6, None)        # (n_lineage,) mean over all needed cols
    log.info("  background frac/lineage: %s",
             {L: round(float(b), 4) for L, b in zip(LINEAGES, bg)})

    # accessibility at each variant tile: (n_variant_tile x n_lineage)
    fa = np.array([[frac[lin_codes[L], col_pos[t]] for L in LINEAGES] for t in tiles])
    tile2row = {t: i for i, t in enumerate(tiles)}

    # normalized accessibility = lineage frac / lineage background (fold-over-own-depth)
    rows = []
    for _, r in v.iterrows():
        raw = fa[tile2row[r["tile"]]]                 # (n_lineage,)
        norm = raw / bg
        order = np.argsort(norm)[::-1]
        top, second = order[0], order[1]
        max_lin = LINEAGES[top]
        fold = norm[top] / max(norm[second], 1e-9)    # top vs 2nd lineage
        fold_rest = norm[top] / max(np.mean(np.delete(norm, top)), 1e-9)
        d = dict(variant_id=r["variant_id"], gene=r.get("linked_gene", r.get("nearest_gene")),
                 max_pip=r.get("max_pip", np.nan),
                 specific_lineage=max_lin, fold_vs_2nd=round(float(fold), 2),
                 fold_vs_rest=round(float(fold_rest), 2))
        for L, rw, nm in zip(LINEAGES, raw, norm):
            d[f"frac_{L}"] = round(float(rw), 4)
            d[f"norm_{L}"] = round(float(nm), 3)
        rows.append(d)
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)
    log.info("wrote %s (%d variants)", args.out, len(out))

    # quick summary at a couple of fold thresholds
    for fold in (1.5, 2.0, 3.0):
        spec = out[out["fold_vs_2nd"] >= fold]
        by = spec["specific_lineage"].value_counts().to_dict()
        log.info("fold_vs_2nd>=%.1f: %d/%d variants lineage-specific; by lineage: %s",
                 fold, len(spec), len(out), by)
    log.info("done in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
