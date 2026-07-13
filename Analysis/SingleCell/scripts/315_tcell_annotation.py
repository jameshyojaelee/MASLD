#!/usr/bin/env python3
"""
315: T-cell subcluster annotation + GLP1R query (secondary)

Clone of 314 for the T-cell compartment. Scores CD8 T-cell exhaustion
(PDCD1, HAVCR2, TOX, NR4A2, CD8A + LAG3/TIGIT/CTLA4/ENTPD1) per Leiden
subcluster and queries the incretin/glucagon receptor family. GLP1R is not
expected on human T cells; this is a specificity check for the LSEC headline
and a look for any exhausted-CD8 GLP1R signal.

Inputs:
    {outdir}/tcell_subcluster.h5ad   (from 313)
    {outdir}/tcell_leiden_resolution_scores.csv

Outputs (to RNA-seq/results/glp1ra/lsec/):
    tcell_subcluster_glp1r.csv     per-subcluster CD8/exhaustion means + receptor %/CP10k/pseudobulk

Environment: spatial (CPU)
"""

import argparse
import logging
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

PANELS = {
    "CD8_exhaustion": ["PDCD1", "HAVCR2", "TOX", "NR4A2", "CD8A",
                       "LAG3", "TIGIT", "CTLA4", "ENTPD1"],
    "CD8_cytotoxic": ["CD8A", "CD8B", "GZMB", "GZMK", "NKG7", "PRF1", "GNLY"],
    "CD4_helper": ["CD4", "IL7R", "CCR7", "TCF7"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
}
RECEPTORS = ["GLP1R", "GIPR", "GCGR", "DPP4", "GLP2R"]
CONTROLS = ["CD3D", "CD8A", "PTPRC", "ALB"]  # T-cell high; ALB low = ambient check


def choose_resolution(adata, sil_df, min_frac=0.005, max_frac=0.70, n_lo=4, n_hi=15):
    n_total = adata.n_obs
    best_res, best_sil = None, -np.inf
    for _, row in sil_df.iterrows():
        res = row["resolution"]
        sizes = adata.obs[f"leiden_{res}"].value_counts()
        mn, mx = sizes.min() / n_total, sizes.max() / n_total
        ok = (mn >= min_frac) and (mx <= max_frac) and (n_lo <= row["n_clusters"] <= n_hi)
        log.info("  res=%.1f: %d clusters min_frac=%.3f max_frac=%.3f sil=%.4f -> %s",
                 res, row["n_clusters"], mn, mx, row["silhouette"], "PASS" if ok else "fail")
        if ok and row["silhouette"] > best_sil:
            best_sil, best_res = row["silhouette"], res
    if best_res is None:
        best_res = float(sil_df.loc[sil_df["silhouette"].idxmax(), "resolution"])
        log.warning("No resolution passed; best-silhouette fallback res=%.1f", best_res)
    return best_res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(
        BASE, "Analysis/SingleCell/results_gpu_v2/tcell_subtypes"))
    ap.add_argument("--tag", default="tcell")
    ap.add_argument("--results-dir", default=os.path.join(BASE, "RNA-seq/results/glp1ra/lsec"))
    ap.add_argument("--leiden-res", type=float, default=None)
    args = ap.parse_args()
    os.makedirs(args.results_dir, exist_ok=True)

    log.info("=" * 60)
    log.info("315: T-cell annotation + GLP1R query")
    log.info("=" * 60)

    h5ad = os.path.join(args.outdir, f"{args.tag}_subcluster.h5ad")
    adata = sc.read_h5ad(h5ad)
    log.info("Shape: %s", adata.shape)
    sil_df = pd.read_csv(os.path.join(args.outdir, f"{args.tag}_leiden_resolution_scores.csv"))

    res = args.leiden_res if args.leiden_res is not None else choose_resolution(adata, sil_df)
    adata.obs["subcluster"] = adata.obs[f"leiden_{res}"].astype(str)
    subclusters = sorted(adata.obs["subcluster"].unique(), key=lambda x: int(x))
    log.info("Resolution %.1f -> %d subclusters", res, len(subclusters))

    # Score panels (log1p-CP10k in .X)
    for name, genes in PANELS.items():
        present = [g for g in genes if g in adata.var_names]
        log.info("  %s: %d/%d present", name, len(present), len(genes))
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{name}")

    from scipy import sparse
    counts = adata.layers["counts"].tocsr()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0

    def gene_cp10k_pos(gene):
        if gene not in adata.var_names:
            return None, None
        j = adata.var_names.get_loc(gene)
        col = counts[:, j].toarray().ravel()
        return col / lib * 1e4, (col > 0)

    score_cols = [f"score_{n}" for n in PANELS]
    sc_means = adata.obs.groupby("subcluster")[score_cols].mean()

    rows = []
    for st in subclusters:
        m = (adata.obs["subcluster"] == st).values
        n = int(m.sum())
        rec = {"subcluster": st, "n_cells": n}
        for scn in score_cols:
            rec[scn] = float(sc_means.loc[st, scn])
        for g in RECEPTORS + CONTROLS:
            cp, pos = gene_cp10k_pos(g)
            if cp is None:
                rec[f"{g}_pct_pos"] = np.nan
                rec[f"{g}_mean_cp10k"] = np.nan
                rec[f"{g}_pseudobulk"] = np.nan
                rec[f"{g}_n_pos"] = np.nan
                continue
            npos = int(pos[m].sum())
            rec[f"{g}_n_pos"] = npos
            rec[f"{g}_pct_pos"] = npos / n if n else np.nan
            rec[f"{g}_mean_cp10k"] = float(cp[m].mean())
            j = adata.var_names.get_loc(g)
            rec[f"{g}_pseudobulk"] = float(counts[:, j].toarray().ravel()[m].sum())
        rows.append(rec)

    out_df = pd.DataFrame(rows)
    # crude subtype label = argmax panel score per subcluster
    lab = sc_means.idxmax(axis=1).str.replace("score_", "")
    out_df["dominant_panel"] = out_df["subcluster"].map(lab.to_dict())
    out_path = os.path.join(args.results_dir, "tcell_subcluster_glp1r.csv")
    out_df.to_csv(out_path, index=False)
    log.info("Saved %s", out_path)
    log.info("GLP1R by T subcluster:\n%s",
             out_df[["subcluster", "dominant_panel", "n_cells",
                     "GLP1R_n_pos", "GLP1R_pct_pos", "GLP1R_pseudobulk"]].to_string(index=False))

    # exhausted-CD8 vs rest GLP1R (whichever subcluster is most exhausted)
    exh_sub = sc_means["score_CD8_exhaustion"].idxmax()
    glp_cp, glp_pos = gene_cp10k_pos("GLP1R")
    if glp_pos is not None:
        m = (adata.obs["subcluster"] == exh_sub).values
        npos, n = int(glp_pos[m].sum()), int(m.sum())
        lo, hi = proportion_confint(npos, n, alpha=0.05, method="wilson") if n else (np.nan, np.nan)
        log.info("Most-exhausted-CD8 subcluster=%s: GLP1R+ %d/%d (%.3f%%, 95%% CI %.3f-%.3f%%)",
                 exh_sub, npos, n, 100 * npos / n, 100 * lo, 100 * hi)

    log.info("Done.")


if __name__ == "__main__":
    main()
