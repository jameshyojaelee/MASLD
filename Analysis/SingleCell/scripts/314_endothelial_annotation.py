#!/usr/bin/env python3
"""
314: Endothelial LSEC subtype annotation + incretin-receptor query (HEADLINE)

Reads the 313 endothelial subcluster h5ad (full gene set), assigns each Leiden
subcluster to pericentral-LSEC / periportal-LSEC / vascular-or-other-EC using the
curated Guilliams-2022 LSEC zonation panels, and queries the incretin/glucagon
receptor family (GLP1R, GIPR, GCGR, DPP4, GLP2R) per subcluster.

Motivation (Drucker, Cell Metab 2026): in MOUSE, Glp1r localizes to PERICENTRAL
LSECs. This script tests the HUMAN snRNA/scRNA replication — with pre-registered
honesty about the detection floor (human droplet data is far less sensitive than
the GEM-X Flex chemistry Drucker used; GLP1R is a low-abundance GPCR).

Pre-registered verdict logic (per teammate brief):
  * If GLP1R is enriched in pericentral-LSEC vs other EC (Fisher exact on GLP1R+
    counts, one-sided greater) -> HUMAN REPLICATION.
  * If GLP1R sits at the detection floor even in pericentral LSEC (a handful of
    positive cells, CI overlapping other-EC) -> INCONCLUSIVE / underpowered,
    NOT a refutation of Drucker.

Inputs:
    {outdir}/endothelial_subcluster.h5ad          (from 313)
    {outdir}/endothelial_leiden_resolution_scores.csv

Outputs (to RNA-seq/results/glp1ra/lsec/):
    endothelial_subcluster_receptors.csv   per-subcluster receptor %/CP10k/pseudobulk
    endothelial_lsec_sublabel.csv          per-cell barcode -> subcluster -> LSEC_subtype
    endothelial_glp1r_pericentral_vs_other.csv  the Fisher test + Wilson CIs
    endothelial_subcluster_marker_means.csv  marker validation table
    endothelial_umap_markers.pdf           UMAP + marker/receptor dotplot (optional)

Environment: spatial (CPU)
"""

import argparse
import logging
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

# ---------------------------------------------------------------------------
# Curated Guilliams-2022 LSEC zonation panels (from
# Analysis/SingleCell/scripts/build_liver_celltype_markers.py). CD32B == FCGR2B.
# ---------------------------------------------------------------------------
PANELS = {
    # pericentral / central-vein-proximal LSEC (scavenger, fenestrated)
    "LSEC_central": ["CLEC4G", "CLEC4M", "STAB2", "STAB1", "FCN3", "OIT3",
                     "F8", "PLPP1", "FCGR2B", "DNASE1L3", "CTSL", "LYVE1", "FCN2"],
    # periportal LSEC (more continuous-EC-like)
    "LSEC_portal": ["VWF", "PECAM1", "CD34", "CDH5", "ENG", "SELP",
                    "ICAM2", "TM4SF1", "RAMP3", "NTS", "MGP", "RBP7",
                    "DLL4", "GJA5", "EFNB2"],
    # pan-LSEC scavenger IDENTITY (separates true LSEC from macrovascular EC)
    "LSEC_identity": ["CLEC4G", "CLEC4M", "STAB2", "STAB1", "FCN2", "FCN3",
                      "OIT3", "FCGR2B", "DNASE1L3", "LYVE1"],
    # macrovascular / non-sinusoidal EC (portal vein, central vein, arterial, lymphatic)
    "vascular_general": ["WNT2", "PECAM1", "CDH5", "VWF", "CD34", "DLL4",
                         "GJA5", "EFNB2", "SELE", "VCAM1"],
}

RECEPTORS = ["GLP1R", "GIPR", "GCGR", "DPP4", "GLP2R"]
# PECAM1/VWF/STAB2 high; ALB low = ambient/hepatocyte-contamination check
CONTROLS = ["PECAM1", "VWF", "STAB2", "ALB"]


def choose_resolution(adata, sil_df, min_frac=0.005, max_frac=0.70, n_lo=4, n_hi=15):
    """Highest-silhouette resolution with a sane cluster structure."""
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
        log.warning("No resolution passed; using best-silhouette fallback res=%.1f", best_res)
    return best_res


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else x - x.mean()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(
        BASE, "Analysis/SingleCell/results_gpu_v2/endothelial_subtypes"))
    ap.add_argument("--tag", default="endothelial")
    ap.add_argument("--results-dir", default=os.path.join(BASE, "RNA-seq/results/glp1ra/lsec"))
    ap.add_argument("--leiden-res", type=float, default=None,
                    help="force a resolution instead of auto-select")
    args = ap.parse_args()
    os.makedirs(args.results_dir, exist_ok=True)

    log.info("=" * 60)
    log.info("314: Endothelial LSEC annotation + receptor query")
    log.info("=" * 60)

    h5ad = os.path.join(args.outdir, f"{args.tag}_subcluster.h5ad")
    log.info("Loading %s ...", h5ad)
    adata = sc.read_h5ad(h5ad)
    log.info("Shape: %s", adata.shape)

    sil_df = pd.read_csv(os.path.join(args.outdir, f"{args.tag}_leiden_resolution_scores.csv"))
    log.info("Resolution scores:\n%s", sil_df.to_string(index=False))

    # ── 1. Choose Leiden resolution ──────────────────────────────────────
    res = args.leiden_res if args.leiden_res is not None else choose_resolution(adata, sil_df)
    leiden_key = f"leiden_{res}"
    adata.obs["subcluster"] = adata.obs[leiden_key].astype(str)
    subclusters = sorted(adata.obs["subcluster"].unique(), key=lambda x: int(x))
    log.info("Selected resolution %.1f -> %d subclusters: %s", res, len(subclusters), subclusters)

    # ── 2. Score LSEC zonation panels (log1p-CP10k in .X) ────────────────
    log.info("Scoring panels ...")
    for name, genes in PANELS.items():
        present = [g for g in genes if g in adata.var_names]
        log.info("  %s: %d/%d genes present: %s", name, len(present), len(genes), present)
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{name}")

    # z-score each panel across cells so cross-panel comparison is scale-free
    for name in PANELS:
        adata.obs[f"z_{name}"] = zscore(adata.obs[f"score_{name}"].values)
    # continuous pericentral zonation axis (positive = pericentral)
    adata.obs["pericentral_zonation"] = adata.obs["z_LSEC_central"] - adata.obs["z_LSEC_portal"]

    # ── 3. Assign subcluster LSEC subtype ────────────────────────────────
    # LSEC vs macrovascular: scavenger-LSEC identity vs vascular-general identity.
    # Among LSEC subclusters: pericentral if z_central > z_portal else periportal.
    sc_mean = adata.obs.groupby("subcluster")[
        ["z_LSEC_identity", "z_vascular_general", "z_LSEC_central", "z_LSEC_portal",
         "pericentral_zonation"]
    ].mean()

    assign = {}
    for st in subclusters:
        row = sc_mean.loc[st]
        is_lsec = row["z_LSEC_identity"] > row["z_vascular_general"]
        if not is_lsec:
            assign[st] = "vascular_other_EC"
        else:
            assign[st] = "pericentral_LSEC" if row["z_LSEC_central"] > row["z_LSEC_portal"] else "periportal_LSEC"
        log.info("  subcluster %s: z_id=%.2f z_vasc=%.2f z_cen=%.2f z_por=%.2f -> %s",
                 st, row["z_LSEC_identity"], row["z_vascular_general"],
                 row["z_LSEC_central"], row["z_LSEC_portal"], assign[st])
    adata.obs["lsec_subtype"] = adata.obs["subcluster"].map(assign).astype(str)

    log.info("LSEC subtype cell counts:\n%s",
             adata.obs["lsec_subtype"].value_counts().to_string())

    # ── 4. Marker-mean validation table (raw-count-derived CP10k means) ──
    counts = adata.layers["counts"].tocsr()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0

    def gene_cp10k_and_pos(gene):
        """Return (per-cell CP10k vector, per-cell positive boolean) for a gene."""
        if gene not in adata.var_names:
            return None, None
        j = adata.var_names.get_loc(gene)
        col = counts[:, j].toarray().ravel()
        cp10k = col / lib * 1e4
        return cp10k, (col > 0)

    marker_genes = sorted(set(
        [g for gs in PANELS.values() for g in gs] + RECEPTORS + CONTROLS
    ))
    mm_rows = []
    for st in subclusters:
        m = (adata.obs["subcluster"] == st).values
        rec = {"subcluster": st, "lsec_subtype": assign[st], "n_cells": int(m.sum())}
        for g in marker_genes:
            cp, _ = gene_cp10k_and_pos(g)
            rec[g] = float(cp[m].mean()) if cp is not None else np.nan
        mm_rows.append(rec)
    mm_df = pd.DataFrame(mm_rows)
    mm_path = os.path.join(args.results_dir, "endothelial_subcluster_marker_means.csv")
    mm_df.to_csv(mm_path, index=False)
    log.info("Saved marker means to %s", mm_path)

    # ── 5. Receptor query per subcluster ─────────────────────────────────
    log.info("Receptor query ...")
    query_genes = RECEPTORS + CONTROLS
    rq_rows = []
    for st in subclusters:
        m = (adata.obs["subcluster"] == st).values
        n = int(m.sum())
        for g in query_genes:
            cp, pos = gene_cp10k_and_pos(g)
            if cp is None:
                rq_rows.append({"subcluster": st, "lsec_subtype": assign[st], "gene": g,
                                "n_cells": n, "n_pos": np.nan, "pct_pos": np.nan,
                                "pct_ci_lo": np.nan, "pct_ci_hi": np.nan,
                                "mean_cp10k": np.nan, "pseudobulk_sum_counts": np.nan})
                continue
            npos = int(pos[m].sum())
            pct = npos / n if n else np.nan
            lo, hi = proportion_confint(npos, n, alpha=0.05, method="wilson") if n else (np.nan, np.nan)
            j = adata.var_names.get_loc(g)
            pb = float(counts[:, j].toarray().ravel()[m].sum())
            rq_rows.append({
                "subcluster": st, "lsec_subtype": assign[st], "gene": g,
                "n_cells": n, "n_pos": npos, "pct_pos": pct,
                "pct_ci_lo": lo, "pct_ci_hi": hi,
                "mean_cp10k": float(cp[m].mean()), "pseudobulk_sum_counts": pb,
            })
    rq_df = pd.DataFrame(rq_rows)
    rq_path = os.path.join(args.results_dir, "endothelial_subcluster_receptors.csv")
    rq_df.to_csv(rq_path, index=False)
    log.info("Saved receptor table to %s", rq_path)
    log.info("GLP1R by subcluster:\n%s",
             rq_df[rq_df.gene == "GLP1R"][
                 ["subcluster", "lsec_subtype", "n_cells", "n_pos", "pct_pos",
                  "mean_cp10k", "pseudobulk_sum_counts"]].to_string(index=False))

    # ── 6. HEADLINE: GLP1R pericentral-LSEC vs other-EC ──────────────────
    log.info("GLP1R pericentral-LSEC vs other-EC ...")
    glp_cp, glp_pos = gene_cp10k_and_pos("GLP1R")
    peri = (adata.obs["lsec_subtype"] == "pericentral_LSEC").values
    other = ~peri  # every non-pericentral EC (periportal + vascular/other)

    def group_stat(mask, label):
        n = int(mask.sum())
        npos = int(glp_pos[mask].sum())
        lo, hi = proportion_confint(npos, n, alpha=0.05, method="wilson") if n else (np.nan, np.nan)
        pb = float(counts[:, adata.var_names.get_loc("GLP1R")].toarray().ravel()[mask].sum())
        # donor-level: how many donors contribute >=1 GLP1R+ cell
        donors = adata.obs.loc[mask, "sample"].astype(str)
        pos_donors = adata.obs.loc[mask & glp_pos, "sample"].astype(str).nunique() if npos else 0
        return {"group": label, "n_cells": n, "n_glp1r_pos": npos,
                "pct_pos": npos / n if n else np.nan, "pct_ci_lo": lo, "pct_ci_hi": hi,
                "mean_cp10k": float(glp_cp[mask].mean()) if n else np.nan,
                "pseudobulk_sum_counts": pb,
                "n_donors_total": int(donors.nunique()), "n_donors_glp1r_pos": int(pos_donors)}

    g_peri = group_stat(peri, "pericentral_LSEC")
    g_other = group_stat(other, "other_EC")
    # 2x2 Fisher (one-sided greater: pericentral enriched for GLP1R+)
    a, b = g_peri["n_glp1r_pos"], g_peri["n_cells"] - g_peri["n_glp1r_pos"]
    c, d = g_other["n_glp1r_pos"], g_other["n_cells"] - g_other["n_glp1r_pos"]
    orr, p_two = fisher_exact([[a, b], [c, d]], alternative="two-sided")
    _, p_greater = fisher_exact([[a, b], [c, d]], alternative="greater")

    verdict_df = pd.DataFrame([g_peri, g_other])
    verdict_df["fisher_or"] = orr
    verdict_df["fisher_p_two_sided"] = p_two
    verdict_df["fisher_p_greater"] = p_greater
    v_path = os.path.join(args.results_dir, "endothelial_glp1r_pericentral_vs_other.csv")
    verdict_df.to_csv(v_path, index=False)
    log.info("GLP1R pericentral vs other:\n%s", verdict_df.to_string(index=False))
    log.info("Fisher OR=%.3f  p_two=%.3g  p_greater=%.3g", orr, p_two, p_greater)

    # ── 7. Per-cell sublabel for downstream CCC relabel ──────────────────
    sub = pd.DataFrame({
        "cell_barcode": adata.obs_names.astype(str),
        "dataset": adata.obs["dataset"].astype(str).values,
        "sample": adata.obs["sample"].astype(str).values,
        "condition": adata.obs["condition"].astype(str).values,
        "subcluster": adata.obs["subcluster"].values,
        "lsec_subtype": adata.obs["lsec_subtype"].values,
        "pericentral_zonation": adata.obs["pericentral_zonation"].values,
    })
    sub_path = os.path.join(args.results_dir, "endothelial_lsec_sublabel.csv")
    sub.to_csv(sub_path, index=False)
    log.info("Saved sublabel to %s (%d cells)", sub_path, len(sub))

    # ── 8. Optional UMAP + dotplot PDF ───────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.umap(adata, color="lsec_subtype", ax=axes[0], show=False, size=3, frameon=False)
        sc.pl.umap(adata, color="pericentral_zonation", ax=axes[1], show=False, size=3,
                   frameon=False, cmap="RdBu_r")
        fig.tight_layout()
        pdf1 = os.path.join(args.results_dir, "endothelial_umap_subtypes.pdf")
        fig.savefig(pdf1, bbox_inches="tight"); plt.close(fig)

        dot_genes = ["CLEC4G", "STAB2", "FCN2", "OIT3", "VWF", "CD34", "PECAM1",
                     "WNT2", "ALB"] + RECEPTORS
        dot_genes = [g for g in dot_genes if g in adata.var_names]
        ax = sc.pl.dotplot(adata, dot_genes, groupby="lsec_subtype", show=False, return_fig=True)
        pdf2 = os.path.join(args.results_dir, "endothelial_markers_dotplot.pdf")
        ax.savefig(pdf2, bbox_inches="tight")
        log.info("Saved UMAP + dotplot PDFs")
    except Exception as e:
        log.warning("Figure step skipped: %s", e)

    log.info("Done.")


if __name__ == "__main__":
    main()
