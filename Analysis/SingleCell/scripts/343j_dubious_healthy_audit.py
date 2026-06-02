#!/usr/bin/env python
# ============================================================================
# 343j_dubious_healthy_audit.py
#
# Audit donors labeled `disease_stage_coarse == "Healthy"` whose F-stage
# inference (scVI argmax, fib-signature, HSC counting, scVI-kNN) places them
# at F3/F4. Convergent failure across 3+ orthogonal methods is taken as
# evidence the donor is NOT biopsy-grade healthy tissue (peri-tumoral,
# transplant residual, etc.).
#
# Outputs:
#   - dubious_healthy_donors.tsv
#   - dubious_healthy_consensus.tsv
#   - dubious_vs_clean_healthy_expression.tsv
#   - figS_dubious_healthy_investigation.pdf  (4 panels)
#
# Per-cell mode "dominant_fstage" used to define dubious-healthy:
#   disease_stage_coarse == "Healthy" AND
#   (dominant_fstage >= 3 OR F4_frac > 0.4)
#
# Clean-healthy reference:
#   disease_stage_coarse == "Healthy" AND
#   dominant_fstage <= 1 AND F4_frac < 0.1
#
# Submitted via run_343j_dubious_healthy.sbatch:
#   cpu --qos=interactive --time=01:00:00 --mem=64G --cpus-per-task=8
# ============================================================================

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

ST_DIR     = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
ATLAS_H5AD = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"

DONOR_FRAC    = ST_DIR / "donor_fstage_cellfrac.tsv"
DONOR_META    = ST_DIR / "donor_metadata_extended.tsv"
DOC_F         = ST_DIR / "donor_fstage_documented.tsv"
SCVI_PRED     = ST_DIR / "donor_fstage_scvi_predicted.tsv"
FIBSIG_PRED   = ST_DIR / "donor_fstage_fibsignature_predicted.tsv"
HSC_PRED      = ST_DIR / "donor_fstage_hsc_predicted.tsv"
ALT_PRED      = ST_DIR / "donor_fstage_alt_predictions.tsv"

OUT_DUB       = ST_DIR / "dubious_healthy_donors.tsv"
OUT_CONS      = ST_DIR / "dubious_healthy_consensus.tsv"
OUT_EXPR      = ST_DIR / "dubious_vs_clean_healthy_expression.tsv"
FIG_DIR       = PROJECT_ROOT / "figures/supplementary/stage_ccc"
OUT_PDF       = FIG_DIR / "figS_dubious_healthy_investigation.pdf"
FIG_DIR.mkdir(parents=True, exist_ok=True)

FIBROSIS_MARKERS = ["ACTA2", "COL1A1", "COL3A1", "TIMP1", "LOX", "PDGFRB"]

# Cell types pooled for expression characterization (fibroblasts + hepatocytes)
EXPR_LINEAGES = ["Hepatocytes", "Fibroblasts"]

CONTROL_GRAY = "#9E9E9E"
DUB_COLOR    = "#D32F2F"   # crimson - alerts to "dubious"
F4_COLOR     = "#7B1FA2"   # purple - genuine F4
F0_COLOR     = "#1976D2"   # blue - clean F0 reference


# ---------------------------------------------------------------------------
def load_inputs():
    print(f"[load] donor_fstage_cellfrac: {DONOR_FRAC}")
    frac = pd.read_csv(DONOR_FRAC, sep="\t")
    print(f"        {len(frac)} donors")

    meta = pd.read_csv(DONOR_META, sep="\t", low_memory=False)
    print(f"[load] donor_metadata_extended: {len(meta)} donors")

    scvi = pd.read_csv(SCVI_PRED, sep="\t")
    fib  = pd.read_csv(FIBSIG_PRED, sep="\t")
    hsc  = pd.read_csv(HSC_PRED, sep="\t")
    alt  = pd.read_csv(ALT_PRED, sep="\t")
    return frac, meta, scvi, fib, hsc, alt


def flag_dubious(frac: pd.DataFrame) -> pd.DataFrame:
    """Healthy donors with cell-mode F3/F4 or F4_frac > 0.4."""
    h = frac[frac["disease_stage_coarse"].astype(str) == "Healthy"].copy()
    h["dominant_fstage"] = pd.to_numeric(h["dominant_fstage"], errors="coerce")
    h["F4_frac"] = pd.to_numeric(h["F4_frac"], errors="coerce")
    h["dubious_healthy"] = (
        (h["dominant_fstage"] >= 3) | (h["F4_frac"] > 0.4)
    )
    h["clean_healthy"] = (h["dominant_fstage"] <= 1) & (h["F4_frac"] < 0.1)
    return h


def consensus_table(dub_ids, scvi, fib, hsc, alt):
    """Per-dubious donor, count how many of 4 methods flag F=F4."""
    # scVI argmax (Andrews-trained, projected): flag if predicted = 4
    s = scvi.set_index("sample")["F_stage_predicted_argmax"]
    # Fib-sig classifier: flag if F_stage_fibsig = 4
    f = fib.set_index("sample")["F_stage_fibsig"]
    # HSC rule-based: flag if F_stage_hsc = 4
    h = hsc.set_index("sample")["F_stage_hsc"]
    # scVI kNN (alt predictions): flag if F_stage_knn_scvi = 4
    k = alt.set_index("sample")["F_stage_knn_scvi"]

    rows = []
    for sid in dub_ids:
        s_v = pd.to_numeric(s.get(sid, np.nan), errors="coerce")
        f_v = pd.to_numeric(f.get(sid, np.nan), errors="coerce")
        h_v = pd.to_numeric(h.get(sid, np.nan), errors="coerce")
        k_v = pd.to_numeric(k.get(sid, np.nan), errors="coerce")
        flags = {
            "flag_scvi_argmax": int(s_v == 4) if pd.notna(s_v) else np.nan,
            "flag_fibsig":      int(f_v == 4) if pd.notna(f_v) else np.nan,
            "flag_hsc":         int(h_v == 4) if pd.notna(h_v) else np.nan,
            "flag_scvi_knn":    int(k_v == 4) if pd.notna(k_v) else np.nan,
        }
        n_flag = int(np.nansum(list(flags.values())))
        n_meth = int(np.sum([pd.notna(v) for v in flags.values()]))
        rows.append({"sample": sid,
                     "scvi_argmax": s_v, "fibsig": f_v,
                     "hsc": h_v, "scvi_knn": k_v,
                     **flags,
                     "n_methods_flagging_F4": n_flag,
                     "n_methods_evaluated": n_meth})
    return pd.DataFrame(rows)


def donor_mean_expression(adata, samples, genes, lineages):
    """Per donor, mean log1p expression of `genes` pooled over `lineages`."""
    obs = adata.obs
    mask_ct = obs["cell_type"].astype(str).isin(lineages)
    mask_smp = obs["sample"].astype(str).isin(set(map(str, samples)))
    mask = mask_ct & mask_smp
    idx = np.where(mask.values)[0]
    print(f"  [expr] selected {len(idx):,} cells for {len(samples)} donors in {lineages}")
    if len(idx) == 0:
        return pd.DataFrame(columns=["sample"] + genes)

    var = adata.var
    gene_to_col = {}
    for g in genes:
        if g in var.index:
            gene_to_col[g] = var.index.get_loc(g)
        elif "gene_name" in var.columns and (var["gene_name"] == g).any():
            gene_to_col[g] = int(np.where(var["gene_name"].values == g)[0][0])
    print(f"  [expr] resolved {len(gene_to_col)}/{len(genes)} marker genes")

    # Read in chunks to avoid loading the full 40 GB matrix.
    chunk = 50_000
    samples_arr = obs["sample"].astype(str).values
    accum_sum = {g: {} for g in gene_to_col}
    accum_n   = {}
    for start in range(0, len(idx), chunk):
        sel = idx[start:start + chunk]
        sub = adata.X[sel, :] if not adata.isbacked else adata.X[sel, :]
        # Convert to dense for the marker columns only
        if hasattr(sub, "toarray"):
            for g, col in gene_to_col.items():
                vals = np.asarray(sub[:, col].todense()).ravel()
                ss = samples_arr[sel]
                for s_id, v in zip(ss, vals):
                    accum_sum[g].setdefault(s_id, 0.0)
                    accum_sum[g][s_id] += float(v)
        else:
            arr = np.asarray(sub)
            for g, col in gene_to_col.items():
                vals = arr[:, col]
                ss = samples_arr[sel]
                for s_id, v in zip(ss, vals):
                    accum_sum[g].setdefault(s_id, 0.0)
                    accum_sum[g][s_id] += float(v)
        ss = samples_arr[sel]
        for s_id in ss:
            accum_n[s_id] = accum_n.get(s_id, 0) + 1

    rows = []
    for s_id, n in accum_n.items():
        row = {"sample": s_id, "n_cells_pooled": n}
        for g in gene_to_col:
            row[g] = accum_sum[g].get(s_id, 0.0) / max(n, 1)
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    print("=" * 70)
    print("343j: dubious-healthy donor audit")
    print("=" * 70)

    frac, meta, scvi, fib, hsc, alt = load_inputs()

    # ---- 1) Define dubious-healthy ---------------------------------------
    h = flag_dubious(frac)
    print(f"\n[1] Healthy donors total: {len(h)}")
    print(f"    dubious-healthy: {int(h['dubious_healthy'].sum())}")
    print(f"    clean-healthy:   {int(h['clean_healthy'].sum())}")

    # ---- 2) Per-dataset tally -------------------------------------------
    tally = (h.groupby("dataset")
               .agg(n_healthy=("sample", "size"),
                    n_dubious=("dubious_healthy", "sum"),
                    n_clean=("clean_healthy", "sum"))
               .reset_index())
    tally["pct_dubious"] = 100 * tally["n_dubious"] / tally["n_healthy"].clip(lower=1)
    tally = tally.sort_values("n_dubious", ascending=False)
    print("\n[2] Per-dataset dubious-healthy tally:")
    print(tally.to_string(index=False))

    dub = h[h["dubious_healthy"]][
        ["sample", "dataset", "n_cells", "dominant_fstage",
         "F0_frac", "F1_frac", "F2_frac", "F3_frac", "F4_frac", "entropy"]
    ].sort_values(["dataset", "F4_frac"], ascending=[True, False])
    dub.to_csv(OUT_DUB, sep="\t", index=False)
    print(f"\n[2] wrote {OUT_DUB}  ({len(dub)} rows)")

    # ---- 3) Cross-method consensus --------------------------------------
    cons = consensus_table(dub["sample"].tolist(), scvi, fib, hsc, alt)
    cons = cons.merge(dub[["sample", "dataset", "dominant_fstage",
                           "F4_frac"]], on="sample", how="left")
    cons = cons[["sample", "dataset", "dominant_fstage", "F4_frac",
                 "scvi_argmax", "fibsig", "hsc", "scvi_knn",
                 "flag_scvi_argmax", "flag_fibsig",
                 "flag_hsc", "flag_scvi_knn",
                 "n_methods_flagging_F4", "n_methods_evaluated"]]
    cons.to_csv(OUT_CONS, sep="\t", index=False)
    n_ge3 = int((cons["n_methods_flagging_F4"] >= 3).sum())
    n_ge2 = int((cons["n_methods_flagging_F4"] >= 2).sum())
    print(f"\n[3] cross-method consensus written to {OUT_CONS}")
    print(f"    flagged F4 by >=3 methods: {n_ge3} / {len(cons)}")
    print(f"    flagged F4 by >=2 methods: {n_ge2} / {len(cons)}")

    # ---- 4) Biological characterization ---------------------------------
    print("\n[4] loading atlas (backed) for fibrosis marker pooling ...")
    adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
    print(f"     {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    if "cell_type" not in adata.obs.columns:
        sys.exit("FATAL: 'cell_type' missing from obs")

    # Donor groups for expression comparison
    dub_ids   = dub["sample"].tolist()
    clean_ids = h.loc[h["clean_healthy"], "sample"].tolist()

    # Andrews documented F0 / F4 (Cirrhosis ground truth)
    docf = pd.read_csv(DOC_F, sep="\t")
    docf["F_stage_documented"] = pd.to_numeric(
        docf["F_stage_documented"], errors="coerce")
    f4_ids = docf.loc[docf["F_stage_documented"] == 4, "sample"].tolist()
    f0_ids = docf.loc[docf["F_stage_documented"] == 0, "sample"].tolist()
    print(f"     n dubious={len(dub_ids)}  clean={len(clean_ids)}  "
          f"Andrews-F4={len(f4_ids)}  Andrews-F0={len(f0_ids)}")

    all_target_ids = set(dub_ids) | set(clean_ids) | set(f4_ids) | set(f0_ids)
    expr = donor_mean_expression(adata, all_target_ids,
                                 FIBROSIS_MARKERS, EXPR_LINEAGES)

    def tag(sid):
        if sid in set(dub_ids):   return "dubious_healthy"
        if sid in set(clean_ids): return "clean_healthy"
        if sid in set(f4_ids):    return "andrews_F4"
        if sid in set(f0_ids):    return "andrews_F0"
        return "other"
    expr["group"] = expr["sample"].map(tag)
    expr = expr[expr["group"] != "other"]
    print(f"     expression matrix: {len(expr)} donors x "
          f"{len(FIBROSIS_MARKERS)} markers")

    # Per-marker Wilcoxon dubious vs clean
    wilcox = []
    for g in FIBROSIS_MARKERS:
        if g not in expr.columns:
            continue
        a_vals = expr.loc[expr["group"] == "dubious_healthy", g].dropna().values
        b_vals = expr.loc[expr["group"] == "clean_healthy", g].dropna().values
        if len(a_vals) >= 3 and len(b_vals) >= 3:
            stat, pval = stats.mannwhitneyu(a_vals, b_vals,
                                            alternative="two-sided")
            wilcox.append({"gene": g,
                           "mean_dubious": float(np.mean(a_vals)),
                           "mean_clean":   float(np.mean(b_vals)),
                           "n_dubious": len(a_vals),
                           "n_clean":   len(b_vals),
                           "wilcoxon_U": float(stat),
                           "wilcoxon_p": float(pval)})
    wdf = pd.DataFrame(wilcox)
    print("\n[4] dubious-vs-clean fibrosis marker Wilcoxon:")
    if not wdf.empty:
        print(wdf.to_string(index=False))

    # Cell-type composition: fibroblast / macrophage / hepatocyte fractions
    comp_cols = [c for c in meta.columns
                 if c.startswith("frac_") or c == "n_cells"]
    comp = meta[["sample", "dataset"] + comp_cols].copy()
    comp["group"] = comp["sample"].map(tag)
    comp = comp[comp["group"].isin(
        ["dubious_healthy", "clean_healthy", "andrews_F4", "andrews_F0"])]

    composition_summary = (comp.groupby("group")[comp_cols].mean()
                                .reset_index())
    print("\n[4] mean lineage fraction by group:")
    print(composition_summary.to_string(index=False))

    # Wilcoxon on fibroblast fraction
    fib_col = "frac_Fibroblasts"
    fib_w = {}
    if fib_col in comp.columns:
        a_v = comp.loc[comp["group"] == "dubious_healthy",
                       fib_col].dropna().values
        b_v = comp.loc[comp["group"] == "clean_healthy",
                       fib_col].dropna().values
        if len(a_v) >= 3 and len(b_v) >= 3:
            U, p = stats.mannwhitneyu(a_v, b_v, alternative="two-sided")
            fib_w = {"n_dubious": len(a_v), "n_clean": len(b_v),
                     "wilcoxon_U": float(U), "wilcoxon_p": float(p),
                     "mean_dubious": float(np.mean(a_v)),
                     "mean_clean":   float(np.mean(b_v))}
            print(f"\n[4] fibroblast fraction Wilcoxon: U={U:.1f} p={p:.3g}  "
                  f"(dub mean={np.mean(a_v):.4f}, clean mean={np.mean(b_v):.4f})")

    # ---- 5) Comparison to Andrews F4 vs F0 ------------------------------
    # For each dubious donor: Spearman rho across the 6-marker profile vs
    # mean F4 profile and mean F0 profile.  Paired test: is rho_F4 > rho_F0?
    print("\n[5] correlation contrast: dubious vs Andrews F4 / F0 reference")
    f4_profile = expr.loc[expr["group"] == "andrews_F4",
                          FIBROSIS_MARKERS].mean(axis=0).values
    f0_profile = expr.loc[expr["group"] == "andrews_F0",
                          FIBROSIS_MARKERS].mean(axis=0).values

    corr_rows = []
    for sid in dub_ids:
        row = expr.loc[expr["sample"] == sid, FIBROSIS_MARKERS]
        if row.empty:
            continue
        v = row.values.ravel().astype(float)
        if np.all(np.isnan(v)) or len(v) < 3:
            continue
        rho_f4 = stats.spearmanr(v, f4_profile, nan_policy="omit").correlation
        rho_f0 = stats.spearmanr(v, f0_profile, nan_policy="omit").correlation
        corr_rows.append({"sample": sid, "rho_andrews_F4": rho_f4,
                          "rho_andrews_F0": rho_f0,
                          "delta_F4_minus_F0": rho_f4 - rho_f0})
    cdf = pd.DataFrame(corr_rows)
    if not cdf.empty:
        wstat, wp = stats.wilcoxon(cdf["rho_andrews_F4"].fillna(0),
                                   cdf["rho_andrews_F0"].fillna(0),
                                   alternative="greater")
        print(f"     n={len(cdf)}  median rho_F4={cdf['rho_andrews_F4'].median():.3f}"
              f"  median rho_F0={cdf['rho_andrews_F0'].median():.3f}")
        print(f"     Wilcoxon (rho_F4 > rho_F0): W={wstat:.1f}, p={wp:.3g}")
        corr_summary = {"n_dubious": len(cdf),
                        "median_rho_F4": float(cdf['rho_andrews_F4'].median()),
                        "median_rho_F0": float(cdf['rho_andrews_F0'].median()),
                        "wilcoxon_W": float(wstat),
                        "wilcoxon_p_F4_gt_F0": float(wp)}
    else:
        corr_summary = {}
        print("     no dubious donors with sufficient expression data")

    # ---- Combined expression export -------------------------------------
    expr_out = expr.merge(cdf, on="sample", how="left")
    expr_out.to_csv(OUT_EXPR, sep="\t", index=False)
    print(f"\n[4-5] wrote {OUT_EXPR}")

    # ---- 6) Figure: 4 panels --------------------------------------------
    print(f"\n[6] writing PDF: {OUT_PDF}")
    with PdfPages(OUT_PDF) as pdf:
        # Panel A: per-dataset dubious-healthy bar
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        idx = np.arange(len(tally))
        ax.barh(idx, tally["n_healthy"], color=CONTROL_GRAY,
                label="clean / borderline Healthy", height=0.75)
        ax.barh(idx, tally["n_dubious"], color=DUB_COLOR,
                label="dubious-healthy (dom>=3 or F4>0.4)", height=0.75)
        for i, (n, pct) in enumerate(zip(tally["n_dubious"],
                                         tally["pct_dubious"])):
            if n > 0:
                ax.text(n + 0.4, i, f"{int(n)} ({pct:.0f}%)",
                        va="center", fontsize=8)
        ax.set_yticks(idx)
        ax.set_yticklabels(tally["dataset"])
        ax.set_xlabel("n Healthy-labeled donors")
        ax.set_title("A. Dubious-healthy donors per dataset\n"
                     "(Healthy label + per-cell F-mode >= F3 OR F4_frac > 0.4)")
        ax.legend(loc="lower right", fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Panel B: fibrosis markers, 4 groups
        fig, ax = plt.subplots(figsize=(7.5, 4.4))
        groups_order = ["andrews_F0", "clean_healthy",
                        "dubious_healthy", "andrews_F4"]
        colors_g = {"andrews_F0": F0_COLOR,
                    "clean_healthy": CONTROL_GRAY,
                    "dubious_healthy": DUB_COLOR,
                    "andrews_F4": F4_COLOR}
        positions = np.arange(len(FIBROSIS_MARKERS))
        width = 0.20
        for j, grp in enumerate(groups_order):
            vals = [expr.loc[expr["group"] == grp, g].mean()
                    if g in expr.columns else np.nan
                    for g in FIBROSIS_MARKERS]
            ax.bar(positions + (j - 1.5) * width, vals, width,
                   label=grp, color=colors_g[grp])
        ax.set_xticks(positions)
        ax.set_xticklabels(FIBROSIS_MARKERS, rotation=30, ha="right")
        ax.set_ylabel("mean log1p expression\n(hepatocytes + fibroblasts pooled)")
        ax.set_title("B. Fibrosis markers across donor groups")
        ax.legend(fontsize=8, frameon=False, ncol=2)
        ax.spines[["top", "right"]].set_visible(False)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Panel C: fibroblast fraction box
        fig, ax = plt.subplots(figsize=(5.5, 4.0))
        comp_plot = comp[comp["group"].isin(groups_order)].copy()
        data_box, labels_box, colors_box = [], [], []
        for grp in groups_order:
            v = comp_plot.loc[comp_plot["group"] == grp,
                              "frac_Fibroblasts"].dropna().values
            if len(v) > 0:
                data_box.append(v)
                labels_box.append(f"{grp}\n(n={len(v)})")
                colors_box.append(colors_g[grp])
        bp = ax.boxplot(data_box, patch_artist=True, widths=0.6,
                        showfliers=False)
        for patch, c in zip(bp["boxes"], colors_box):
            patch.set_facecolor(c); patch.set_alpha(0.7)
        for i, v in enumerate(data_box):
            ax.scatter(np.full_like(v, i + 1, dtype=float)
                       + np.random.uniform(-0.08, 0.08, len(v)),
                       v, color="black", s=6, alpha=0.5)
        ax.set_xticklabels(labels_box, fontsize=8)
        ax.set_ylabel("Fibroblast fraction")
        title_c = "C. Fibroblast fraction by group"
        if fib_w:
            title_c += f"\n(dub vs clean Wilcoxon p={fib_w['wilcoxon_p']:.2g})"
        ax.set_title(title_c)
        ax.spines[["top", "right"]].set_visible(False)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Panel D: rho-F4 vs rho-F0 for each dubious donor
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        if not cdf.empty:
            x = cdf["rho_andrews_F0"].values
            y = cdf["rho_andrews_F4"].values
            ax.scatter(x, y, color=DUB_COLOR, s=22, alpha=0.7,
                       edgecolors="black", linewidths=0.4)
            lo, hi = -1.05, 1.05
            ax.plot([lo, hi], [lo, hi], color="gray", lw=0.8, ls="--")
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
            ax.set_xlabel("Spearman rho vs Andrews-F0 profile")
            ax.set_ylabel("Spearman rho vs Andrews-F4 profile")
            ax.set_title("D. Dubious-healthy donor profile correlation\n"
                         f"(n={len(cdf)}; one-sided Wilcoxon "
                         f"rho_F4>rho_F0 p={corr_summary.get('wilcoxon_p_F4_gt_F0', np.nan):.2g})")
            ax.spines[["top", "right"]].set_visible(False)
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"[6] PDF written")

    # ---- 7) Patch donor_metadata_extended.tsv with clean_healthy_flag ----
    # Logic mirrors the definition above so Script 344 can reproduce it.
    print("\n[7] adding clean_healthy_flag to donor_metadata_extended.tsv")
    meta = meta.merge(
        h[["sample", "dataset", "dominant_fstage", "F4_frac",
           "dubious_healthy", "clean_healthy"]],
        on=["sample", "dataset"], how="left",
        suffixes=("", "_fstage"))
    meta["clean_healthy_flag"] = meta["clean_healthy"].fillna(False).astype(bool)
    meta["dubious_healthy_flag"] = meta["dubious_healthy"].fillna(False).astype(bool)
    meta = meta.drop(columns=["dubious_healthy", "clean_healthy"])
    meta.to_csv(DONOR_META, sep="\t", index=False)
    print(f"     n clean_healthy_flag = TRUE:    "
          f"{int(meta['clean_healthy_flag'].sum())}")
    print(f"     n dubious_healthy_flag = TRUE:  "
          f"{int(meta['dubious_healthy_flag'].sum())}")
    print(f"     wrote {DONOR_META}")

    print("\n[done]")


if __name__ == "__main__":
    main()
