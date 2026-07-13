#!/usr/bin/env python3
"""
42b_govaere2026_cosmx_de.py — Per-cell-type differential expression on the
Govaere 2026 CosMx SMI 1000-plex deposit (GSE312698, Leuven_1..4 = 522,145
cells × 968 target genes after QC).

Inputs:
  Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad  (built by 41_*)

Disease-stage assignment (slide-level; documented in W2.B loader):
  Leuven_1 (GSM9351704)  "Explant_x2"          → 2 end-stage MASH explants  → MASH
  Leuven_2 (GSM9351705)  "F3_MASL_Normal_F3"   → 1 Normal + 1 MASL + 2 F3   → no_MASH (mixed-with-normal)
  Leuven_3 (GSM9351706)  "F1_F4"               → F1 + F4                    → MASH
  Leuven_4 (GSM9351707)  "F3_F3_MASL"          → 2 F3 + 1 MASL              → MASH

The GEO Sample_description does not ship a FOV-to-tissue map, so we cannot
isolate the Normal cells inside Leuven_2. We therefore label Leuven_2 as the
no-MASH comparator slide (mixed-with-normal). This is conservative — it dilutes
MASH-up genes (Leuven_2 has 3 MASH tissues mixed in) and is documented here as
a noisy MASH-enrichment contrast, NOT a clean MASH-vs-normal contrast. The
primary value lies in the per-cell-type direction + recovery of paper markers.

Outputs (Analysis/Spatial/results/govaere2026/):
  cosmx_de_<cell_type>_MASH_vs_noMASH.csv   — SLIDE-level direction per cell type.
      A6 FIX (2026-06-20): the cell-level Wilcoxon over ~297K cells in a 3-vs-1
      SLIDE design is pseudoreplicated (cells within a slide are not independent),
      and produced 779-883 anti-conservative "sig" hep/KC genes. With a 3-vs-1
      slide design NO inferential p-value is valid. We therefore NO LONGER write a
      cell-level p-value column. The CSV now carries:
        • logfoldchange  = SLIDE-level direction (mean MASH slides − mean no_MASH
                           slides, log1p-CPM scale) — the experimental-unit effect.
        • pval_adj       = all-NaN sentinel. The atlas consumer (06_integration.py
                           → spatial_govaere2026_cosmx_<ct>_mash_padj) maps this
                           column, so the atlas now carries the slide DIRECTION but
                           never a threshold-passable significance value.
        • slide_logfc / slide_direction_concordant / n_slides_*  = slide evidence.
        • cell_logfoldchange / cell_score  = DESCRIPTIVE cell-level ranking only
                           (drives the MetMac marker-recovery check); NOT significance.
  cosmx_de_summary.tsv                       — per-contrast n_genes + slide-direction
                                               up/down + cell/slide concordance
  cosmx_celltype_composition.tsv             — sample × cell_type cell counts + proportions

Method:
  • sc.tl.rank_genes_groups(method="wilcoxon", reference="no_MASH") inside each
    cell-type subset with ≥30 cells per group.
  • The h5ad ships log-normalised counts at .layers["lognorm"]; we use that for
    Wilcoxon (target_sum=1e4, log1p) — matches the loader's settings.
  • LFC is the scanpy `logfoldchanges` (natural-log of mean+1 ratio).

Validation:
  • Print MetMac/KC composition shift across slides.
  • Print top 10 macrophage-cluster MASH-up genes (here all macrophage lineage
    cells fall into the "KC" cluster — see loader notes; MetMac/TransMac were
    not separated by leiden res=0.5 on this panel).
  • Check whether canonical MetMac markers (GPNMB, LPL, FABP5, HS3ST2) sit in
    the top 50 KC-cluster MASH-up genes — this is the published lipid-associated
    macrophage signature regardless of cluster label.

Environment:
  micromamba activate spatial

SLURM:
  --partition=io --qos=interactive --cpus-per-task=4 --mem=64G --time=4:00:00
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Dict, List

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Donor/slide-aware helpers (F148): slide-level pseudobulk for the MASH contrast.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from spatial_stats import pseudobulk_by_donor  # noqa: E402

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ADATA_PATH   = PROJECT_ROOT / "Analysis/Spatial/results/preprocessed/cosmx_govaere2026.h5ad"
OUT_DIR      = PROJECT_ROOT / "Analysis/Spatial/results/govaere2026"

# ── Slide-level disease assignment ───────────────────────────────────────────
# Anchor: paper Methods (10 tissues total: 9 MASLD + 1 Normal) + GEO sample
# descriptions retained in .obs["sample_disease"].
SAMPLE_DISEASE = {
    "Leuven_1": {"sample_disease_str": "Explant_x2",        "n_tissues": 2, "n_mash": 2, "n_normal": 0, "stage": "MASH"},
    "Leuven_2": {"sample_disease_str": "F3_MASL_Normal_F3", "n_tissues": 4, "n_mash": 3, "n_normal": 1, "stage": "no_MASH"},
    "Leuven_3": {"sample_disease_str": "F1_F4",             "n_tissues": 2, "n_mash": 2, "n_normal": 0, "stage": "MASH"},
    "Leuven_4": {"sample_disease_str": "F3_F3_MASL",        "n_tissues": 3, "n_mash": 3, "n_normal": 0, "stage": "MASH"},
}

# Cell types we want to test if present (n>=30 per group). Paper's
# macrophage subclusters MetMac/TransMac/preMac are listed; if absent in the
# h5ad (loader's leiden + argmax merged them), we run DE for the cell types
# we do have and document the absence.
TARGET_CELLTYPES = ["MetMac", "KC", "Hepatocyte", "TransMac"]

# Canonical MetMac markers (Govaere 2026 Fig 1d + Sup Table 11)
METMAC_MARKERS = ["GPNMB", "LPL", "FABP5", "HS3ST2", "TREM2", "CD9"]

MIN_CELLS_PER_GROUP = 30


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def assign_disease_stage(adata: ad.AnnData) -> ad.AnnData:
    """Map .obs.sample_id → MASH / no_MASH using the slide-level table above."""
    if "sample_id" not in adata.obs.columns:
        sys.exit("ERROR: .obs.sample_id missing from h5ad")

    stage_map = {sid: meta["stage"] for sid, meta in SAMPLE_DISEASE.items()}
    adata.obs["disease_stage"] = adata.obs["sample_id"].map(stage_map).astype("category")
    if adata.obs["disease_stage"].isna().any():
        bad = adata.obs.loc[adata.obs["disease_stage"].isna(), "sample_id"].unique()
        sys.exit(f"ERROR: unmapped sample_ids: {bad}")
    _log("Disease-stage assignment (slide-level):")
    for sid, meta in SAMPLE_DISEASE.items():
        n = int((adata.obs["sample_id"] == sid).sum())
        _log(f"  {sid:>9} ({meta['sample_disease_str']:>20}) → {meta['stage']:>8}   "
             f"({n:>7,} cells; {meta['n_mash']} MASH / {meta['n_normal']} normal tissues per slide)")
    return adata


def per_celltype_de(adata: ad.AnnData, cell_type: str,
                    out_dir: pathlib.Path) -> Dict[str, object]:
    """Wilcoxon DE for one cell type: MASH vs no_MASH."""
    mask = adata.obs["cell_type"] == cell_type
    n_total = int(mask.sum())
    if n_total == 0:
        _log(f"  [{cell_type}] SKIP — 0 cells in atlas")
        return {"cell_type": cell_type, "skipped": True, "reason": "absent",
                "n_cells_mash": 0, "n_cells_nomash": 0,
                "n_slide_up": 0, "n_slide_down": 0, "n_direction_concordant": 0}
    sub = adata[mask].copy()
    n_mash   = int((sub.obs["disease_stage"] == "MASH").sum())
    n_nomash = int((sub.obs["disease_stage"] == "no_MASH").sum())
    _log(f"  [{cell_type}] subset {n_total:,} cells "
         f"(MASH={n_mash:,}; no_MASH={n_nomash:,})")
    if n_mash < MIN_CELLS_PER_GROUP or n_nomash < MIN_CELLS_PER_GROUP:
        _log(f"  [{cell_type}] SKIP — group size below {MIN_CELLS_PER_GROUP}")
        return {"cell_type": cell_type, "skipped": True, "reason": "low_n",
                "n_cells_mash": n_mash, "n_cells_nomash": n_nomash,
                "n_slide_up": 0, "n_slide_down": 0, "n_direction_concordant": 0}

    # Use the log-normalised layer (per loader: target_sum=1e4 + log1p).
    if "lognorm" not in sub.layers:
        # Fall back: log-normalise on the fly from counts at .X
        _log(f"  [{cell_type}] WARN — .layers['lognorm'] missing; re-computing on subset")
        sub.layers["counts"] = sub.X.copy()
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sub.layers["lognorm"] = sub.X.copy()
    else:
        sub.X = sub.layers["lognorm"]

    sub.obs["disease_stage"] = sub.obs["disease_stage"].astype("category")
    # Drop unused categories so rank_genes_groups doesn't choke
    sub.obs["disease_stage"] = sub.obs["disease_stage"].cat.remove_unused_categories()

    sc.tl.rank_genes_groups(
        sub,
        groupby="disease_stage",
        groups=["MASH"],
        reference="no_MASH",
        method="wilcoxon",
        pts=True,            # populate pct_nz_group, pct_nz_reference
        use_raw=False,
        rankby_abs=False,
    )

    rg = sub.uns["rank_genes_groups"]
    # ── A6 pseudoreplication fix (2026-06-20) ────────────────────────────────
    # The cell-level Wilcoxon is retained ONLY as a descriptive ranking signal
    # (`cell_score` drives the marker-recovery validation block below). It is
    # NOT written as a `pval`/`pval_adj` significance column, because disease
    # status is assigned at the SLIDE level (3 MASH vs 1 no_MASH slide). With
    # ~297K cells nested in 4 slides, cells within a slide are not independent;
    # the cell-level Wilcoxon p-values are pseudoreplicated and anti-conservative
    # (they returned 779-883 "sig" hep/KC genes, e.g. SERPINA1 padj≈0 while its
    # SLIDE-level direction is the OPPOSITE sign). With a 3-vs-1 slide design no
    # inferential test is valid at all. The honest, experimental-unit evidence is
    # the slide-level direction (`logfoldchange` = slide_logfc /
    # `slide_direction_concordant` / `n_slides_*`) computed below. The atlas-facing
    # `pval_adj` column is therefore emitted as all-NaN (see below) so the atlas
    # carries the slide-level DIRECTION but never a cell-level significance level.
    de_df = pd.DataFrame({
        "gene":             rg["names"]["MASH"],
        # Descriptive cell-level ranking — kept under explicit `cell_*` names so
        # it can never be mistaken for donor-level significance.
        "cell_logfoldchange": rg["logfoldchanges"]["MASH"],
        "cell_score":         rg["scores"]["MASH"],
    })
    if "pts" in rg:
        pts = rg["pts"]
        de_df["pct_nz_MASH"]    = pts["MASH"].reindex(de_df["gene"]).values
        de_df["pct_nz_no_MASH"] = pts["no_MASH"].reindex(de_df["gene"]).values
    de_df["n_cells_MASH"]    = n_mash
    de_df["n_cells_no_MASH"] = n_nomash
    de_df["cell_type"]       = cell_type

    # ── Slide-level pseudobulk direction (F148) ──────────────────────────────
    # Aggregate to one lognorm-mean profile per slide, then contrast MASH vs
    # no_MASH at the slide level. slide_logfc = mean(MASH slides) − mean(no_MASH
    # slides) on the log1p-CPM scale (a log-scale difference, NOT scanpy's
    # natural-log fold change — they are on different scales but agree in sign).
    # This is the actual experimental-unit evidence; with 3 MASH vs 1 no_MASH
    # slide it supports only a direction, not a p-value.
    try:
        pb = pseudobulk_by_donor(
            sub, donor_col="sample_id", layer="lognorm",
            agg="mean", obs_cols=["disease_stage"])
        mash_slides = pb.index[pb["disease_stage"].astype(str) == "MASH"]
        nomash_slides = pb.index[pb["disease_stage"].astype(str) == "no_MASH"]
        gene_cols = [c for c in pb.columns if c != "disease_stage"]
        mash_mean = pb.loc[mash_slides, gene_cols].mean(axis=0)
        nomash_mean = pb.loc[nomash_slides, gene_cols].mean(axis=0)
        slide_logfc = (mash_mean - nomash_mean)
        de_df["slide_logfc"] = slide_logfc.reindex(de_df["gene"]).values
        de_df["n_slides_MASH"]    = int(len(mash_slides))
        de_df["n_slides_no_MASH"] = int(len(nomash_slides))
        # Does the slide-level direction agree with the descriptive cell-level
        # logFC sign? (A concordance QC of the cell ranking, not a test.)
        de_df["slide_direction_concordant"] = (
            np.sign(de_df["slide_logfc"]) == np.sign(de_df["cell_logfoldchange"]))
        _log(f"  [{cell_type}] slide-level pseudobulk: "
             f"{len(mash_slides)} MASH vs {len(nomash_slides)} no_MASH slide(s); "
             f"direction-concordant genes = "
             f"{int(de_df['slide_direction_concordant'].sum())}/{de_df.shape[0]}")
    except Exception as e:  # defensive: never let pseudobulk break the DE writer
        _log(f"  [{cell_type}] WARN — slide pseudobulk failed ({str(e)[:100]}); "
             f"slide_logfc set NaN")
        de_df["slide_logfc"] = np.nan
        de_df["n_slides_MASH"]    = np.nan
        de_df["n_slides_no_MASH"] = np.nan
        de_df["slide_direction_concordant"] = np.nan

    # ── Atlas-facing columns (A6 fix) ────────────────────────────────────────
    # `06_integration.py` reads `logfoldchange` + `pval_adj` from this CSV into
    # spatial_govaere2026_cosmx_<ct>_mash_{logfc,padj}. We make BOTH honest:
    #   • `logfoldchange` := the SLIDE-level direction (experimental-unit effect),
    #     NOT scanpy's per-cell natural-log fold change.
    #   • `pval_adj` := all-NaN. A 3-vs-1 slide design supports a direction only,
    #     never a p-value, so we refuse to ship a significance number. The atlas
    #     `*_mash_padj` column therefore carries no threshold-passable value; any
    #     downstream significance gate on it is a guaranteed no-op (as intended).
    de_df["logfoldchange"] = de_df["slide_logfc"]
    de_df["pval_adj"]      = np.nan

    # No valid cell-level significance exists; report the slide-level DIRECTION
    # summary instead (n genes whose slide direction is up / down). This replaces
    # the retired pseudoreplicated `n_sig padj<0.05` count.
    n_slide_up   = int((de_df["slide_logfc"] > 0).sum())
    n_slide_down = int((de_df["slide_logfc"] < 0).sum())
    n_concordant = int(de_df["slide_direction_concordant"].fillna(False).sum())
    _log(f"  [{cell_type}] NO cell-level p-value reported (3-vs-1 slide design). "
         f"slide-direction up={n_slide_up} down={n_slide_down}; "
         f"cell/slide direction-concordant={n_concordant}/{de_df.shape[0]}")

    out_path = out_dir / f"cosmx_de_{cell_type}_MASH_vs_noMASH.csv"
    de_df.to_csv(out_path, index=False)
    _log(f"  [{cell_type}] → {out_path.name}  ({out_path.stat().st_size/1e3:.1f} KB)")

    return {
        "cell_type": cell_type,
        "skipped": False,
        "n_cells_mash": n_mash,
        "n_cells_nomash": n_nomash,
        "n_genes_tested": int(de_df.shape[0]),
        # A6 fix: no valid cell-level significance under a 3-vs-1 slide design.
        # Report slide-level DIRECTION counts + cell/slide direction concordance.
        "n_slide_up": n_slide_up,
        "n_slide_down": n_slide_down,
        "n_direction_concordant": n_concordant,
        "out_path": str(out_path),
        "de_df": de_df,
    }


def composition_table(adata: ad.AnnData, out_dir: pathlib.Path) -> pd.DataFrame:
    """Per-sample × per-cell-type counts + proportions."""
    comp = (adata.obs
            .groupby(["sample_id", "cell_type"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index())
    sample_totals = comp.groupby("sample_id", observed=True)["n_cells"].transform("sum")
    comp["proportion"] = comp["n_cells"] / sample_totals
    stage_map = {sid: meta["stage"] for sid, meta in SAMPLE_DISEASE.items()}
    comp["disease_stage"] = comp["sample_id"].map(stage_map)
    comp["sample_disease_str"] = comp["sample_id"].map(
        {sid: meta["sample_disease_str"] for sid, meta in SAMPLE_DISEASE.items()})
    out_path = out_dir / "cosmx_celltype_composition.tsv"
    comp.to_csv(out_path, sep="\t", index=False)
    _log(f"Wrote composition table → {out_path.name} ({comp.shape[0]} rows)")
    return comp


def report_validation(de_results: List[Dict], comp: pd.DataFrame) -> None:
    """Echo the validation checks called out in the task spec."""
    _log("=" * 70)
    _log("VALIDATION")
    _log("=" * 70)

    # 1. Composition shifts
    _log("Per-cell-type proportion by slide (MASH composition descending):")
    pivot = comp.pivot_table(index="cell_type", columns="sample_id",
                             values="proportion", fill_value=0.0)
    for ct in pivot.index:
        row = " ".join(f"{sid}={pivot.loc[ct, sid]:.4f}" for sid in pivot.columns)
        _log(f"  {ct:>14}   {row}")

    # MetMac/KC trend across MASH-only slides vs the mixed-normal slide.
    if "KC" in pivot.index:
        kc_mash_only = pivot.loc["KC", ["Leuven_1", "Leuven_3", "Leuven_4"]].mean()
        kc_mixed     = pivot.loc["KC", "Leuven_2"] if "Leuven_2" in pivot.columns else float("nan")
        _log(f"  KC mean proportion (MASH-only slides) = {kc_mash_only:.4f}")
        _log(f"  KC proportion (Leuven_2 mixed-with-normal) = {kc_mixed:.4f}")
        _log("  (Paper: MetMac 5.7%→12% in MASH; KC decrease. KC here lumps all macrophage lineage.)")

    # 2. Top macrophage-cluster MASH-up genes — KC cluster is the only macrophage
    #    lineage cluster present.
    mac_result = next((r for r in de_results if r["cell_type"] == "KC" and not r.get("skipped")), None)
    if mac_result is None:
        _log("KC DE not available — cannot run MetMac marker recovery check.")
        return
    de_df = mac_result["de_df"].copy()
    # Marker recovery ranks by the DESCRIPTIVE cell-level score/logFC (NOT a
    # significance call); `pval_adj` is intentionally NaN (3-vs-1 slide design),
    # so we report the slide-level direction-concordance flag instead.
    de_df_up = de_df[de_df["cell_logfoldchange"] > 0].sort_values("cell_score", ascending=False)
    top10 = de_df_up.head(10)
    _log("Top 10 MASH-up genes in KC cluster (macrophage lineage proxy; cell-level ranking, descriptive only):")
    for _, row in top10.iterrows():
        _log(f"  {row['gene']:>10}  cell_lfc={row['cell_logfoldchange']:+.3f}  "
             f"slide_lfc={row.get('slide_logfc', float('nan')):+.3f}  "
             f"slide_concordant={row.get('slide_direction_concordant', float('nan'))}  "
             f"pct_MASH={row.get('pct_nz_MASH', float('nan')):.3f}")
    # MetMac marker recovery
    de_df_up_50 = de_df_up.head(50)["gene"].tolist()
    found = [g for g in METMAC_MARKERS if g in de_df_up_50]
    missing = [g for g in METMAC_MARKERS if g not in de_df_up_50]
    _log(f"MetMac markers in top 50 KC MASH-up: {found}")
    _log(f"MetMac markers absent from top 50:    {missing}")


def main():
    parser = argparse.ArgumentParser(description="Per-cell-type DE for Govaere 2026 CosMx")
    parser.add_argument("--adata", type=str, default=str(ADATA_PATH))
    parser.add_argument("--out-dir", type=str, default=str(OUT_DIR))
    args = parser.parse_args()

    adata_path = pathlib.Path(args.adata)
    out_dir    = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _log("=" * 70)
    _log("42b_govaere2026_cosmx_de.py")
    _log(f"  AnnData : {adata_path}")
    _log(f"  Out dir : {out_dir}")
    _log("=" * 70)

    _log(f"Reading AnnData ({adata_path.stat().st_size/1e6:.1f} MB)")
    adata = sc.read_h5ad(adata_path)
    _log(f"  shape: {adata.shape}")
    _log(f"  .obs columns: {list(adata.obs.columns)}")
    _log(f"  .layers: {list(adata.layers.keys())}")
    _log(f"  sample_id values: {sorted(adata.obs['sample_id'].unique().tolist())}")
    _log(f"  cell_type values: {sorted(adata.obs['cell_type'].unique().tolist())}")

    adata = assign_disease_stage(adata)

    # Composition table — independent of DE
    comp = composition_table(adata, out_dir)

    # Per-cell-type DE
    de_results: List[Dict] = []
    for ct in TARGET_CELLTYPES:
        res = per_celltype_de(adata, ct, out_dir)
        de_results.append(res)

    # Summary
    summary_rows = []
    for r in de_results:
        summary_rows.append({
            "cell_type":         r["cell_type"],
            "contrast":          "MASH_vs_no_MASH",
            "n_cells_MASH":      r.get("n_cells_mash", 0),
            "n_cells_no_MASH":   r.get("n_cells_nomash", 0),
            "n_genes_tested":    r.get("n_genes_tested", 0),
            # A6 fix: slide-level DIRECTION summary replaces the retired
            # pseudoreplicated cell-level n_sig_padj05 counts (3-vs-1 slide design
            # = direction only, no valid p-value).
            "n_slide_up":            r.get("n_slide_up", 0),
            "n_slide_down":          r.get("n_slide_down", 0),
            "n_direction_concordant": r.get("n_direction_concordant", 0),
            "skipped":           r.get("skipped", False),
            "skip_reason":       r.get("reason", ""),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_dir / "cosmx_de_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    _log(f"Wrote DE summary → {summary_path.name}")

    # Validation
    report_validation(de_results, comp)
    _log("DONE")


if __name__ == "__main__":
    main()
