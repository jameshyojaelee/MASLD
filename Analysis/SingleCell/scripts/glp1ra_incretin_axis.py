#!/usr/bin/env python
"""
GLP-1RA incretin-axis per-cell-type expression map (Analysis A, Fig 5 panel).

Reads per-cell-type pseudobulk (raw count SUMS, genes x donors) from the scRNA
atlas, CPM-normalizes per donor, and reports per-cell-type mean/median CPM and
donor-detection rate for the incretin/target axis:
    GLP1R, GIPR, GCGR, GLP2R, GCG, DPP4
plus positive-control lineage markers so we can show the pipeline detects genes
that ARE expressed.

Headline claim to verify: GLP1R ~ 0 in Hepatocytes (and stellate/Fibroblasts),
while DPP4 is high in Hepatocytes -> the approved drug's receptor is NOT on the
parenchymal cells that execute MASH; its liver benefit is indirect.

Output: RNA-seq/results/glp1ra/scrna_incretin_axis/
  - incretin_axis_celltype_cpm_long.csv     (tidy long)
  - incretin_axis_meancpm_matrix.csv        (genes x cell types, mean CPM)
  - incretin_axis_detection_matrix.csv      (genes x cell types, frac donors detected)
"""
import os
import glob
import sys
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PB_DIR = os.path.join(ROOT, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
OUT = os.path.join(ROOT, "RNA-seq/results/glp1ra/scrna_incretin_axis")
os.makedirs(OUT, exist_ok=True)

AXIS_GENES = ["GLP1R", "GIPR", "GCGR", "GLP2R", "GCG", "DPP4"]
# positive-control lineage markers (should be high in the matching cell type)
CTRL_GENES = ["ALB", "APOB", "PECAM1", "PTPRC", "COL1A1", "DCN", "EPCAM", "CD68"]
GENES = AXIS_GENES + CTRL_GENES


def pick_files():
    """One file per cell type; prefer underscore variant (canonical refresh) over
    the older space-named duplicate."""
    raw = glob.glob(os.path.join(PB_DIR, "*_pseudobulk.csv"))
    by_norm = {}
    for f in raw:
        ct = os.path.basename(f).replace("_pseudobulk.csv", "")
        norm = ct.replace(" ", "_")
        if norm not in by_norm:
            by_norm[norm] = f
        else:
            prev_has_space = " " in os.path.basename(by_norm[norm])
            if prev_has_space and " " not in os.path.basename(f):
                by_norm[norm] = f
    return dict(sorted(by_norm.items()))


def main():
    files = pick_files()
    print(f"[info] {len(files)} cell-type pseudobulk files:")
    for ct in files:
        print("   ", ct)

    rows = []
    for ct, f in files.items():
        df = pd.read_csv(f, index_col=0)
        # numeric only
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        libsize = df.sum(axis=0)              # per-donor library size
        libsize = libsize.replace(0, np.nan)
        n_donors = df.shape[1]
        for g in GENES:
            if g not in df.index:
                rows.append(dict(cell_type=ct, gene=g, present=False,
                                 mean_cpm=np.nan, median_cpm=np.nan,
                                 frac_donors_detected=np.nan, n_donors=n_donors))
                continue
            counts = df.loc[g]
            if isinstance(counts, pd.DataFrame):   # duplicate gene rows -> sum
                counts = counts.sum(axis=0)
            cpm = counts / libsize * 1e6
            rows.append(dict(
                cell_type=ct, gene=g, present=True,
                mean_cpm=float(np.nanmean(cpm)),
                median_cpm=float(np.nanmedian(cpm)),
                frac_donors_detected=float((counts > 0).mean()),
                n_donors=n_donors,
            ))

    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(OUT, "incretin_axis_celltype_cpm_long.csv"), index=False)

    mean_mat = long.pivot(index="gene", columns="cell_type", values="mean_cpm").reindex(GENES)
    det_mat = long.pivot(index="gene", columns="cell_type", values="frac_donors_detected").reindex(GENES)
    mean_mat.to_csv(os.path.join(OUT, "incretin_axis_meancpm_matrix.csv"))
    det_mat.to_csv(os.path.join(OUT, "incretin_axis_detection_matrix.csv"))

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.float_format", lambda x: f"{x:8.2f}")

    print("\n===== MEAN CPM (incretin axis + controls) x cell type =====")
    print(mean_mat.round(2).to_string())

    # Headline check
    print("\n===== HEADLINE CHECK (Hepatocytes) =====")
    hep_cols = [c for c in mean_mat.columns if c.lower().startswith("hepatocyte")]
    if hep_cols:
        hc = hep_cols[0]
        for g in AXIS_GENES + ["ALB", "DPP4"]:
            if g in mean_mat.index:
                print(f"  {g:6s}  mean_cpm={mean_mat.loc[g, hc]:10.3f}   "
                      f"detect={det_mat.loc[g, hc]:.2f}")
        glp1r = mean_mat.loc["GLP1R", hc] if "GLP1R" in mean_mat.index else np.nan
        dpp4 = mean_mat.loc["DPP4", hc] if "DPP4" in mean_mat.index else np.nan
        print(f"\n  => GLP1R hepatocyte mean CPM = {glp1r:.3f} ; DPP4 = {dpp4:.3f} ; "
              f"ratio DPP4/GLP1R = {dpp4/glp1r if glp1r else float('inf'):.1f}")

    # Where IS GLP1R highest?
    if "GLP1R" in mean_mat.index:
        print("\n===== GLP1R mean CPM ranked across cell types =====")
        print(mean_mat.loc["GLP1R"].sort_values(ascending=False).round(3).to_string())
    print("\n[done] outputs in", OUT)


if __name__ == "__main__":
    main()
