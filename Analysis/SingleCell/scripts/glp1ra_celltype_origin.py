#!/usr/bin/env python
"""
C1 — Genome-wide cell-type-of-origin axis (for the GLP-1RA "reach" test).

Per gene, the BASELINE (not DE) cell-type of expression: which liver cell type
expresses it, and how hepatocyte-specific is it. This is the axis needed to test
whether semaglutide's reversal is predicted by where a gene is expressed — the
refractory-metabolic genes (BHMT/FADS2/ACSL5...) are NOT DE-resolvable in snRNA,
so we use baseline expression, not DE attribution.

Method: for each per-cell-type pseudobulk expression matrix (genes x donors, raw
count sums), CPM-normalize each donor, average across donors -> per-gene mean CPM
in that cell type. Assemble a genes x cell-types mean-CPM matrix, then per gene:
  hep_fraction  = Hepatocyte CPM / sum(all cell-type CPM)
  nonparen_frac = 1 - hep_fraction
  tau           = sum(1 - x_i/x_max)/(n-1)   (tissue-specificity index, 0=ubiquitous,1=specific)
  dominant_celltype, is_hep_dominant

Output: RNA-seq/results/glp1ra/celltype_origin_axis.csv
Env: spatial (pandas/numpy).
"""
import os
import glob
import numpy as np
import pandas as pd

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PB_DIR = os.path.join(ROOT, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
OUT = os.path.join(ROOT, "RNA-seq/results/glp1ra/celltype_origin_axis.csv")

HEPATOCYTE = "Hepatocytes"
MIN_TOTAL_CPM = 1.0   # genes below this summed CPM across cell types -> low-confidence origin


def pick_files():
    """One file per cell type; prefer underscore variant over the space-named dup."""
    raw = glob.glob(os.path.join(PB_DIR, "*_pseudobulk.csv"))
    by_norm = {}
    for f in raw:
        ct = os.path.basename(f).replace("_pseudobulk.csv", "")
        norm = ct.replace(" ", "_")
        if norm not in by_norm or (" " in os.path.basename(by_norm[norm]) and " " not in os.path.basename(f)):
            by_norm[norm] = f
    return dict(sorted(by_norm.items()))


def main():
    files = pick_files()
    print(f"[info] {len(files)} cell types: {list(files)}")

    # per-cell-type mean CPM vector (Series indexed by gene)
    means = {}
    for ct, f in files.items():
        df = pd.read_csv(f, index_col=0)
        df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        # collapse duplicate gene rows if any
        if df.index.duplicated().any():
            df = df.groupby(level=0).sum()
        libsize = df.sum(axis=0).replace(0, np.nan)      # per-donor library size
        cpm = df.divide(libsize, axis=1) * 1e6            # genes x donors CPM
        means[ct] = cpm.mean(axis=1)                       # mean CPM across donors
        print(f"   {ct:28s} genes={df.shape[0]} donors={df.shape[1]}")

    mat = pd.DataFrame(means).fillna(0.0)                  # genes x cell types
    cts = list(mat.columns)
    if HEPATOCYTE not in cts:
        raise SystemExit(f"[FATAL] '{HEPATOCYTE}' column not found in {cts}")

    row_sum = mat.sum(axis=1)
    row_max = mat.max(axis=1)
    n_ct = mat.shape[1]

    out = pd.DataFrame(index=mat.index)
    out.index.name = "gene"
    out["total_cpm"] = row_sum
    out["hep_cpm"] = mat[HEPATOCYTE]
    out["hep_fraction"] = np.where(row_sum > 0, mat[HEPATOCYTE] / row_sum, np.nan)
    out["nonparen_fraction"] = 1 - out["hep_fraction"]
    # tau specificity index
    tau = ((1 - mat.divide(row_max.replace(0, np.nan), axis=0)).sum(axis=1)) / (n_ct - 1)
    out["tau"] = tau
    out["dominant_celltype"] = mat.idxmax(axis=1)
    out["is_hep_dominant"] = out["dominant_celltype"] == HEPATOCYTE
    out["origin_confident"] = row_sum >= MIN_TOTAL_CPM
    # attach the per-cell-type mean CPM for transparency
    for ct in cts:
        out[f"cpm__{ct.replace(' ', '_')}"] = mat[ct]

    out.reset_index().to_csv(OUT, index=False)
    print(f"\n[done] wrote {OUT}  ({out.shape[0]} genes)")

    # ---- sanity checks ----
    def show(g):
        if g in out.index:
            r = out.loc[g]
            print(f"  {g:8s} hep_frac={r.hep_fraction:5.2f}  tau={r.tau:4.2f}  "
                  f"dominant={r.dominant_celltype:20s} total_cpm={r.total_cpm:9.1f}")
        else:
            print(f"  {g:8s} (absent)")
    print("\n=== sanity (expect: ALB/APOB hep-dominant ~1; COL1A1/DCN fibroblast; "
          "GCGR hep-high; PTPRC immune) ===")
    for g in ["ALB", "APOB", "GCGR", "GLP1R", "DPP4", "COL1A1", "DCN", "PTPRC", "PECAM1",
              "BHMT", "FADS2", "ACSL5", "TIMP1"]:
        show(g)
    print("\nhep_fraction distribution (confident genes):")
    print(out.loc[out.origin_confident, "hep_fraction"].describe().round(3).to_string())


if __name__ == "__main__":
    main()
