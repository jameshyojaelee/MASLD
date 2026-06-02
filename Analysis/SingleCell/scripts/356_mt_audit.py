#!/usr/bin/env python
"""
356_mt_audit.py  (S5)

Per-cohort × per-cell-type mt% audit + sensitivity reruns.

Audit:
    - load atlas; compute per-cohort, per-cell-type pct_counts_mt distribution
    - emit histograms (matplotlib PDF)
    - sensitivity: for hepatocytes, recompute the Progressor / Moderate /
      Stable proportion (5 meta-subtype label) at mt% cutoffs {15, 20, 25}
    - compare against the production mt% cutoff currently in use
      (default 20 if obs column 'pct_counts_mt' present and any cell retained)

Output: Analysis/SingleCell/results_gpu_v2/mt_audit/REPORT.md
        + mt_audit_summary.tsv + per-cohort histograms PDF
"""

from __future__ import annotations
import os
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
H5AD = PROJECT_ROOT / "Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad"
HEP_LABELS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/meta_subtype_labels.tsv"
OUT_ROOT = Path(os.environ.get("S5_OUT_ROOT", str(PROJECT_ROOT)))
OUT = OUT_ROOT / "Analysis/SingleCell/results_gpu_v2/mt_audit"
OUT.mkdir(parents=True, exist_ok=True)

MT_CUTOFFS = [15, 20, 25]


def main():
    print(f"[load] {H5AD}", flush=True)
    A = ad.read_h5ad(H5AD, backed="r")
    obs = A.obs.copy()
    print(f"[load] {A.shape}; cols={obs.columns.tolist()[:15]}", flush=True)

    mt_col = None
    for cand in ("pct_counts_mt", "percent_mt", "pct_mt", "mt_pct"):
        if cand in obs.columns:
            mt_col = cand
            break
    if mt_col is None:
        print("[fatal] no mt% column on adata.obs", flush=True)
        (OUT / "REPORT.md").write_text(
            "# mt audit (S5)\n\n"
            "No mitochondrial percentage column found on integrated atlas.\n"
            "Expected one of: pct_counts_mt / percent_mt / pct_mt / mt_pct.\n")
        return

    ds_col = "dataset" if "dataset" in obs.columns else "sample"
    ct_col = "cell_type"

    summary = (
        obs.groupby([ds_col, ct_col])[mt_col]
           .agg(["count", "mean", "median",
                 lambda x: float((x < 15).mean()),
                 lambda x: float((x < 20).mean()),
                 lambda x: float((x < 25).mean())])
           .rename(columns={"<lambda_0>": "frac_lt_15",
                            "<lambda_1>": "frac_lt_20",
                            "<lambda_2>": "frac_lt_25"})
           .reset_index()
    )
    summary.to_csv(OUT / "mt_audit_summary.tsv", sep="\t", index=False)
    print(summary.head().to_string(index=False), flush=True)

    # Histograms: per-cohort, faceted by major cell types
    major_cts = ["Hepatocytes", "Macrophages", "Fibroblasts",
                 "Endothelial cells", "Cholangiocytes"]
    datasets = sorted(obs[ds_col].unique())
    n_ds = len(datasets)
    fig, axes = plt.subplots(n_ds, 1, figsize=(8, 2.0 * max(n_ds, 1)),
                             squeeze=False)
    for i, ds in enumerate(datasets):
        ax = axes[i, 0]
        for ct in major_cts:
            d = obs[(obs[ds_col] == ds) & (obs[ct_col] == ct)][mt_col]
            if len(d) == 0:
                continue
            ax.hist(d.clip(0, 60), bins=40, alpha=0.4, label=ct)
        ax.set_title(f"{ds}  ({mt_col})")
        ax.set_xlabel("mt %"); ax.set_ylabel("cells")
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT / "mt_histograms_by_cohort.pdf")
    plt.close(fig)
    print(f"[write] {OUT / 'mt_histograms_by_cohort.pdf'}", flush=True)

    # Sensitivity: Progressor proportion in hepatocytes at mt < {15,20,25}
    sens_rows = []
    if HEP_LABELS.exists():
        lab = pd.read_csv(HEP_LABELS, sep="\t")
        # Assume column 'meta_subtype' and join on cell barcode index
        join_col = "cell" if "cell" in lab.columns else lab.columns[0]
        lab = lab.set_index(join_col)
        common = obs.index.intersection(lab.index)
        hep_obs = obs.loc[common].copy()
        hep_obs["meta_subtype"] = lab.loc[common, "meta_subtype"]
        hep_obs = hep_obs[hep_obs[ct_col] == "Hepatocytes"]
        for cut in MT_CUTOFFS:
            d = hep_obs[hep_obs[mt_col] < cut]
            vc = d["meta_subtype"].value_counts(normalize=True)
            row = {"mt_cutoff": cut, "n_cells": int(len(d))}
            for k, v in vc.items():
                row[k] = float(v)
            sens_rows.append(row)
    else:
        print(f"[skip sensitivity] no {HEP_LABELS}", flush=True)
    sens = pd.DataFrame(sens_rows)
    sens.to_csv(OUT / "progressor_proportion_sensitivity.tsv",
                sep="\t", index=False)

    # REPORT.md
    lines = []
    lines.append("# mt audit (S5)\n")
    lines.append(f"Source h5ad: `{H5AD}`\n")
    lines.append(f"mt column: `{mt_col}`\n\n")
    lines.append("## Per-cohort × cell-type mt summary (top 30 rows)\n\n")
    lines.append(summary.head(30).to_markdown(index=False))
    lines.append("\n\n## Hepatocyte Progressor sensitivity to mt cutoff\n\n")
    if not sens.empty:
        lines.append(sens.to_markdown(index=False))
        lines.append("\n\nInterpretation: stability of meta_subtype proportions "
                     "across mt% ∈ {15, 20, 25} indicates the hepatocyte "
                     "subclustering is not driven by mitochondrial-content "
                     "filtering threshold.\n")
    else:
        lines.append("Hepatocyte meta_subtype labels not on disk; sensitivity "
                     "skipped.\n")
    (OUT / "REPORT.md").write_text("\n".join(lines))
    print(f"[write] {OUT / 'REPORT.md'}", flush=True)


if __name__ == "__main__":
    main()
