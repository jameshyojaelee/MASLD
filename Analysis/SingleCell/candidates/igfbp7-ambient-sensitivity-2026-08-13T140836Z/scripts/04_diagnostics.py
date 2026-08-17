#!/usr/bin/env python
"""Ambient-contamination diagnostics for the IGFBP7 SENSITIVITY ANALYSIS.

SCOPE (binding): sensitivity analysis only. Frozen program memberships and L1
weights are read-only; nothing is refit, reweighted or re-selected.

Question: within HEPATOCYTE barcodes, how much of each frozen program's
member-gene signal is ambient-attributable, and does the Stromal ECM (IGFBP7)
program behave like fibroblast/stellate transcripts that hepatocytes should not
express, or like bona fide hepatocyte transcripts?
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
WORK = CAND / "work"
RES = CAND / "results"

PROGRAM_ROOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/hotspot"
)
HEP_LABEL = "Hepatocytes"
FIB_LABEL = "Fibroblasts"


def main() -> None:
    # decontX runs as one shard per sequencing dataset; concatenate the shards.
    shards = sorted(RES.glob("02_ambient_by_gene_celltype*.tsv"))
    if not shards:
        raise RuntimeError("no decontX ambient tables found")
    gene_tab = pd.concat([pd.read_csv(p, sep="\t") for p in shards], ignore_index=True)
    seen = sorted(gene_tab["dataset"].unique())
    print(f"[diagnostics] {len(shards)} decontX shards covering datasets: {seen}", flush=True)
    # R writes logicals as the strings TRUE/FALSE; astype(bool) would make both True.
    is_corrected = gene_tab["corrected"].map(
        lambda v: str(v).strip().upper() in {"TRUE", "T", "1"}
    )
    uncorrected = sorted(gene_tab.loc[~is_corrected, "dataset"].unique())
    if uncorrected:
        print(f"[diagnostics] UNCORRECTED datasets carried through raw: {uncorrected}", flush=True)
    manifest = json.loads((WORK / "manifest.json").read_text(encoding="utf-8"))
    membership = pd.read_csv(PROGRAM_ROOT / "program_membership_v2.tsv", sep="\t")
    registry = pd.read_csv(PROGRAM_ROOT / "program_registry_v2.tsv", sep="\t")

    # Datasets that could not be corrected pass their raw counts through, which
    # forces their ambient fraction to exactly 0. Including them would bias the
    # pooled diagnostic toward zero, so they are excluded here and reported.
    gene_tab = gene_tab[is_corrected].copy()

    # Pool over datasets by SUMMING counts, so the ambient fraction is a count
    # ratio and not an average of ratios.
    pooled = (
        gene_tab.groupby(["cell_type", "gene"], as_index=False)[
            ["n_cells", "sum_raw", "sum_decont", "n_pos_raw", "n_pos_decont", "n_pos_decont_strict"]
        ].sum()
    )
    pooled["ambient_fraction"] = np.where(
        pooled["sum_raw"] > 0, 1 - pooled["sum_decont"] / pooled["sum_raw"], np.nan
    )
    pooled["pct_pos_raw"] = 100 * pooled["n_pos_raw"] / pooled["n_cells"]
    pooled["pct_pos_decont"] = 100 * pooled["n_pos_decont"] / pooled["n_cells"]
    pooled.to_csv(RES / "04_ambient_pooled_by_gene_celltype.tsv", sep="\t", index=False)

    hep = pooled[pooled["cell_type"] == HEP_LABEL].set_index("gene")
    fib = pooled[pooled["cell_type"] == FIB_LABEL].set_index("gene")

    # ---- per-gene panel for every frozen hepatocyte program -------------
    mem = membership[membership["cell_type"] == "hepatocytes"].copy()
    mem["ambient_fraction_hep"] = mem["source_gene"].map(hep["ambient_fraction"])
    mem["ambient_fraction_fib"] = mem["source_gene"].map(fib["ambient_fraction"])
    mem["pct_pos_raw_hep"] = mem["source_gene"].map(hep["pct_pos_raw"])
    mem["pct_pos_decont_hep"] = mem["source_gene"].map(hep["pct_pos_decont"])
    mem["sum_raw_hep"] = mem["source_gene"].map(hep["sum_raw"])
    mem["sum_raw_fib"] = mem["source_gene"].map(fib["sum_raw"])
    mem = mem.merge(
        registry[["cell_type", "module", "module_name", "robust_display", "primary_selected"]],
        on=["cell_type", "module"],
        how="left",
    )
    keep = [
        "cell_type", "module", "module_name", "source_gene", "original_l1_weight",
        "ambient_fraction_hep", "ambient_fraction_fib", "pct_pos_raw_hep",
        "pct_pos_decont_hep", "sum_raw_hep", "sum_raw_fib",
        "robust_display", "primary_selected",
    ]
    mem[keep].sort_values(["module", "original_l1_weight"], ascending=[True, False]).to_csv(
        RES / "04_program_member_ambient.tsv", sep="\t", index=False
    )

    # ---- L1-weighted program-level ambient burden -----------------------
    prog_rows = []
    for module, grp in mem.groupby("module"):
        w = grp["original_l1_weight"].to_numpy(float)
        af = grp["ambient_fraction_hep"].to_numpy(float)
        ok = np.isfinite(af)
        prog_rows.append(
            {
                "cell_type": "hepatocytes",
                "module": int(module),
                "module_name": grp["module_name"].iloc[0],
                "robust_display": grp["robust_display"].iloc[0],
                "primary_selected": grp["primary_selected"].iloc[0],
                "n_members": len(grp),
                "n_members_with_ambient_estimate": int(ok.sum()),
                "l1_weight_covered": float(w[ok].sum()),
                "l1_weighted_ambient_fraction_hep": float((w[ok] * af[ok]).sum() / w[ok].sum()),
                "median_ambient_fraction_hep": float(np.nanmedian(af)),
                "frac_members_ambient_gt_0.5": float(np.nanmean(af > 0.5)),
                "frac_members_ambient_gt_0.8": float(np.nanmean(af > 0.8)),
                # aggregate count-ratio version: total ambient counts / total counts
                "program_count_ambient_fraction_hep": float(
                    1
                    - grp["source_gene"].map(hep["sum_decont"]).sum()
                    / grp["source_gene"].map(hep["sum_raw"]).sum()
                ),
                "fibroblast_to_hepatocyte_count_ratio": float(
                    grp["sum_raw_fib"].sum() / max(grp["sum_raw_hep"].sum(), 1.0)
                ),
            }
        )
    prog = pd.DataFrame(prog_rows).sort_values("l1_weighted_ambient_fraction_hep", ascending=False)
    prog.to_csv(RES / "04_program_ambient_burden.tsv", sep="\t", index=False)

    # ---- control panels -------------------------------------------------
    ctrl_rows = []
    panels = {
        "fibroblast_positive_control": manifest["fibro_pos_controls"],
        "hepatocyte_negative_control": manifest["hep_neg_controls"],
    }
    for panel, genes in panels.items():
        for g in genes:
            if g not in hep.index:
                ctrl_rows.append({"panel": panel, "gene": g, "note": "absent from atlas"})
                continue
            ctrl_rows.append(
                {
                    "panel": panel,
                    "gene": g,
                    "ambient_fraction_hep": float(hep.loc[g, "ambient_fraction"]),
                    "ambient_fraction_fib": float(fib.loc[g, "ambient_fraction"]) if g in fib.index else np.nan,
                    "pct_pos_raw_hep": float(hep.loc[g, "pct_pos_raw"]),
                    "pct_pos_decont_hep": float(hep.loc[g, "pct_pos_decont"]),
                    "sum_raw_hep": float(hep.loc[g, "sum_raw"]),
                    "sum_raw_fib": float(fib.loc[g, "sum_raw"]) if g in fib.index else np.nan,
                    "note": "",
                }
            )
    ctrl = pd.DataFrame(ctrl_rows)
    ctrl.to_csv(RES / "04_control_panels.tsv", sep="\t", index=False)

    # per-dataset replication of the headline contrast, so the verdict does not
    # rest on any single cohort
    rep_rows = []
    panels_named = {
        "IGFBP7 program (hep module 8)": mem.loc[mem["module"] == 8, "source_gene"].tolist(),
        "BICC1 program (hep module 20)": mem.loc[mem["module"] == 20, "source_gene"].tolist(),
        "STAB2 endothelial program (hep module 2)": mem.loc[mem["module"] == 2, "source_gene"].tolist(),
        "DOCK2 leukocyte program (hep module 1)": mem.loc[mem["module"] == 1, "source_gene"].tolist(),
        "ALB program (hep module 12)": mem.loc[mem["module"] == 12, "source_gene"].tolist(),
        "CYP program (hep module 21)": mem.loc[mem["module"] == 21, "source_gene"].tolist(),
        "fibroblast positive controls": manifest["fibro_pos_controls"],
        "hepatocyte negative controls": manifest["hep_neg_controls"],
    }
    for dataset, grp in gene_tab[gene_tab["cell_type"] == HEP_LABEL].groupby("dataset"):
        gi = grp.set_index("gene")
        af = np.where(gi["sum_raw"] > 0, 1 - gi["sum_decont"] / gi["sum_raw"], np.nan)
        af = pd.Series(af, index=gi.index)
        for name, genes in panels_named.items():
            v = af.reindex(genes).dropna()
            rep_rows.append(
                {
                    "dataset": dataset,
                    "n_hepatocytes": int(gi["n_cells"].iloc[0]),
                    "set": name,
                    "n_genes": int(len(v)),
                    "median_ambient_fraction_hep": float(v.median()) if len(v) else np.nan,
                }
            )
    pd.DataFrame(rep_rows).to_csv(RES / "04_per_dataset_replication.tsv", sep="\t", index=False)

    # ---- headline comparison -------------------------------------------
    def panel_summary(name: str, af: np.ndarray) -> dict:
        af = af[np.isfinite(af)]
        return {
            "set": name,
            "n_genes": int(len(af)),
            "median_ambient_fraction_hep": float(np.median(af)) if len(af) else np.nan,
            "mean_ambient_fraction_hep": float(np.mean(af)) if len(af) else np.nan,
            "q25": float(np.percentile(af, 25)) if len(af) else np.nan,
            "q75": float(np.percentile(af, 75)) if len(af) else np.nan,
            "frac_gt_0.5": float(np.mean(af > 0.5)) if len(af) else np.nan,
            "frac_gt_0.8": float(np.mean(af > 0.8)) if len(af) else np.nan,
        }

    comp = [
        panel_summary(
            "Stromal ECM (IGFBP7) members [hep module 8]",
            mem.loc[mem["module"] == 8, "ambient_fraction_hep"].to_numpy(float),
        ),
        panel_summary(
            "Ductular injury (BICC1) members [hep module 20]",
            mem.loc[mem["module"] == 20, "ambient_fraction_hep"].to_numpy(float),
        ),
        panel_summary(
            "Secretory plasma protein (ALB) members [hep module 12]",
            mem.loc[mem["module"] == 12, "ambient_fraction_hep"].to_numpy(float),
        ),
        panel_summary(
            "Xenobiotic / drug metab (CYP) members [hep module 21]",
            mem.loc[mem["module"] == 21, "ambient_fraction_hep"].to_numpy(float),
        ),
        panel_summary(
            "Sinusoidal endothelial (STAB2) members [hep module 2]",
            mem.loc[mem["module"] == 2, "ambient_fraction_hep"].to_numpy(float),
        ),
        panel_summary(
            "fibroblast positive controls",
            ctrl.loc[ctrl["panel"] == "fibroblast_positive_control", "ambient_fraction_hep"].to_numpy(float)
            if "ambient_fraction_hep" in ctrl
            else np.array([]),
        ),
        panel_summary(
            "hepatocyte negative controls",
            ctrl.loc[ctrl["panel"] == "hepatocyte_negative_control", "ambient_fraction_hep"].to_numpy(float)
            if "ambient_fraction_hep" in ctrl
            else np.array([]),
        ),
    ]
    pd.DataFrame(comp).to_csv(RES / "04_headline_comparison.tsv", sep="\t", index=False)

    print(pd.DataFrame(comp).to_string(index=False), flush=True)
    print("\n[diagnostics] program ambient burden (top 10 and bottom 5):", flush=True)
    cols = [
        "module", "module_name", "robust_display",
        "l1_weighted_ambient_fraction_hep", "program_count_ambient_fraction_hep",
        "fibroblast_to_hepatocyte_count_ratio",
    ]
    print(prog[cols].head(10).to_string(index=False), flush=True)
    print(prog[cols].tail(5).to_string(index=False), flush=True)
    print("[diagnostics] DONE", flush=True)


if __name__ == "__main__":
    main()
