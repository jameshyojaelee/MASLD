#!/usr/bin/env python3
"""Single-cell layer web data (data contract §5.1).

Sources (under Analysis/SingleCell/results_gpu_v2/):
  cell_type_proportions.csv                      -> composition (275 samples)
  disease_signatures/celltype_disease_activity.csv -> per-cell-type disease activity
  disease_signatures/fgsea_celltype_F_transitions.csv -> fgsea_transitions
  pseudobulk_de/allcell_pseudobulk_de.csv        -> sc_pseudobulk_de.parquet (ALL-CELL)
  hepatocyte_subtypes/{subtype_axis_scores,subtype_disease_enrichment}.csv -> hep_subtypes
  hepatocyte_subtypes/subtype_markers.csv        -> sc_hep_markers.parquet (dotplot)

Outputs (-> --output-dir):
  singlecell_summary.json
  sc_pseudobulk_de.parquet     (all-cell, long, sorted by symbol)
  sc_hep_markers.parquet       (top-N markers per hep subtype, sorted by symbol)

Run:
  micromamba run -n spatial python scripts/portal/generate_singlecell_data.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _portal_io import (dump_json, load_ensembl_symbol_map, project_root, r,
                        strip_version, write_parquet)

# Given by the atlas headline (composition CSV holds fractions, not counts).
N_CELLS = 1_232_318
N_DATASETS = 7
N_SAMPLES = 275

CELL_TYPES = [
    "B cells", "Basophils", "Cholangiocytes", "Circulating NK/NKT",
    "Endothelial cells", "Fibroblasts", "Hepatocytes", "Macrophages",
    "Mono+mono derived cells", "Neutrophils", "Plasma cells", "Resident NK",
    "T cells", "cDC1s", "cDC2s", "pDCs",
]

# celltype_disease_activity.csv: rows = cell types, cols = ct_<CellType>_MASLD
# signatures. Per-cell-type "activity" = the diagonal (a cell type scored on its
# own MASLD disease signature). The 'monorived' typo in the Mono column is why
# this map is explicit rather than string-normalized. B cells / Macrophages cols
# are all-NaN in the source, so they drop out (finite-diagonal only).
ACT_COL_TO_CT = {
    "ct_B_cells_MASLD": "B cells",
    "ct_Cholangiocytes_MASLD": "Cholangiocytes",
    "ct_Circulating_NK_NKT_MASLD": "Circulating NK/NKT",
    "ct_Endothelial_cells_MASLD": "Endothelial cells",
    "ct_Fibroblasts_MASLD": "Fibroblasts",
    "ct_Hepatocytes_MASLD": "Hepatocytes",
    "ct_Macrophages_MASLD": "Macrophages",
    "ct_Mono+monorived_cells_MASLD": "Mono+mono derived cells",
    "ct_Plasma_cells_MASLD": "Plasma cells",
    "ct_Resident_NK_MASLD": "Resident NK",
    "ct_T_cells_MASLD": "T cells",
}

TOP_MARKERS_PER_SUBTYPE = 30


def build_composition(sc_dir: Path):
    df = pd.read_csv(sc_dir / "cell_type_proportions.csv")
    present = [c for c in CELL_TYPES if c in df.columns]
    comp = []
    for _, x in df.iterrows():
        comp.append({
            "sample": str(x["sample"]),
            "dataset": str(x["dataset"]),
            "condition": str(x["condition"]),
            "condition_harmonized": str(x["condition_harmonized"]),
            "fractions": {ct: r(x[ct], 4) for ct in present},
        })
    return comp, present


def build_activity(sc_dir: Path):
    act = pd.read_csv(sc_dir / "disease_signatures/celltype_disease_activity.csv",
                      index_col=0)
    out = []
    for col, ct in ACT_COL_TO_CT.items():
        if col in act.columns and ct in act.index:
            v = r(act.loc[ct, col], 4)
            if v is not None:
                out.append({"cell_type": ct, "activity": v})
    out.sort(key=lambda d: d["cell_type"])
    return out


def build_fgsea(sc_dir: Path):
    f = pd.read_csv(sc_dir / "disease_signatures/fgsea_celltype_F_transitions.csv")
    rows = []
    for _, x in f.iterrows():
        padj = r(x["padj"], 4)
        rows.append({
            "cell_type": str(x["cell_type"]),
            "contrast": str(x["contrast"]),
            "nes": r(x["NES"], 3),
            "padj": padj,
            "size": int(x["size"]),
            "leading_edge_n": int(x["leadingEdge_n"]),
            "significant": bool(padj is not None and padj < 0.05),
        })
    return rows


def build_hep_subtypes(sc_dir: Path):
    enr = pd.read_csv(sc_dir / "hepatocyte_subtypes/subtype_disease_enrichment.csv")
    axes = pd.read_csv(sc_dir / "hepatocyte_subtypes/subtype_axis_scores.csv")
    axes = axes.set_index("subtype")
    rows = []
    for _, x in enr.iterrows():
        st = x["subtype"]
        ax = axes.loc[st] if st in axes.index else None
        rows.append({
            "subtype": str(st),
            "n_cells": int(x["n_cells"]),
            "odds_ratio": r(x["odds_ratio"], 4),
            "fisher_padj": r(x["fisher_padj"], 6),
            "enrichment_class": str(x["enrichment_class"]),
            "stage_spearman_rho": r(x["stage_spearman_rho"], 4),
            "axes": {
                "periportal": r(ax["axis_periportal"], 4) if ax is not None else None,
                "pericentral": r(ax["axis_pericentral"], 4) if ax is not None else None,
                "lipid_accumulation": r(ax["axis_lipid_accumulation"], 4) if ax is not None else None,
                "inflammatory": r(ax["axis_inflammatory"], 4) if ax is not None else None,
                "ferroptosis": r(ax["axis_ferroptosis"], 4) if ax is not None else None,
            },
        })
    return rows


def build_pseudobulk_parquet(sc_dir: Path, out: Path, emap: dict):
    de = pd.read_csv(sc_dir / "pseudobulk_de/allcell_pseudobulk_de.csv")
    # gene col is (mostly) Ensembl; resolve to human_symbol, fallback to id.
    sym = de["gene"].map(lambda g: emap.get(strip_version(g), str(g)))
    pq = pd.DataFrame({
        "symbol": sym,
        "logfc": de["lfc"].round(3),
        "ave_expr": de["ave"].round(3),
        "tstat": de["t_stat"].round(3),
        "pval": de["pvalue"].round(6),
        "padj": de["padj"].round(6),
    })
    write_parquet(pq, out / "sc_pseudobulk_de.parquet")
    return len(pq)


def build_hep_markers_parquet(sc_dir: Path, out: Path):
    mk = pd.read_csv(sc_dir / "hepatocyte_subtypes/subtype_markers.csv")
    mk = mk.sort_values(["subtype", "scores"], ascending=[True, False])
    top = mk.groupby("subtype", sort=True).head(TOP_MARKERS_PER_SUBTYPE)
    pq = pd.DataFrame({
        "symbol": top["names"].astype(str),
        "subtype": top["subtype"].astype(str),
        "logfoldchanges": top["logfoldchanges"].round(3),
        "pvals_adj": top["pvals_adj"].round(6),
        "pct_nz_group": top["pct_nz_group"].round(4),
    })
    write_parquet(pq, out / "sc_hep_markers.parquet")
    return len(pq)


def main():
    root = project_root()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    args = ap.parse_args()
    out = Path(args.output_dir)
    sc_dir = root / "Analysis/SingleCell/results_gpu_v2"

    comp, present_cts = build_composition(sc_dir)
    activity = build_activity(sc_dir)
    fgsea = build_fgsea(sc_dir)
    hep = build_hep_subtypes(sc_dir)

    summary = {
        "n_cells": N_CELLS,
        "n_datasets": N_DATASETS,
        "n_samples": N_SAMPLES,
        "cell_types": present_cts,
        "composition": comp,
        "celltype_disease_activity": activity,
        "fgsea_transitions": fgsea,
        "hep_subtypes": hep,
    }
    kb = dump_json(summary, out / "singlecell_summary.json")
    print(f"  -> singlecell_summary.json: {kb:.1f} KB "
          f"(composition={len(comp)}, activity={len(activity)}, "
          f"fgsea={len(fgsea)}, hep_subtypes={len(hep)})")

    emap = load_ensembl_symbol_map(root)
    n_pb = build_pseudobulk_parquet(sc_dir, out, emap)
    print(f"  -> sc_pseudobulk_de.parquet: {n_pb} rows")
    n_mk = build_hep_markers_parquet(sc_dir, out)
    print(f"  -> sc_hep_markers.parquet: {n_mk} rows")


if __name__ == "__main__":
    main()
