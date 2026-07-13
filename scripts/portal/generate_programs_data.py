#!/usr/bin/env python3
"""Programs layer web data (data contract §5.4).

Framing: k=6 NMF programs = CONTINUOUS interpretive axes P1-P6, NOT discrete
subtypes (2026-07-04 audit). Not patient stratification.

Sources:
  RNA-seq/results/subtypes/program_labels.csv    -> nmf_programs (P1-P6 labels)
  RNA-seq/results/subtypes/nmf_assignments.csv    -> mean_activity_by_stage
  RNA-seq/results/subtypes/subtype_markers.csv    -> program_gene_membership.parquet
  .../integration/metadata/unified_metadata.csv   -> fibrosis_stage per sample
  Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv -> hotspot_modules
  Analysis/SingleCell/results_gpu_v2/mcp/crn_changepoints_clean.csv  -> changepoints

Outputs (-> --output-dir):
  programs_summary.json
  program_gene_membership.parquet   (long, sorted by symbol)

Run:
  micromamba run -n spatial python scripts/portal/generate_programs_data.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _portal_io import (dump_json, load_ensembl_symbol_map, project_root, r,
                        strip_version, write_parquet)

PROGRAMS = ["P1", "P2", "P3", "P4", "P5", "P6"]
STAGE_LABELS = {0.0: "F0", 1.0: "F1", 2.0: "F2", 3.0: "F3", 4.0: "F4"}


def _tf(v):
    return bool(str(v).strip().upper() == "TRUE")


def build_nmf_programs(root: Path):
    labels = pd.read_csv(root / "RNA-seq/results/subtypes/program_labels.csv")
    nmf = pd.read_csv(root / "RNA-seq/results/subtypes/nmf_assignments.csv")
    meta = pd.read_csv(
        root / "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/"
        "unified_metadata.csv", usecols=["sample_id", "fibrosis_stage"],
        low_memory=False)
    stage = dict(zip(meta["sample_id"], meta["fibrosis_stage"]))
    nmf = nmf.copy()
    nmf["stage"] = nmf["sample_id"].map(stage)
    nmf = nmf.dropna(subset=["stage"])

    # mean program activity per fibrosis stage
    by_stage = {}
    for st in sorted(nmf["stage"].unique()):
        g = nmf[nmf["stage"] == st]
        by_stage[st] = {p: g[p].mean() for p in PROGRAMS if p in g.columns}

    programs = []
    for _, x in labels.iterrows():
        code = str(x["program_code"])
        genes = [t.strip() for t in str(x["top_genes"]).split(",") if t.strip()]
        mabs = [{"stage": STAGE_LABELS.get(st, str(st)),
                 "mean_score": r(by_stage[st].get(code), 4)}
                for st in sorted(by_stage.keys())]
        programs.append({
            "code": code,
            "label": str(x["biological_label"]),
            "label_category": str(x["label_category"]),
            "top_pathway": str(x["top_pathway"]),
            "top_genes": genes,
            "mean_activity_by_stage": mabs,
        })
    return programs


def build_hotspot(root: Path):
    am = pd.read_csv(
        root / "Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv",
        sep="\t")
    mods = []
    for _, x in am.iterrows():
        mods.append({
            "cell_type": str(x["cell_type"]),
            "module": str(x["module"]),
            "module_name": str(x["module_name"]),
            "module_hallmark": str(x["module_hallmark"]),
            "is_novel": _tf(x["is_novel"]),
            "stability_score": r(x["stability_score"], 4),
            "bulk_replicated": _tf(x["bulk_replicated"]),
            "n_genes": int(x["n_genes_overlapping"]) if pd.notna(x["n_genes_overlapping"]) else None,
            "mean_bulk_logfc": r(x["mean_bulk_logFC"], 4),
            "progression_module": _tf(x["progression_module"]),
        })
    n_novel = sum(1 for m in mods if m["is_novel"])
    return mods, n_novel


def build_changepoints(root: Path):
    cp = pd.read_csv(
        root / "Analysis/SingleCell/results_gpu_v2/mcp/crn_changepoints_clean.csv")
    return [{
        "program": str(x["program"]),
        "best_model": str(x["best_model"]),
        "breakpoint_stage": r(x["breakpoint_stage"], 3),
        "delta_aic": r(x["delta_aic"], 3),
    } for _, x in cp.iterrows()]


def build_membership_parquet(root: Path, out: Path, emap: dict):
    labels = pd.read_csv(root / "RNA-seq/results/subtypes/program_labels.csv")
    code2label = dict(zip(labels["program_code"], labels["biological_label"]))
    mk = pd.read_csv(root / "RNA-seq/results/subtypes/subtype_markers.csv")
    sym = mk["gene"].map(lambda g: emap.get(strip_version(g), str(g)))
    pq = pd.DataFrame({
        "symbol": sym,
        "program_code": mk["program_code"].astype(str),
        "program_label": mk["program_code"].map(code2label).astype(str),
        "logfc": mk["logFC"].round(3),
        "padj": mk["adj.P.Val"].round(6),
        "direction": mk["direction"].astype(str),
    })
    write_parquet(pq, out / "program_gene_membership.parquet")
    return len(pq)


def main():
    root = project_root()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    args = ap.parse_args()
    out = Path(args.output_dir)

    programs = build_nmf_programs(root)
    hotspot, n_novel = build_hotspot(root)
    changepoints = build_changepoints(root)

    summary = {
        "nmf_programs": programs,
        "hotspot_modules": hotspot,
        "n_hotspot_modules": len(hotspot),
        "n_novel": n_novel,
        "changepoints": changepoints,
    }
    kb = dump_json(summary, out / "programs_summary.json")
    print(f"  -> programs_summary.json: {kb:.1f} KB "
          f"(programs={len(programs)}, hotspot={len(hotspot)}, novel={n_novel}, "
          f"changepoints={len(changepoints)})")

    emap = load_ensembl_symbol_map(root)
    n = build_membership_parquet(root, out, emap)
    print(f"  -> program_gene_membership.parquet: {n} rows")


if __name__ == "__main__":
    main()
