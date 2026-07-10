#!/usr/bin/env python3
"""Spatial layer web data (data contract §5.2).

Sources (under Analysis/Spatial/results/):
  zonation/zonation_scores.csv               -> n_samples, zonation_dist, by_condition
  zonation/deg_zonation_classification.csv   -> spatial_zonation.parquet (gene table)
  zonation/zonation_disruption_scores.csv    -> top_disrupted
  svg/differential_svgs.csv                  -> svg_categories, top SVGs, joined into parquet

Outputs (-> --output-dir):
  spatial_summary.json
  spatial_zonation.parquet   (long, sorted by symbol)

Run:
  micromamba run -n spatial python scripts/portal/generate_spatial_data.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from _portal_io import dump_json, project_root, r, write_parquet

TOP_N = 30


def main():
    root = project_root()
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", default=str(root / "masld-atlas-v2/public/data"))
    args = ap.parse_args()
    out = Path(args.output_dir)
    sp = root / "Analysis/Spatial/results"

    # --- zonation spot scores ---
    zs = pd.read_csv(sp / "zonation/zonation_scores.csv")
    n_samples = int(zs["sample_id"].nunique())
    datasets = sorted(str(d) for d in zs["dataset"].dropna().unique())
    zonation_dist = [
        {"zonation_bin": str(b), "n": int(n)}
        for b, n in zs["zonation_bin"].value_counts().sort_index().items()
    ]
    by_cond = []
    for cond, g in zs.groupby("condition"):
        by_cond.append({
            "condition": str(cond),
            "mean_zonation_score": r(g["zonation_score"].mean(), 4),
            "n": int(len(g)),
        })
    by_cond.sort(key=lambda d: d["condition"])

    # --- zonation disruption (top disrupted genes) ---
    dis = pd.read_csv(sp / "zonation/zonation_disruption_scores.csv")
    dis = dis.sort_values("disruption_score", ascending=False).head(TOP_N)
    top_disrupted = [{
        "gene": str(x["gene"]),
        "rho_healthy": r(x["rho_healthy"], 4),
        "rho_masld": r(x["rho_masld"], 4),
        "disruption_score": r(x["disruption_score"], 4),
        "underpowered": bool(str(x["underpowered"]).strip().upper() == "TRUE"),
    } for _, x in dis.iterrows()]

    # --- differential SVGs ---
    svg = pd.read_csv(sp / "svg/differential_svgs.csv", index_col=0)
    svg.index.name = "gene"
    svg = svg.reset_index()
    svg_categories = [
        {"category": str(c), "n": int(n)}
        for c, n in svg["category"].value_counts().items()
    ]
    top_svg = svg.reindex(svg["delta_I"].abs().sort_values(ascending=False).index).head(TOP_N)
    top_differential_svgs = [{
        "gene": str(x["gene"]),
        "morans_i_healthy": r(x["morans_I_healthy"], 4),
        "morans_i_masld": r(x["morans_I_masld"], 4),
        "delta_i": r(x["delta_I"], 4),
        "category": str(x["category"]),
    } for _, x in top_svg.iterrows()]

    summary = {
        "n_samples": n_samples,
        "datasets": datasets,
        "zonation_dist": zonation_dist,
        "zonation_by_condition": by_cond,
        "top_disrupted": top_disrupted,
        "svg_categories": svg_categories,
        "top_differential_svgs": top_differential_svgs,
    }
    kb = dump_json(summary, out / "spatial_summary.json")
    print(f"  -> spatial_summary.json: {kb:.1f} KB "
          f"(n_samples={n_samples}, datasets={len(datasets)}, "
          f"bins={len(zonation_dist)}, svg_cats={len(svg_categories)})")

    # --- gene-level zonation parquet (join differential SVGs by gene) ---
    deg = pd.read_csv(sp / "zonation/deg_zonation_classification.csv", index_col=0)
    svg_join = svg[["gene", "morans_I_healthy", "morans_I_masld", "delta_I", "category"]].rename(
        columns={"morans_I_healthy": "morans_i_healthy",
                 "morans_I_masld": "morans_i_masld",
                 "delta_I": "delta_i",
                 "category": "svg_category"})
    m = deg.merge(svg_join, on="gene", how="left")
    pq = pd.DataFrame({
        "symbol": m["gene"].astype(str),
        "zonation_class": m["zonation_class"].astype(str),
        "spearman_rho": m["spearman_rho"].round(4),
        "kruskal_pval": m["kruskal_pval"].round(6),
        "n_donors": m["n_donors"].astype("Int64"),
        "mean_PP1": m["mean_PP1"].round(4),
        "mean_PP2": m["mean_PP2"].round(4),
        "mean_Mid": m["mean_Mid"].round(4),
        "mean_PC2": m["mean_PC2"].round(4),
        "mean_PC1": m["mean_PC1"].round(4),
        "morans_i_healthy": m["morans_i_healthy"].round(4),
        "morans_i_masld": m["morans_i_masld"].round(4),
        "delta_i": m["delta_i"].round(4),
        "svg_category": m["svg_category"],
    })
    write_parquet(pq, out / "spatial_zonation.parquet")
    print(f"  -> spatial_zonation.parquet: {len(pq)} rows "
          f"({pq['svg_category'].notna().sum()} with SVG join)")


if __name__ == "__main__":
    main()
