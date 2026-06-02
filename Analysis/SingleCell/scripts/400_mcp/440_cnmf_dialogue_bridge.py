#!/usr/bin/env python
"""
440_cnmf_dialogue_bridge.py — cross-tool concordance between cNMF GEPs and DIALOGUE meta-MCPs.

Inputs:
  results_gpu_v2/mcp/cnmf_annot/{name}/program_topgenes.k{K}.tsv (top-100 gene list per cNMF GEP)
  results_gpu_v2/mcp/cnmf_runs/{name}/{name}/{name}.gene_spectra_score.k_{K}.dt_0_03.txt (full spectra)
  results_gpu_v2/mcp/dialogue/meta_mcp_gene_loadings.tsv (aggregated top gene lists per meta-MCP x celltype)

Outputs:
  results_gpu_v2/mcp/integration/cnmf_dialogue_concordance.csv
    cols: cnmf_run, cnmf_k, cnmf_program, dialogue_meta_mcp, celltype,
          jaccard_top100, hypergeom_p, hypergeom_q,
          pearson_r (spectra x loading, on gene universe overlap), sample_spearman
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"


def hypergeom(overlap: int, set_a: int, set_b: int, universe: int) -> float:
    return float(stats.hypergeom.sf(overlap - 1, universe, set_a, set_b))


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty(n, float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        j = order[i]
        prev = min(prev, p[j] * n / (i + 1))
        ranked[j] = prev
    return ranked


def main(args: argparse.Namespace) -> None:
    cnmf_top = pd.read_csv(MCP_ROOT / "cnmf_annot" / args.name / f"program_topgenes.k{args.k}.tsv", sep="\t")
    dia_load = pd.read_csv(MCP_ROOT / "dialogue" / "meta_mcp_gene_loadings.tsv", sep="\t")
    # Take top-100 DIALOGUE per (meta_mcp, celltype) if not already
    dia_top = dia_load.groupby(["meta_mcp_id", "celltype"]).apply(lambda d: d.nlargest(100, "n_occur")).reset_index(drop=True)

    rows = []
    cnmf_programs = cnmf_top["program"].unique()
    for p in cnmf_programs:
        ca = cnmf_top.loc[cnmf_top["program"] == p, "gene_name"].tolist()
        set_a = set(ca)
        for (mm, ct), dsub in dia_top.groupby(["meta_mcp_id", "celltype"]):
            cb = dsub["gene_name"].tolist()
            set_b = set(cb)
            universe = 15000  # approx gene universe
            overlap = len(set_a & set_b)
            jac = len(set_a & set_b) / max(1, len(set_a | set_b))
            hp = hypergeom(overlap, len(set_a), len(set_b), universe)
            rows.append({
                "cnmf_program": p,
                "dialogue_meta_mcp": mm,
                "celltype": ct,
                "overlap_top100": overlap,
                "jaccard_top100": jac,
                "hypergeom_p": hp,
            })
    out = pd.DataFrame(rows)
    out["hypergeom_q"] = bh(out["hypergeom_p"].values)
    out = out.sort_values("hypergeom_q")
    outdir = MCP_ROOT / "integration"
    outdir.mkdir(exist_ok=True, parents=True)
    out.to_csv(outdir / "cnmf_dialogue_concordance.csv", index=False)
    print(f"[440] wrote {len(out)} rows to cnmf_dialogue_concordance.csv")
    print(f"[440] top hits (q<0.05): {int((out['hypergeom_q']<0.05).sum())}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="global")
    ap.add_argument("--k", type=int, required=True)
    main(ap.parse_args())
