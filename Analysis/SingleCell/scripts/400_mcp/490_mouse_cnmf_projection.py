#!/usr/bin/env python
"""
490_mouse_cnmf_projection.py — Project human cNMF programs onto mouse scRNA via orthologs.

Approach:
  Human cNMF top-100 genes (by spectra_score) → convert to human Ensembl via gencode_v49 metadata →
  map to mouse Ensembl via ortholog table (ortholog_one2one only) → find mouse gene symbols →
  sc.tl.score_genes per program per cell.

Outputs:
  mouse_projection_{name}_k{k}_mean_by_celltype_condition.tsv
  mouse_projection_{name}_k{k}_gene_mapping.tsv  (# of human top genes mapped to mouse)
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
MOUSE_ATLAS = PROJECT_ROOT / "Analysis/Deconvolution/reference/reference_mouse.h5ad"
ORTHO = PROJECT_ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
GENCODE_HUMAN = PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
OUT = MCP / "validation" / "mouse"
OUT.mkdir(parents=True, exist_ok=True)


def build_human_to_mouse_symbol():
    # human symbol → human ensembl
    gm = pd.read_csv(GENCODE_HUMAN, sep="\t")
    sym2ens = dict(zip(gm["gene_name"], gm["ensembl_base"]))

    # human ensembl → mouse ensembl (1:1 orthologs)
    o = pd.read_csv(ORTHO, sep="\t")
    o = o[o["orthology_type"] == "ortholog_one2one"]
    h2m_ens = dict(zip(o["human_ensembl_gene_id"], o["mouse_ensembl_gene_id"]))

    return sym2ens, h2m_ens


def build_mouse_ens_to_symbol(mouse_adata):
    # Some atlases use gene_ids column; others put ensembl in var index. Fall back to symbol=index.
    if "gene_ids" in mouse_adata.var.columns:
        return dict(zip(mouse_adata.var["gene_ids"], mouse_adata.var.index))
    return {s: s for s in mouse_adata.var.index}


def main(args) -> None:
    name = args.name
    k = args.k
    topg = pd.read_csv(MCP / "cnmf_annot" / name / f"program_topgenes.k{k}.tsv", sep="\t")
    print(f"[490] cNMF programs: {topg['program'].nunique()}; top genes total: {len(topg)}")

    sym2ens, h2m_ens = build_human_to_mouse_symbol()

    print(f"[490] Load mouse atlas {MOUSE_ATLAS}")
    a = sc.read_h5ad(MOUSE_ATLAS)
    print(f"[490] Mouse atlas: {a.shape}")
    a.obs_names_make_unique()
    m_ens2sym = build_mouse_ens_to_symbol(a)

    # Normalize if needed
    if "log1p" not in a.uns:
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)

    # Build case-conversion fallback: human symbol XXX -> mouse symbol Xxx
    mouse_var_set = set(a.var_names)

    def human_to_mouse_symbol(hs: str):
        # Primary: ortholog chain
        e = sym2ens.get(hs)
        if e:
            m = h2m_ens.get(e)
            if m and m in m_ens2sym:
                ms = m_ens2sym[m]
                if ms in mouse_var_set:
                    return ms
        # Fallback: case conversion
        cand = hs[0].upper() + hs[1:].lower() if len(hs) > 1 else hs
        if cand in mouse_var_set:
            return cand
        return None

    mapping_rows = []
    for p, g in topg.groupby("program"):
        h_syms = g["gene_name"].tolist()[: args.top_n]
        m_syms = [human_to_mouse_symbol(s) for s in h_syms]
        m_syms = [s for s in m_syms if s]
        mapping_rows.append({
            "program": p,
            "n_human_top": len(h_syms),
            "n_mapped_to_mouse": len(m_syms),
        })
        if len(m_syms) >= 5:
            sc.tl.score_genes(a, m_syms, score_name=f"prog_{p}")

    map_df = pd.DataFrame(mapping_rows)
    map_df.to_csv(OUT / f"mouse_projection_{name}_k{k}_gene_mapping.tsv", sep="\t", index=False)
    print(f"[490] mapping: mean genes per-program mapped to mouse = {map_df['n_mapped_to_mouse'].mean():.1f}")

    prog_cols = [c for c in a.obs.columns if c.startswith("prog_")]
    if prog_cols:
        # Aggregate by sample; also compute per-cell-type if present
        group_cols = [c for c in ["sample", "cell_type", "condition", "leiden"] if c in a.obs.columns]
        for gc in group_cols:
            try:
                agg = a.obs.groupby(gc, observed=True)[prog_cols].mean()
                agg.to_csv(OUT / f"mouse_projection_{name}_k{k}_mean_by_{gc}.tsv", sep="\t")
                print(f"[490] wrote aggregation by {gc}: {agg.shape}")
            except Exception as e:
                print(f"[490] groupby {gc} failed: {e}")
        # Global program score stats
        score_stats = a.obs[prog_cols].agg(["mean", "std", "min", "max"]).T
        score_stats.to_csv(OUT / f"mouse_projection_{name}_k{k}_score_stats.tsv", sep="\t")
        print(f"[490] scored {len(prog_cols)} programs on {a.shape[0]:,} mouse cells")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="global")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--top-n", type=int, default=100)
    main(ap.parse_args())
