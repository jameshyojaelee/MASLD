#!/usr/bin/env python
"""
480_spatial_validation.py — Score cNMF GEPs on Visium spots and measure spatial coherence.

Inputs:
  cNMF gene-spectra (program top genes)
  Visium samples: GSE192741 (6,546 spots) + Vu et al. (10 arrays, 17.5K spots)

Outputs per program:
  Moran's I (spatial autocorrelation)
  zone enrichment (periportal / midzonal / pericentral)
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from libpysal.weights import KNN
from esda.moran import Moran

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
MCP = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
SPATIAL = PROJECT_ROOT / "Analysis/Spatial/results"
OUT = MCP / "validation" / "spatial"
OUT.mkdir(parents=True, exist_ok=True)


def load_programs(name: str, k: int, top_n: int = 100):
    topg = pd.read_csv(MCP / "cnmf_annot" / name / f"program_topgenes.k{k}.tsv", sep="\t")
    programs = {}
    for p, g in topg.groupby("program"):
        programs[str(p)] = g.nlargest(top_n, "spectra_score")["gene_name"].tolist()
    return programs


def score_visium(adata, programs):
    """Score each program via scanpy tl.score_genes."""
    for p, genes in programs.items():
        avail = [g for g in genes if g in adata.var_names]
        if len(avail) < 5:
            continue
        try:
            sc.tl.score_genes(adata, avail, score_name=f"program_{p}", use_raw=False)
        except Exception as e:
            print(f"[480] score_genes {p} failed: {e}")
    return adata


def moran_autocorr(adata, prog_cols, permutations: int = 999):
    coords = np.column_stack([adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]])
    w = KNN.from_array(coords, k=8)
    w.transform = "r"
    rows = []
    for c in prog_cols:
        y = adata.obs[c].values.astype(float)
        y = y - y.mean()
        try:
            m = Moran(y, w, permutations=permutations)
            rows.append({
                "program": c,
                "morans_I": float(m.I),
                "p_sim": float(m.p_sim),
                "p_z_sim": float(getattr(m, "p_z_sim", np.nan)),
                "EI_sim": float(getattr(m, "EI_sim", np.nan)),
                "VI_sim": float(getattr(m, "VI_sim", np.nan)),
                "n_perm": int(permutations),
            })
        except Exception as e:
            rows.append({
                "program": c,
                "morans_I": np.nan,
                "p_sim": np.nan,
                "p_z_sim": np.nan,
                "EI_sim": np.nan,
                "VI_sim": np.nan,
                "n_perm": int(permutations),
            })
    return pd.DataFrame(rows)


def main(args: argparse.Namespace) -> None:
    programs = load_programs(args.name, args.k, args.top_n)
    print(f"[480] {len(programs)} programs loaded.")

    # Find Visium AnnData candidates
    candidates = list((SPATIAL / "preprocessed").glob("*.h5ad"))
    candidates += list((PROJECT_ROOT / "data/vu_et_al_2025/raw/Visium").glob("*/outs/*.h5ad"))
    if not candidates:
        print("[480] no Visium h5ad found under Analysis/Spatial/results/preprocessed/ or Vu raw")
        return

    agg_rows = []
    for vf in candidates:
        print(f"[480] Load {vf}")
        try:
            ad_sp = sc.read_h5ad(vf)
        except Exception as e:
            print(f"[480] skip {vf}: {e}")
            continue
        ad_sp = score_visium(ad_sp, programs)
        prog_cols = [c for c in ad_sp.obs.columns if c.startswith("program_")]
        if not prog_cols:
            continue
        moran = moran_autocorr(ad_sp, prog_cols, permutations=args.permutations)
        moran["sample"] = vf.stem
        agg_rows.append(moran)

    if agg_rows:
        out = pd.concat(agg_rows, ignore_index=True)
        # BH-FDR over all (program, sample) tests using p_sim
        try:
            from scipy.stats import false_discovery_control
            mask = out["p_sim"].notna()
            if mask.any():
                out.loc[mask, "q_bh"] = false_discovery_control(
                    out.loc[mask, "p_sim"].values, method="bh"
                )
        except Exception as e:
            print(f"[480] BH-FDR fallback: {e}")
            mask = out["p_sim"].notna()
            if mask.any():
                pv = out.loc[mask, "p_sim"].values
                order = np.argsort(pv)
                ranks = np.empty_like(order)
                ranks[order] = np.arange(1, len(pv) + 1)
                qv = pv * len(pv) / ranks
                qv = np.minimum.accumulate(qv[order[::-1]])[::-1]
                qfull = np.empty_like(pv)
                qfull[order] = qv
                out.loc[mask, "q_bh"] = qfull

        suffix = f"_{args.permutations}perms" if args.permutations != 999 else ""
        out_f = OUT / f"spatial_moran_{args.name}_k{args.k}{suffix}.tsv"
        out.to_csv(out_f, sep="\t", index=False)
        print(f"[480] wrote {len(out)} rows -> {out_f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="global")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--permutations", type=int, default=999,
                    help="Moran permutations (Echo A16: 9999)")
    main(ap.parse_args())
