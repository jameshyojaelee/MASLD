#!/usr/bin/env python
"""Reproduction gate for the IGFBP7 SENSITIVITY ANALYSIS.

Rebuilds the frozen hepatocyte Hotspot run and checks that re-scoring the frozen
programs on RAW counts reproduces the STORED per-cell scores. This must pass
before any ambient-corrected score is interpreted. Read-only; scores nothing new.
"""
from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import hotspot
import numpy as np
import pandas as pd
from hotspot import modules as hotspot_modules

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
RES = CAND / "results"
RES.mkdir(parents=True, exist_ok=True)
HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
SCORE_501 = ROOT / "Analysis/SingleCell/scripts/hotspot_modules/501_run_hotspot.py"
MEMBERSHIP = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
)

spec = importlib.util.spec_from_file_location("hotspot_run_501", SCORE_501)
run501 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run501)

rm = json.loads((HS_ROOT / "hepatocytes/run_metadata.json").read_text(encoding="utf-8"))
adata = run501.load_atlas("hepatocytes", smoke=False)
excluded = rm.get("exclude_datasets") or []
if excluded:
    adata = adata[~adata.obs["dataset"].isin(set(excluded))].copy()
adata = run501.strip_confounders(adata)
adata = run501.filter_detected(adata, min_frac=0.01)
assert adata.n_obs == int(rm["n_cells"]), "cell census drift"
assert adata.n_vars == int(rm["n_genes_kept"]), "gene census drift"
run501.ensure_raw_layer(adata)
latent = run501.resolve_latent(adata)
hs = hotspot.Hotspot(adata, layer_key="counts", model="danb",
                     latent_obsm_key=latent, umi_counts_obs_key="n_counts")
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
print(f"[gate] rebuilt {adata.n_obs:,} cells x {adata.n_vars:,} genes, latent={latent}", flush=True)

stored = pd.read_parquet(HS_ROOT / "hepatocytes/cell_scores.parquet")
stored["module"] = stored["module"].astype(int)
mem = pd.read_csv(MEMBERSHIP, sep="\t")
mem = mem[mem["cell_type"] == "hepatocytes"]
cell_index = pd.Index(adata.obs_names.astype(str))

rows = []
for module in sorted(mem["module"].astype(int).unique()):
    genes = mem.loc[mem["module"].astype(int) == module, "source_gene"].astype(str).tolist()
    counts = hs._counts_from_anndata(adata[:, genes], hs.layer_key, dense=True)
    s = hotspot_modules.compute_scores(
        counts, hs.model, hs.umi_counts.values, hs.neighbors.values, hs.weights.values
    )
    v = stored[stored["module"] == module].set_index("cell_id")["score"].reindex(cell_index).to_numpy(float)
    r = float(np.corrcoef(s, v)[0, 1])
    e = float(np.max(np.abs(s - v)))
    rows.append({"module": module, "n_genes": len(genes), "reproduction_r": r, "max_abs_error": e})
    print(f"[gate] module {module:>2}  n={len(genes):>4}  r={r:.10f}  max_abs_err={e:.6g}", flush=True)

out = pd.DataFrame(rows)
out.to_csv(RES / "02b_reproduction_gate.tsv", sep="\t", index=False)
print(f"\n[gate] min r = {out['reproduction_r'].min():.10f}", flush=True)
print(f"[gate] modules below 0.995: {out.loc[out['reproduction_r'] < 0.995, 'module'].tolist()}", flush=True)
