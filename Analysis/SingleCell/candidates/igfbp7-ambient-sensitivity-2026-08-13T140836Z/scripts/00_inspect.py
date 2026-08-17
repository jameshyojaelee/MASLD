#!/usr/bin/env python
"""Read-only inspection of the substrate for the IGFBP7 ambient-RNA sensitivity.

SENSITIVITY ANALYSIS ONLY. Nothing here rediscovers, refits, reweights, renames
or re-selects any Hotspot program. It only reports the structure of the frozen
inputs so the correction and re-scoring steps can be written correctly.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(os.environ["MASLD_PROJECT_ROOT"])
CAND = Path(os.environ["CAND_ROOT"])
OUT = CAND / "results"
OUT.mkdir(parents=True, exist_ok=True)

GLOBAL_H5 = ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
HEP_H5 = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_hepatocytes.h5ad"
HS_ROOT = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"

report: dict = {}

print("[inspect] reading global atlas obs (backed)", flush=True)
g = ad.read_h5ad(GLOBAL_H5, backed="r")
report["global"] = {
    "path": str(GLOBAL_H5),
    "n_obs": int(g.n_obs),
    "n_vars": int(g.n_vars),
    "obs_columns": list(map(str, g.obs.columns)),
    "has_raw": g.raw is not None,
    "n_raw_vars": int(g.raw.n_vars) if g.raw is not None else None,
    "obs_names_head": list(map(str, g.obs_names[:3])),
    "var_names_head": list(map(str, g.var_names[:5])),
    "layers": list(map(str, g.layers.keys())),
    "obsm": list(map(str, g.obsm.keys())),
}
if "cell_type" in g.obs:
    report["global"]["cell_type_counts"] = {
        str(k): int(v) for k, v in g.obs["cell_type"].astype(str).value_counts().items()
    }
if "dataset" in g.obs:
    report["global"]["dataset_counts"] = {
        str(k): int(v) for k, v in g.obs["dataset"].astype(str).value_counts().items()
    }
if "cell_type" in g.obs and "dataset" in g.obs:
    ct_ds = (
        pd.crosstab(g.obs["dataset"].astype(str), g.obs["cell_type"].astype(str))
    )
    ct_ds.to_csv(OUT / "00_global_celltype_by_dataset.tsv", sep="\t")
    report["global"]["celltype_by_dataset_written"] = True

# integer-ness of the raw layer, sampled
if g.raw is not None:
    idx = np.sort(np.random.default_rng(20260813).choice(g.n_obs, size=min(2000, g.n_obs), replace=False))
    sub = g.raw.X[idx, :]
    d = np.asarray(sub[sub != 0]).ravel()[:200000]
    report["global"]["raw_is_integer_sampled"] = bool(np.allclose(d, np.round(d)))
    report["global"]["raw_sample_max"] = float(d.max()) if d.size else None
del g

print("[inspect] reading hepatocyte scoring atlas obs (backed)", flush=True)
h = ad.read_h5ad(HEP_H5, backed="r")
report["hepatocyte_atlas"] = {
    "path": str(HEP_H5),
    "n_obs": int(h.n_obs),
    "n_vars": int(h.n_vars),
    "obs_columns": list(map(str, h.obs.columns)),
    "layers": list(map(str, h.layers.keys())),
    "obsm": list(map(str, h.obsm.keys())),
    "has_raw": h.raw is not None,
    "obs_names_head": list(map(str, h.obs_names[:3])),
}
hep_names = pd.Index(h.obs_names.astype(str))
if "sample" in h.obs:
    report["hepatocyte_atlas"]["n_samples"] = int(h.obs["sample"].astype(str).nunique())
if "dataset" in h.obs:
    report["hepatocyte_atlas"]["dataset_counts"] = {
        str(k): int(v) for k, v in h.obs["dataset"].astype(str).value_counts().items()
    }
del h

print("[inspect] checking barcode compatibility global <-> hepatocyte", flush=True)
g = ad.read_h5ad(GLOBAL_H5, backed="r")
gnames = pd.Index(g.obs_names.astype(str))
inter = hep_names.intersection(gnames)
report["barcode_join"] = {
    "n_hepatocyte_cells": int(len(hep_names)),
    "n_global_cells": int(len(gnames)),
    "n_intersect": int(len(inter)),
    "hepatocyte_fully_contained": bool(len(inter) == len(hep_names)),
}
del g

print("[inspect] stored cell scores", flush=True)
cs = pd.read_parquet(HS_ROOT / "hepatocytes/cell_scores.parquet")
report["stored_cell_scores"] = {
    "path": str(HS_ROOT / "hepatocytes/cell_scores.parquet"),
    "n_rows": int(len(cs)),
    "n_cells": int(cs["cell_id"].nunique()),
    "modules": sorted(int(m) for m in cs["module"].unique()),
    "cell_id_head": list(map(str, cs["cell_id"].head(3))),
    "cell_ids_match_hep_atlas": bool(
        set(cs["cell_id"].astype(str).unique()) == set(hep_names)
    ),
}

(OUT / "00_inspection.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2), flush=True)
print("[inspect] DONE", flush=True)
