#!/usr/bin/env python3
"""Extract a compact UMAP+metadata CSV for fig2's all-cell scRNA bridge panel.

Source: Analysis/SingleCell/results_gpu_v2/disease_signatures/scored_atlas.h5ad
        (150MB, obsm/X_umap + obs/cell_type/condition/dataset/sample for 1.23M cells)

Stage mapping reuses the logic in aggregate_nmf_celltype_stage.py: condition ->
{Healthy:0, NAFLD/MASL:1, NASH/MASH:2, Cirrhotic:3}; GSE244832 unlabeled cells
are recovered from donor_pairing.csv.

Output: Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz
        columns: umap_1, umap_2, cell_type, disease_status, disease_stage_coarse
"""

import os
import gzip
import h5py
import numpy as np
import pandas as pd

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SCORED = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/disease_signatures/scored_atlas.h5ad"
)
DONOR_PAIRING = os.path.join(BASE, "data/GSE244832/metadata/donor_pairing.csv")
PREP_METHOD = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/cell_type_proportions.csv"
)
OUT = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/atlas_umap_for_fig2.csv.gz"
)

STAGE_LABELS = {0: "Healthy", 1: "Steatosis", 2: "Steatohepatitis", 3: "Cirrhosis"}


def read_obs_categorical(f, key):
    obj = f[f"obs/{key}"]
    if "codes" in obj:
        cats = obj["categories"][:]
        if cats.dtype.kind == "O":
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c)
                             for c in cats])
        elif cats.dtype.kind == "S":
            cats = cats.astype(str)
        return cats[obj["codes"][:]]
    vals = obj[:]
    if vals.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else str(v)
                         for v in vals])
    return vals


def build_stage_vector(samples, datasets, conditions):
    cond_to_stage = {
        "Healthy": 0, "Mixed": 0,
        "NAFLD": 1, "MASL": 1,
        "NASH": 2, "MASH": 2,
        "Cirrhotic": 3,
    }
    stage = np.array([cond_to_stage.get(c, np.nan) for c in conditions], dtype=float)

    if os.path.exists(DONOR_PAIRING):
        dp = pd.read_csv(DONOR_PAIRING)
        srr_to_stage = {}
        g244_map = {"NORMAL": 0, "MASL": 1, "MASH": 2}
        for _, row in dp.iterrows():
            cond = g244_map.get(row["condition"])
            if cond is None or pd.isna(row["rna_srrs"]):
                continue
            for srr in str(row["rna_srrs"]).split(";"):
                srr = srr.strip()
                if srr:
                    srr_to_stage[srr] = cond
        g244_mask = (datasets == "GSE244832") & np.isnan(stage)
        for i in np.where(g244_mask)[0]:
            s = srr_to_stage.get(samples[i])
            if s is not None:
                stage[i] = s
    return stage


print(f"Reading {SCORED}")
with h5py.File(SCORED, "r") as f:
    umap = f["obsm/X_umap"][:]            # (N, 2)
    cell_type = read_obs_categorical(f, "cell_type")
    sample = read_obs_categorical(f, "sample")
    dataset = read_obs_categorical(f, "dataset")
    condition = read_obs_categorical(f, "condition")

n = umap.shape[0]
print(f"  {n:,} cells; UMAP shape={umap.shape}")

stage_vec = build_stage_vector(sample, dataset, condition)
stage_lab = np.array(["Unknown"] * n, dtype=object)
mask = ~np.isnan(stage_vec)
stage_lab[mask] = [STAGE_LABELS[int(s)] for s in stage_vec[mask]]

# Binary disease_status mirrors figS_singlecell.R fallback mapping
disease_status = np.where(
    np.isin(condition, ["Healthy", "Mixed"]), "Control",
    np.where(np.isin(condition, ["MASLD", "NAFLD", "MASL", "NASH", "MASH", "Cirrhotic"]),
             "MASLD", "Unknown"),
)

# Join preparation_method from cell_type_proportions.csv (sample-level lookup).
# Used downstream to filter to balanced prep methods (nuclei / cd45_negative /
# unsorted) so cell-type density-ratio plots are not confounded by sorted
# CD45+/CD45- enrichment.
prep_lookup = pd.read_csv(PREP_METHOD, usecols=["sample", "preparation_method"])
prep_lookup = prep_lookup.drop_duplicates("sample").set_index("sample")
prep_method = prep_lookup["preparation_method"].reindex(sample).fillna("unknown").to_numpy()

df = pd.DataFrame({
    "umap_1": umap[:, 0].astype(np.float32),
    "umap_2": umap[:, 1].astype(np.float32),
    "cell_type": cell_type,
    "sample": sample,
    "dataset": dataset,
    "preparation_method": prep_method,
    "disease_status": disease_status,
    "disease_stage_coarse": stage_lab,
})

n_staged = (df["disease_stage_coarse"] != "Unknown").sum()
print(f"  staged cells: {n_staged:,} ({100*n_staged/n:.1f}%)")
print("  cell-type counts:")
print(df["cell_type"].value_counts().to_string())
print("  disease_stage_coarse counts:")
print(df["disease_stage_coarse"].value_counts().to_string())
print("  preparation_method counts:")
print(df["preparation_method"].value_counts().to_string())
print("  preparation_method x stage cross-tab:")
print(pd.crosstab(df["preparation_method"], df["disease_stage_coarse"]).to_string())

print(f"Writing {OUT}")
df.to_csv(OUT, index=False, compression="gzip")
print(f"  wrote {os.path.getsize(OUT)/1e6:.1f} MB")
