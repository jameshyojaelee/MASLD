"""Create corrected obs beside immutable expression data using HDF5 links."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

parser = argparse.ArgumentParser()
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--metadata", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
source = args.source.resolve()
assert source != args.output.resolve() and not args.output.exists()
labels = pd.read_csv(args.metadata, sep="\t").set_index("sample")
assert not labels.index.has_duplicates
with h5py.File(source, "r") as original, h5py.File(args.output, "x") as target:
    for key, value in original.attrs.items():
        target.attrs[key] = value
    for key in original:
        if key != "obs":
            target[key] = h5py.ExternalLink(str(source), "/" + key)
    obs = read_elem(original["obs"])
    assert set(obs["sample"].astype(str).unique()) == set(labels.index)
    matched = labels.loc[obs["sample"].astype(str)]
    assert np.array_equal(obs["dataset"].astype(str), matched["dataset"].astype(str))
    for key in ("condition", "condition_harmonized"):
        if key in obs:
            obs[key + "_original"] = obs[key].copy()
    for key in ("donor", "condition", "condition_harmonized", "disease_stage_coarse",
                "source_diagnosis", "source_histology", "source_etiology", "qualification_reason"):
        obs[key] = pd.Categorical(matched[key].fillna("").astype(str).to_numpy())
    for key in ("binary_eligible", "exclude_stage_analysis", "stage_primary_eligible"):
        obs[key] = matched[key].astype(str).str.lower().eq("true").to_numpy()
    obs["diagnosis_release"] = pd.Categorical(["source-diagnosis-correction-20261026"] * len(obs))
    write_elem(target, "obs", obs)
with h5py.File(args.output, "r") as checked, h5py.File(source, "r") as original:
    assert checked["raw/X/data"].id == original["raw/X/data"].id
    assert checked["X"].id == original["X"].id
    assert checked["raw/X/data"].shape == original["raw/X/data"].shape
    assert len(read_elem(checked["obs"])) == len(obs)
args.output.with_suffix(".json").write_text(json.dumps({
    "source": str(source), "metadata": str(args.metadata.resolve()),
    "n_cells": len(obs), "expression": "HDF5 external links; original data unchanged",
    "source_bytes": source.stat().st_size,
}, indent=2))
print("Corrected obs written; raw counts and X resolve to the original HDF5 datasets.")
