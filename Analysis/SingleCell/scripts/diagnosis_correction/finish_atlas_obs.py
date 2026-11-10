"""Validate the corrected AnnData obs and add the source-measured F field."""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
labels = pd.read_csv(task / "metadata/library_metadata_source_F.tsv", sep="\t").set_index("sample")
path = task / "integrated_atlas_source_qualified_v2.h5ad"
with h5py.File(path, "r+") as f:
    assert isinstance(f.get("obs", getlink=True), h5py.HardLink)
    for col in f["obs"]:
        obj = f["obs"][col]
        if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
            obj.attrs.setdefault("ordered", False)
    obs = read_elem(f["obs"])
    matched = labels.loc[obs["sample"].astype(str)]
    obs["F_stage_documented"] = pd.to_numeric(matched["F_stage_documented"], errors="coerce").to_numpy()
    obs["F_stage_source"] = pd.Categorical(np.where(np.isfinite(obs["F_stage_documented"]), "deposited_numeric_SAF", "not_numerically_documented"))
    assert len(obs) == 1232318
    assert np.array_equal(obs.condition.astype(str), matched.condition.astype(str))
    write_elem(f, "obs", obs)
with h5py.File(path, "r") as f:
    checked = read_elem(f["obs"])
    assert len(checked) == 1232318 and "preparation_method" in checked
    assert "F_stage_documented" in checked
    assert checked.loc[checked.donor.eq("GSE202379_P98"), "F_stage_documented"].eq(1).all()
    assert checked.loc[checked.donor.isin(["GSE202379_PHL1", "GSE202379_PHL2"]), "F_stage_documented"].isna().all()
(task / "ATLAS_OBS_VALIDATED").touch()
print("AnnData obs roundtrip passed; expression links unchanged.")
