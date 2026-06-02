#!/usr/bin/env python
"""
410b_patch_stage_metadata.py — rewrite the obs table of existing prep h5ads after fixing
the stage-mapping bug (`.any()` falsy-zero bug in 410). Only touches obs; keeps X/layers.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
HEP_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"
INPUTS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"


def stage_map() -> pd.DataFrame:
    hep = pd.read_csv(HEP_META)

    def mode_or_nan(s):
        d = s.dropna()
        if len(d) == 0:
            return np.nan
        return d.mode().iloc[0]

    return hep.groupby("sample").agg(
        disease_stage_coarse=("disease_stage_coarse", mode_or_nan),
        disease_stage_numeric=("disease_stage_numeric", mode_or_nan),
        condition_binary=("condition_binary", "first"),
    ).reset_index()


def patch_file(fp: Path, sm: pd.DataFrame) -> None:
    print(f"[410b] {fp.name}...", flush=True)
    a = ad.read_h5ad(fp)
    obs = a.obs.copy()
    # Drop previous (buggy) columns
    for c in ("disease_stage_coarse", "disease_stage_numeric", "condition_binary"):
        if c in obs.columns and c != "condition_binary":
            obs = obs.drop(columns=[c])
    # Merge fresh
    merged = obs.reset_index().merge(sm, how="left", on="sample")
    merged = merged.set_index(obs.index.name or "cell_barcode") if obs.index.name else merged.set_index(merged.columns[0])
    # Make sure same order
    a.obs = merged.loc[a.obs.index.astype(str).intersection(merged.index.astype(str))] if False else merged
    # A simpler approach: reindex merged by current obs index
    merged2 = obs.reset_index().rename(columns={obs.index.name or 'index': 'cell_barcode'}).merge(
        sm[["sample", "disease_stage_coarse", "disease_stage_numeric"]],
        how="left", on="sample"
    ).set_index("cell_barcode")
    merged2 = merged2.loc[a.obs.index.astype(str)]
    a.obs = merged2
    # Write back
    a.write_h5ad(fp, compression="gzip")
    # Summary
    stages = merged2["disease_stage_numeric"].value_counts(dropna=False).to_dict()
    print(f"[410b]   done. stage dist: {stages}", flush=True)


def rebuild_donor_metadata(sm: pd.DataFrame) -> None:
    """Rebuild donor_metadata.tsv using fresh stage map and recomputed cell-type fractions."""
    glob_f = INPUTS / "atlas_cnmf_global.h5ad"
    if not glob_f.exists():
        print("[410b] global h5ad missing; skip donor_metadata rebuild")
        return
    a = ad.read_h5ad(glob_f, backed="r")
    donor = a.obs.groupby("sample", observed=True).agg(
        dataset=("dataset", "first"),
        condition=("condition", "first"),
        n_cells=("cell_type", "size"),
    ).reset_index()
    ct_wide = a.obs.groupby(["sample", "cell_type"], observed=True).size().unstack(fill_value=0)
    ct_wide = ct_wide.divide(ct_wide.sum(axis=1), axis=0).add_prefix("frac_")
    donor = donor.merge(sm, on="sample", how="left").merge(ct_wide.reset_index(), on="sample", how="left")
    out = INPUTS / "donor_metadata.tsv"
    donor.to_csv(out, sep="\t", index=False)
    a.file.close()
    print(f"[410b] rebuilt donor_metadata.tsv ({len(donor)} rows)")
    print(donor["disease_stage_numeric"].value_counts(dropna=False).to_dict())


def main() -> None:
    sm = stage_map()
    print(f"[410b] stage_map size={len(sm)}; stage=0 samples={(sm['disease_stage_numeric']==0).sum()}")
    targets = list(INPUTS.glob("atlas_cnmf_*.h5ad"))
    print(f"[410b] {len(targets)} h5ads to patch")
    for fp in targets:
        try:
            patch_file(fp, sm)
        except Exception as e:
            print(f"[410b] ERROR patching {fp.name}: {e}")
    rebuild_donor_metadata(sm)


if __name__ == "__main__":
    main()
