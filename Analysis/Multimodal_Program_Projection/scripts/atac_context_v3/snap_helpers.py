#!/usr/bin/env python3
"""SnapATAC2-backed dataset helpers shared by peak and count producers."""

from __future__ import annotations

import csv
from pathlib import Path

import anndata as ad
import numpy as np
import snapatac2 as snap

from atac_context_v3_lib import LABEL_MAP, ContractError, project_root


ROOT = project_root()
COHORT_SPEC = {
    "GSE244832": {
        "h5ad": ROOT / "Analysis/ATAC/Human_Multiome/results/label_transfer/snapatac2_relabeled.h5ad",
        "metadata": ROOT / "data/GSE244832/metadata/donor_pairing.csv",
        "per_donor": ROOT / "Analysis/ATAC/Human_Multiome/results/snapatac2/per_donor",
    },
    "GSE281367": {
        "h5ad": ROOT / "Analysis/ATAC/Human_External/snapatac2_fast/snapatac2_processed.h5ad",
        "metadata": ROOT / "data/GSE281367/metadata/donor_pairing.csv",
        "per_donor": ROOT / "Analysis/ATAC/Human_External/snapatac2_fast/per_donor",
    },
}


def metadata_conditions(cohort: str) -> dict[str, str]:
    with COHORT_SPEC[cohort]["metadata"].open("r", encoding="utf-8", newline="") as handle:
        return {
            str(row["donor_id"]): str(row["condition"]).upper()
            for row in csv.DictReader(handle)
        }


def label_lookup(cohort: str) -> dict[tuple[str, str], str]:
    source = ad.read_h5ad(COHORT_SPEC[cohort]["h5ad"], backed="r")
    obs = source.obs[["donor_id", "cell_type"]].copy()
    obs["lineage"] = obs["cell_type"].astype(str).map(LABEL_MAP[cohort]).fillna("excluded")
    result = {
        (str(donor), str(cell)): str(lineage)
        for cell, donor, lineage in zip(obs.index, obs["donor_id"], obs["lineage"])
    }
    source.file.close()
    return result


def build_dataset(cohort: str, filename: Path):
    if filename.exists():
        raise ContractError(f"refusing to overwrite AnnDataSet: {filename}")
    conditions = metadata_conditions(cohort)
    adatas = []
    for donor in sorted(conditions):
        path = COHORT_SPEC[cohort]["per_donor"] / f"{donor}.h5ad"
        if not path.is_file():
            raise ContractError(f"missing donor h5ad: {path}")
        adatas.append((donor, snap.read(str(path), backed="r")))
    dataset = snap.AnnDataSet(adatas=adatas, filename=str(filename), add_key="donor_id")
    labels = label_lookup(cohort)
    donor_values = np.asarray(dataset.obs["donor_id"]).astype(str).tolist()
    cell_values = np.asarray(dataset.obs_names).astype(str).tolist()
    lineage_values = [
        labels.get((donor, cell), "excluded")
        for donor, cell in zip(donor_values, cell_values)
    ]
    condition_values = [conditions[donor] for donor in donor_values]
    dataset.obs["lineage_v3"] = np.asarray(lineage_values, dtype=str)
    dataset.obs["condition_v3"] = np.asarray(condition_values, dtype=str)
    return dataset, adatas


def close_dataset(dataset, adatas) -> None:
    try:
        dataset.close()
    finally:
        for _, value in adatas:
            try:
                value.close()
            except Exception:
                pass
