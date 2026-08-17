"""Export the read-only incumbent Harmony representation as an audited bundle."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from .config import repo_path, write_json_exclusive
from .contracts import (
    ContractError, DeterministicGzipTextWriter, sha256_path, verify_contract_lock,
)
from .data import read_library_manifest


def export_harmony(
    config: dict[str, Any], contract_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    contract_dir = Path(contract_lock).resolve().parent
    libraries = read_library_manifest(contract_dir / "library_manifest.tsv")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    cells_path = output / "embedding_cells.tsv.gz"
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    with h5py.File(atlas, "r") as handle:
        if "X_pca_harmony" not in handle["obsm"]:
            raise ContractError("canonical atlas lacks obsm/X_pca_harmony")
        latent = handle["obsm/X_pca_harmony"][:]
    if latent.shape[0] != lock["observed_contract"]["descriptive"]["cells"]:
        raise ContractError("Harmony cell census differs from the contract")
    np.save(latent_path, latent)
    cell_manifest = contract_dir / "cell_manifest.tsv.gz"
    with gzip.open(cell_manifest, "rt", newline="") as source, DeterministicGzipTextWriter(
        cells_path
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        fieldnames = [
            "row_index", "cell_id", "library_id", "donor_id", "dataset",
            "audit_cell_type", "predicted_cell_type", "preparation", "technical_batch",
            "strict_reference", "primary_query", "query_control", "analysis_eligible",
        ]
        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for count, row in enumerate(reader, start=1):
            library = libraries[row["library_id"]]
            writer.writerow({
                "row_index": count - 1,
                "cell_id": row["cell_id"],
                "library_id": row["library_id"],
                "donor_id": row["donor_id"],
                "dataset": row["dataset"],
                "audit_cell_type": row["cell_type"],
                "predicted_cell_type": row["cell_type"],
                "preparation": library["preparation"],
                "technical_batch": library["technical_batch"],
                "strict_reference": row["strict_reference"],
                "primary_query": row["primary_query"],
                "query_control": row["query_control"],
                "analysis_eligible": row["analysis_eligible"],
            })
    if count != latent.shape[0]:
        raise ContractError("Harmony latent and cell manifest row counts differ")
    manifest = {
        "schema_version": "masld-cl-embedding-v1",
        "method": "incumbent_scalesc_pca_harmony",
        "model_kind": "all_lineage",
        "n_cells": int(latent.shape[0]),
        "n_latent": int(latent.shape[1]),
        "latent_file": latent_path.name,
        "cells_file": cells_path.name,
        "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "config_sha256": config["_config_sha256"],
        "contract_lock_sha256": lock["lock_sha256"],
    }
    write_json_exclusive(output / "embedding_manifest.json", manifest)
    return manifest


def export_uncorrected_pca(
    config: dict[str, Any], contract_lock: str | Path, output: str | Path,
    dimensions: int = 30,
) -> dict[str, Any]:
    """Export the uncorrected PCA with the same audited cell roster as Harmony."""
    if int(dimensions) != int(config["architecture"]["n_latent"]):
        raise ContractError("uncorrected PCA export must match the configured latent dimension")
    lock = verify_contract_lock(config, contract_lock, full_hash=False)
    contract_dir = Path(contract_lock).resolve().parent
    libraries = read_library_manifest(contract_dir / "library_manifest.tsv")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    latent_path = output / "embedding_latent.npy"
    cells_path = output / "embedding_cells.tsv.gz"
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    with h5py.File(atlas, "r") as handle:
        if "X_pca" not in handle["obsm"] or handle["obsm/X_pca"].shape[1] < dimensions:
            raise ContractError("canonical atlas lacks the required uncorrected PCA")
        latent = handle["obsm/X_pca"][:, :dimensions]
    if latent.shape[0] != lock["observed_contract"]["descriptive"]["cells"]:
        raise ContractError("uncorrected PCA cell census differs from the contract")
    np.save(latent_path, latent)
    cell_manifest = contract_dir / "cell_manifest.tsv.gz"
    with gzip.open(cell_manifest, "rt", newline="") as source, DeterministicGzipTextWriter(
        cells_path
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        fieldnames = [
            "row_index", "cell_id", "library_id", "donor_id", "dataset",
            "audit_cell_type", "predicted_cell_type", "preparation", "technical_batch",
            "strict_reference", "primary_query", "query_control", "analysis_eligible",
        ]
        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for count, row in enumerate(reader, start=1):
            library = libraries[row["library_id"]]
            writer.writerow({
                "row_index": count - 1, "cell_id": row["cell_id"],
                "library_id": row["library_id"], "donor_id": row["donor_id"],
                "dataset": row["dataset"], "audit_cell_type": row["cell_type"],
                "predicted_cell_type": row["cell_type"],
                "preparation": library["preparation"],
                "technical_batch": library["technical_batch"],
                "strict_reference": row["strict_reference"],
                "primary_query": row["primary_query"],
                "query_control": row["query_control"],
                "analysis_eligible": row["analysis_eligible"],
            })
    if count != latent.shape[0]:
        raise ContractError("uncorrected PCA and cell manifest row counts differ")
    manifest = {
        "schema_version": "masld-cl-embedding-v1",
        "method": "uncorrected_scalesc_pca_30d",
        "model_kind": "all_lineage", "n_cells": int(latent.shape[0]),
        "n_latent": int(latent.shape[1]), "latent_file": latent_path.name,
        "cells_file": cells_path.name, "latent_sha256": sha256_path(latent_path),
        "cells_sha256": sha256_path(cells_path),
        "config_sha256": config["_config_sha256"],
        "contract_lock_sha256": lock["lock_sha256"],
    }
    write_json_exclusive(output / "embedding_manifest.json", manifest)
    return manifest
