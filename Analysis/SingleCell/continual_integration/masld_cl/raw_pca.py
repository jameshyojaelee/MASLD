"""Donor-pseudobulk raw-count PCA authority for post-lock preservation tests."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy import sparse

from .config import repo_path, write_json_exclusive
from .contracts import ContractError, sha256_path, verify_contract_lock
from .firewall import load_outcome_selection_lock


def log_cpm_pca(counts: np.ndarray, components: int) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim != 2 or np.any(counts < 0) or not np.isfinite(counts).all():
        raise ContractError("pseudobulk counts must be a finite non-negative matrix")
    totals = counts.sum(axis=1)
    if np.any(totals <= 0):
        raise ContractError("pseudobulk donor has zero total counts")
    values = np.log1p(counts * (1_000_000.0 / totals[:, None]))
    values -= values.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(values, full_matrices=False)
    n_components = min(int(components), max(1, len(values) - 1), len(singular))
    return u[:, :n_components] * singular[:n_components]


def _read_cell_groups(config: dict[str, Any], contract_dir: Path):
    powered = set(config["evaluation"]["powered_query_studies"])
    lineages = set(config["lineages"])
    group_keys: list[tuple[str, str, str]] = []
    group_index: dict[tuple[str, str, str], int] = {}
    cell_group: list[int] = []
    donor_control: dict[tuple[str, str], bool] = {}
    with gzip.open(contract_dir / "cell_manifest.tsv.gz", "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            eligible = row["analysis_eligible"].lower() == "true"
            study, lineage, donor = row["dataset"], row["cell_type"], row["donor_id"]
            if eligible and study in powered and lineage in lineages:
                key = (study, lineage, donor)
                if key not in group_index:
                    group_index[key] = len(group_keys)
                    group_keys.append(key)
                cell_group.append(group_index[key])
                control = row["query_control"].lower() == "true"
                donor_key = (study, donor)
                if donor_key in donor_control and donor_control[donor_key] != control:
                    raise ContractError(f"query-control status varies within donor: {donor}")
                donor_control[donor_key] = control
            else:
                cell_group.append(-1)
    return np.asarray(cell_group, dtype=np.int32), group_keys, donor_control


def _aggregate_raw_counts(
    atlas: Path, cell_group: np.ndarray, n_groups: int, chunk_rows: int = 8192,
) -> np.ndarray:
    with h5py.File(atlas, "r") as handle:
        raw = handle["raw/X"]
        n_cells, n_genes = map(int, raw.attrs["shape"])
        if len(cell_group) != n_cells:
            raise ContractError("cell manifest order differs from raw/X")
        aggregate = np.zeros((n_groups, n_genes), dtype=np.float64)
        for start in range(0, n_cells, chunk_rows):
            stop = min(start + chunk_rows, n_cells)
            codes = cell_group[start:stop]
            selected = np.flatnonzero(codes >= 0)
            if not len(selected):
                continue
            pointers = raw["indptr"][start : stop + 1].astype(np.int64)
            data_start, data_stop = int(pointers[0]), int(pointers[-1])
            block = sparse.csr_matrix(
                (
                    raw["data"][data_start:data_stop],
                    raw["indices"][data_start:data_stop],
                    pointers - data_start,
                ),
                shape=(stop - start, n_genes),
            )
            aggregator = sparse.csr_matrix(
                (
                    np.ones(len(selected), dtype=np.float64),
                    (codes[selected], selected),
                ),
                shape=(n_groups, stop - start),
            )
            partial = (aggregator @ block).tocsr()
            for group in np.unique(codes[selected]):
                row = partial.getrow(int(group))
                aggregate[int(group), row.indices] += row.data
    if np.any(aggregate < 0) or not np.all(aggregate == np.floor(aggregate)):
        raise ContractError("raw pseudobulk aggregation lost count integrity")
    return aggregate


def _write_scope(
    output: Path, scope_id: str, keys: list[tuple[str, str, str]],
    counts: np.ndarray, control: dict[tuple[str, str], bool], components: int,
) -> dict[str, Any]:
    studies = np.asarray([x[0] for x in keys], dtype=object)
    donors = np.asarray([x[2] for x in keys], dtype=object)
    controls = np.asarray([control[(x[0], x[2])] for x in keys], dtype=bool)
    normalized = counts.astype(np.float64)
    totals = normalized.sum(axis=1)
    if np.any(totals <= 0):
        raise ContractError(f"zero-count pseudobulk in {scope_id}")
    normalized = np.log1p(normalized * (1_000_000.0 / totals[:, None]))
    if scope_id.startswith("POOLED_PRIMARY|"):
        for study in sorted(set(studies)):
            normalized[studies == study] -= normalized[studies == study].mean(axis=0)
        centered = normalized
        u, singular, _ = np.linalg.svd(centered, full_matrices=False)
        count = min(components, max(1, len(centered) - 1))
        latent = u[:, :count] * singular[:count]
    else:
        latent = log_cpm_pca(counts, components)
    path = output / f"{scope_id.replace('|', '__')}.npz"
    np.savez_compressed(
        path, donor_id=donors.astype(str), study=studies.astype(str),
        query_control=controls, latent=latent,
    )
    return {
        "scope": scope_id,
        "file": path.name,
        "sha256": sha256_path(path),
        "n_donors": len(donors),
        "n_controls": int(controls.sum()),
        "n_cases": int((~controls).sum()),
        "testable": bool(controls.sum() >= 3 and (~controls).sum() >= 3),
    }


def build_raw_pca_authority(
    config: dict[str, Any], contract_lock: str | Path,
    selection_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    contract = verify_contract_lock(config, contract_lock, full_hash=False)
    selection = load_outcome_selection_lock(selection_lock, config)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    contract_dir = Path(contract_lock).resolve().parent
    cell_group, group_keys, controls = _read_cell_groups(config, contract_dir)
    aggregate = _aggregate_raw_counts(
        repo_path(config, config["input"]["atlas_h5ad"]), cell_group, len(group_keys)
    )
    components = config["evaluation"]["raw_pca"]["components"]
    scopes = []
    for lineage in config["lineages"]:
        pooled_indices = [i for i, key in enumerate(group_keys) if key[1] == lineage]
        pooled_keys = [group_keys[i] for i in pooled_indices]
        scopes.append(_write_scope(
            output, f"POOLED_PRIMARY|{lineage}", pooled_keys,
            aggregate[pooled_indices], controls, components,
        ))
        for study in config["evaluation"]["powered_query_studies"]:
            indices = [
                i for i, key in enumerate(group_keys)
                if key[0] == study and key[1] == lineage
            ]
            keys = [group_keys[i] for i in indices]
            scopes.append(_write_scope(
                output, f"{study}|{lineage}", keys, aggregate[indices], controls, components,
            ))
    manifest = {
        "schema_version": "masld-cl-raw-pca-v1",
        "selection_lock_sha256": selection["lock_sha256"],
        "contract_lock_sha256": contract["lock_sha256"],
        "config_sha256": config["_config_sha256"],
        "count_space": config["evaluation"]["raw_pca"]["count_space"],
        "unit": "biological_donor_pseudobulk",
        "scopes": scopes,
    }
    write_json_exclusive(output / "raw_pca_manifest.json", manifest)
    return manifest
