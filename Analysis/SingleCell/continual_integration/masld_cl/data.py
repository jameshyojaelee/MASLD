"""Prepare the fixed-count, fixed-cell HVG object for model fitting."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import (
    ContractError, read_h5ad_column, sha256_path, verify_contract_lock,
)
from .firewall import validate_program_firewall


def read_library_manifest(path: str | Path) -> dict[str, dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ContractError("library manifest is empty")
    result = {row["library_id"]: row for row in rows}
    if len(result) != len(rows):
        raise ContractError("library manifest contains duplicate IDs")
    return result


def compound_cell_ids(libraries: np.ndarray, source_ids: np.ndarray) -> np.ndarray:
    result = np.asarray([f"{library}|{cell}" for library, cell in zip(libraries, source_ids)], dtype=object)
    if len(set(result)) != len(result):
        raise ContractError("compound cell IDs are not unique")
    return result


def _hash_strings(values: Any) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _hash_sparse(matrix: Any) -> str:
    matrix = matrix.tocsr()
    digest = hashlib.sha256()
    digest.update(json.dumps(list(matrix.shape)).encode("ascii"))
    for name, values in (
        ("indptr", matrix.indptr), ("indices", matrix.indices), ("data", matrix.data)
    ):
        digest.update(name.encode("ascii"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(np.ascontiguousarray(values).view(np.uint8))
    return digest.hexdigest()


def _hash_prepared_obs(obs: Any) -> str:
    fields = [
        "library_id", "assay_id", "donor_id", "dataset", "preparation_method",
        "native_condition_authoritative",
        "harmonized_stage_authoritative", "technical_batch", "audit_cell_type",
        "strict_reference", "primary_query", "query_control", "analysis_eligible",
    ]
    digest = hashlib.sha256()
    digest.update(_hash_strings(obs.index).encode("ascii"))
    for field in fields:
        if field not in obs:
            raise ContractError(f"prepared metadata field is missing: {field}")
        digest.update(field.encode("utf-8"))
        digest.update(_hash_strings(obs[field].astype(str)).encode("ascii"))
    return digest.hexdigest()


def prepare_hvg_anndata(
    config: dict[str, Any],
    contract_lock: str | Path,
    library_manifest: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    """Materialize all canonical cells over reference-selected HVGs.

    This is an atlas-scale operation and must run on a compute node. It verifies
    the complete raw-count digest before extracting data.
    """
    import anndata as ad
    import h5py
    import pandas as pd
    import scanpy as sc
    from anndata.io import sparse_dataset

    lock = verify_contract_lock(config, contract_lock, full_hash=True)
    validate_program_firewall(config)
    atlas = Path(lock["atlas_realpath"])
    expected_library_hash = lock["manifest_sha256"]["library_manifest.tsv"]
    if sha256_path(library_manifest) != expected_library_hash:
        raise ContractError("library manifest does not match the production contract")
    manifest = read_library_manifest(library_manifest)
    # The canonical file contains legacy nullable ``uns`` metadata that is not
    # required here and is not readable by every pinned AnnData release. Read
    # only the two source obs columns needed to build the authoritative fields.
    source = h5py.File(atlas, "r")
    try:
        sample_categories, sample_codes = read_h5ad_column(source["obs"], "sample")
        type_categories, type_codes = read_h5ad_column(source["obs"], "cell_type")
        source_ids = source["obs/_index"].asstr()[:]
        obs = pd.DataFrame(
            {
                "sample": sample_categories[sample_codes].astype(str),
                "cell_type": type_categories[type_codes].astype(str),
            },
            index=pd.Index(source_ids.astype(str), name="source_cell_id"),
        )
        libraries = obs["sample"].astype(str).to_numpy()
        missing = sorted(set(libraries) - set(manifest))
        if missing:
            raise ContractError(f"prepared data has libraries absent from manifest: {missing}")
        mapped = [manifest[x] for x in libraries]
        obs["library_id"] = libraries
        obs["assay_id"] = [x["assay_id"] for x in mapped]
        obs["donor_id"] = [x["donor_id"] for x in mapped]
        obs["dataset"] = [x["dataset"] for x in mapped]
        obs["preparation_method"] = [x["preparation"] for x in mapped]
        obs["native_condition_authoritative"] = [x["native_condition"] for x in mapped]
        obs["harmonized_stage_authoritative"] = [x["harmonized_stage"] for x in mapped]
        obs["technical_batch"] = [x["technical_batch"] for x in mapped]
        for key in ("strict_reference", "primary_query", "query_control", "analysis_eligible"):
            obs[key] = np.asarray([x[key].lower() == "true" for x in mapped], dtype=bool)
        obs["audit_cell_type"] = obs["cell_type"].astype(str)
        obs.index = compound_cell_ids(libraries, source_ids.astype(str))

        reference_indices = np.flatnonzero(obs["strict_reference"].to_numpy())
        if len(reference_indices) != config["expected_contract"]["strict_reference"]["cells"]:
            raise ContractError("strict-reference cell count changed before HVG selection")
        with h5py.File(atlas, "r") as handle:
            raw = sparse_dataset(handle["raw/X"])
            reference_counts = raw[reference_indices, :].tocsr()
            genes = handle["raw/var/_index"].asstr()[:]
            donor_values = obs.iloc[reference_indices]["donor_id"].astype(str).to_numpy()
            unique_donors = sorted(set(donor_values))
            detected = np.zeros(reference_counts.shape[1], dtype=np.int16)
            for donor in unique_donors:
                detected += (reference_counts[donor_values == donor].getnnz(axis=0) > 0)
            eligible = detected >= config["features"]["min_reference_donors_detected"]
            if int(eligible.sum()) < config["features"]["n_hvg"]:
                raise ContractError("fewer reference-detected genes than requested HVGs")
            reference = ad.AnnData(
                X=reference_counts[:, eligible],
                obs=obs.iloc[reference_indices].copy(),
                var=pd.DataFrame(index=pd.Index(genes[eligible], name="gene")),
            )
            reference.layers["counts"] = reference.X
            sc.pp.highly_variable_genes(
                reference,
                layer="counts",
                flavor="seurat_v3",
                n_top_genes=config["features"]["n_hvg"],
                batch_key="donor_id",
                subset=False,
            )
            selected_names = reference.var_names[reference.var["highly_variable"]].astype(str)
            if len(selected_names) != config["features"]["n_hvg"]:
                raise ContractError("HVG selection did not return exactly the configured count")
            gene_index = {gene: i for i, gene in enumerate(genes)}
            selected_indices = np.asarray([gene_index[x] for x in selected_names], dtype=np.int64)
            counts = raw[:, selected_indices].tocsr()
        if np.any(counts.data < 0) or not np.all(counts.data == np.floor(counts.data)):
            raise ContractError("prepared model matrix is not non-negative integer counts")
        counts.data = counts.data.astype(np.int32, copy=False)
        prepared = ad.AnnData(
            X=counts,
            obs=obs,
            var=pd.DataFrame(index=pd.Index(selected_names, name="gene")),
        )
        prepared.layers["counts"] = prepared.X
        prepared.uns["masld_cl_contract"] = {
            "schema_version": lock["schema_version"],
            "config_sha256": config["_config_sha256"],
            "gene_order_sha256": lock["gene_order_sha256"],
            "raw_counts_sha256": lock["raw_counts_sha256"],
            "descriptive_cells": lock["observed_contract"]["descriptive"]["cells"],
            "model_gene_count": len(selected_names),
            "model_gene_order_sha256": _hash_strings(selected_names),
            "prepared_obs_sha256": _hash_prepared_obs(obs),
            "prepared_counts_sha256": _hash_sparse(counts),
        }
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            raise FileExistsError(output)
        prepared.write_h5ad(output, compression="gzip")
        prepared_file_hash = sha256_path(output)
        prepared_lock = {
            "schema_version": "masld-cl-prepared-v1",
            "config_sha256": config["_config_sha256"],
            "contract_lock_sha256": lock["lock_sha256"],
            "raw_counts_sha256": lock["raw_counts_sha256"],
            "prepared_realpath": str(output.resolve()),
            "prepared_size_bytes": output.stat().st_size,
            "prepared_file_sha256": prepared_file_hash,
            "model_gene_order_sha256": prepared.uns["masld_cl_contract"]["model_gene_order_sha256"],
            "prepared_obs_sha256": prepared.uns["masld_cl_contract"]["prepared_obs_sha256"],
            "prepared_counts_sha256": prepared.uns["masld_cl_contract"]["prepared_counts_sha256"],
        }
        prepared_lock["lock_sha256"] = hashlib.sha256(
            canonical_json_bytes(prepared_lock)
        ).hexdigest()
        prepared_lock_path = output.with_suffix(output.suffix + ".lock.json")
        write_json_exclusive(prepared_lock_path, prepared_lock)
        summary = {
            "output": str(output.resolve()),
            "n_cells": prepared.n_obs,
            "n_genes": prepared.n_vars,
            "n_reference_cells": len(reference_indices),
            "n_reference_donors": len(unique_donors),
            "config_sha256": config["_config_sha256"],
            "raw_counts_sha256": lock["raw_counts_sha256"],
            "prepared_file_sha256": prepared_file_hash,
            "prepared_lock": str(prepared_lock_path.resolve()),
        }
        with output.with_suffix(".summary.json").open("x") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return summary
    finally:
        source.close()


def verify_prepared_file(
    prepared_path: str | Path, prepared_lock_path: str | Path,
    config: dict[str, Any], contract_lock: dict[str, Any],
) -> dict[str, Any]:
    path = Path(prepared_path).resolve()
    with Path(prepared_lock_path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-prepared-v1":
        raise ContractError("unsupported prepared-data lock schema")
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if hashlib.sha256(canonical_json_bytes(payload)).hexdigest() != lock.get("lock_sha256"):
        raise ContractError("prepared-data lock content hash mismatch")
    expected = {
        "config_sha256": config["_config_sha256"],
        "contract_lock_sha256": contract_lock["lock_sha256"],
        "raw_counts_sha256": contract_lock["raw_counts_sha256"],
        "prepared_realpath": str(path),
        "prepared_size_bytes": path.stat().st_size,
    }
    for key, value in expected.items():
        if lock.get(key) != value:
            raise ContractError(f"prepared-data lock mismatch: {key}")
    if sha256_path(path) != lock.get("prepared_file_sha256"):
        raise ContractError("prepared H5AD content hash mismatch")
    return lock


def verify_prepared(
    adata: Any, config: dict[str, Any], contract_lock: dict[str, Any],
    prepared_lock: dict[str, Any],
) -> None:
    contract = adata.uns.get("masld_cl_contract", {})
    expected = contract_lock["observed_contract"]["descriptive"]["cells"]
    if adata.n_obs != expected or contract.get("descriptive_cells") != expected:
        raise ContractError("prepared object cell census does not match the contract")
    if adata.n_vars != config["features"]["n_hvg"]:
        raise ContractError("prepared object does not contain exactly the locked HVGs")
    if contract.get("config_sha256") != config["_config_sha256"]:
        raise ContractError("prepared object config hash mismatch")
    if contract.get("raw_counts_sha256") != contract_lock.get("raw_counts_sha256"):
        raise ContractError("prepared object raw-count hash mismatch")
    if contract.get("model_gene_order_sha256") != prepared_lock.get("model_gene_order_sha256"):
        raise ContractError("prepared model gene order changed")
    if _hash_strings(adata.var_names.astype(str)) != prepared_lock.get("model_gene_order_sha256"):
        raise ContractError("prepared H5AD genes are missing or permuted")
    if contract.get("prepared_obs_sha256") != prepared_lock.get("prepared_obs_sha256"):
        raise ContractError("prepared metadata hash mismatch")
    if contract.get("prepared_counts_sha256") != prepared_lock.get("prepared_counts_sha256"):
        raise ContractError("prepared count-matrix hash mismatch")
    if "counts" not in adata.layers:
        raise ContractError("prepared object is missing the integer count layer")
    count_data = adata.layers["counts"].data
    if np.any(~np.isfinite(count_data)) or np.any(count_data < 0) or not np.all(count_data == np.floor(count_data)):
        raise ContractError("prepared count layer is noninteger, negative, or non-finite")
    forbidden = set(config["input"]["forbidden_conditioning_fields"])
    if config["features"]["selection_batch_key"] in forbidden:
        raise ContractError("configured batch key is a forbidden biological field")
