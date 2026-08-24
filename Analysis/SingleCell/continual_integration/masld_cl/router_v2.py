"""Frozen CellTypist broad-lineage router, independent of integration coordinates."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .external_gse212837 import _load_csr, _read_frame_column
from .firewall import validate_program_firewall


POLICY_SCHEMA = "masld-cl-router-v2-policy-v43"


class RouterV2Error(ContractError):
    """Raised when the locked router or its mapping inputs differ."""


def load_router_v2_policy(
    config: dict[str, Any], value: str | Path, *, verify_diagnostic_hash: bool = True,
):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "router_v2_policy_v43.json"
    )
    if path != expected:
        raise RouterV2Error("Router V2 requires its source-controlled V43 policy")
    with path.open() as handle:
        policy = json.load(handle)
    expected_labels = {
        "B cells", "Basophils", "Cholangiocytes", "Circulating NK/NKT",
        "Endothelial cells", "Fibroblasts", "Hepatocytes", "Macrophages",
        "Mig.cDCs", "Mono+mono derived cells", "Neutrophils", "Plasma cells",
        "Resident NK", "T cells", "cDC1s", "cDC2s", "pDCs",
    }
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("selection_status", {}).get("gse212837_is_diagnostic_only") is not True
        or policy.get("selection_status", {}).get("gse212837_may_change_method_or_thresholds") is not False
        or policy.get("mapping", {}).get("query_author_labels_or_conditions_used") is not False
        or policy.get("mapping", {}).get("integration_coordinates_used") is not False
        or set(policy.get("collapsed_labels", {})) != expected_labels
        or policy.get("firewall", {}).get("external_labels_cannot_select_or_revise_the_router") is not True
    ):
        raise RouterV2Error("Router V2 policy identity or firewall differs")
    model = Path(policy["model"]["path"]).resolve()
    if sha256_path(model) != policy["model"]["sha256"]:
        raise RouterV2Error("frozen CellTypist model changed")
    diagnostic = (path.parent / policy["diagnostic_gse212837"]["h5ad"]["path"]).resolve()
    if (
        verify_diagnostic_hash
        and sha256_path(diagnostic) != policy["diagnostic_gse212837"]["h5ad"]["sha256"]
    ):
        raise RouterV2Error("GSE212837 diagnostic input changed")
    return path, policy, model, diagnostic


def collapse_router_labels(labels, mapping: dict[str, str]) -> np.ndarray:
    labels = np.asarray(labels).astype(str)
    unknown = sorted(set(labels) - set(mapping))
    if unknown:
        raise RouterV2Error(f"unknown CellTypist labels: {unknown}")
    return np.asarray([mapping[label] for label in labels], dtype=str)


def _validate_counts(counts) -> None:
    if counts.shape[0] == 0 or counts.shape[1] == 0:
        raise RouterV2Error("router count matrix is empty")
    if counts.nnz and (
        not np.isfinite(counts.data).all()
        or np.any(counts.data < 0)
        or not np.all(counts.data == np.floor(counts.data))
    ):
        raise RouterV2Error("router counts must be finite nonnegative integers")
    totals = np.asarray(counts.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise RouterV2Error("router input contains a zero-count cell")


def _route_counts(
    counts, gene_names, cell_ids, policy: dict[str, Any], model_path: Path,
) -> dict[str, Any]:
    import anndata as ad
    import celltypist
    import scanpy as sc
    from celltypist import models

    _validate_counts(counts)
    genes = np.asarray(gene_names).astype(str)
    if len(set(genes)) != len(genes):
        raise RouterV2Error("router gene symbols are duplicated")
    model = models.Model.load(model=str(model_path))
    if (
        celltypist.__version__ != policy["model"]["celltypist_version"]
        or len(model.features) != policy["model"]["features"]
        or len(model.cell_types) != policy["model"]["cell_types"]
    ):
        raise RouterV2Error("CellTypist runtime or frozen-model vocabulary differs")
    lookup = {gene: index for index, gene in enumerate(genes)}
    present = [gene for gene in model.features if gene in lookup]
    if len(present) < policy["mapping"]["minimum_model_features_present"]:
        raise RouterV2Error("too few frozen CellTypist model features are present")
    positions = np.asarray([lookup[gene] for gene in present], dtype=np.int64)
    router = ad.AnnData(
        X=counts[:, positions].tocsr().astype(np.float32),
        obs={"cell_id": np.asarray(cell_ids).astype(str)},
        var={"gene": present},
    )
    router.obs_names = router.obs["cell_id"].astype(str)
    router.var_names = np.asarray(present, dtype=str)
    if not router.obs_names.is_unique:
        raise RouterV2Error("router cell identifiers are duplicated")
    np.random.seed(int(policy["mapping"]["seed"]))
    sc.settings.n_jobs = 1
    sc.pp.normalize_total(router, target_sum=1e4)
    sc.pp.log1p(router)
    predictions = celltypist.annotate(
        router, model=model, mode=policy["mapping"]["prediction_mode"],
        majority_voting=True,
    )
    labels = predictions.predicted_labels
    required = {"predicted_labels", "majority_voting"}
    if not required.issubset(labels.columns) or len(labels) != router.n_obs:
        raise RouterV2Error("CellTypist did not return raw and majority-vote labels")
    raw = labels["predicted_labels"].astype(str).to_numpy()
    majority = labels["majority_voting"].astype(str).to_numpy()
    confidence = predictions.probability_matrix.max(axis=1).to_numpy(dtype=float)
    if not np.isfinite(confidence).all() or np.any((confidence < 0) | (confidence > 1)):
        raise RouterV2Error("CellTypist confidence values differ")
    return {
        "raw_labels": raw,
        "majority_labels": majority,
        "raw_collapsed": collapse_router_labels(raw, policy["collapsed_labels"]),
        "majority_collapsed": collapse_router_labels(majority, policy["collapsed_labels"]),
        "confidence": confidence,
        "model_features_present": len(present),
    }


def route_gse212837_diagnostic(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Route GSE212837 without opening author labels or conditions."""
    import h5py

    validate_program_firewall(config)
    policy_path, policy, model_path, source = load_router_v2_policy(config, policy_value)
    spec = policy["diagnostic_gse212837"]
    with h5py.File(source, "r") as handle:
        cell_ids = _read_frame_column(handle["obs"], "_index")
        donors = _read_frame_column(handle["obs"], "Donor")
        libraries = _read_frame_column(handle["obs"], "batch_name")
        genes = _read_frame_column(handle[spec["gene_group"]], "_index")
        counts = _load_csr(handle[spec["count_group"]])
    observed = {
        "cells": len(cell_ids), "donors": len(set(donors)),
        "libraries": len(set(libraries)),
    }
    if observed != {key: spec[key] for key in observed}:
        raise RouterV2Error(f"GSE212837 diagnostic census differs: {observed}")
    routed = _route_counts(counts, genes, cell_ids, policy, model_path)
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    cells_path = output / "routing_cells.tsv.gz"
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "row_index", "cell_id", "donor_id", "library_id",
            "raw_label", "majority_label", "raw_collapsed",
            "majority_collapsed", "confidence",
        ])
        for index, row in enumerate(zip(
            cell_ids, donors, libraries, routed["raw_labels"],
            routed["majority_labels"], routed["raw_collapsed"],
            routed["majority_collapsed"], routed["confidence"],
        )):
            writer.writerow([index, *row])
    result = {
        "schema_version": "masld-cl-router-v2-gse212837-diagnostic-v43",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "source": {"path": str(source), "sha256": sha256_path(source)},
        "model": {"path": str(model_path), "sha256": sha256_path(model_path)},
        "observed_census": observed,
        "model_features_present": routed["model_features_present"],
        "cells_file": cells_path.name,
        "cells_sha256": sha256_path(cells_path),
        "author_labels_or_conditions_opened": False,
        "integration_coordinates_opened_or_changed": False,
        "diagnostic_only_no_method_or_threshold_revision": True,
    }
    write_json_exclusive(output / "routing_manifest.json", result)
    return result
