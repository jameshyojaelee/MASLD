"""Condition-blind mapping of the independent GSE212837 liver cohort."""

from __future__ import annotations

import csv
import gc
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path
from .embedding import load_embedding
from .external_reference_common import common_gene_mapping
from .firewall import validate_program_firewall
from .orthogonal_bridge import (
    apply_scaled_orthogonal_bridge,
    fit_scaled_orthogonal_bridge,
)


POLICY_SCHEMA = "masld-cl-external-gse212837-policy-v40"


class ExternalGSE212837Error(ContractError):
    """Raised when the independent external mapping changes identity."""


def load_external_gse212837_policy(
    config: dict[str, Any], value: str | Path, *, verify_source_hashes: bool = True,
) -> tuple[Path, dict[str, Any], dict[str, Path]]:
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "external_gse212837_policy_v40.json"
    )
    if path != expected:
        raise ExternalGSE212837Error("GSE212837 requires its source-controlled V40 policy")
    with path.open() as handle:
        policy = json.load(handle)
    forbidden = set(policy.get("mapping", {}).get("forbidden_mapping_inputs", []))
    required_forbidden = {
        "celltype_pred", "Dataset Filter Name", "control_or_case", "disease",
        "stage", "fibrosis", "Sex", "Age", "Cas13",
    }
    required_models = {
        "all_lineage", "hepatocytes", "macrophages", "fibroblasts",
        "cholangiocytes", "t_cells",
    }
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("status", {}).get("cohort_used_during_method_development") is not False
        or policy.get("status", {}).get("expression_outcomes_read_before_policy_lock") is not False
        or forbidden != required_forbidden
        or set(policy.get("models", {})) != required_models
        or policy.get("firewall", {}).get("query_labels_and_conditions_unavailable_to_mapping") is not True
        or policy.get("firewall", {}).get("failure_is_reported_without_method_revision_on_this_cohort") is not True
    ):
        raise ExternalGSE212837Error("GSE212837 policy identity or firewall differs")
    sources: dict[str, Path] = {}
    for name, spec in policy["sources"].items():
        source = (path.parent / spec["path"]).resolve()
        if verify_source_hashes and sha256_path(source) != spec["sha256"]:
            raise ExternalGSE212837Error(f"GSE212837 source changed: {name}")
        sources[name] = source
    for model, spec in policy["models"].items():
        source = (path.parent / spec["reference_embedding"]).resolve()
        if verify_source_hashes and sha256_path(source) != spec["reference_embedding_sha256"]:
            raise ExternalGSE212837Error(
                f"GSE212837 reference embedding changed: {model}"
            )
        sources[f"reference.{model}"] = source
    return path, policy, sources


def _read_frame_column(group, name: str) -> np.ndarray:
    value = group[name]
    if hasattr(value, "shape"):
        result = value.asstr()[:] if value.dtype.kind in {"O", "S", "U"} else value[:]
    elif set(value) >= {"categories", "codes"}:
        categories = value["categories"]
        categories = (
            categories.asstr()[:]
            if categories.dtype.kind in {"O", "S", "U"}
            else categories[:].astype(str)
        )
        codes = value["codes"][:]
        if np.any(codes < 0):
            raise ExternalGSE212837Error(f"missing categorical values in {name}")
        result = categories[codes]
    else:
        raise ExternalGSE212837Error(f"unsupported H5AD column encoding: {name}")
    return np.asarray(result).astype(str)


def _load_csr(group):
    from scipy.sparse import csr_matrix

    shape = tuple(map(int, group.attrs["shape"]))
    result = csr_matrix(
        (group["data"][:], group["indices"][:], group["indptr"][:]),
        shape=shape,
    )
    if result.shape != shape:
        raise ExternalGSE212837Error("sparse count matrix shape differs")
    return result


def _normalize_log_common(counts, target_sum: float = 1e4):
    counts = counts.tocsr().astype(np.float32, copy=True)
    if counts.nnz and (
        np.any(counts.data < 0)
        or not np.isfinite(counts.data).all()
        or not np.all(counts.data == np.floor(counts.data))
    ):
        raise ExternalGSE212837Error("common-universe counts are not finite integers")
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    if np.any(totals <= 0):
        raise ExternalGSE212837Error("a cell has zero common-universe library size")
    scales = (float(target_sum) / totals).astype(np.float32)
    counts.data *= np.repeat(scales, np.diff(counts.indptr))
    np.log1p(counts.data, out=counts.data)
    return counts


def _balanced_positions(keys, cap: int, seed: int) -> np.ndarray:
    keys = np.asarray([
        str(key) if isinstance(key, str) else "\x1f".join(map(str, key))
        for key in keys
    ], dtype=str)
    rng = np.random.default_rng(int(seed))
    keep: list[int] = []
    for key in sorted(set(keys.tolist())):
        positions = np.flatnonzero(keys == key)
        if len(positions) > cap:
            positions = np.sort(rng.choice(positions, cap, replace=False))
        keep.extend(map(int, positions))
    result = np.asarray(sorted(keep), dtype=np.int64)
    if not len(result):
        raise ExternalGSE212837Error("balanced fit roster is empty")
    return result


def _transform_sparse(pca, matrix, chunk_size: int = 4096) -> np.ndarray:
    output = np.empty((matrix.shape[0], int(pca.n_components_)), dtype=np.float32)
    for start in range(0, matrix.shape[0], chunk_size):
        end = min(start + chunk_size, matrix.shape[0])
        output[start:end] = pca.transform(matrix[start:end].toarray()).astype(np.float32)
    if not np.isfinite(output).all():
        raise ExternalGSE212837Error("PCA projection contains non-finite values")
    return output


def _donor_balanced_centroid(values: np.ndarray, donors: np.ndarray):
    donors = np.asarray(donors).astype(str)
    unique = sorted(set(donors))
    if not unique:
        raise ExternalGSE212837Error("donor-balanced centroid has no donors")
    centroids = np.vstack([np.asarray(values)[donors == donor].mean(axis=0) for donor in unique])
    return centroids.mean(axis=0), len(unique)


def _donor_centroids(values: np.ndarray, donors: np.ndarray):
    donors = np.asarray(donors).astype(str)
    unique = np.asarray(sorted(set(donors)), dtype=str)
    centroids = np.vstack([np.asarray(values)[donors == donor].mean(axis=0) for donor in unique])
    return unique, centroids


def _write_cells(path: Path, cell_ids, donors, libraries, routing):
    with DeterministicGzipTextWriter(path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["row_index", "cell_id", "donor_id", "library_id", "routing_label"])
        for index, row in enumerate(zip(cell_ids, donors, libraries, routing)):
            writer.writerow([index, *row])


def build_external_gse212837_mapping(
    config: dict[str, Any], policy_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Build a query-label-free mapping and do not read evaluation-only fields."""
    import h5py
    import pandas as pd
    from scipy import sparse
    from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    validate_program_firewall(config)
    policy_path, policy, sources = load_external_gse212837_policy(config, policy_value)
    expected = policy["expected_external"]
    mapping_policy = policy["mapping"]
    pca_policy = mapping_policy["pca"]

    with sources["prepared_lock"].open() as handle:
        prepared_lock = json.load(handle)
    if (
        prepared_lock.get("prepared_file_sha256") != policy["sources"]["prepared_counts"]["sha256"]
        or prepared_lock.get("model_gene_order_sha256") != "000cdaaf8f2fc85d381f127aba2b694dc86f87ac9d2d04cfb154099ffb1a439e"
    ):
        raise ExternalGSE212837Error("prepared-count lock differs")

    reference_info, reference_target, reference_cells = load_embedding(
        sources["reference.all_lineage"]
    )
    if (
        reference_info.get("model_kind") != "all_lineage"
        or len(reference_cells) != 216957
        or not reference_cells["strict_reference"].to_numpy(dtype=bool).all()
    ):
        raise ExternalGSE212837Error("strict-reference all-lineage embedding differs")

    # Read only identifiers needed for mapping. Author labels and conditions are not opened.
    with h5py.File(sources["external_h5ad"], "r") as external_handle:
        external_cell_ids = _read_frame_column(external_handle["obs"], "_index")
        external_donors = _read_frame_column(external_handle["obs"], "Donor")
        external_libraries = _read_frame_column(external_handle["obs"], "batch_name")
        external_gene_ids = _read_frame_column(external_handle["raw/var"], "gene_ids")
        external_gene_names = _read_frame_column(external_handle["raw/var"], "_index")
        observed_external = {
            "cells": len(external_cell_ids),
            "biological_donors": len(set(external_donors)),
            "libraries": len(set(external_libraries)),
            "raw_genes": len(external_gene_ids),
        }
        if observed_external != {key: expected[key] for key in observed_external}:
            raise ExternalGSE212837Error(
                f"GSE212837 census differs: {observed_external}"
            )

        with h5py.File(sources["prepared_counts"], "r") as prepared_handle:
            prepared_cell_ids = _read_frame_column(prepared_handle["obs"], "_index")
            prepared_genes = _read_frame_column(prepared_handle["var"], "gene")
            gene_map = common_gene_mapping(
                prepared_genes.tolist(), external_gene_ids.tolist(),
                external_gene_names.tolist(),
            )
            if (
                gene_map["mapped_genes"] != expected["mapped_common_genes"]
                or gene_map["mapped_genes"] < expected["minimum_common_genes"]
            ):
                raise ExternalGSE212837Error("GSE212837 common-gene universe differs")
            prepared_index = pd.Index(prepared_cell_ids)
            reference_rows = prepared_index.get_indexer(
                reference_cells["cell_id"].astype(str).to_numpy()
            )
            if np.any(reference_rows < 0) or len(set(reference_rows.tolist())) != len(reference_rows):
                raise ExternalGSE212837Error("reference cells do not map uniquely to prepared counts")
            prepared_full = _load_csr(prepared_handle["X"])
            reference_counts = prepared_full[reference_rows][:, gene_map["canonical_indices"]]
            del prepared_full
            gc.collect()

        external_full = _load_csr(external_handle["raw/X"])
        external_counts = external_full[:, gene_map["external_indices"]]
        del external_full
        gc.collect()

    reference_norm = _normalize_log_common(reference_counts)
    external_norm = _normalize_log_common(external_counts)
    del reference_counts, external_counts
    gc.collect()

    seed = int(pca_policy["seed"])
    reference_fit = _balanced_positions(
        list(zip(
            reference_cells["donor_id"].astype(str),
            reference_cells["audit_cell_type"].astype(str),
        )),
        int(pca_policy["reference_cap_per_donor_label"]), seed,
    )
    query_fit = _balanced_positions(
        external_donors, int(pca_policy["query_cap_per_donor"]), seed + 1,
    )
    fit_matrix = sparse.vstack(
        [reference_norm[reference_fit], external_norm[query_fit]], format="csr"
    ).toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=int(pca_policy["dimensions"]), svd_solver="randomized",
        random_state=seed, iterated_power=int(pca_policy["iterated_power"]),
        whiten=False,
    )
    pca.fit(fit_matrix)
    del fit_matrix
    reference_source = _transform_sparse(pca, reference_norm)
    external_source = _transform_sparse(pca, external_norm)

    # Independent unintegrated comparator: a query-only, condition-blind PCA.
    query_pca = PCA(
        n_components=int(pca_policy["dimensions"]), svd_solver="randomized",
        random_state=seed + 2, iterated_power=int(pca_policy["iterated_power"]),
        whiten=False,
    )
    query_pca.fit(external_norm[query_fit].toarray().astype(np.float32, copy=False))
    external_raw_pca = _transform_sparse(query_pca, external_norm)
    del reference_norm, external_norm
    gc.collect()

    source_mean, target_mean, rotation, scale = fit_scaled_orthogonal_bridge(
        reference_source[reference_fit], np.asarray(reference_target)[reference_fit]
    )
    external_all_unadapted = apply_scaled_orthogonal_bridge(
        external_source, source_mean, target_mean, rotation, scale
    )

    # Routing is trained only on frozen strict-reference labels, before translation.
    lda = LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")
    lda.fit(
        np.asarray(reference_target)[reference_fit],
        reference_cells.iloc[reference_fit]["audit_cell_type"].astype(str).to_numpy(),
    )
    routing_labels = lda.predict(external_all_unadapted).astype(str)

    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        output / "projection_basis.npz",
        common_canonical_indices=gene_map["canonical_indices"],
        common_external_indices=gene_map["external_indices"],
        pca_mean=pca.mean_, pca_components=pca.components_,
        pca_explained_variance=pca.explained_variance_,
        query_pca_mean=query_pca.mean_, query_pca_components=query_pca.components_,
        reference_fit_positions=reference_fit, query_fit_positions=query_fit,
    )

    all_reference_lookup = pd.Index(reference_cells["cell_id"].astype(str))
    model_outputs: dict[str, Any] = {}
    for model_kind, spec in policy["models"].items():
        ref_info, ref_target, ref_cells = load_embedding(sources[f"reference.{model_kind}"])
        expected_reference_kind = (
            "all_lineage" if spec["lineage_label"] is None else spec["lineage_label"]
        )
        if ref_info.get("model_kind") != expected_reference_kind:
            raise ExternalGSE212837Error(f"reference model kind differs: {model_kind}")
        source_rows = all_reference_lookup.get_indexer(ref_cells["cell_id"].astype(str))
        if np.any(source_rows < 0) or len(set(source_rows.tolist())) != len(source_rows):
            raise ExternalGSE212837Error(f"reference source roster differs: {model_kind}")
        if spec["lineage_label"] is None:
            external_rows = np.arange(len(external_cell_ids), dtype=np.int64)
            model_reference_source = reference_source
        else:
            external_rows = np.flatnonzero(routing_labels == spec["lineage_label"])
            model_reference_source = reference_source[source_rows]
        if not len(external_rows):
            raise ExternalGSE212837Error(f"reference-only routing produced no {model_kind} cells")
        if model_kind == "all_lineage":
            fit_rows = reference_fit
            unadapted = external_all_unadapted
        else:
            cap = int(config["sampling"]["lineage_cap_per_donor"])
            fit_rows = _balanced_positions(
                ref_cells["donor_id"].astype(str).to_numpy(), cap, seed,
            )
            bridge = fit_scaled_orthogonal_bridge(
                model_reference_source[fit_rows], np.asarray(ref_target)[fit_rows]
            )
            unadapted = apply_scaled_orthogonal_bridge(
                external_source[external_rows], *bridge
            )
        query_centroid, n_query_donors = _donor_balanced_centroid(
            unadapted, external_donors[external_rows]
        )
        reference_centroid, n_reference_donors = _donor_balanced_centroid(
            np.asarray(ref_target), ref_cells["donor_id"].astype(str).to_numpy()
        )
        if min(n_query_donors, n_reference_donors) < int(mapping_policy["minimum_donors"]):
            raise ExternalGSE212837Error(f"underpowered stage-blind offset: {model_kind}")
        offset = np.asarray(reference_centroid - query_centroid, dtype=np.float64)
        adapted = (np.asarray(unadapted, dtype=np.float64) + offset).astype(np.float32)
        centered_error = float(np.max(np.abs(
            (adapted.astype(np.float64) - adapted.mean(axis=0, dtype=np.float64))
            - (np.asarray(unadapted, dtype=np.float64)
               - np.asarray(unadapted, dtype=np.float64).mean(axis=0))
        ), initial=0.0))
        model_dir = output / model_kind
        model_dir.mkdir()
        np.savez_compressed(
            model_dir / "external_coordinates.npz",
            unadapted=np.asarray(unadapted, dtype=np.float32),
            adapted=adapted,
            raw_study_pca=external_raw_pca[external_rows],
            offset=offset,
        )
        _write_cells(
            model_dir / "external_cells.tsv.gz",
            external_cell_ids[external_rows], external_donors[external_rows],
            external_libraries[external_rows], routing_labels[external_rows],
        )
        model_manifest = {
            "schema_version": "masld-cl-external-gse212837-model-v40",
            "config_sha256": config["_config_sha256"],
            "model_kind": model_kind,
            "lineage_label": spec["lineage_label"],
            "n_cells": int(len(external_rows)),
            "n_query_donors": n_query_donors,
            "n_reference_donors": n_reference_donors,
            "reference_embedding": {
                "path": str(sources[f"reference.{model_kind}"]),
                "sha256": sha256_path(sources[f"reference.{model_kind}"]),
            },
            "coordinates_file": "external_coordinates.npz",
            "coordinates_sha256": sha256_path(model_dir / "external_coordinates.npz"),
            "cells_file": "external_cells.tsv.gz",
            "cells_sha256": sha256_path(model_dir / "external_cells.tsv.gz"),
            "translation_offset": list(map(float, offset)),
            "maximum_centered_coordinate_error": centered_error,
            "mapping_inputs": ["raw/X", "Donor", "batch_name", "cell_id"],
            "query_labels_or_conditions_read": False,
            "author_labels_read": False,
            "reference_coordinates_written_or_changed": False,
            "independently_evaluable": spec["independently_evaluable"],
        }
        if not spec["independently_evaluable"]:
            model_manifest["untestable_reason"] = spec["untestable_reason"]
        write_json_exclusive(model_dir / "mapping_manifest.json", model_manifest)
        model_outputs[model_kind] = {
            "manifest": str((model_dir / "mapping_manifest.json").resolve()),
            "manifest_sha256": sha256_path(model_dir / "mapping_manifest.json"),
        }

    routing_counts = {
        str(label): int(count)
        for label, count in zip(*np.unique(routing_labels, return_counts=True))
    }
    basis_path = output / "projection_basis.npz"
    result = {
        "schema_version": "masld-cl-external-gse212837-bundle-v40",
        "config_sha256": config["_config_sha256"],
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "external_source": {
            "path": str(sources["external_h5ad"]),
            "sha256": policy["sources"]["external_h5ad"]["sha256"],
        },
        "cells": len(external_cell_ids),
        "donors": len(set(external_donors)),
        "libraries": len(set(external_libraries)),
        "common_genes": gene_map["mapped_genes"],
        "common_mapping_methods": gene_map["mapping_methods"],
        "joint_pca_reference_fit_cells": int(len(reference_fit)),
        "joint_pca_query_fit_cells": int(len(query_fit)),
        "routing_counts": routing_counts,
        "projection_basis": {"path": str(basis_path), "sha256": sha256_path(basis_path)},
        "models": model_outputs,
        "query_labels_or_conditions_read": False,
        "evaluation_locked_until_bundle_complete": True,
    }
    write_json_exclusive(output / "bundle_manifest.json", result)
    return result
