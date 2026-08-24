"""Label-blind preparation and validation helpers for GSE296875."""

from __future__ import annotations

import hashlib
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, DeterministicGzipTextWriter, sha256_path


SOURCE_SCHEMA = "masld-cl-gse296875-source-lock-v44"
MAPPING_SCHEMA = "masld-cl-gse296875-mapping-policy-v45"


class GSE296875RouterError(ContractError):
    """Raised when an external-router input differs from its locked contract."""


def load_gse296875_source_lock(
    config: dict[str, Any], value: str | Path, *, verify_large_hashes: bool = False,
):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "gse296875_source_lock_v44.json"
    )
    if path != expected:
        raise GSE296875RouterError("GSE296875 requires its source-controlled V44 lock")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != SOURCE_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or set(policy.get("raw_rna_sources", {}))
        != {f"well{index}" for index in range(1, 9)}
        or policy.get("structural_inspection", {}).get("metadata_values_allowed") is not False
        or policy.get("mapping_contract", {}).get("processed_expression_used") is not False
        or policy.get("mapping_contract", {}).get("author_labels_or_azimuth_labels_used") is not False
        or policy.get("firewall", {}).get("source_and_method_locked_before_any_author_label_values") is not True
        or policy.get("firewall", {}).get("failure_is_reported_without_tuning_on_gse296875") is not True
    ):
        raise GSE296875RouterError("GSE296875 source lock identity or firewall differs")

    parent = (path.parent / policy["parent_router_policy"]["path"]).resolve()
    manifest = (path.parent / policy["raw_rna_manifest"]["path"]).resolve()
    script = (path.parent / policy["structural_inspection"]["script"]).resolve()
    for source, expected_hash, name in (
        (parent, policy["parent_router_policy"]["sha256"], "parent router policy"),
        (manifest, policy["raw_rna_manifest"]["sha256"], "raw RNA manifest"),
        (script, policy["structural_inspection"]["script_sha256"], "inspection script"),
    ):
        if sha256_path(source) != expected_hash:
            raise GSE296875RouterError(f"GSE296875 {name} changed")

    processed = (path.parent / policy["processed_source"]["path"]).resolve()
    raw_sources: dict[str, Path] = {}
    if processed.stat().st_size != int(policy["processed_source"]["bytes"]):
        raise GSE296875RouterError("GSE296875 processed source size changed")
    for well, spec in policy["raw_rna_sources"].items():
        source = (path.parent / spec["path"]).resolve()
        if source.stat().st_size != int(spec["bytes"]):
            raise GSE296875RouterError(f"GSE296875 raw source size changed: {well}")
        raw_sources[well] = source
    if verify_large_hashes:
        if sha256_path(processed) != policy["processed_source"]["sha256"]:
            raise GSE296875RouterError("GSE296875 processed source hash changed")
        for well, source in raw_sources.items():
            if sha256_path(source) != policy["raw_rna_sources"][well]["sha256"]:
                raise GSE296875RouterError(f"GSE296875 raw source hash changed: {well}")
    return path, policy, processed, raw_sources


def load_gse296875_mapping_policy(
    config: dict[str, Any], value: str | Path,
):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "gse296875_mapping_policy_v45.json"
    )
    if path != expected:
        raise GSE296875RouterError("GSE296875 requires its source-controlled V45 mapping policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != MAPPING_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("structure", {}).get("metadata_values_opened") is not False
        or policy.get("mapping_metadata", {}).get("all_other_metadata_values_forbidden_during_mapping") is not True
        or policy.get("router", {}).get("method_or_threshold_revision_from_gse296875") is not False
        or set(policy.get("author_label_evaluation", {}).get("expected_classes_and_aliases", {}))
        != {"B cells", "Cholangiocytes", "Hepatocytes", "Kupffer cells", "LSECs", "Mesenchymal cells", "NK-T cells"}
        or policy.get("firewall", {}).get("author_and_azimuth_label_values_not_opened_before_this_lock") is not True
        or policy.get("firewall", {}).get("failure_is_reported_without_tuning_on_gse296875") is not True
    ):
        raise GSE296875RouterError("GSE296875 mapping policy identity or firewall differs")
    source_lock_path = (path.parent / policy["source_lock"]["path"]).resolve()
    if sha256_path(source_lock_path) != policy["source_lock"]["sha256"]:
        raise GSE296875RouterError("GSE296875 source lock changed after structural inspection")
    source_lock = load_gse296875_source_lock(config, source_lock_path)
    structure = (path.parent / policy["structure"]["path"]).resolve()
    export_script = (path.parent / policy["mapping_metadata"]["export_script"]).resolve()
    router_policy = (path.parent / policy["router"]["policy_path"]).resolve()
    v40_policy = (path.parent / policy["v40_baseline"]["policy_path"]).resolve()
    checks = (
        (structure, policy["structure"]["sha256"], "structure"),
        (export_script, policy["mapping_metadata"]["export_script_sha256"], "metadata exporter"),
        (router_policy, policy["router"]["policy_sha256"], "router policy"),
        (v40_policy, policy["v40_baseline"]["policy_sha256"], "V40 policy"),
    )
    for source, expected_hash, name in checks:
        if sha256_path(source) != expected_hash:
            raise GSE296875RouterError(f"GSE296875 locked {name} changed")
    return path, policy, source_lock, router_policy, v40_policy


def prepare_gse296875_query_counts(
    config: dict[str, Any], mapping_policy_value: str | Path,
    metadata_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Prepare exact raw query counts without opening any annotation values."""
    import h5py
    import pandas as pd
    from celltypist import models
    from scipy import sparse

    from .external_gse212837 import (
        _read_frame_column,
        load_external_gse212837_policy,
    )
    from .firewall import validate_program_firewall
    from .router_v2 import _validate_counts, load_router_v2_policy

    validate_program_firewall(config)
    (
        mapping_policy_path, mapping_policy, source_lock_result,
        router_policy_path, v40_policy_path,
    ) = load_gse296875_mapping_policy(config, mapping_policy_value)
    source_lock_path, source_lock, _, raw_sources = source_lock_result
    _, router_policy, model_path, _ = load_router_v2_policy(
        config, router_policy_path, verify_diagnostic_hash=False,
    )
    _, _, v40_sources = load_external_gse212837_policy(
        config, v40_policy_path, verify_source_hashes=False,
    )

    metadata_path = Path(metadata_value).resolve()
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    expected_columns = ["cell_id", "donor_id", "well_id", "raw_barcode"]
    metadata_policy = mapping_policy["mapping_metadata"]
    if list(metadata.columns) != expected_columns:
        raise GSE296875RouterError("GSE296875 mapping metadata columns differ")
    if (
        len(metadata) != int(metadata_policy["expected_cells"])
        or metadata["cell_id"].duplicated().any()
        or metadata[["well_id", "raw_barcode"]].duplicated().any()
        or metadata.isna().any().any()
        or metadata["donor_id"].nunique() != int(metadata_policy["expected_biological_donors"])
        or set(metadata["well_id"]) != set(metadata_policy["expected_well_values"])
    ):
        raise GSE296875RouterError("GSE296875 mapping metadata census differs")

    model = models.Model.load(model=str(model_path))
    if (
        len(model.features) != router_policy["model"]["features"]
        or len(set(model.features)) != len(model.features)
    ):
        raise GSE296875RouterError("frozen router feature registry differs")
    with h5py.File(v40_sources["prepared_counts"], "r") as handle:
        canonical_genes = _read_frame_column(handle["var"], "gene")
    desired_genes = list(map(str, model.features))

    counts_by_well = []
    baseline_counts_by_well = []
    metadata_by_well = []
    present_genes = None
    baseline_genes = None
    well_diagnostics: dict[str, Any] = {}
    for well in metadata_policy["expected_well_values"]:
        cells = metadata.loc[metadata["well_id"] == well].copy()
        if cells.empty:
            raise GSE296875RouterError(f"GSE296875 mapping well is empty: {well}")
        prefix = f"{well}_"
        if not cells["raw_barcode"].str.startswith(prefix).all():
            raise GSE296875RouterError(f"GSE296875 source barcode prefix differs: {well}")
        cells["raw_barcode"] = cells["raw_barcode"].str[len(prefix):]
        if not cells["raw_barcode"].str.fullmatch(r"[ACGT]+-[0-9]+").all():
            raise GSE296875RouterError(f"GSE296875 raw barcode encoding differs: {well}")
        subset = load_10x_raw_subset(
            raw_sources[well], cells["raw_barcode"].to_numpy(), desired_genes,
            match_feature_ids=False,
        )
        baseline_subset = load_10x_canonical_subset(
            raw_sources[well], cells["raw_barcode"].to_numpy(), canonical_genes,
        )
        if (
            subset["rna_feature_registry_sha256"]
            != source_lock["raw_rna_registry"]["rna_feature_id_and_symbol_sha256"]
        ):
            raise GSE296875RouterError(f"GSE296875 RNA feature registry differs: {well}")
        if present_genes is None:
            present_genes = subset["genes"]
        elif not np.array_equal(present_genes, subset["genes"]):
            raise GSE296875RouterError("GSE296875 present-gene order differs across wells")
        if baseline_genes is None:
            baseline_genes = baseline_subset["genes"]
        elif not np.array_equal(baseline_genes, baseline_subset["genes"]):
            raise GSE296875RouterError("GSE296875 baseline-gene order differs across wells")
        counts_by_well.append(subset["counts"])
        baseline_counts_by_well.append(baseline_subset["counts"])
        metadata_by_well.append(cells)
        well_diagnostics[well] = {
            "cells": int(len(cells)),
            "source_barcodes": subset["source_barcodes"],
            "source_features": subset["source_features"],
            "source_rna_features": subset["source_rna_features"],
            "v40_common_genes": int(len(baseline_subset["genes"])),
            "v40_mapping_methods": baseline_subset["mapping_methods"],
            "v40_colliding_external_targets_excluded": baseline_subset["colliding_external_targets_excluded"],
        }
    assert present_genes is not None and baseline_genes is not None
    query_counts = sparse.vstack(counts_by_well, format="csr")
    baseline_counts = sparse.vstack(baseline_counts_by_well, format="csr")
    query_cells = pd.concat(metadata_by_well, ignore_index=True)
    _validate_counts(query_counts)
    _validate_counts(baseline_counts)
    present_gene_set = set(present_genes)
    model_features_present = int(sum(gene in present_gene_set for gene in model.features))
    if model_features_present != int(source_lock["raw_rna_registry"]["celltypist_model_features_present"]):
        raise GSE296875RouterError("GSE296875 frozen-model feature overlap differs")

    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    counts_path = output / "query_counts.npz"
    genes_path = output / "query_genes.npy"
    baseline_counts_path = output / "v40_query_counts.npz"
    baseline_genes_path = output / "v40_query_genes.npy"
    cells_path = output / "query_cells.tsv.gz"
    sparse.save_npz(counts_path, query_counts, compressed=True)
    np.save(genes_path, np.asarray(present_genes, dtype=str), allow_pickle=False)
    sparse.save_npz(baseline_counts_path, baseline_counts, compressed=True)
    np.save(baseline_genes_path, np.asarray(baseline_genes, dtype=str), allow_pickle=False)
    with DeterministicGzipTextWriter(cells_path) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(expected_columns)
        writer.writerows(query_cells[expected_columns].itertuples(index=False, name=None))
    result = {
        "schema_version": "masld-cl-gse296875-query-counts-v45",
        "config_sha256": config["_config_sha256"],
        "mapping_policy": {
            "path": str(mapping_policy_path), "sha256": sha256_path(mapping_policy_path),
        },
        "source_lock": {"path": str(source_lock_path), "sha256": sha256_path(source_lock_path)},
        "mapping_metadata": {"path": str(metadata_path), "sha256": sha256_path(metadata_path)},
        "cells": int(query_counts.shape[0]),
        "genes": int(query_counts.shape[1]),
        "v40_genes": int(baseline_counts.shape[1]),
        "donors": int(query_cells["donor_id"].nunique()),
        "wells": int(query_cells["well_id"].nunique()),
        "model_features_present": model_features_present,
        "well_diagnostics": well_diagnostics,
        "counts_file": counts_path.name,
        "counts_sha256": sha256_path(counts_path),
        "genes_file": genes_path.name,
        "genes_sha256": sha256_path(genes_path),
        "v40_counts_file": baseline_counts_path.name,
        "v40_counts_sha256": sha256_path(baseline_counts_path),
        "v40_genes_file": baseline_genes_path.name,
        "v40_genes_sha256": sha256_path(baseline_genes_path),
        "cells_file": cells_path.name,
        "cells_sha256": sha256_path(cells_path),
        "author_or_azimuth_labels_opened": False,
        "integration_coordinates_opened": False,
        "processed_expression_used": False,
    }
    write_json_exclusive(output / "query_counts_manifest.json", result)
    return result


def _read_strings(dataset) -> np.ndarray:
    if dataset.dtype.kind not in {"O", "S", "U"}:
        raise GSE296875RouterError("expected a string dataset in the 10x HDF5 file")
    return np.asarray(dataset.asstr()[:]).astype(str)


def load_10x_raw_subset(
    value: str | Path, desired_barcodes, desired_genes, *, match_feature_ids: bool = False,
):
    """Load selected raw barcodes and aggregate duplicate RNA gene symbols."""
    import h5py
    import pandas as pd
    from scipy import sparse

    path = Path(value).resolve()
    barcodes_requested = np.asarray(desired_barcodes).astype(str)
    genes_requested = np.asarray(desired_genes).astype(str)
    if (
        not len(barcodes_requested)
        or len(set(barcodes_requested)) != len(barcodes_requested)
        or not len(genes_requested)
        or len(set(genes_requested)) != len(genes_requested)
    ):
        raise GSE296875RouterError("requested barcodes and genes must be nonempty and unique")

    with h5py.File(path, "r") as handle:
        if "matrix" not in handle:
            raise GSE296875RouterError("10x HDF5 matrix group is absent")
        matrix_group = handle["matrix"]
        required = {"barcodes", "data", "indices", "indptr", "shape", "features"}
        if not required.issubset(matrix_group):
            raise GSE296875RouterError("10x HDF5 matrix fields differ")
        features_group = matrix_group["features"]
        if not {"id", "name", "feature_type"}.issubset(features_group):
            raise GSE296875RouterError("10x feature fields differ")

        source_barcodes = _read_strings(matrix_group["barcodes"])
        barcode_index = pd.Index(source_barcodes)
        if not barcode_index.is_unique:
            raise GSE296875RouterError("10x source barcodes are duplicated")
        positions = barcode_index.get_indexer(barcodes_requested)
        if np.any(positions < 0):
            missing = barcodes_requested[positions < 0]
            raise GSE296875RouterError(
                f"processed barcodes are missing from the raw matrix: {missing[:5].tolist()}"
            )

        feature_ids = _read_strings(features_group["id"])
        feature_names = _read_strings(features_group["name"])
        feature_types = _read_strings(features_group["feature_type"])
        shape = tuple(map(int, matrix_group["shape"][:]))
        if (
            len(shape) != 2
            or shape != (len(feature_names), len(source_barcodes))
            or len(feature_ids) != shape[0]
            or len(feature_types) != shape[0]
        ):
            raise GSE296875RouterError("10x matrix and feature dimensions differ")

        rna_rows = np.flatnonzero(feature_types == "Gene Expression")
        requested_lookup = {gene: index for index, gene in enumerate(genes_requested)}
        rows_by_gene: dict[str, list[int]] = {}
        for row in rna_rows:
            name = feature_names[row]
            identifier = feature_ids[row]
            name_match = name if name in requested_lookup else None
            id_match = identifier if match_feature_ids and identifier in requested_lookup else None
            if name_match is not None and id_match is not None and name_match != id_match:
                raise GSE296875RouterError("a raw feature ambiguously matches two requested genes")
            match = name_match if name_match is not None else id_match
            if match is not None:
                rows_by_gene.setdefault(match, []).append(int(row))
        present_genes = [gene for gene in genes_requested if gene in rows_by_gene]
        if not present_genes:
            raise GSE296875RouterError("no requested genes occur in the raw RNA matrix")

        selected_rows: list[int] = []
        aggregate_columns: list[int] = []
        for output_column, gene in enumerate(present_genes):
            rows = rows_by_gene[gene]
            selected_rows.extend(rows)
            aggregate_columns.extend([output_column] * len(rows))

        data = matrix_group["data"][:]
        indices = matrix_group["indices"][:]
        indptr = matrix_group["indptr"][:]
        if (
            np.any(data < 0)
            or not np.all(data == np.floor(data))
            or len(indptr) != shape[1] + 1
        ):
            raise GSE296875RouterError("10x raw counts are not nonnegative integers")
        full = sparse.csc_matrix((data, indices, indptr), shape=shape)
        selected = full[:, positions][selected_rows, :].T.tocsr()
        aggregation = sparse.csr_matrix(
            (
                np.ones(len(selected_rows), dtype=np.int8),
                (np.arange(len(selected_rows)), np.asarray(aggregate_columns)),
            ),
            shape=(len(selected_rows), len(present_genes)),
        )
        counts = (selected @ aggregation).tocsr()

    registry = hashlib.sha256()
    for row in rna_rows:
        registry.update(feature_ids[row].encode())
        registry.update(b"\0")
        registry.update(feature_names[row].encode())
        registry.update(b"\n")
    return {
        "counts": counts,
        "genes": np.asarray(present_genes, dtype=str),
        "rna_feature_registry_sha256": registry.hexdigest(),
        "source_barcodes": int(len(source_barcodes)),
        "source_features": int(shape[0]),
        "source_rna_features": int(len(rna_rows)),
    }


def load_10x_canonical_subset(
    value: str | Path, desired_barcodes, canonical_genes,
):
    """Load the exact one-to-one V40 canonical feature map from a raw 10x file."""
    import h5py
    import pandas as pd
    from scipy import sparse

    from .external_reference_common import common_gene_mapping

    path = Path(value).resolve()
    barcodes_requested = np.asarray(desired_barcodes).astype(str)
    genes_requested = np.asarray(canonical_genes).astype(str)
    if (
        not len(barcodes_requested)
        or len(set(barcodes_requested)) != len(barcodes_requested)
        or not len(genes_requested)
        or len(set(genes_requested)) != len(genes_requested)
    ):
        raise GSE296875RouterError("V40 requested barcodes and genes must be nonempty and unique")
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        source_barcodes = _read_strings(group["barcodes"])
        barcode_index = pd.Index(source_barcodes)
        if not barcode_index.is_unique:
            raise GSE296875RouterError("10x source barcodes are duplicated")
        positions = barcode_index.get_indexer(barcodes_requested)
        if np.any(positions < 0):
            raise GSE296875RouterError("processed barcodes are missing from the V40 raw matrix")
        feature_ids = _read_strings(group["features/id"])
        feature_names = _read_strings(group["features/name"])
        feature_types = _read_strings(group["features/feature_type"])
        rna_rows = np.flatnonzero(feature_types == "Gene Expression")
        mapping = common_gene_mapping(
            genes_requested.tolist(), feature_ids[rna_rows].tolist(), feature_names[rna_rows].tolist(),
        )
        raw_rows = rna_rows[mapping["external_indices"]]
        present_genes = genes_requested[mapping["canonical_indices"]]
        shape = tuple(map(int, group["shape"][:]))
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        if np.any(data < 0) or not np.all(data == np.floor(data)):
            raise GSE296875RouterError("10x V40 raw counts are not nonnegative integers")
        full = sparse.csc_matrix((data, indices, indptr), shape=shape)
        counts = full[:, positions][raw_rows, :].T.tocsr()
    return {
        "counts": counts,
        "genes": np.asarray(present_genes, dtype=str),
        "mapping_methods": mapping["mapping_methods"],
        "colliding_external_targets_excluded": int(mapping["colliding_external_targets_excluded"]),
    }


def build_v40_lda_baseline_routing(
    config: dict[str, Any], v40_policy_value: str | Path,
    query_counts, query_genes, query_donors,
):
    """Apply the frozen V40 reference-only LDA algorithm to a new raw query."""
    import gc

    import h5py
    import pandas as pd
    from scipy import sparse
    from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

    from .embedding import load_embedding
    from .external_gse212837 import (
        _balanced_positions,
        _load_csr,
        _normalize_log_common,
        _read_frame_column,
        _transform_sparse,
        load_external_gse212837_policy,
    )
    from .orthogonal_bridge import (
        apply_scaled_orthogonal_bridge,
        fit_scaled_orthogonal_bridge,
    )

    _, policy, sources = load_external_gse212837_policy(
        config, v40_policy_value, verify_source_hashes=False,
    )
    with sources["prepared_lock"].open() as handle:
        prepared_lock = json.load(handle)
    if prepared_lock.get("prepared_file_sha256") != policy["sources"]["prepared_counts"]["sha256"]:
        raise GSE296875RouterError("V40 prepared-count lock differs")
    reference_info, reference_target, reference_cells = load_embedding(
        sources["reference.all_lineage"]
    )
    if (
        reference_info.get("model_kind") != "all_lineage"
        or len(reference_cells) != 216957
        or not reference_cells["strict_reference"].to_numpy(dtype=bool).all()
    ):
        raise GSE296875RouterError("V40 strict-reference embedding differs")

    query_genes = np.asarray(query_genes).astype(str)
    query_donors = np.asarray(query_donors).astype(str)
    if (
        len(query_genes) != query_counts.shape[1]
        or len(set(query_genes)) != len(query_genes)
        or len(query_donors) != query_counts.shape[0]
        or len(set(query_donors)) < 3
    ):
        raise GSE296875RouterError("V40 query count metadata differs")
    query_lookup = {gene: index for index, gene in enumerate(query_genes)}

    with h5py.File(sources["prepared_counts"], "r") as prepared_handle:
        prepared_cell_ids = _read_frame_column(prepared_handle["obs"], "_index")
        prepared_genes = _read_frame_column(prepared_handle["var"], "gene")
        canonical_positions = np.asarray([
            index for index, gene in enumerate(prepared_genes) if gene in query_lookup
        ], dtype=np.int64)
        query_positions = np.asarray([
            query_lookup[prepared_genes[index]] for index in canonical_positions
        ], dtype=np.int64)
        minimum = int(policy["expected_external"]["minimum_common_genes"])
        if len(canonical_positions) < minimum:
            raise GSE296875RouterError(
                f"V40 baseline has only {len(canonical_positions)} common genes"
            )
        prepared_index = pd.Index(prepared_cell_ids)
        reference_rows = prepared_index.get_indexer(
            reference_cells["cell_id"].astype(str).to_numpy()
        )
        if np.any(reference_rows < 0) or len(set(reference_rows.tolist())) != len(reference_rows):
            raise GSE296875RouterError("V40 reference cells do not map uniquely")
        prepared_full = _load_csr(prepared_handle["X"])
        reference_counts = prepared_full[reference_rows][:, canonical_positions]
        del prepared_full
        gc.collect()

    query_common = query_counts[:, query_positions]
    reference_norm = _normalize_log_common(reference_counts)
    query_norm = _normalize_log_common(query_common)
    del reference_counts, query_common
    gc.collect()

    pca_policy = policy["mapping"]["pca"]
    seed = int(pca_policy["seed"])
    reference_fit = _balanced_positions(
        list(zip(
            reference_cells["donor_id"].astype(str),
            reference_cells["audit_cell_type"].astype(str),
        )),
        int(pca_policy["reference_cap_per_donor_label"]), seed,
    )
    query_fit = _balanced_positions(
        query_donors, int(pca_policy["query_cap_per_donor"]), seed + 1,
    )
    fit_matrix = sparse.vstack(
        [reference_norm[reference_fit], query_norm[query_fit]], format="csr"
    ).toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=int(pca_policy["dimensions"]),
        svd_solver=pca_policy["solver"],
        random_state=seed,
        iterated_power=int(pca_policy["iterated_power"]),
        whiten=False,
    )
    pca.fit(fit_matrix)
    del fit_matrix
    reference_source = _transform_sparse(pca, reference_norm)
    query_source = _transform_sparse(pca, query_norm)
    del reference_norm, query_norm
    gc.collect()

    bridge = fit_scaled_orthogonal_bridge(
        reference_source[reference_fit], np.asarray(reference_target)[reference_fit]
    )
    query_unadapted = apply_scaled_orthogonal_bridge(query_source, *bridge)
    lda = LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")
    lda.fit(
        np.asarray(reference_target)[reference_fit],
        reference_cells.iloc[reference_fit]["audit_cell_type"].astype(str).to_numpy(),
    )
    labels = lda.predict(query_unadapted).astype(str)
    return {
        "labels": labels,
        "common_genes": int(len(canonical_positions)),
        "reference_fit_cells": int(len(reference_fit)),
        "query_fit_cells": int(len(query_fit)),
        "query_donors": int(len(set(query_donors))),
        "algorithm": "V40 joint-common-gene PCA, strict-reference bridge, and reference-only shrinkage LDA",
    }


def route_gse296875_query(
    config: dict[str, Any], mapping_policy_value: str | Path,
    query_manifest_value: str | Path, output_value: str | Path,
) -> dict[str, Any]:
    """Run the locked Router V2 and V40 baseline before author labels are opened."""
    import gc

    import pandas as pd
    from scipy import sparse

    from .firewall import validate_program_firewall
    from .router_v2 import (
        _route_counts,
        collapse_router_labels,
        load_router_v2_policy,
    )

    validate_program_firewall(config)
    mapping_policy_path, mapping_policy, _, router_policy_path, v40_policy_path = (
        load_gse296875_mapping_policy(config, mapping_policy_value)
    )
    _, router_policy, model_path, _ = load_router_v2_policy(
        config, router_policy_path, verify_diagnostic_hash=False,
    )
    query_manifest_path = Path(query_manifest_value).resolve()
    with query_manifest_path.open() as handle:
        query_manifest = json.load(handle)
    if (
        query_manifest.get("schema_version") != "masld-cl-gse296875-query-counts-v45"
        or query_manifest.get("mapping_policy", {}).get("sha256") != sha256_path(mapping_policy_path)
        or query_manifest.get("author_or_azimuth_labels_opened") is not False
        or query_manifest.get("integration_coordinates_opened") is not False
        or query_manifest.get("processed_expression_used") is not False
    ):
        raise GSE296875RouterError("GSE296875 prepared-query identity or firewall differs")
    query_directory = query_manifest_path.parent
    counts_path = query_directory / query_manifest["counts_file"]
    genes_path = query_directory / query_manifest["genes_file"]
    baseline_counts_path = query_directory / query_manifest["v40_counts_file"]
    baseline_genes_path = query_directory / query_manifest["v40_genes_file"]
    cells_path = query_directory / query_manifest["cells_file"]
    for source, expected_hash, name in (
        (counts_path, query_manifest["counts_sha256"], "counts"),
        (genes_path, query_manifest["genes_sha256"], "genes"),
        (baseline_counts_path, query_manifest["v40_counts_sha256"], "V40 counts"),
        (baseline_genes_path, query_manifest["v40_genes_sha256"], "V40 genes"),
        (cells_path, query_manifest["cells_sha256"], "cells"),
    ):
        if sha256_path(source) != expected_hash:
            raise GSE296875RouterError(f"GSE296875 prepared-query {name} changed")
    counts = sparse.load_npz(counts_path).tocsr()
    genes = np.load(genes_path, allow_pickle=False).astype(str)
    baseline_counts = sparse.load_npz(baseline_counts_path).tocsr()
    baseline_genes = np.load(baseline_genes_path, allow_pickle=False).astype(str)
    cells = pd.read_csv(cells_path, sep="\t", compression="gzip", dtype=str)
    if (
        counts.shape != (len(cells), len(genes))
        or baseline_counts.shape != (len(cells), len(baseline_genes))
        or len(cells) != int(mapping_policy["mapping_metadata"]["expected_cells"])
        or cells["cell_id"].duplicated().any()
    ):
        raise GSE296875RouterError("GSE296875 prepared-query dimensions differ")

    routed = _route_counts(
        counts, genes, cells["cell_id"].to_numpy(), router_policy, model_path,
    )
    gc.collect()
    baseline = build_v40_lda_baseline_routing(
        config, v40_policy_path, baseline_counts, baseline_genes,
        cells["donor_id"].to_numpy(),
    )
    baseline_collapsed = collapse_router_labels(
        baseline["labels"], router_policy["collapsed_labels"],
    )

    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    cells_output = output / "routing_cells.tsv.gz"
    with DeterministicGzipTextWriter(cells_output) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "row_index", "cell_id", "donor_id", "well_id", "raw_barcode",
            "raw_label", "majority_label", "raw_collapsed", "majority_collapsed",
            "confidence", "v40_label", "v40_collapsed",
        ])
        for index, row in enumerate(zip(
            cells["cell_id"], cells["donor_id"], cells["well_id"], cells["raw_barcode"],
            routed["raw_labels"], routed["majority_labels"], routed["raw_collapsed"],
            routed["majority_collapsed"], routed["confidence"], baseline["labels"],
            baseline_collapsed,
        )):
            writer.writerow([index, *row])
    result = {
        "schema_version": "masld-cl-gse296875-routing-v45",
        "config_sha256": config["_config_sha256"],
        "mapping_policy": {
            "path": str(mapping_policy_path), "sha256": sha256_path(mapping_policy_path),
        },
        "query_counts_manifest": {
            "path": str(query_manifest_path), "sha256": sha256_path(query_manifest_path),
        },
        "cells": int(len(cells)),
        "donors": int(cells["donor_id"].nunique()),
        "wells": int(cells["well_id"].nunique()),
        "model_features_present": routed["model_features_present"],
        "v40_baseline": {key: value for key, value in baseline.items() if key != "labels"},
        "cells_file": cells_output.name,
        "cells_sha256": sha256_path(cells_output),
        "author_or_azimuth_labels_opened": False,
        "integration_coordinates_opened_or_changed": False,
        "case_stage_program_hero_gene_umap_or_cas13_opened": False,
        "mapping_complete": True,
    }
    manifest_path = output / "routing_manifest.json"
    write_json_exclusive(manifest_path, result)
    manifest_sha256 = sha256_path(manifest_path)
    token_path = output / "routing_complete.sha256"
    with token_path.open("x") as handle:
        handle.write(manifest_sha256 + "\n")
    return {**result, "manifest_sha256": manifest_sha256, "routing_token": str(token_path)}


def _donor_balanced_class_metrics(frame, prediction: str, policy: dict[str, Any]):
    validation = policy["validation"]
    rows = []
    for label in sorted(set(frame["truth"])):
        truth_positive = frame["truth"] == label
        positive_cells = int(truth_positive.sum())
        positive_donors = int(frame.loc[truth_positive, "donor_id"].nunique())
        donor_scores = []
        for _, group in frame.groupby("donor_id", sort=True):
            truth = group["truth"].to_numpy() == label
            predicted = group[prediction].to_numpy() == label
            tp = int(np.sum(truth & predicted))
            fp = int(np.sum(~truth & predicted))
            fn = int(np.sum(truth & ~predicted))
            denominator = 2 * tp + fp + fn
            if denominator:
                donor_scores.append(2 * tp / denominator)
        evaluable = bool(
            positive_cells >= int(validation["minimum_cells_for_class_evaluation"])
            and positive_donors >= int(validation["minimum_donors_for_class_evaluation"])
            and donor_scores
        )
        value = float(np.mean(donor_scores)) if donor_scores else 0.0
        rows.append({
            "label": str(label),
            "positive_cells": positive_cells,
            "positive_donors": positive_donors,
            "donors_in_score": int(len(donor_scores)),
            "donor_balanced_f1": value,
            "evaluable": evaluable,
            "pass": bool(
                not evaluable
                or value >= float(validation["minimum_evaluable_class_f1"])
            ),
        })
    return rows


def evaluate_gse296875_router(
    config: dict[str, Any], mapping_policy_value: str | Path,
    routing_manifest_value: str | Path, author_labels_value: str | Path,
    output_value: str | Path,
) -> dict[str, Any]:
    """Evaluate the untouched external labels only after mapping completion."""
    import pandas as pd

    from .metrics import macro_f1
    from .router_v2_evaluation import paired_donor_bootstrap

    mapping_policy_path, policy, _, _, _ = load_gse296875_mapping_policy(
        config, mapping_policy_value,
    )
    manifest_path = Path(routing_manifest_value).resolve()
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    manifest_hash = sha256_path(manifest_path)
    token_path = manifest_path.parent / "routing_complete.sha256"
    token = token_path.read_text().splitlines()
    if (
        manifest.get("schema_version") != "masld-cl-gse296875-routing-v45"
        or manifest.get("mapping_policy", {}).get("sha256") != sha256_path(mapping_policy_path)
        or manifest.get("author_or_azimuth_labels_opened") is not False
        or manifest.get("mapping_complete") is not True
        or token != [manifest_hash]
    ):
        raise GSE296875RouterError("GSE296875 routing completion identity differs")
    routing_cells_path = manifest_path.parent / manifest["cells_file"]
    if sha256_path(routing_cells_path) != manifest["cells_sha256"]:
        raise GSE296875RouterError("GSE296875 routing cells changed")
    routed = pd.read_csv(routing_cells_path, sep="\t", compression="gzip", dtype=str)
    labels_path = Path(author_labels_value).resolve()
    labels = pd.read_csv(labels_path, sep="\t", dtype=str)
    if (
        list(labels.columns) != ["cell_id", "author_label"]
        or len(labels) != int(policy["mapping_metadata"]["expected_cells"])
        or labels["cell_id"].duplicated().any()
        or labels.isna().any().any()
    ):
        raise GSE296875RouterError("GSE296875 author-label export differs")
    frame = routed.merge(labels, on="cell_id", validate="one_to_one")
    aliases = policy["author_label_evaluation"]["expected_classes_and_aliases"]
    reverse = {}
    for canonical, values in aliases.items():
        for value in values:
            if value in reverse and reverse[value] != canonical:
                raise GSE296875RouterError("GSE296875 author-label aliases overlap")
            reverse[value] = canonical
    frame["author_canonical"] = frame["author_label"].map(reverse)
    if frame["author_canonical"].isna().any():
        unknown = sorted(frame.loc[frame["author_canonical"].isna(), "author_label"].unique())
        raise GSE296875RouterError(f"unlocked GSE296875 author labels: {unknown}")
    frame["truth"] = frame["author_canonical"].map(
        policy["author_label_evaluation"]["author_to_router_broad"]
    )
    if frame["truth"].isna().any() or frame["truth"].nunique() != 7:
        raise GSE296875RouterError("GSE296875 canonical author classes differ")

    predictions = ["majority_collapsed", "raw_collapsed", "v40_collapsed"]
    minimum_classes = int(policy["validation"]["minimum_classes_per_donor"])
    donor_scores: dict[str, dict[str, float]] = {}
    for prediction in predictions:
        values = {}
        for donor, group in frame.groupby("donor_id", sort=True):
            if group["truth"].nunique() >= minimum_classes:
                values[str(donor)] = macro_f1(group["truth"], group[prediction])
        donor_scores[prediction] = values
    rosters = [set(value) for value in donor_scores.values()]
    if not rosters or any(roster != rosters[0] for roster in rosters[1:]) or len(rosters[0]) < 3:
        raise GSE296875RouterError("GSE296875 donor metric rosters differ")
    means = {
        prediction: float(np.mean(list(values.values())))
        for prediction, values in donor_scores.items()
    }
    bootstrap = paired_donor_bootstrap(
        donor_scores["majority_collapsed"], donor_scores["v40_collapsed"],
        int(policy["validation"]["bootstrap_replicates"]),
        int(policy["validation"]["bootstrap_seed"]),
    )
    class_metrics = _donor_balanced_class_metrics(frame, "majority_collapsed", policy)
    gate_pass = bool(
        means["majority_collapsed"]
        >= float(policy["validation"]["minimum_donor_balanced_macro_f1"])
        and bootstrap["ci_low"]
        > float(policy["validation"]["paired_improvement_lower_95_ci_vs_v40_must_exceed"])
        and all(row["pass"] for row in class_metrics)
    )
    result = {
        "schema_version": "masld-cl-gse296875-router-validation-v46",
        "config_sha256": config["_config_sha256"],
        "scientific_role": "untouched_external_validation_no_tuning",
        "mapping_policy": {
            "path": str(mapping_policy_path), "sha256": sha256_path(mapping_policy_path),
        },
        "routing_manifest": {"path": str(manifest_path), "sha256": manifest_hash},
        "author_labels": {"path": str(labels_path), "sha256": sha256_path(labels_path)},
        "cells": int(len(frame)),
        "donors_evaluated": int(len(rosters[0])),
        "author_classes": sorted(frame["author_canonical"].unique()),
        "donor_balanced_macro_f1": {
            "majority_voting": means["majority_collapsed"],
            "raw_prediction": means["raw_collapsed"],
            "v40_latent_lda": means["v40_collapsed"],
            "locked_minimum": policy["validation"]["minimum_donor_balanced_macro_f1"],
        },
        "paired_improvement_majority_vs_v40": bootstrap,
        "per_class": class_metrics,
        "external_gate_pass": gate_pass,
        "author_labels_opened_only_after_mapping_complete": True,
        "result_may_change_router_baseline_or_thresholds": False,
        "failure_consequence": "reject Router V2 and leave Harmony/V42 unchanged",
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "router_external_validation.json", result)
    with (output / "per_donor_metrics.tsv").open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["donor_id", "majority_macro_f1", "raw_macro_f1", "v40_macro_f1"])
        for donor in sorted(rosters[0]):
            writer.writerow([
                donor,
                donor_scores["majority_collapsed"][donor],
                donor_scores["raw_collapsed"][donor],
                donor_scores["v40_collapsed"][donor],
            ])
    return result
