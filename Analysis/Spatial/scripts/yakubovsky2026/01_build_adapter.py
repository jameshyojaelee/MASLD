#!/usr/bin/env python3
"""Build the source-native Yakubovsky human Visium adapter after gate freeze.

The deposited matrix is streamed from MATLAB v7.3 HDF5 in physical-chunk-
aligned tiles.  It is
the source authors' nonnegative, mean-background-subtracted abundance—not raw
UMI counts—and is preserved without rounding.  The source-defined binary lipid
label is preserved exactly as deposited; no image classification, continuous
lipid field, or threshold reconstruction is used.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from yakubovsky_common import (
    ADAPTER_STAGE_ARTIFACTS,
    ALLOWED_LIPID_CLASSES,
    DATASET,
    MIN_DONOR_SPOTS,
    RELEASE_ID,
    ContractError,
    atomic_write_frame,
    bulk_decode_matlab_cellstr,
    bool_value,
    decode_matlab_char,
    default_paths,
    gene_axis_sha256,
    numeric_vector,
    sha256_file,
    validate_gene_axis_manifest,
    validate_registry_contract,
    validate_source_gate_summary,
    validate_stage_chain,
    write_stage_seal,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def h5ad_lipid_zone_value(lipid_class: str, inference_eligible: bool) -> str:
    """Return a serialization-safe three-state observation value."""
    if not inference_eligible:
        return "not_applicable"
    if lipid_class not in ALLOWED_LIPID_CLASSES:
        raise ContractError(f"Invalid source-defined binary lipid class: {lipid_class}")
    return "true" if lipid_class == "lipid_zone" else "false"


def source_spot_analysis_eligible(
    structural_join: bool,
    background_corrected_spot_sum: float,
) -> bool:
    """Retain reference-only spots but prohibit them from any inference."""

    return bool(
        structural_join
        and np.isfinite(background_corrected_spot_sum)
        and background_corrected_spot_sum > 0
    )


def source_spot_binary_lipid_state(
    donor_passes_gate: bool,
    structural_join: bool,
    authoritative_ordinal_join: bool,
    deposited_class: str,
) -> tuple[bool, str]:
    """Apply the donor gate at spot resolution without fabricating labels."""

    if not donor_passes_gate or not structural_join or not authoritative_ordinal_join:
        return False, "not_applicable"
    if deposited_class not in ALLOWED_LIPID_CLASSES:
        raise ContractError(f"Invalid source-defined binary lipid class: {deposited_class}")
    return True, deposited_class


def freeze_file_map(path: Path) -> dict[str, dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"record_type", "name", "value", "sha256", "lipid_outcome_read"}
    if not required.issubset(frame.columns):
        raise ContractError(f"Freeze manifest lacks {sorted(required)}")
    files = frame[frame["record_type"] == "file_input"]
    if files["name"].duplicated().any():
        raise ContractError("Freeze manifest contains duplicate file-input roles")
    if any(bool_value(value) for value in files["lipid_outcome_read"]):
        raise ContractError("Freeze manifest does not certify outcome-blind freezing")
    return files.set_index("name").to_dict(orient="index")


def verify_frozen_files(frozen: dict[str, dict[str, str]]) -> None:
    for role, record in frozen.items():
        path = Path(record["value"])
        if not path.is_file():
            raise ContractError(f"Frozen input disappeared ({role}): {path}")
        observed = sha256_file(path)
        if observed != record["sha256"]:
            raise ContractError(
                f"Frozen input drift ({role}): observed {observed}, expected {record['sha256']}"
            )


def map_source_genes(
    donor_gene_lists: dict[str, list[str]], gencode_path: Path
) -> tuple[dict[str, str | None], pd.DataFrame, list[str], dict[str, list[str]]]:
    gencode = pd.read_csv(gencode_path, sep="\t", dtype=str)
    base_to_symbols = (
        gencode.groupby("ensembl_base")["gene_name"].agg(lambda values: sorted(set(values)))
    )
    symbol_to_bases = (
        gencode.groupby("gene_name")["ensembl_base"].agg(lambda values: sorted(set(values)))
    )
    appearances: Counter[str] = Counter()
    source_order: list[str] = []
    for genes in donor_gene_lists.values():
        for gene in genes:
            appearances[gene] += 1
            if gene not in appearances or appearances[gene] == 1:
                source_order.append(gene)

    mapping: dict[str, str | None] = {}
    rows: list[dict[str, Any]] = []
    symbol_sources: dict[str, list[str]] = defaultdict(list)
    for source_gene in dict.fromkeys(source_order):
        stripped = source_gene.strip()
        canonical: str | None = None
        status = ""
        if not stripped:
            status = "excluded_empty_identifier"
        elif stripped.startswith("ENSG"):
            base = stripped.split(".", 1)[0]
            candidates = base_to_symbols.get(base, [])
            if isinstance(candidates, float) and np.isnan(candidates):
                candidates = []
            if len(candidates) == 1:
                canonical = candidates[0]
                status = "ensembl_version_stripped_unique_gencode_v49"
            elif len(candidates) == 0:
                status = "excluded_ensembl_not_in_gencode_v49"
            else:
                status = "excluded_ambiguous_ensembl_to_symbol"
        else:
            candidates = symbol_to_bases.get(stripped, [])
            if isinstance(candidates, float) and np.isnan(candidates):
                candidates = []
            if len(candidates) == 1:
                canonical = stripped
                status = "source_symbol_unique_gencode_v49"
            elif len(candidates) == 0:
                status = "excluded_symbol_not_in_gencode_v49"
            else:
                status = "excluded_ambiguous_symbol_in_gencode_v49"
        mapping[source_gene] = canonical
        if canonical is not None:
            symbol_sources[canonical].append(source_gene)
        rows.append(
            {
                "source_gene": source_gene,
                "canonical_symbol": canonical or "",
                "mapping_status": status,
                "n_donor_gene_lists": appearances[source_gene],
                "ensembl_version_stripped": source_gene.startswith("ENSG") and "." in source_gene,
            }
        )
    canonical_order: list[str] = []
    for source_gene in source_order:
        canonical = mapping[source_gene]
        if canonical is not None and canonical not in canonical_order:
            canonical_order.append(canonical)
    return mapping, pd.DataFrame(rows), canonical_order, symbol_sources


def stream_donor_abundance(
    matrix: h5py.Dataset,
    n_spots: int,
    local_genes: list[str],
    source_mapping: dict[str, str | None],
    global_index: dict[str, int],
    zonation_landmark_local_indices: np.ndarray,
    dense_tile_mib: float,
    orientation_hint: str | None = None,
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    plan = build_matrix_read_plan(
        matrix,
        n_spots,
        len(local_genes),
        dense_tile_mib=dense_tile_mib,
        orientation_hint=orientation_hint,
    )
    local_to_global = np.asarray(
        [global_index.get(source_mapping[gene] or "", -1) for gene in local_genes],
        dtype=int,
    )
    row_blocks: list[sparse.csr_matrix] = []
    zonation_landmark_abundance = np.zeros(n_spots, dtype=float)
    source_spot_sum = np.zeros(n_spots, dtype=float)
    source_detected = np.zeros(n_spots, dtype=np.int64)
    tiles_by_spot_slab: dict[tuple[int, int], list[tuple[int, int, int, int]]] = (
        defaultdict(list)
    )
    for tile in plan.tiles:
        tiles_by_spot_slab[(tile[0], tile[1])].append(tile)
    for (spot_start, spot_stop), tiles in tiles_by_spot_slab.items():
        n_slab_spots = spot_stop - spot_start
        # Long-double partial accumulation makes the result insensitive to the
        # bounded gene-slab grouping to substantially tighter than float64.
        sum_accumulator = np.zeros(n_slab_spots, dtype=np.longdouble)
        landmark_accumulator = np.zeros(n_slab_spots, dtype=np.longdouble)
        detected_accumulator = np.zeros(n_slab_spots, dtype=np.int64)
        coo_rows: list[np.ndarray] = []
        coo_columns: list[np.ndarray] = []
        coo_values: list[np.ndarray] = []
        for tile in tiles:
            _, _, gene_start, gene_stop = tile
            block = read_logical_tile(matrix, plan.orientation, tile)
            expected_shape = (n_slab_spots, gene_stop - gene_start)
            if block.shape != expected_shape:
                raise ContractError(
                    f"Unexpected chunk-aligned matrix tile {block.shape}; expected {expected_shape}"
                )
            if not np.isfinite(block).all() or (block < 0).any():
                raise ContractError(
                    "Source background-corrected abundance has non-finite or negative values"
                )
            sum_accumulator += np.sum(block, axis=1, dtype=np.longdouble)
            detected_accumulator += np.count_nonzero(block, axis=1)
            in_tile_landmarks = zonation_landmark_local_indices[
                (zonation_landmark_local_indices >= gene_start)
                & (zonation_landmark_local_indices < gene_stop)
            ]
            if len(in_tile_landmarks):
                landmark_accumulator += np.sum(
                    block[:, in_tile_landmarks - gene_start],
                    axis=1,
                    dtype=np.longdouble,
                )
            mapped_positions = np.flatnonzero(local_to_global[gene_start:gene_stop] >= 0)
            if len(mapped_positions):
                retained = block[:, mapped_positions]
                row, column = np.nonzero(retained)
                if len(row):
                    coo_rows.append(row.astype(np.int64, copy=False))
                    coo_columns.append(
                        local_to_global[gene_start + mapped_positions[column]].astype(
                            np.int64, copy=False
                        )
                    )
                    coo_values.append(retained[row, column])
        if coo_values:
            sparse_block = sparse.coo_matrix(
                (
                    np.concatenate(coo_values),
                    (np.concatenate(coo_rows), np.concatenate(coo_columns)),
                ),
                shape=(n_slab_spots, len(global_index)),
            ).tocsr()
            sparse_block.sum_duplicates()
            sparse_block.eliminate_zeros()
        else:
            sparse_block = sparse.csr_matrix(
                (n_slab_spots, len(global_index)), dtype=matrix.dtype
            )
        row_blocks.append(sparse_block)
        source_spot_sum[spot_start:spot_stop] = np.asarray(
            sum_accumulator, dtype=float
        )
        source_detected[spot_start:spot_stop] = detected_accumulator
        zonation_landmark_abundance[spot_start:spot_stop] = np.asarray(
            landmark_accumulator, dtype=float
        )
    abundance = sparse.vstack(row_blocks, format="csr")
    audit = {
        "source_matrix_shape": "x".join(map(str, matrix.shape)),
        "source_matrix_dtype": str(matrix.dtype),
        "source_matrix_semantics": "mean_background_subtracted_abundance_truncated_at_zero",
        "matrix_orientation": plan.orientation,
        "matrix_traversal": plan.traversal,
        "source_matrix_hdf5_chunks": plan.physical_chunk_shape,
        "logical_chunk_shape_spot_by_gene": plan.logical_chunk_shape,
        "dense_tile_target_mib": dense_tile_mib,
        "n_matrix_read_calls": len(plan.tiles),
        "n_physical_hdf5_chunks": plan.n_physical_chunks,
        "optimized_estimated_chunk_touches": plan.optimized_chunk_touches,
        "legacy_256spot_estimated_chunk_touches": plan.legacy_chunk_touches,
        "estimated_chunk_touch_reduction_factor": plan.estimated_reduction_factor,
        "estimated_chunk_touch_reduction_percent": plan.estimated_reduction_percent,
        "optimized_estimated_logical_decompressed_gib": plan.optimized_logical_gib,
        "legacy_256spot_estimated_logical_decompressed_gib": plan.legacy_logical_gib,
        "physical_chunk_coverage_min": plan.chunk_coverage_min,
        "physical_chunk_coverage_max": plan.chunk_coverage_max,
        "minimum_atomic_chunk_promoted_over_budget": plan.atomic_chunk_promoted,
        "maximum_dense_tile_mib": plan.maximum_tile_bytes / (1024**2),
        "n_local_genes": len(local_genes),
        "n_retained_local_gene_columns": int(np.sum(local_to_global >= 0)),
        "n_output_genes": len(global_index),
        "n_output_nonzero": abundance.nnz,
    }
    return (
        abundance,
        zonation_landmark_abundance,
        source_spot_sum,
        source_detected,
        audit,
    )


@dataclass(frozen=True)
class MatrixReadPlan:
    orientation: str
    traversal: str
    tiles: tuple[tuple[int, int, int, int], ...]
    physical_chunk_shape: str
    logical_chunk_shape: str
    n_physical_chunks: int
    optimized_chunk_touches: int
    legacy_chunk_touches: int
    optimized_logical_gib: float
    legacy_logical_gib: float
    estimated_reduction_factor: float
    estimated_reduction_percent: float
    chunk_coverage_min: int | str
    chunk_coverage_max: int | str
    atomic_chunk_promoted: bool
    maximum_tile_bytes: int


def infer_matrix_orientation(
    shape: tuple[int, ...],
    n_spots: int,
    n_genes: int,
    orientation_hint: str | None = None,
) -> str:
    allowed = {"h5py_spot_by_gene", "h5py_gene_by_spot_transposed_on_read"}
    if orientation_hint is not None and orientation_hint not in allowed:
        raise ContractError(f"Invalid matrix orientation hint: {orientation_hint}")
    spot_match = tuple(shape) == (n_spots, n_genes)
    gene_match = tuple(shape) == (n_genes, n_spots)
    if spot_match and gene_match and orientation_hint is None:
        raise ContractError(
            "Square source matrix has ambiguous spot/gene axes; an explicit orientation is required"
        )
    if orientation_hint == "h5py_spot_by_gene" and not spot_match:
        raise ContractError("Spot-by-gene orientation hint contradicts source matrix shape")
    if orientation_hint == "h5py_gene_by_spot_transposed_on_read" and not gene_match:
        raise ContractError("Gene-by-spot orientation hint contradicts source matrix shape")
    if orientation_hint is not None:
        return orientation_hint
    if spot_match:
        return "h5py_spot_by_gene"
    if gene_match:
        return "h5py_gene_by_spot_transposed_on_read"
    raise ContractError(
        f"Source abundance matrix shape {shape} does not match "
        f"{n_spots} spots and {n_genes} genes"
    )


def build_matrix_read_plan(
    matrix: h5py.Dataset,
    n_spots: int,
    n_genes: int,
    *,
    dense_tile_mib: float,
    orientation_hint: str | None = None,
) -> MatrixReadPlan:
    if not np.isfinite(dense_tile_mib) or dense_tile_mib <= 0:
        raise ContractError("Dense matrix tile budget must be finite and positive")
    orientation = infer_matrix_orientation(
        tuple(matrix.shape), n_spots, n_genes, orientation_hint
    )
    target_bytes = max(1, int(dense_tile_mib * 1024**2))
    itemsize = int(np.dtype(matrix.dtype).itemsize)
    tiles: list[tuple[int, int, int, int]] = []
    if matrix.chunks is None:
        if orientation == "h5py_spot_by_gene":
            spots_per_tile = max(1, target_bytes // max(1, n_genes * itemsize))
            for spot_start in range(0, n_spots, spots_per_tile):
                tiles.append(
                    (spot_start, min(spot_start + spots_per_tile, n_spots), 0, n_genes)
                )
        else:
            genes_per_tile = max(1, target_bytes // max(1, n_spots * itemsize))
            for gene_start in range(0, n_genes, genes_per_tile):
                tiles.append(
                    (0, n_spots, gene_start, min(gene_start + genes_per_tile, n_genes))
                )
        maximum = max(
            (s1 - s0) * (g1 - g0) * itemsize for s0, s1, g0, g1 in tiles
        )
        return MatrixReadPlan(
            orientation=orientation,
            traversal="contiguous_bounded_fallback",
            tiles=tuple(tiles),
            physical_chunk_shape="contiguous",
            logical_chunk_shape="not_applicable",
            n_physical_chunks=0,
            optimized_chunk_touches=0,
            legacy_chunk_touches=0,
            optimized_logical_gib=float("nan"),
            legacy_logical_gib=float("nan"),
            estimated_reduction_factor=float("nan"),
            estimated_reduction_percent=float("nan"),
            chunk_coverage_min="not_applicable",
            chunk_coverage_max="not_applicable",
            atomic_chunk_promoted=False,
            maximum_tile_bytes=maximum,
        )

    physical_chunks = tuple(int(value) for value in matrix.chunks)
    if orientation == "h5py_spot_by_gene":
        spot_chunk, gene_chunk = physical_chunks
    else:
        gene_chunk, spot_chunk = physical_chunks
    n_spot_chunks = int(np.ceil(n_spots / spot_chunk))
    n_gene_chunks = int(np.ceil(n_genes / gene_chunk))
    atomic_chunk_bytes = spot_chunk * gene_chunk * itemsize
    chunks_per_tile = max(1, target_bytes // max(1, atomic_chunk_bytes))
    if n_spot_chunks <= chunks_per_tile:
        spot_chunks_per_tile = n_spot_chunks
        gene_chunks_per_tile = max(1, chunks_per_tile // n_spot_chunks)
    else:
        spot_chunks_per_tile = chunks_per_tile
        gene_chunks_per_tile = 1
    spot_chunks_per_tile = min(spot_chunks_per_tile, n_spot_chunks)
    gene_chunks_per_tile = min(gene_chunks_per_tile, n_gene_chunks)
    coverage = np.zeros((n_spot_chunks, n_gene_chunks), dtype=np.uint8)
    for spot_chunk_start in range(0, n_spot_chunks, spot_chunks_per_tile):
        spot_chunk_stop = min(
            spot_chunk_start + spot_chunks_per_tile, n_spot_chunks
        )
        for gene_chunk_start in range(0, n_gene_chunks, gene_chunks_per_tile):
            gene_chunk_stop = min(
                gene_chunk_start + gene_chunks_per_tile, n_gene_chunks
            )
            coverage[
                spot_chunk_start:spot_chunk_stop,
                gene_chunk_start:gene_chunk_stop,
            ] += 1
            tiles.append(
                (
                    spot_chunk_start * spot_chunk,
                    min(spot_chunk_stop * spot_chunk, n_spots),
                    gene_chunk_start * gene_chunk,
                    min(gene_chunk_stop * gene_chunk, n_genes),
                )
            )
    if coverage.min() != 1 or coverage.max() != 1:
        raise ContractError("Chunk-aligned HDF5 plan does not cover each physical chunk once")
    legacy_touches = 0
    legacy_bytes = 0
    for spot_start in range(0, n_spots, 256):
        spot_stop = min(spot_start + 256, n_spots)
        first_chunk = spot_start // spot_chunk
        last_chunk = (spot_stop - 1) // spot_chunk
        for spot_chunk_index in range(first_chunk, last_chunk + 1):
            chunk_height = min(
                spot_chunk, n_spots - spot_chunk_index * spot_chunk
            )
            legacy_touches += n_gene_chunks
            legacy_bytes += chunk_height * n_genes * itemsize
    optimized_bytes = n_spots * n_genes * itemsize
    reduction = legacy_bytes / optimized_bytes
    maximum = max(
        (s1 - s0) * (g1 - g0) * itemsize for s0, s1, g0, g1 in tiles
    )
    return MatrixReadPlan(
        orientation=orientation,
        traversal="hdf5_physical_chunk_aligned_exact_once",
        tiles=tuple(tiles),
        physical_chunk_shape="x".join(map(str, physical_chunks)),
        logical_chunk_shape=f"{spot_chunk}x{gene_chunk}",
        n_physical_chunks=n_spot_chunks * n_gene_chunks,
        optimized_chunk_touches=n_spot_chunks * n_gene_chunks,
        legacy_chunk_touches=legacy_touches,
        optimized_logical_gib=optimized_bytes / (1024**3),
        legacy_logical_gib=legacy_bytes / (1024**3),
        estimated_reduction_factor=reduction,
        estimated_reduction_percent=100 * (1 - 1 / reduction),
        chunk_coverage_min=int(coverage.min()),
        chunk_coverage_max=int(coverage.max()),
        atomic_chunk_promoted=atomic_chunk_bytes > target_bytes,
        maximum_tile_bytes=maximum,
    )


def read_logical_tile(
    matrix: h5py.Dataset,
    orientation: str,
    tile: tuple[int, int, int, int],
) -> np.ndarray:
    spot_start, spot_stop, gene_start, gene_stop = tile
    if orientation == "h5py_spot_by_gene":
        destination = np.empty(
            (spot_stop - spot_start, gene_stop - gene_start), dtype=matrix.dtype
        )
        matrix.read_direct(
            destination,
            source_sel=np.s_[spot_start:spot_stop, gene_start:gene_stop],
        )
        return destination
    if orientation == "h5py_gene_by_spot_transposed_on_read":
        destination = np.empty(
            (gene_stop - gene_start, spot_stop - spot_start), dtype=matrix.dtype
        )
        matrix.read_direct(
            destination,
            source_sel=np.s_[gene_start:gene_stop, spot_start:spot_stop],
        )
        return destination.T
    raise ContractError(f"Unexpected matrix orientation: {orientation}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--dense-tile-mib", type=float, default=64.0)
    args = parser.parse_args()
    paths = default_paths(args.base)
    output = (args.output_dir or paths["candidate"]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    final_h5ad = output / "yakubovsky_human_adapter.h5ad"
    final_tables = [
        output / "sample_manifest.tsv",
        output / "gene_mapping_audit.tsv",
        output / "design_audit.tsv",
        output / "source_manifest.tsv",
        output / "adapter_execution_manifest.tsv",
    ]
    existing = [str(path) for path in [final_h5ad, *final_tables] if path.exists()]
    if existing:
        raise ContractError("Refusing to overwrite adapter output(s): " + ", ".join(existing))

    freeze_path = output / "analysis_freeze_manifest.tsv"
    if not freeze_path.is_file():
        raise ContractError("Outcome-blind analysis freeze must exist before adapter construction")
    validate_stage_chain(output, "freeze")
    frozen = freeze_file_map(freeze_path)
    verify_frozen_files(frozen)
    registry = validate_registry_contract(paths)
    gate = validate_source_gate_summary(paths["source_gate"])
    if frozen["plan10_source_gate"]["sha256"] != gate.source_gate_sha256:
        raise ContractError("Frozen Plan 10 gate hash does not match the current terminal gate")
    if frozen["plan20_registry"]["sha256"] != registry.registry_sha256:
        raise ContractError("Frozen Plan 20 registry hash does not match the current registry")
    if "yakubovsky_gene_axis_equivalence" not in frozen:
        raise ContractError("Pre-outcome freeze lacks the exact gene-axis certificate")
    gene_axis = validate_gene_axis_manifest(
        paths["gene_axis_manifest"], frozen["yakubovsky_v_mat"]["sha256"]
    )
    if (
        frozen["yakubovsky_gene_axis_equivalence"]["sha256"]
        != gene_axis.manifest_sha256
    ):
        raise ContractError("Frozen Yakubovsky gene-axis certificate hash drift")

    required_join = [
        "patient",
        "main_feature",
        "spot_index_1based",
        "barcode",
        "composite_spot_id",
        "x_coordinate",
        "y_coordinate",
        "eta",
        "zone_index",
        "structural_expression_coordinate_zonation_join",
        "source_defined_binary_lipid_class",
        "authoritative_ordinal_lipid_join",
    ]
    join_header = pd.read_csv(paths["source_join"], sep="\t", nrows=0).columns
    if not set(required_join).issubset(join_header):
        raise ContractError(f"Plan 10 join audit lacks {sorted(required_join)}")
    # ``usecols`` is a load-bearing firewall: the same source audit also carries
    # an unauthorised continuous lipid field, which this binary adapter never
    # materializes.
    join = pd.read_csv(
        paths["source_join"],
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=required_join,
    )
    if join["composite_spot_id"].duplicated().any():
        raise ContractError("Plan 10 join audit contains duplicate composite spot IDs")
    join = join.set_index("composite_spot_id", drop=False)

    with h5py.File(paths["v_mat"], "r") as handle:
        references = np.asarray(handle["v"][...]).reshape(-1, order="F")
        if len(references) != len(gene_axis.rows):
            raise ContractError(
                f"v.mat donor count disagrees with gene-axis certificate: "
                f"{len(references)} != {len(gene_axis.rows)}"
            )
        donor_info: list[dict[str, Any]] = []
        cell_axes: dict[str, h5py.Dataset] = {}
        seen_patients: set[str] = set()
        for reference, expected_axis in zip(
            references, gene_axis.rows.itertuples(index=False), strict=True
        ):
            group = handle[reference]
            patient = decode_matlab_char(group["patient"])
            if patient in seen_patients:
                raise ContractError(f"Duplicate patient struct in v.mat: {patient}")
            seen_patients.add(patient)
            if patient != expected_axis.donor or group.name != expected_axis.group_path:
                raise ContractError(
                    "v.mat donor/group order disagrees with exact gene-axis certificate: "
                    f"{patient}@{group.name} != "
                    f"{expected_axis.donor}@{expected_axis.group_path}"
                )
            n_spots = int(group["spot_name"].size)
            matrix_shape = tuple(int(value) for value in group["mat"].shape)
            if matrix_shape not in {
                (n_spots, gene_axis.n_genes),
                (gene_axis.n_genes, n_spots),
            }:
                raise ContractError(
                    f"Source matrix axes disagree with certified common gene axis for {patient}"
                )
            cell_axes[f"barcodes::{patient}"] = group["spot_name"]
            if patient == gene_axis.master_donor:
                cell_axes["master_gene_axis"] = group["gene_name"]
            donor_info.append(
                {
                    "patient": patient,
                    "group_path": group.name,
                    "main_feature": decode_matlab_char(group["main_feature"]),
                }
            )
        resolved_cellstr = bulk_decode_matlab_cellstr(
            cell_axes, handle, allow_null=False
        )
        common_genes = resolved_cellstr.values["master_gene_axis"]
        if len(common_genes) != gene_axis.n_genes:
            raise ContractError("Decoded master gene-axis length disagrees with certificate")
        if gene_axis_sha256(common_genes) != gene_axis.common_sha256:
            raise ContractError("Decoded master gene-axis content disagrees with certificate")
        for info in donor_info:
            info["barcodes"] = resolved_cellstr.values[
                f"barcodes::{info['patient']}"
            ]

    donor_gene_lists = {
        info["patient"]: common_genes for info in donor_info
    }

    expected_join_rows = sum(len(info["barcodes"]) for info in donor_info)
    if len(join) != expected_join_rows:
        raise ContractError(
            "Plan 10 join row universe does not equal the complete v.mat spot universe: "
            f"{len(join)} != {expected_join_rows}"
        )

    mapping, gene_audit, canonical_genes, symbol_sources = map_source_genes(
        donor_gene_lists, paths["gencode"]
    )
    if not canonical_genes:
        raise ContractError("No source genes mapped unambiguously")
    global_index = {gene: index for index, gene in enumerate(canonical_genes)}

    abundance_parts: list[sparse.csr_matrix] = []
    obs_rows: list[dict[str, Any]] = []
    donor_audits: list[dict[str, Any]] = []
    with h5py.File(paths["v_mat"], "r") as handle:
        for info in donor_info:
            patient = info["patient"]
            group = handle[info["group_path"]]
            barcodes = info["barcodes"]
            n_spots = len(barcodes)
            coordinates = np.asarray(group["coor"][...], dtype=float).T
            eta = numeric_vector(group["eta"])
            zone_index = numeric_vector(group["zon_struct/zone_index"])
            if coordinates.shape != (n_spots, 2) or len(eta) != n_spots or len(zone_index) != n_spots:
                raise ContractError(f"Source axes do not align for {patient}")
            zonation_landmark_indices: list[int] = []
            for field in ("LM_pc_ind", "LM_pp_ind"):
                if field in group:
                    raw = numeric_vector(group[field])
                    if not np.allclose(raw, np.round(raw)):
                        raise ContractError(f"Non-integer source landmark index in {patient}/{field}")
                    zonation_landmark_indices.extend((raw.astype(int) - 1).tolist())
            zonation_landmark_indices_array = np.asarray(
                sorted(
                    set(index for index in zonation_landmark_indices if index >= 0)
                ),
                dtype=int,
            )
            if (
                len(zonation_landmark_indices_array)
                and zonation_landmark_indices_array.max()
                >= len(donor_gene_lists[patient])
            ):
                raise ContractError(f"Out-of-range source landmark index for {patient}")
            (
                abundance,
                zonation_landmark_raw,
                spot_sum,
                detected,
                matrix_audit,
            ) = stream_donor_abundance(
                group["mat"],
                n_spots,
                donor_gene_lists[patient],
                mapping,
                global_index,
                zonation_landmark_indices_array,
                args.dense_tile_mib,
            )
            zonation_landmark_per_million = np.divide(
                1e6 * zonation_landmark_raw,
                spot_sum,
                out=np.zeros_like(spot_sum),
                where=spot_sum > 0,
            )
            zonation_landmark_expression = np.log1p(
                zonation_landmark_per_million
            )
            mapped_spot_sum = np.asarray(abundance.sum(axis=1)).ravel().astype(float)
            abundance_parts.append(abundance)

            donor_class_counts: Counter[str] = Counter()
            n_structural_join = 0
            n_analysis_eligible = 0
            for index, barcode in enumerate(barcodes, start=1):
                composite = f"{patient}:{barcode}"
                if composite not in join.index:
                    raise ContractError(f"v.mat spot absent from Plan 10 join audit: {composite}")
                source = join.loc[composite]
                if isinstance(source, pd.DataFrame):
                    raise ContractError(f"Duplicate Plan 10 join row: {composite}")
                if (
                    source["patient"] != patient
                    or source["barcode"] != barcode
                    or source["main_feature"] != info["main_feature"]
                ):
                    raise ContractError(f"Source identity drift for {composite}")
                structural_join = bool_value(
                    source["structural_expression_coordinate_zonation_join"]
                )
                n_structural_join += int(structural_join)
                if int(source["spot_index_1based"]) != index:
                    raise ContractError(f"Spot-index drift for {composite}")
                if structural_join:
                    if not np.allclose(
                        [float(source["x_coordinate"]), float(source["y_coordinate"])],
                        coordinates[index - 1],
                        atol=1e-10,
                        rtol=0,
                    ):
                        raise ContractError(f"Coordinate drift for {composite}")
                    if not np.isclose(
                        float(source["eta"]), eta[index - 1], atol=1e-12, rtol=0
                    ):
                        raise ContractError(f"Zonation drift for {composite}")
                    if not np.isclose(
                        float(source["zone_index"]),
                        zone_index[index - 1],
                        atol=1e-12,
                        rtol=0,
                    ):
                        raise ContractError(f"Categorical zonation drift for {composite}")
                analysis_eligible = source_spot_analysis_eligible(
                    structural_join,
                    spot_sum[index - 1],
                )
                n_analysis_eligible += int(analysis_eligible)
                binary_inference_eligible, lipid_class = (
                    source_spot_binary_lipid_state(
                        patient in gate.passing_donors,
                        structural_join,
                        bool_value(source["authoritative_ordinal_lipid_join"]),
                        source["source_defined_binary_lipid_class"],
                    )
                )
                donor_class_counts[lipid_class] += 1
                obs_rows.append(
                    {
                        "spot_id": composite,
                        "dataset": DATASET,
                        "donor": patient,
                        "section": patient,
                        "barcode": barcode,
                        "main_feature": info["main_feature"],
                        "x_coordinate": coordinates[index - 1, 0],
                        "y_coordinate": coordinates[index - 1, 1],
                        "zonation_eta": eta[index - 1],
                        "source_zone_index": zone_index[index - 1],
                        "source_structural_expression_coordinate_zonation_join": structural_join,
                        "source_defined_binary_lipid_class": lipid_class,
                        "lipid_zone": h5ad_lipid_zone_value(
                            lipid_class, binary_inference_eligible
                        ),
                        "binary_lipid_inference_eligible": binary_inference_eligible,
                        "analysis_eligible": analysis_eligible,
                        "spot_eligibility_basis": "authoritative_structural_join_and_positive_source_abundance_no_hepatocyte_filter",
                        "background_corrected_spot_sum": spot_sum[index - 1],
                        "mapped_background_corrected_spot_sum": mapped_spot_sum[index - 1],
                        "detected_genes": detected[index - 1],
                        "source_zonation_landmark_expression": zonation_landmark_expression[
                            index - 1
                        ],
                    }
                )
            n_binary_lipid_eligible = sum(
                donor_class_counts[label] for label in ALLOWED_LIPID_CLASSES
            )
            if patient in gate.passing_donors:
                if (
                    n_binary_lipid_eligible < MIN_DONOR_SPOTS
                    or n_binary_lipid_eligible / n_spots < 0.90
                    or not all(
                        donor_class_counts[label] > 0
                        for label in ALLOWED_LIPID_CLASSES
                    )
                ):
                    raise ContractError(
                        f"Adapter cannot rederive the donor-level ordinal gate for {patient}"
                    )
            donor_audits.append(
                {
                    "audit_stage": "adapter",
                    "dataset": DATASET,
                    "donor": patient,
                    "section": patient,
                    "n_spots": n_spots,
                    "n_source_structural_join_spots": n_structural_join,
                    "n_source_structural_join_failed_reference_spots": n_spots
                    - n_structural_join,
                    "n_analysis_eligible_spots": n_analysis_eligible,
                    "n_lipid_zone_spots": donor_class_counts["lipid_zone"],
                    "n_non_lipid_zone_spots": donor_class_counts["non_lipid_zone"],
                    "n_binary_lipid_eligible_spots": n_binary_lipid_eligible,
                    "binary_lipid_eligible_spot_fraction": n_binary_lipid_eligible
                    / n_spots,
                    "binary_lipid_inference_eligible": patient in gate.passing_donors,
                    "binary_label_nondegenerate": (
                        all(donor_class_counts[label] > 0 for label in ALLOWED_LIPID_CLASSES)
                        if patient in gate.passing_donors
                        else "not_applicable"
                    ),
                    "authoritative_join_fraction": n_structural_join / n_spots,
                    "n_source_zonation_landmark_genes": len(
                        zonation_landmark_indices_array
                    ),
                    "common_gene_axis_reused": True,
                    "common_gene_axis_sha256": gene_axis.common_sha256,
                    "gene_axis_manifest_sha256": gene_axis.manifest_sha256,
                    "bulk_cellstr_n_references": resolved_cellstr.n_references,
                    "bulk_cellstr_n_unique_objects": resolved_cellstr.n_unique_objects,
                    "bulk_cellstr_n_objects_visited": resolved_cellstr.n_objects_visited,
                    **matrix_audit,
                }
            )

    abundance_all = sparse.vstack(abundance_parts, format="csr")
    obs = pd.DataFrame(obs_rows).set_index("spot_id", drop=False)
    if len(obs) != abundance_all.shape[0] or abundance_all.shape[1] != len(canonical_genes):
        raise ContractError("Adapter abundance and annotation dimensions disagree")
    if obs.index.duplicated().any():
        raise ContractError("Adapter spot IDs are not unique")
    if abundance_all.data.size and (
        (abundance_all.data < 0).any() or not np.isfinite(abundance_all.data).all()
    ):
        raise ContractError("Collapsed source background-corrected abundance is invalid")
    spot_sum_all = obs["background_corrected_spot_sum"].to_numpy(float)
    scale = np.divide(
        1e6,
        spot_sum_all,
        out=np.zeros_like(spot_sum_all),
        where=spot_sum_all > 0,
    )
    normalized = (sparse.diags(scale) @ abundance_all).tocsr().astype(np.float32)
    normalized.data = np.log1p(normalized.data)

    var = pd.DataFrame(index=pd.Index(canonical_genes, name="gene_symbol"))
    var["source_gene_ids"] = [
        ";".join(sorted(set(symbol_sources[gene]))) for gene in canonical_genes
    ]
    var["n_source_gene_ids"] = [len(set(symbol_sources[gene])) for gene in canonical_genes]
    adata = ad.AnnData(
        X=None,
        obs=obs,
        var=var,
        shape=abundance_all.shape,
    )
    adata.layers["source_background_corrected_abundance"] = abundance_all
    adata.layers["log1p_source_abundance_per_million"] = normalized
    adata.obsm["spatial"] = obs[["x_coordinate", "y_coordinate"]].to_numpy(float)
    adata.uns.update(
        {
            "release_id": RELEASE_ID,
            "dataset": DATASET,
            "source_matrix_location": "layers/source_background_corrected_abundance",
            "source_matrix_semantics": "mean_background_subtracted_abundance_truncated_at_zero_not_raw_counts",
            "normalization_denominator": "background_corrected_spot_sum_before_ambiguous_gene_exclusion",
            "analysis_layer": "log1p_source_abundance_per_million",
            "analysis_transform": "log1p(1e6*abundance/background_corrected_spot_sum)",
            "count_likelihood_used": False,
            "lipid_exposure_type": "source_defined_binary_lipid_zone_vs_non_lipid_zone",
            "non_lipid_zone_interpretation": "complement_among_retained_spots_not_lipid_free",
            "continuous_lipid_fields_read": False,
            "cell2location_used": False,
            "authoritative_hepatocyte_eligibility_available": False,
            "spot_eligibility_basis": "authoritative_structural_join_and_positive_source_abundance_no_hepatocyte_filter",
            "full_source_spot_universe_preserved": True,
            "structural_join_required_for_analysis": True,
            "source_zonation_landmark_interpretation": "periportal_and_pericentral_zonation_landmark_expression_not_hepatocyte_identity_or_purity",
            "registry_sha256": registry.registry_sha256,
            "membership_sha256": registry.membership_sha256,
            "source_gate_sha256": gate.source_gate_sha256,
            "source_join_sha256": sha256_file(paths["source_join"]),
            "analysis_freeze_manifest_sha256": sha256_file(freeze_path),
            "common_gene_axis_sha256": gene_axis.common_sha256,
            "gene_axis_manifest_sha256": gene_axis.manifest_sha256,
            "common_gene_axis_reused_after_exact_all_donor_verification": True,
            "matrix_traversal": "hdf5_physical_chunk_aligned_exact_once_or_contiguous_bounded_fallback",
            "dense_tile_target_mib": args.dense_tile_mib,
        }
    )
    with tempfile.NamedTemporaryFile(
        suffix=".h5ad", dir=output, prefix=".yakubovsky_human_adapter.", delete=False
    ) as handle:
        temporary_h5ad = Path(handle.name)
    try:
        adata.write_h5ad(temporary_h5ad, compression="gzip")
        os.replace(temporary_h5ad, final_h5ad)
    finally:
        if temporary_h5ad.exists():
            temporary_h5ad.unlink()

    sample_manifest = obs.reset_index(drop=True)
    atomic_write_frame(output / "sample_manifest.tsv", sample_manifest)
    atomic_write_frame(output / "gene_mapping_audit.tsv", gene_audit)
    atomic_write_frame(output / "design_audit.tsv", pd.DataFrame(donor_audits))
    source_manifest = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "role": role,
                "path": record["value"],
                "sha256": record["sha256"],
            }
            for role, record in frozen.items()
        ]
    )
    atomic_write_frame(output / "source_manifest.tsv", source_manifest)
    execution = pd.DataFrame(
        [
            {
                "release_id": RELEASE_ID,
                "stage": "source_native_adapter",
                "producer": str(Path(__file__).resolve()),
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
                "completed_utc": utc_now(),
                "n_donors": obs["donor"].nunique(),
                "n_sections": obs[["donor", "section"]].drop_duplicates().shape[0],
                "n_spots": len(obs),
                "n_genes": len(canonical_genes),
                "source_matrix_semantics": "mean_background_subtracted_abundance_truncated_at_zero",
                "binary_labels_opened": True,
                "continuous_lipid_fields_read": False,
                "image_reclassification_used": False,
                "cell2location_used": False,
                "count_likelihood_used": False,
                "authoritative_hepatocyte_eligibility_available": False,
                "spot_eligibility_basis": "authoritative_structural_join_and_positive_source_abundance_no_hepatocyte_filter",
                "full_source_spot_universe_preserved": True,
                "structural_join_required_for_analysis": True,
                "source_zonation_landmarks_used_only_as_zonation_overadjustment_sensitivity": True,
                "common_gene_axis_reused_after_exact_all_donor_verification": True,
                "common_gene_axis_sha256": gene_axis.common_sha256,
                "gene_axis_manifest_sha256": gene_axis.manifest_sha256,
                "bulk_cellstr_n_references": resolved_cellstr.n_references,
                "bulk_cellstr_n_unique_objects": resolved_cellstr.n_unique_objects,
                "bulk_cellstr_n_objects_visited": resolved_cellstr.n_objects_visited,
                "dense_tile_target_mib": args.dense_tile_mib,
                "exit_state": "pass",
            }
        ]
    )
    atomic_write_frame(output / "adapter_execution_manifest.tsv", execution)
    write_stage_seal(output, "adapter", ADAPTER_STAGE_ARTIFACTS)
    print(
        json.dumps(
            {
                "release_id": RELEASE_ID,
                "status": "adapter_complete",
                "n_donors": int(obs["donor"].nunique()),
                "n_passing_binary_lipid_donors": len(gate.passing_donors),
                "n_spots": len(obs),
                "n_genes": len(canonical_genes),
                "h5ad_sha256": sha256_file(final_h5ad),
                "continuous_lipid_fields_read": False,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
