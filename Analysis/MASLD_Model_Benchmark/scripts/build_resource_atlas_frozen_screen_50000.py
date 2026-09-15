#!/usr/bin/env python3
"""Build the deterministic 50,000-cell Resource Atlas frozen-screen fixture."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
from pathlib import Path
from typing import Any, Iterable, Iterator


ONTOLOGY: dict[str, tuple[str, ...]] = {
    "cholangiocyte": ("Cholangiocytes",),
    "endothelial": ("Endothelial cells",),
    "hepatocyte": ("Hepatocytes",),
    "immune": (
        "B cells",
        "Basophils",
        "Circulating NK/NKT",
        "Macrophages",
        "Mono+mono derived cells",
        "Neutrophils",
        "Plasma cells",
        "Resident NK",
        "T cells",
        "cDC1s",
        "cDC2s",
        "pDCs",
    ),
    "mesenchymal_stromal": ("Fibroblasts",),
}
ONTOLOGY_ID = "broad_liver_five_v1"
SUBSET_ID = "resource_atlas_frozen_screen_50000_v1"
CELL_BUDGET = 50_000
CELLS_PER_CLASS = 10_000
SELECTION_SEED = 20260824
RAW_HASH_CHUNK = 4_000_000
REQUIRED_MANIFEST_FIELDS = {
    "cell_id",
    "source_cell_id",
    "library_id",
    "donor_id",
    "dataset",
    "cell_type",
    "analysis_eligible",
}
OBS_FIELDS = (
    "row_id",
    "source_row_index",
    "source_cell_id",
    "library_id",
    "donor_id",
    "dataset",
    "source_cell_type",
    "broad_label",
    "selection_class_rank",
    "donor_round_index",
    "selection_priority_sha256",
    "n_counts",
)


class FrozenScreenError(ValueError):
    """Raised when an input or output does not meet the frozen-screen requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def write_json_exclusive(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(
            value,
            handle,
            sort_keys=True,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")


def _manifest_rows(path: Path) -> Iterator[tuple[int, dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not REQUIRED_MANIFEST_FIELDS.issubset(
            reader.fieldnames
        ):
            raise FrozenScreenError("cell manifest lacks required source fields")
        for row_index, row in enumerate(reader):
            yield row_index, row


def _raw_to_broad() -> dict[str, str]:
    result: dict[str, str] = {}
    for broad, source_labels in ONTOLOGY.items():
        for source_label in source_labels:
            if source_label in result:
                raise FrozenScreenError(
                    f"source label occurs twice in ontology: {source_label}"
                )
            result[source_label] = broad
    return result


def priority_digest(*parts: object, seed: int = SELECTION_SEED) -> bytes:
    payload = "\x1f".join(map(str, (seed, *parts))).encode("utf-8")
    return hashlib.sha256(payload).digest()


def priority(*parts: object, seed: int = SELECTION_SEED) -> int:
    return int.from_bytes(priority_digest(*parts, seed=seed), "big")


def load_contract(
    path: Path,
    *,
    expected_file_sha256: str,
    expected_lock_sha256: str,
    expected_raw_counts_sha256: str,
    expected_cell_order_sha256: str,
    expected_gene_order_sha256: str,
) -> dict[str, Any]:
    if sha256_file(path) != expected_file_sha256:
        raise FrozenScreenError("production contract-lock file SHA-256 differs")
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema_version") != "masld-cl-contract-v1":
        raise FrozenScreenError("unsupported production contract schema")
    payload = {key: value for key, value in contract.items() if key != "lock_sha256"}
    if canonical_sha256(payload) != contract.get("lock_sha256"):
        raise FrozenScreenError("production contract-lock content hash differs")
    expected = {
        "lock_sha256": expected_lock_sha256,
        "raw_counts_sha256": expected_raw_counts_sha256,
        "cell_order_sha256": expected_cell_order_sha256,
        "gene_order_sha256": expected_gene_order_sha256,
    }
    for key, wanted in expected.items():
        if contract.get(key) != wanted:
            raise FrozenScreenError(f"production contract {key} differs")
    if not contract.get("production_ready"):
        raise FrozenScreenError("production contract is not release-ready")
    required_manifests = {
        "assay_manifest.tsv",
        "cell_manifest.tsv.gz",
        "donor_manifest.tsv",
        "library_manifest.tsv",
    }
    if set(contract.get("manifest_sha256", {})) != required_manifests:
        raise FrozenScreenError("production contract manifest roster differs")
    return contract


def scan_capacities(
    cell_manifest: Path,
) -> tuple[dict[str, dict[str, int]], dict[str, Any]]:
    raw_to_broad = _raw_to_broad()
    capacities: dict[str, dict[str, int]] = {broad: {} for broad in ONTOLOGY}
    cell_order = hashlib.sha256()
    rows = 0
    eligible_rows = 0
    for row_index, row in _manifest_rows(cell_manifest):
        rows = row_index + 1
        row_id = row["cell_id"]
        if row_id != f"{row['library_id']}|{row['source_cell_id']}":
            raise FrozenScreenError(f"compound row identity differs at row {row_index}")
        cell_order.update(row_id.encode("utf-8"))
        cell_order.update(b"\0")
        if row["analysis_eligible"] != "True":
            continue
        broad = raw_to_broad.get(row["cell_type"])
        if broad is None:
            raise FrozenScreenError(
                f"analysis-eligible cell has unmapped type {row['cell_type']!r}"
            )
        donor = row["donor_id"]
        capacities[broad][donor] = capacities[broad].get(donor, 0) + 1
        eligible_rows += 1
    if not rows:
        raise FrozenScreenError("cell manifest is empty")
    return capacities, {
        "manifest_rows": rows,
        "analysis_eligible_rows": eligible_rows,
        "cell_order_sha256": cell_order.hexdigest(),
    }


def allocate_round_robin(
    capacities: dict[str, int], *, target: int, donor_order: Iterable[str]
) -> dict[str, int]:
    order = list(donor_order)
    if not order or len(order) != len(set(order)) or set(order) != set(capacities):
        raise FrozenScreenError("donor order does not match capacity roster")
    if any(value < 0 for value in capacities.values()):
        raise FrozenScreenError("donor capacity cannot be negative")
    if sum(capacities.values()) < target:
        raise FrozenScreenError(
            f"donor capacities provide {sum(capacities.values())}, require {target}"
        )
    quotas = {donor: 0 for donor in order}
    selected = 0
    while selected < target:
        progressed = False
        for donor in order:
            if quotas[donor] >= capacities[donor]:
                continue
            quotas[donor] += 1
            selected += 1
            progressed = True
            if selected == target:
                break
        if not progressed:
            raise FrozenScreenError("round-robin quota allocation stalled")
    if sum(quotas.values()) != target or any(
        quotas[donor] > capacities[donor] for donor in order
    ):
        raise FrozenScreenError("round-robin quota allocation is inconsistent")
    return quotas


def select_cells(
    cell_manifest: Path, *, cells_per_class: int = CELLS_PER_CLASS
) -> tuple[list[dict[str, str | int]], dict[str, Any]]:
    capacities, scan = scan_capacities(cell_manifest)
    donor_orders: dict[str, list[str]] = {}
    quotas: dict[str, dict[str, int]] = {}
    for broad in sorted(ONTOLOGY):
        donor_orders[broad] = sorted(
            capacities[broad],
            key=lambda donor: (priority("donor", broad, donor), donor),
        )
        quotas[broad] = allocate_round_robin(
            capacities[broad],
            target=cells_per_class,
            donor_order=donor_orders[broad],
        )

    pools: dict[
        str,
        dict[str, list[tuple[int, int, str, dict[str, str | int]]]],
    ] = {
        broad: {donor: [] for donor, quota in quotas[broad].items() if quota > 0}
        for broad in ONTOLOGY
    }
    raw_to_broad = _raw_to_broad()
    for row_index, row in _manifest_rows(cell_manifest):
        if row["analysis_eligible"] != "True":
            continue
        broad = raw_to_broad[row["cell_type"]]
        donor = row["donor_id"]
        quota = quotas[broad].get(donor, 0)
        if quota == 0:
            continue
        row_id = row["cell_id"]
        score_bytes = priority_digest("cell", broad, donor, row_id)
        score = int.from_bytes(score_bytes, "big")
        record: dict[str, str | int] = {
            "source_row_index": row_index,
            "row_id": row_id,
            "source_cell_id": row["source_cell_id"],
            "library_id": row["library_id"],
            "donor_id": donor,
            "dataset": row["dataset"],
            "source_cell_type": row["cell_type"],
            "broad_label": broad,
            "selection_priority_sha256": score_bytes.hex(),
        }
        heap = pools[broad][donor]
        item = (-score, -row_index, row_id, record)
        if len(heap) < quota:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    selected: list[dict[str, str | int]] = []
    capacity_rows: list[dict[str, str | int]] = []
    for broad in sorted(ONTOLOGY):
        donor_rows: dict[str, list[dict[str, str | int]]] = {}
        for donor in donor_orders[broad]:
            expected = quotas[broad][donor]
            observed = len(pools[broad].get(donor, []))
            if observed != expected:
                raise FrozenScreenError(
                    f"selected pool differs for {broad}/{donor}: {observed} != {expected}"
                )
            donor_rows[donor] = [
                item[3]
                for item in sorted(
                    pools[broad].get(donor, []),
                    key=lambda item: (-item[0], -item[1], item[2]),
                )
            ]
            capacity_rows.append(
                {
                    "broad_label": broad,
                    "donor_order_rank": donor_orders[broad].index(donor),
                    "donor_id": donor,
                    "source_capacity": capacities[broad][donor],
                    "selected_quota": expected,
                }
            )
        class_rows: list[dict[str, str | int]] = []
        round_index = 0
        while len(class_rows) < cells_per_class:
            progressed = False
            for donor in donor_orders[broad]:
                values = donor_rows[donor]
                if round_index >= len(values):
                    continue
                row = values[round_index]
                row["selection_class_rank"] = len(class_rows)
                row["donor_round_index"] = round_index
                class_rows.append(row)
                progressed = True
                if len(class_rows) == cells_per_class:
                    break
            if not progressed:
                raise FrozenScreenError(f"selection stalled for class {broad}")
            round_index += 1
        selected.extend(class_rows)

    expected_rows = cells_per_class * len(ONTOLOGY)
    if len(selected) != expected_rows:
        raise FrozenScreenError(
            f"selection produced {len(selected)} cells, require {expected_rows}"
        )
    row_ids = [str(row["row_id"]) for row in selected]
    if len(set(row_ids)) != len(row_ids):
        raise FrozenScreenError("selection repeats an immutable row ID")
    selected.sort(key=lambda row: int(row["source_row_index"]))
    capacity_rows.sort(
        key=lambda row: (str(row["broad_label"]), int(row["donor_order_rank"]))
    )
    return selected, {
        **scan,
        "capacity_rows": capacity_rows,
        "donor_orders": donor_orders,
        "capacities": capacities,
        "quotas": quotas,
    }


def _decode_array(dataset: Any) -> Any:
    values = dataset[:]
    if values.dtype.kind in {"S", "O", "U"}:
        return dataset.asstr()[:]
    return values


def _read_h5ad_column(group: Any, key: str) -> tuple[Any, Any]:
    import h5py
    import numpy as np

    node = group[key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = _decode_array(node["categories"])
        codes = node["codes"][:].astype(np.int64)
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise FrozenScreenError(
                f"source observation column {key!r} has invalid codes"
            )
        return categories, codes
    values = _decode_array(node)
    categories, codes = np.unique(values, return_inverse=True)
    return categories, codes.astype(np.int64)


def _hash_strings(values: Iterable[object]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _copy_selected_intervals(
    *,
    dataset: Any,
    output: Any,
    source_starts: Any,
    source_ends: Any,
    target_starts: Any,
    digest: Any,
    component: str,
    n_features: int,
) -> None:
    import numpy as np

    interval = 0
    n_intervals = len(source_starts)
    for chunk_start in range(0, len(dataset), RAW_HASH_CHUNK):
        chunk_end = min(chunk_start + RAW_HASH_CHUNK, len(dataset))
        values = dataset[chunk_start:chunk_end]
        if component == "data":
            if not np.all(np.isfinite(values)):
                raise FrozenScreenError("raw/X contains non-finite counts")
            if np.any(values < 0) or not np.all(values == np.floor(values)):
                raise FrozenScreenError("raw/X is not nonnegative integer UMI data")
        elif component == "indices":
            if np.any(values < 0) or np.any(values >= n_features):
                raise FrozenScreenError("raw/X contains an out-of-range feature index")
        digest.update(np.ascontiguousarray(values).view(np.uint8))
        while interval < n_intervals:
            source_start = int(source_starts[interval])
            source_end = int(source_ends[interval])
            if source_end <= source_start:
                interval += 1
                continue
            if source_start >= chunk_end:
                break
            overlap_start = max(source_start, chunk_start)
            overlap_end = min(source_end, chunk_end)
            if overlap_end > overlap_start:
                target_start = int(target_starts[interval]) + overlap_start - source_start
                target_end = target_start + overlap_end - overlap_start
                output[target_start:target_end] = values[
                    overlap_start - chunk_start : overlap_end - chunk_start
                ]
            if source_end <= chunk_end:
                interval += 1
            else:
                break
    while interval < n_intervals and source_ends[interval] <= source_starts[interval]:
        interval += 1
    if interval != n_intervals:
        raise FrozenScreenError(f"selected {component} extraction did not finish")


def materialize_subset(
    *,
    atlas: Path,
    selected: list[dict[str, str | int]],
    output_h5ad: Path,
    contract: dict[str, Any],
) -> dict[str, Any]:
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if len(selected) != CELL_BUDGET:
        raise FrozenScreenError(
            f"materialization requires exactly {CELL_BUDGET} selected cells"
        )
    row_indices = np.asarray(
        [int(row["source_row_index"]) for row in selected], dtype=np.int64
    )
    if len(row_indices) == 0 or np.any(row_indices[1:] <= row_indices[:-1]):
        raise FrozenScreenError("selected source row indices are not strictly ordered")

    with h5py.File(atlas, "r") as handle:
        raw_x = handle["raw/X"]
        source_shape = tuple(int(value) for value in raw_x.attrs["shape"])
        expected_shape = tuple(int(value) for value in contract["raw_counts"]["shape"])
        if source_shape != expected_shape:
            raise FrozenScreenError("source raw-count shape differs from contract")
        if len(raw_x["data"]) != int(contract["raw_counts"]["nnz"]):
            raise FrozenScreenError("source raw-count nnz differs from contract")

        source_features = list(map(str, _decode_array(handle["raw/var/_index"])))
        if _hash_strings(source_features) != contract["gene_order_sha256"]:
            raise FrozenScreenError("source gene order differs from contract")
        source_gene_ids = list(map(str, _decode_array(handle["raw/var/gene_ids"])))
        if len(source_features) != len(source_gene_ids):
            raise FrozenScreenError("source raw-variable axes differ")
        ensembl_ids = [value.split(".", 1)[0] for value in source_gene_ids]
        if (
            any(not value.startswith("ENSG") for value in ensembl_ids)
            or len(set(ensembl_ids)) != len(ensembl_ids)
        ):
            raise FrozenScreenError("source gene_ids cannot form a unique Ensembl axis")

        obs = handle["obs"]
        sample_categories, sample_codes = _read_h5ad_column(obs, "sample")
        cell_categories, cell_codes = _read_h5ad_column(obs, "_index")
        if len(sample_codes) != source_shape[0] or len(cell_codes) != source_shape[0]:
            raise FrozenScreenError("source observation axis differs from raw counts")
        cell_digest = hashlib.sha256()
        selected_cursor = 0
        for source_row, (sample_code, cell_code) in enumerate(
            zip(sample_codes, cell_codes, strict=True)
        ):
            library_id = str(sample_categories[sample_code])
            source_cell_id = str(cell_categories[cell_code])
            row_id = f"{library_id}|{source_cell_id}"
            cell_digest.update(row_id.encode("utf-8"))
            cell_digest.update(b"\0")
            if (
                selected_cursor < len(selected)
                and source_row == int(selected[selected_cursor]["source_row_index"])
            ):
                row = selected[selected_cursor]
                if (
                    row_id != str(row["row_id"])
                    or library_id != str(row["library_id"])
                    or source_cell_id != str(row["source_cell_id"])
                ):
                    raise FrozenScreenError(
                        f"selected row identity differs at source row {source_row}"
                    )
                selected_cursor += 1
        if selected_cursor != len(selected):
            raise FrozenScreenError("not all selected rows occur in the source Atlas")
        if cell_digest.hexdigest() != contract["cell_order_sha256"]:
            raise FrozenScreenError("source cell identity/order differs from contract")

        source_indptr_native = raw_x["indptr"][:]
        source_indptr = source_indptr_native.astype(np.int64, copy=False)
        if (
            len(source_indptr) != source_shape[0] + 1
            or source_indptr[0] != 0
            or source_indptr[-1] != len(raw_x["data"])
            or np.any(source_indptr[1:] < source_indptr[:-1])
        ):
            raise FrozenScreenError("source raw CSR indptr is invalid")
        source_starts = source_indptr[row_indices]
        source_ends = source_indptr[row_indices + 1]
        subset_indptr = np.empty(len(selected) + 1, dtype=np.int64)
        subset_indptr[0] = 0
        np.cumsum(source_ends - source_starts, out=subset_indptr[1:])
        selected_nnz = int(subset_indptr[-1])
        selected_indices = np.empty(selected_nnz, dtype=raw_x["indices"].dtype)
        selected_data = np.empty(selected_nnz, dtype=raw_x["data"].dtype)

        digest = hashlib.sha256()
        digest.update(json.dumps(source_shape).encode("ascii"))
        digest.update(b"indptr")
        digest.update(str(raw_x["indptr"].dtype).encode("ascii"))
        digest.update(np.ascontiguousarray(source_indptr_native).view(np.uint8))
        digest.update(b"indices")
        digest.update(str(raw_x["indices"].dtype).encode("ascii"))
        _copy_selected_intervals(
            dataset=raw_x["indices"],
            output=selected_indices,
            source_starts=source_starts,
            source_ends=source_ends,
            target_starts=subset_indptr[:-1],
            digest=digest,
            component="indices",
            n_features=source_shape[1],
        )
        digest.update(b"data")
        digest.update(str(raw_x["data"].dtype).encode("ascii"))
        _copy_selected_intervals(
            dataset=raw_x["data"],
            output=selected_data,
            source_starts=source_starts,
            source_ends=source_ends,
            target_starts=subset_indptr[:-1],
            digest=digest,
            component="data",
            n_features=source_shape[1],
        )
        if digest.hexdigest() != contract["raw_counts_sha256"]:
            raise FrozenScreenError("source raw-count SHA-256 differs from contract")

    matrix = sparse.csr_matrix(
        (selected_data, selected_indices, subset_indptr),
        shape=(len(selected), len(ensembl_ids)),
    )
    n_counts_float = np.asarray(matrix.sum(axis=1)).ravel()
    if (
        np.any(~np.isfinite(n_counts_float))
        or np.any(n_counts_float <= 0)
        or np.any(n_counts_float != np.floor(n_counts_float))
    ):
        raise FrozenScreenError("selected raw UMI totals are invalid")
    n_counts = n_counts_float.astype(np.int64)

    obs_frame = pd.DataFrame(selected)
    obs_frame["n_counts"] = n_counts
    obs_frame = obs_frame.loc[:, list(OBS_FIELDS)].set_index("row_id", drop=False)
    if list(obs_frame.columns) != list(OBS_FIELDS):
        raise FrozenScreenError("output observation field roster differs")
    var_frame = pd.DataFrame(
        {
            "ensembl_id": ensembl_ids,
            "source_feature_id": source_features,
            "source_gene_id": source_gene_ids,
        },
        index=pd.Index(ensembl_ids, name="ensembl_id_index"),
    )
    value = ad.AnnData(X=matrix, obs=obs_frame, var=var_frame)
    value.uns["masld_bench_frozen_screen_contract"] = {
        "subset_id": SUBSET_ID,
        "ontology_id": ONTOLOGY_ID,
        "cell_budget": CELL_BUDGET,
        "cells_per_class": CELLS_PER_CLASS,
        "selection_seed": SELECTION_SEED,
        "selection_policy": "class_then_donor_round_robin_then_sha256_cell_priority",
        "selection_label_use": (
            "source_cell_type_and_broad_label_used_only_for_prespecified_"
            "class_stratified_sampling"
        ),
        "matrix_construction_label_blind": True,
        "matrix": "raw_nonnegative_integer_UMI_counts",
        "sealed_outcomes_used": False,
        "histology_used": False,
    }
    value.write_h5ad(output_h5ad, compression="gzip")

    content_identity = {
        "data_sha256": hashlib.sha256(matrix.data.tobytes()).hexdigest(),
        "indices_sha256": hashlib.sha256(matrix.indices.tobytes()).hexdigest(),
        "indptr_sha256": hashlib.sha256(matrix.indptr.tobytes()).hexdigest(),
        "row_ids_sha256": _hash_strings(value.obs_names),
        "ensembl_ids_sha256": _hash_strings(ensembl_ids),
        "n_counts_sha256": hashlib.sha256(n_counts.tobytes()).hexdigest(),
    }
    with h5py.File(output_h5ad, "r") as written:
        if tuple(int(item) for item in written["X"].attrs["shape"]) != value.shape:
            raise FrozenScreenError("written H5AD matrix shape differs")
        for component in ("data", "indices", "indptr"):
            dataset = written[f"X/{component}"]
            observed = hashlib.sha256()
            for start in range(0, len(dataset), RAW_HASH_CHUNK):
                block = dataset[start : start + RAW_HASH_CHUNK]
                observed.update(np.ascontiguousarray(block).view(np.uint8))
            if observed.hexdigest() != content_identity[f"{component}_sha256"]:
                raise FrozenScreenError(f"written H5AD {component} differs")
    return {
        "shape": list(value.shape),
        "nnz": int(value.X.nnz),
        "content_identity": content_identity,
        "raw_umi_counts_preserved": True,
        "source_raw_counts_sha256_verified": True,
    }


def _write_tsv(path: Path, fields: tuple[str, ...], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def freeze_output_tree(output: Path, *, manifest_sha256: str) -> None:
    records: list[dict[str, Any]] = []
    for path in sorted(output.rglob("*")):
        if path.is_symlink():
            raise FrozenScreenError(f"derived fixture contains a symlink: {path}")
        if path.is_file() and path.name not in {"ARTIFACTS.json", "COMPLETE"}:
            records.append(
                {
                    "path": path.relative_to(output).as_posix(),
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
            )
    artifacts = {
        "schema_version": "masld-bench-artifacts-v1",
        "metadata": {
            "artifact_class": "resource_atlas_frozen_screen_fixture",
            "artifact_id": manifest_sha256,
            "subset_id": SUBSET_ID,
            "row_count": CELL_BUDGET,
            "biological_unit": "donor",
            "sealed_outcomes_read": False,
            "histology_read": False,
            "status": "passed",
        },
        "artifacts": records,
    }
    artifacts_path = output / "ARTIFACTS.json"
    write_json_exclusive(artifacts_path, artifacts)
    complete_path = output / "COMPLETE"
    write_json_exclusive(
        complete_path,
        {
            "schema_version": "masld-bench-complete-v1",
            "manifest_sha256": sha256_file(artifacts_path),
            "artifact_count": len(records),
        },
    )
    for record in records:
        path = output / record["path"]
        if sha256_file(path) != record["sha256"] or path.stat().st_size != record[
            "size_bytes"
        ]:
            raise FrozenScreenError(f"derived artifact changed while freezing: {path}")
        path.chmod(0o440)
    artifacts_path.chmod(0o440)
    complete_path.chmod(0o440)
    output.chmod(0o550)


def build(
    *,
    atlas: Path,
    contract_lock: Path,
    cell_manifest: Path,
    output: Path,
    expected_contract_file_sha256: str,
    expected_contract_lock_sha256: str,
    expected_manifest_sha256: str,
    expected_raw_counts_sha256: str,
    expected_cell_order_sha256: str,
    expected_gene_order_sha256: str,
) -> None:
    if output.exists():
        raise FrozenScreenError(f"output already exists: {output}")
    contract = load_contract(
        contract_lock,
        expected_file_sha256=expected_contract_file_sha256,
        expected_lock_sha256=expected_contract_lock_sha256,
        expected_raw_counts_sha256=expected_raw_counts_sha256,
        expected_cell_order_sha256=expected_cell_order_sha256,
        expected_gene_order_sha256=expected_gene_order_sha256,
    )
    expected_from_contract = contract["manifest_sha256"].get(cell_manifest.name)
    if (
        expected_from_contract != expected_manifest_sha256
        or sha256_file(cell_manifest) != expected_manifest_sha256
    ):
        raise FrozenScreenError("cell manifest differs from production contract")
    if atlas.resolve(strict=True) != Path(contract["atlas_realpath"]).resolve(strict=True):
        raise FrozenScreenError("Atlas path differs from production contract")
    if atlas.stat().st_size != int(contract["atlas_size_bytes"]):
        raise FrozenScreenError("Atlas size differs from production contract")

    selected, selection_audit = select_cells(cell_manifest)
    if selection_audit["cell_order_sha256"] != expected_cell_order_sha256:
        raise FrozenScreenError("cell-manifest identity/order differs from contract")
    output.mkdir(parents=True, exist_ok=False)

    ontology = {
        "schema_version": "masld-bench-cell-ontology-v1",
        "ontology_id": ONTOLOGY_ID,
        "source_field": "cell_type",
        "derived_field": "broad_label",
        "classes": {key: list(value) for key, value in sorted(ONTOLOGY.items())},
        "selection_label_use": (
            "source_cell_type_and_broad_label_used_only_for_prespecified_"
            "class_stratified_sampling"
        ),
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
        "histology_used": False,
    }
    ontology["ontology_sha256"] = canonical_sha256(ontology)
    write_json_exclusive(output / f"{ONTOLOGY_ID}.json", ontology)

    selection_fields = (
        "source_row_index",
        "row_id",
        "source_cell_id",
        "library_id",
        "donor_id",
        "dataset",
        "source_cell_type",
        "broad_label",
        "selection_class_rank",
        "donor_round_index",
        "selection_priority_sha256",
    )
    _write_tsv(output / "selection.tsv", selection_fields, selected)
    capacity_fields = (
        "broad_label",
        "donor_order_rank",
        "donor_id",
        "source_capacity",
        "selected_quota",
    )
    _write_tsv(
        output / "donor_capacity.tsv",
        capacity_fields,
        selection_audit["capacity_rows"],
    )

    output_h5ad = output / "resource_atlas_frozen_screen_50000.h5ad"
    matrix = materialize_subset(
        atlas=atlas,
        selected=selected,
        output_h5ad=output_h5ad,
        contract=contract,
    )
    class_counts: dict[str, int] = {}
    dataset_counts: dict[str, int] = {}
    donors_by_class: dict[str, set[str]] = {}
    for row in selected:
        broad = str(row["broad_label"])
        dataset = str(row["dataset"])
        donor = str(row["donor_id"])
        class_counts[broad] = class_counts.get(broad, 0) + 1
        dataset_counts[dataset] = dataset_counts.get(dataset, 0) + 1
        donors_by_class.setdefault(broad, set()).add(donor)
    if class_counts != {broad: CELLS_PER_CLASS for broad in sorted(ONTOLOGY)}:
        raise FrozenScreenError("fixture is not exactly balanced across five classes")

    manifest = {
        "schema_version": "masld-bench-resource-atlas-frozen-screen-v1",
        "dataset_id": "resource_atlas_current",
        "subset_id": SUBSET_ID,
        "source": {
            "atlas_path": atlas.resolve(strict=True).as_posix(),
            "atlas_size_bytes": atlas.stat().st_size,
            "contract_lock_path": contract_lock.resolve(strict=True).as_posix(),
            "contract_lock_file_sha256": expected_contract_file_sha256,
            "contract_lock_id": expected_contract_lock_sha256,
            "raw_counts_sha256": expected_raw_counts_sha256,
            "cell_order_sha256": expected_cell_order_sha256,
            "gene_order_sha256": expected_gene_order_sha256,
            "cell_manifest_path": cell_manifest.resolve(strict=True).as_posix(),
            "cell_manifest_sha256": expected_manifest_sha256,
        },
        "selection": {
            "cell_budget": CELL_BUDGET,
            "cells_per_class": CELLS_PER_CLASS,
            "selection_seed": SELECTION_SEED,
            "selection_policy": (
                "class_then_donor_round_robin_then_sha256_cell_priority"
            ),
            "selection_label_use": (
                "source_cell_type_and_broad_label_used_only_for_prespecified_"
                "class_stratified_sampling"
            ),
            "counts_by_class": dict(sorted(class_counts.items())),
            "donors_by_class": {
                key: len(value) for key, value in sorted(donors_by_class.items())
            },
            "counts_by_dataset": dict(sorted(dataset_counts.items())),
            "donor_capacity_verified": True,
            "immutable_row_ids_unique": True,
            "selection_outcomes_used": False,
            "sealed_outcomes_used": False,
            "histology_used": False,
        },
        "matrix": matrix,
        "input_contract": {
            "matrix": "X_raw_nonnegative_integer_UMI_counts",
            "obs_fields": list(OBS_FIELDS),
            "obs_n_counts_rederived": True,
            "var_ensembl_id": True,
            "labels_used_in_matrix_construction": False,
            "histology_fields_present": [],
            "sealed_fields_present": [],
            "ready": True,
        },
        "artifacts": {
            "h5ad": "resource_atlas_frozen_screen_50000.h5ad",
            "selection": "selection.tsv",
            "donor_capacity": "donor_capacity.tsv",
            "ontology": f"{ONTOLOGY_ID}.json",
        },
    }
    manifest["manifest_sha256"] = canonical_sha256(manifest)
    write_json_exclusive(output / "frozen_screen_manifest.json", manifest)
    freeze_output_tree(output, manifest_sha256=manifest["manifest_sha256"])


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--atlas", required=True, type=Path)
    value.add_argument("--contract-lock", required=True, type=Path)
    value.add_argument("--cell-manifest", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-contract-file-sha256", required=True)
    value.add_argument("--expected-contract-lock-sha256", required=True)
    value.add_argument("--expected-manifest-sha256", required=True)
    value.add_argument("--expected-raw-counts-sha256", required=True)
    value.add_argument("--expected-cell-order-sha256", required=True)
    value.add_argument("--expected-gene-order-sha256", required=True)
    return value


def main() -> int:
    args = parser().parse_args()
    build(
        atlas=args.atlas,
        contract_lock=args.contract_lock,
        cell_manifest=args.cell_manifest,
        output=args.output,
        expected_contract_file_sha256=args.expected_contract_file_sha256,
        expected_contract_lock_sha256=args.expected_contract_lock_sha256,
        expected_manifest_sha256=args.expected_manifest_sha256,
        expected_raw_counts_sha256=args.expected_raw_counts_sha256,
        expected_cell_order_sha256=args.expected_cell_order_sha256,
        expected_gene_order_sha256=args.expected_gene_order_sha256,
    )
    print((args.output / "frozen_screen_manifest.json").resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
