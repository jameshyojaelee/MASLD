#!/usr/bin/env python3
"""Build a deterministic same-nucleus GSE296875 RNA-ATAC smoke view."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any, Iterable, Mapping, Sequence


VIEW_ID = "gse296875_rna_atac_smoke_1000_v1"
PARENT_DATASET_ID = "gse296875"
SELECTION_SEED = 20260821
CELLS_PER_CLASS = 200
EXPECTED_CELLS = 68_398
EXPECTED_DONORS = 39
EXPECTED_WELLS = tuple(f"well{index}" for index in range(1, 9))
SOURCE_LOCK_SHA256 = "6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620"
MAPPING_METADATA_SHA256 = "5d3d946607f55b140ffa4a44369cf35f6e9e91f0f17276afc0e3e9f008554781"
AUTHOR_LABELS_SHA256 = "45c077e9beb8607402bd2ed7bce0586596430d401f290ca86978f91e230ea67d"
PARENT_REGISTRY_SHA256 = "9465b012bbf46611641555a5acae4d1f033d053c900bf944dcec02979c735b6f"
PRIMARY_CONTIGS = tuple(
    [f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"]
)
LABEL_MAP = {
    "Cholangiocytes": "cholangiocyte",
    "Hepatocytes": "hepatocyte",
    "Kupffer": "macrophage",
    "Mesenchymal": "fibroblast",
    "NK-T": "t_cell",
}


class GSE296875SmokeError(ValueError):
    """Raised when source data does not meet the frozen smoke-view requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _parse_exact_integer(value: str, *, label: str) -> int:
    try:
        parsed = Decimal(value)
    except InvalidOperation as error:
        raise GSE296875SmokeError(f"{label} is not numeric") from error
    if not parsed.is_finite() or parsed != parsed.to_integral_value():
        raise GSE296875SmokeError(f"{label} is not an exact finite integer")
    return int(parsed)


def _priority(*parts: object) -> int:
    text = "\x1f".join(map(str, (SELECTION_SEED, *parts)))
    return int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest(), "big")


def _read_exact_tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise GSE296875SmokeError(
                f"{path.name} columns differ from the frozen contract"
            )
        return [dict(row) for row in reader]


def _write_tsv_exclusive(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_json_exclusive(path: Path, value: Any) -> None:
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


def _freeze_smoke_tree(root: Path, metadata: Mapping[str, Any]) -> str:
    """Freeze one builder-owned tree without importing the control plane.

    The scientific builder runs in the existing Python 3.10 Scanpy environment,
    whereas ``masld_bench`` intentionally requires Python 3.11 or newer.  This
    writer follows the registered output file schema; the SLURM job subsequently
    re-verifies the result with the Python 3.11 control-plane verifier.
    """

    if root.is_symlink() or not root.is_dir():
        raise GSE296875SmokeError("smoke artifact root must be a real directory")
    base = root.resolve(strict=True)
    if (base / "ARTIFACTS.json").exists() or (base / "COMPLETE").exists():
        raise GSE296875SmokeError("smoke artifact tree is already frozen")
    records: list[dict[str, Any]] = []
    for path in sorted(base.rglob("*"), key=lambda value: value.as_posix()):
        if path.is_symlink():
            raise GSE296875SmokeError("symlinks are forbidden in smoke artifacts")
        if not path.is_file():
            continue
        resolved = path.resolve(strict=True)
        try:
            relative = resolved.relative_to(base).as_posix()
        except ValueError as error:
            raise GSE296875SmokeError("smoke artifact escapes its root") from error
        records.append(
            {
                "path": relative,
                "sha256": sha256_file(resolved),
                "size_bytes": resolved.stat().st_size,
            }
        )
    manifest_path = base / "ARTIFACTS.json"
    _write_json_exclusive(
        manifest_path,
        {
            "schema_version": "masld-bench-artifacts-v1",
            "metadata": dict(metadata),
            "artifacts": records,
        },
    )
    _write_json_exclusive(
        base / "COMPLETE",
        {
            "schema_version": "masld-bench-complete-v1",
            "manifest_sha256": sha256_file(manifest_path),
            "artifact_count": len(records),
        },
    )
    return sha256_file(manifest_path)


def _check_source(path: Path, expected_sha256: str, label: str) -> Path:
    resolved = path.resolve(strict=True)
    if resolved.is_symlink() or not resolved.is_file():
        raise GSE296875SmokeError(f"{label} is not a regular source file")
    if sha256_file(resolved) != expected_sha256:
        raise GSE296875SmokeError(f"{label} differs from its frozen SHA-256")
    return resolved


def _source_lock(path: Path) -> dict[str, Any]:
    source = _check_source(path, SOURCE_LOCK_SHA256, "GSE296875 source lock")
    value = json.loads(source.read_text(encoding="utf-8"))
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != "masld-cl-gse296875-source-lock-v44"
        or value.get("published_expectations", {}).get(
            "biological_donors_passing_qc"
        )
        != EXPECTED_DONORS
        or value.get("published_expectations", {}).get("wells") != len(EXPECTED_WELLS)
    ):
        raise GSE296875SmokeError("GSE296875 source lock has changed semantics")
    return value


def _joined_metadata(
    mapping_metadata: Path, author_labels: Path
) -> list[dict[str, str]]:
    mapping_path = _check_source(
        mapping_metadata, MAPPING_METADATA_SHA256, "GSE296875 mapping metadata"
    )
    label_path = _check_source(
        author_labels, AUTHOR_LABELS_SHA256, "GSE296875 author labels"
    )
    mapping_rows = _read_exact_tsv(
        mapping_path, ("cell_id", "donor_id", "well_id", "raw_barcode")
    )
    label_rows = _read_exact_tsv(label_path, ("cell_id", "author_label"))
    if len(mapping_rows) != EXPECTED_CELLS or len(label_rows) != EXPECTED_CELLS:
        raise GSE296875SmokeError("GSE296875 source census is not 68,398 nuclei")
    metadata = {row["cell_id"]: row for row in mapping_rows}
    if len(metadata) != EXPECTED_CELLS:
        raise GSE296875SmokeError("GSE296875 mapping cell identifiers are not unique")
    if {row["cell_id"] for row in label_rows} != set(metadata):
        raise GSE296875SmokeError("GSE296875 mapping and label cell sets differ")
    joined: list[dict[str, str]] = []
    for label_row in label_rows:
        source = metadata[label_row["cell_id"]]
        well = source["well_id"]
        encoded_barcode = source["raw_barcode"]
        prefix = f"{well}_"
        if (
            well not in EXPECTED_WELLS
            or not encoded_barcode.startswith(prefix)
            or re.fullmatch(r"[ACGT]+-[0-9]+", encoded_barcode[len(prefix) :]) is None
        ):
            raise GSE296875SmokeError(
                "GSE296875 well-namespaced barcode encoding changed"
            )
        joined.append(
            {
                "cell_id": source["cell_id"],
                "donor_id": source["donor_id"],
                "well_id": well,
                "raw_barcode": encoded_barcode[len(prefix) :],
                "source_label": label_row["author_label"],
            }
        )
    if (
        len({row["donor_id"] for row in joined}) != EXPECTED_DONORS
        or tuple(sorted({row["well_id"] for row in joined})) != EXPECTED_WELLS
    ):
        raise GSE296875SmokeError("GSE296875 donor or well census changed")
    return joined


def select_cells(
    *, mapping_metadata: Path, author_labels: Path, output: Path
) -> None:
    joined = _joined_metadata(mapping_metadata, author_labels)
    pools: dict[str, dict[str, list[dict[str, str]]]] = {
        label: {} for label in LABEL_MAP.values()
    }
    for row in joined:
        broad = LABEL_MAP.get(row["source_label"])
        if broad is None:
            continue
        record = {**row, "broad_label": broad}
        pools[broad].setdefault(row["donor_id"], []).append(record)
    selected: list[dict[str, str]] = []
    for broad in sorted(pools):
        donor_rows = pools[broad]
        if len(donor_rows) != EXPECTED_DONORS:
            raise GSE296875SmokeError(
                f"{broad} is not observed in all 39 GSE296875 donors"
            )
        for donor, rows in donor_rows.items():
            rows.sort(key=lambda row: (_priority("cell", broad, donor, row["cell_id"]), row["cell_id"]))
        donor_order = sorted(
            donor_rows,
            key=lambda donor: (_priority("donor", broad, donor), donor),
        )
        class_rows: list[dict[str, str]] = []
        offset = 0
        while len(class_rows) < CELLS_PER_CLASS:
            progressed = False
            for donor in donor_order:
                rows = donor_rows[donor]
                if offset < len(rows):
                    class_rows.append(rows[offset])
                    progressed = True
                    if len(class_rows) == CELLS_PER_CLASS:
                        break
            if not progressed:
                raise GSE296875SmokeError(
                    f"GSE296875 cannot supply {CELLS_PER_CLASS} {broad} nuclei"
                )
            offset += 1
        selected.extend(class_rows)
    if (
        len(selected) != CELLS_PER_CLASS * len(LABEL_MAP)
        or len({row["cell_id"] for row in selected}) != len(selected)
        or len({row["donor_id"] for row in selected}) != EXPECTED_DONORS
    ):
        raise GSE296875SmokeError("GSE296875 smoke selection census is invalid")
    output.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv_exclusive(
        output,
        (
            "cell_id",
            "donor_id",
            "well_id",
            "raw_barcode",
            "source_label",
            "broad_label",
        ),
        selected,
    )
    receipt = {
        "schema_version": "masld-bench-gse296875-selection-v1",
        "view_id": VIEW_ID,
        "selection_seed": SELECTION_SEED,
        "selection_policy": "lineage_then_donor_round_robin_then_sha256_nucleus_priority",
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
        "atac_values_or_qc_used_for_selection": False,
        "cells": len(selected),
        "donors": len({row["donor_id"] for row in selected}),
        "wells": len({row["well_id"] for row in selected}),
        "counts_by_class": {
            label: sum(row["broad_label"] == label for row in selected)
            for label in sorted(pools)
        },
        "selection_sha256": sha256_file(output),
    }
    _write_json_exclusive(output.with_suffix(".json"), receipt)


def _decode(values: Any) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def _rna_matrix(
    *, selected: Sequence[Mapping[str, str]], source_lock: Mapping[str, Any]
) -> tuple[Any, list[str], list[str], list[dict[str, Any]]]:
    import h5py
    import numpy as np
    from scipy import sparse

    lock_root = Path(source_lock["_path"]).parent
    raw_sources = source_lock.get("raw_rna_sources")
    if not isinstance(raw_sources, Mapping) or set(raw_sources) != set(EXPECTED_WELLS):
        raise GSE296875SmokeError("GSE296875 raw source roster changed")
    by_well: dict[str, list[tuple[int, Mapping[str, str]]]] = {
        well: [] for well in EXPECTED_WELLS
    }
    for index, row in enumerate(selected):
        by_well[row["well_id"]].append((index, row))
    output_rows: list[Any] = [None] * len(selected)
    gene_ids: list[str] | None = None
    gene_names: list[str] | None = None
    source_records: list[dict[str, Any]] = []
    for well in EXPECTED_WELLS:
        raw = raw_sources[well]
        source = (lock_root / str(raw["path"])).resolve(strict=True)
        expected_size = int(raw["bytes"])
        if source.stat().st_size != expected_size:
            raise GSE296875SmokeError(f"{well} raw matrix size changed")
        _check_source(source, str(raw["sha256"]), f"{well} raw matrix")
        source_records.append(
            {
                "well_id": well,
                "gsm": str(raw["gsm"]),
                "path": source.as_posix(),
                "sha256": str(raw["sha256"]),
                "size_bytes": expected_size,
            }
        )
        with h5py.File(source, "r") as handle:
            matrix = handle["matrix"]
            features = matrix["features"]
            feature_types = features["feature_type"][:]
            gene_source_indices = np.flatnonzero(feature_types == b"Gene Expression")
            observed_gene_ids = _decode(features["id"][:][gene_source_indices])
            observed_gene_names = _decode(features["name"][:][gene_source_indices])
            if (
                len(observed_gene_ids) != 36_601
                or len(set(observed_gene_ids)) != len(observed_gene_ids)
                or any(not value.startswith("ENSG") for value in observed_gene_ids)
            ):
                raise GSE296875SmokeError(f"{well} RNA feature axis changed")
            if gene_ids is None:
                gene_ids = observed_gene_ids
                gene_names = observed_gene_names
            elif gene_ids != observed_gene_ids or gene_names != observed_gene_names:
                raise GSE296875SmokeError("GSE296875 well RNA axes are not identical")
            source_to_gene = np.full(len(feature_types), -1, dtype=np.int64)
            source_to_gene[gene_source_indices] = np.arange(
                len(gene_source_indices), dtype=np.int64
            )
            barcodes = _decode(matrix["barcodes"][:])
            barcode_to_index = {barcode: index for index, barcode in enumerate(barcodes)}
            if len(barcode_to_index) != len(barcodes):
                raise GSE296875SmokeError(f"{well} raw barcodes are not unique")
            indptr = matrix["indptr"]
            for output_index, row in by_well[well]:
                barcode_index = barcode_to_index.get(row["raw_barcode"])
                if barcode_index is None:
                    raise GSE296875SmokeError(
                        f"selected barcode is absent from {well}"
                    )
                start = int(indptr[barcode_index])
                end = int(indptr[barcode_index + 1])
                source_indices = matrix["indices"][start:end].astype(np.int64)
                values = matrix["data"][start:end]
                positions = source_to_gene[source_indices]
                keep = positions >= 0
                values = values[keep]
                positions = positions[keep]
                if (
                    len(values) == 0
                    or np.any(values < 0)
                    or np.any(values != np.floor(values))
                ):
                    raise GSE296875SmokeError(
                        "selected RNA column is empty or not integer count data"
                    )
                output_rows[output_index] = sparse.csr_matrix(
                    (values, positions, [0, len(values)]),
                    shape=(1, len(gene_source_indices)),
                )
    if gene_ids is None or gene_names is None or any(row is None for row in output_rows):
        raise GSE296875SmokeError("GSE296875 RNA materialization is incomplete")
    combined = sparse.vstack(output_rows, format="csr")
    return combined, gene_ids, gene_names, source_records


def _atac_matrix(
    *, selected: Sequence[Mapping[str, str]], atac_export: Path
) -> tuple[Any, list[dict[str, Any]], str]:
    from scipy import io, sparse

    matrix_path = atac_export / "atac_counts.mtx"
    features_path = atac_export / "atac_features.tsv"
    cells_path = atac_export / "atac_cells.tsv"
    session_path = atac_export / "R_sessionInfo.txt"
    for path in (matrix_path, features_path, cells_path, session_path):
        if path.is_symlink() or not path.is_file():
            raise GSE296875SmokeError(f"ATAC export lacks {path.name}")
    cells = _read_exact_tsv(cells_path, ("cell_id",))
    if [row["cell_id"] for row in cells] != [row["cell_id"] for row in selected]:
        raise GSE296875SmokeError("RNA and ATAC selected nucleus order differs")
    raw_features = _read_exact_tsv(
        features_path,
        (
            "peak_id",
            "chromosome",
            "source_start_1based_closed",
            "source_end_1based_closed",
        ),
    )
    features: list[dict[str, Any]] = []
    for row in raw_features:
        chromosome = row["chromosome"]
        source_start = _parse_exact_integer(
            row["source_start_1based_closed"], label="ATAC peak start"
        )
        source_end = _parse_exact_integer(
            row["source_end_1based_closed"], label="ATAC peak end"
        )
        if (
            chromosome not in PRIMARY_CONTIGS
            or source_start < 1
            or source_end < source_start
        ):
            raise GSE296875SmokeError("ATAC peak coordinate is outside its contract")
        features.append(
            {
                "peak_id": row["peak_id"],
                "chromosome": chromosome,
                "source_start_1based_closed": source_start,
                "source_end_1based_closed": source_end,
                "bed_start_0based": source_start - 1,
                "bed_end_half_open": source_end,
            }
        )
    if len({row["peak_id"] for row in features}) != len(features):
        raise GSE296875SmokeError("ATAC peak identifiers are not unique")
    peak_by_cell = io.mmread(matrix_path)
    if not sparse.issparse(peak_by_cell):
        peak_by_cell = sparse.coo_matrix(peak_by_cell)
    cell_by_peak = peak_by_cell.transpose().tocsr()
    if cell_by_peak.shape != (len(selected), len(features)):
        raise GSE296875SmokeError("ATAC matrix axes differ from exported identities")
    if cell_by_peak.nnz == 0 or any(cell_by_peak.sum(axis=1).A1 <= 0):
        raise GSE296875SmokeError("selected ATAC profiles are empty")
    if any(cell_by_peak.data < 0) or any(
        cell_by_peak.data != cell_by_peak.data.astype("int64")
    ):
        raise GSE296875SmokeError("selected ATAC profiles are not integer counts")
    return cell_by_peak, features, sha256_file(session_path)


def _write_strings(group: Any, name: str, values: Sequence[str]) -> None:
    import h5py

    group.create_dataset(
        name,
        data=list(values),
        dtype=h5py.string_dtype(encoding="utf-8"),
        compression="gzip",
    )


def _write_csr(group: Any, matrix: Any) -> dict[str, Any]:
    import numpy as np

    value = matrix.tocsr()
    group.create_dataset("data", data=value.data, compression="gzip", shuffle=True)
    group.create_dataset(
        "indices", data=value.indices.astype(np.int64), compression="gzip", shuffle=True
    )
    group.create_dataset(
        "indptr", data=value.indptr.astype(np.int64), compression="gzip", shuffle=True
    )
    group.create_dataset("shape", data=np.asarray(value.shape, dtype=np.int64))
    return {
        "shape": list(value.shape),
        "nnz": int(value.nnz),
        "data_sha256": hashlib.sha256(value.data.tobytes()).hexdigest(),
        "indices_sha256": hashlib.sha256(
            value.indices.astype("int64").tobytes()
        ).hexdigest(),
        "indptr_sha256": hashlib.sha256(
            value.indptr.astype("int64").tobytes()
        ).hexdigest(),
    }


def _materialize_h5(
    *,
    path: Path,
    selected: Sequence[Mapping[str, str]],
    rna: Any,
    gene_ids: Sequence[str],
    gene_names: Sequence[str],
    atac: Any,
    peaks: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    import h5py
    import numpy as np

    with h5py.File(path, "x") as handle:
        handle.attrs["schema_version"] = "masld-bench-multimodal-h5-v1"
        handle.attrs["view_id"] = VIEW_ID
        handle.attrs["pairing_level"] = "same_nucleus"
        handle.attrs["bed_coordinate_system"] = "0_based_half_open"
        obs = handle.create_group("obs")
        for field in (
            "cell_id",
            "donor_id",
            "well_id",
            "source_label",
            "broad_label",
        ):
            _write_strings(obs, field, [row[field] for row in selected])
        _write_strings(obs, "rna_status", ["observed"] * len(selected))
        _write_strings(obs, "atac_status", ["observed"] * len(selected))
        obs.create_dataset(
            "rna_observed_mask",
            data=np.ones(len(selected), dtype=np.uint8),
            compression="gzip",
        )
        obs.create_dataset(
            "atac_observed_mask",
            data=np.ones(len(selected), dtype=np.uint8),
            compression="gzip",
        )
        rna_group = handle.create_group("rna")
        rna_identity = _write_csr(rna_group.create_group("counts_csr"), rna)
        _write_strings(rna_group, "ensembl_id", gene_ids)
        _write_strings(rna_group, "gene_name", gene_names)
        rna_group.attrs["source_scale"] = "raw_integer_umi_counts"
        atac_group = handle.create_group("atac")
        atac_identity = _write_csr(atac_group.create_group("counts_csr"), atac)
        for field in (
            "peak_id",
            "chromosome",
        ):
            _write_strings(atac_group, field, [str(row[field]) for row in peaks])
        for field in (
            "source_start_1based_closed",
            "source_end_1based_closed",
            "bed_start_0based",
            "bed_end_half_open",
        ):
            atac_group.create_dataset(
                field,
                data=np.asarray([int(row[field]) for row in peaks], dtype=np.int64),
                compression="gzip",
            )
        atac_group.attrs["source_scale"] = "raw_integer_peak_counts"
        atac_group.attrs["source_coordinate_system"] = "GRanges_1_based_closed"
        atac_group.attrs["analysis_coordinate_system"] = "BED_0_based_half_open"
    with h5py.File(path, "r") as handle:
        if (
            handle.attrs.get("schema_version") != "masld-bench-multimodal-h5-v1"
            or tuple(handle["rna/counts_csr/shape"][:]) != rna.shape
            or tuple(handle["atac/counts_csr/shape"][:]) != atac.shape
            or len(handle["obs/cell_id"]) != len(selected)
        ):
            raise GSE296875SmokeError("written multimodal HDF5 failed readback")
    return {"rna": rna_identity, "atac": atac_identity}


def _artifact(path: Path) -> dict[str, Any]:
    return {
        "path": path.resolve(strict=True).as_posix(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _authority_documents(
    *,
    authority_root: Path,
    source_evidence: Mapping[str, Any],
    matrix_identity: Mapping[str, Any],
) -> dict[str, Path]:
    authority_root.mkdir(parents=True, exist_ok=False)
    evidence = {
        "rights": {
            "access": "public_GEO_download",
            "redistribution": "source_reference_only",
            "automatic_download": False,
            "source_terms": "source_terms_apply",
        },
        "topology": {
            "pairing": "same_nucleus",
            "rna": "raw_integer_UMI_counts",
            "atac": "raw_integer_processed_unified_peak_counts",
            "held_ATAC_available_to_RNA_only_inference": False,
            "well_role": "technical_batch_not_biological_replicate",
        },
        "donor_join": {
            "donors": EXPECTED_DONORS,
            "wells": len(EXPECTED_WELLS),
            "mapping_metadata_sha256": MAPPING_METADATA_SHA256,
            "join_key": "well_id+raw_10x_barcode",
            "all_modalities_from_one_donor_stay_in_one_outer_fold": True,
        },
        "labels": {
            "cell_state_context_only": sorted(LABEL_MAP.values()),
            "selection_outcomes_used": False,
            "ATAC_used_for_selection": False,
            "pathology_labels_included": False,
            "MASLD_or_MASH_adjudication_available": False,
        },
        "qc": {
            "selected_nuclei": CELLS_PER_CLASS * len(LABEL_MAP),
            "selected_donors": EXPECTED_DONORS,
            "lineage_counts": {
                label: CELLS_PER_CLASS for label in sorted(LABEL_MAP.values())
            },
            "missing_as_zero": False,
            "rna_shape": matrix_identity["rna"]["shape"],
            "atac_shape": matrix_identity["atac"]["shape"],
        },
        "reference": {
            "native_reference": "10x_GRCh38_2020-A",
            "native_gene_annotation": "GENCODE_v32_filtered",
            "analysis_contigs": list(PRIMARY_CONTIGS),
            "source_peak_coordinates": "GRanges_1_based_closed",
            "analysis_peak_coordinates": "BED_0_based_half_open",
            "project_genome": "GRCh38.p14",
            "project_genome_sha256": "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215",
            "project_gtf": "GENCODE_v49",
            "project_gtf_sha256": "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9",
            "sequence_extraction_ready": False,
        },
        "exposure_audit": {
            "role": "development_only",
            "exposure_status": "downstream_demo",
            "sealed_claim_eligible": False,
            "source_lock_sha256": SOURCE_LOCK_SHA256,
            "author_labels_sha256": AUTHOR_LABELS_SHA256,
        },
    }
    paths: dict[str, Path] = {}
    for authority, authority_evidence in evidence.items():
        path = authority_root / f"{authority}.json"
        _write_json_exclusive(
            path,
            {
                "schema_version": "masld-bench-dataset-view-authority-v1",
                "view_id": VIEW_ID,
                "parent_dataset_id": PARENT_DATASET_ID,
                "authority": authority,
                "admitted": True,
                "evidence": {**authority_evidence, "sources": dict(source_evidence)},
            },
        )
        paths[authority] = path
    return paths


def assemble(
    *,
    selection: Path,
    source_lock_path: Path,
    atac_export: Path,
    output: Path,
    authority_root: Path,
    generated_contract: Path,
) -> None:
    source_lock = _source_lock(source_lock_path)
    source_lock["_path"] = source_lock_path.resolve(strict=True).as_posix()
    selected = _read_exact_tsv(
        selection,
        (
            "cell_id",
            "donor_id",
            "well_id",
            "raw_barcode",
            "source_label",
            "broad_label",
        ),
    )
    if len(selected) != CELLS_PER_CLASS * len(LABEL_MAP):
        raise GSE296875SmokeError("selection does not contain exactly 1,000 nuclei")
    output.mkdir(parents=True, exist_ok=False)
    rna, gene_ids, gene_names, raw_source_records = _rna_matrix(
        selected=selected, source_lock=source_lock
    )
    atac, peaks, r_session_sha256 = _atac_matrix(
        selected=selected, atac_export=atac_export
    )
    data_path = output / "gse296875_rna_atac_smoke_1000.h5"
    matrix_identity = _materialize_h5(
        path=data_path,
        selected=selected,
        rna=rna,
        gene_ids=gene_ids,
        gene_names=gene_names,
        atac=atac,
        peaks=peaks,
    )
    selection_path = output / "selection.tsv"
    shutil.copyfile(selection, selection_path)
    ontology_path = output / "rna_atac_five_lineages_v1.json"
    ontology = {
        "schema_version": "masld-bench-cell-ontology-v1",
        "ontology_id": "rna_atac_five_lineages_v1",
        "source_field": "author_label",
        "classes": {
            label: sorted(source for source, target in LABEL_MAP.items() if target == label)
            for label in sorted(LABEL_MAP.values())
        },
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
    }
    ontology["ontology_sha256"] = canonical_sha256(ontology)
    _write_json_exclusive(ontology_path, ontology)
    source_evidence = {
        "source_lock_sha256": SOURCE_LOCK_SHA256,
        "mapping_metadata_sha256": MAPPING_METADATA_SHA256,
        "author_labels_sha256": AUTHOR_LABELS_SHA256,
        "parent_registry_sha256": PARENT_REGISTRY_SHA256,
        "processed_rds_sha256": source_lock["processed_source"]["sha256"],
        "raw_matrix_roster_sha256": canonical_sha256(raw_source_records),
        "R_sessionInfo_sha256": r_session_sha256,
    }
    derivation_path = output / "smoke_subset_manifest.json"
    derivation = {
        "schema_version": "masld-bench-multimodal-smoke-subset-v1",
        "dataset_id": PARENT_DATASET_ID,
        "subset_id": VIEW_ID,
        "source": source_evidence,
        "selection": {
            "cell_budget": len(selected),
            "selection_seed": SELECTION_SEED,
            "selection_policy": "lineage_then_donor_round_robin_then_sha256_nucleus_priority",
            "selection_outcomes_used": False,
            "sealed_outcomes_used": False,
            "atac_values_or_qc_used_for_selection": False,
            "counts_by_class": {
                label: sum(row["broad_label"] == label for row in selected)
                for label in sorted(LABEL_MAP.values())
            },
            "donors": len({row["donor_id"] for row in selected}),
            "wells": len({row["well_id"] for row in selected}),
        },
        "artifacts": {
            "multimodal_h5": _artifact(data_path),
            "selection": _artifact(selection_path),
            "ontology": _artifact(ontology_path),
        },
        "matrix": matrix_identity,
        "multimodal_input_contract": {
            "ready": True,
            "pairing_level": "same_nucleus",
            "rna_scale": "raw_integer_UMI_counts",
            "atac_scale": "raw_integer_peak_counts",
            "held_ATAC_forbidden_at_RNA_only_inference": True,
            "missingness_is_explicit": True,
        },
    }
    derivation["manifest_sha256"] = canonical_sha256(derivation)
    _write_json_exclusive(derivation_path, derivation)
    _freeze_smoke_tree(
        output,
        {
            "artifact_class": "multimodal_smoke_subset",
            "artifact_id": derivation["manifest_sha256"],
            "subset_id": VIEW_ID,
        },
    )
    authority_paths = _authority_documents(
        authority_root=authority_root,
        source_evidence=source_evidence,
        matrix_identity=matrix_identity,
    )
    class_counts = {
        label: CELLS_PER_CLASS for label in sorted(LABEL_MAP.values())
    }
    contract = {
        "schema_version": "masld-bench-dataset-view-v1",
        "view_id": VIEW_ID,
        "parent_dataset_id": PARENT_DATASET_ID,
        "parent_registry_sha256": PARENT_REGISTRY_SHA256,
        "purpose": "compatibility_smoke",
        "allowed_waves": ["smoke"],
        "row_count": len(selected),
        "biological_unit_count": EXPECTED_DONORS,
        "modalities": ["single_nucleus_rna", "single_nucleus_atac"],
        "pairing_levels": ["same_nucleus"],
        "label_visibility": "development_visible",
        "selection_policy": "lineage_then_donor_round_robin_then_sha256_nucleus_priority",
        "selection_seed": SELECTION_SEED,
        "selection_outcomes_used": False,
        "sealed_outcomes_used": False,
        "class_counts": class_counts,
        "artifact_manifest": {
            **_artifact(output / "ARTIFACTS.json"),
            "media_type": "application/json",
            "role": f"dataset_view_manifest:{VIEW_ID}",
        },
        "complete_receipt": {
            **_artifact(output / "COMPLETE"),
            "media_type": "application/json",
            "role": f"dataset_view_complete:{VIEW_ID}",
        },
        "data_artifact": {
            **_artifact(data_path),
            "media_type": "application/x-hdf5",
            "role": f"dataset_view_data:{VIEW_ID}",
        },
        "selection_artifact": {
            **_artifact(selection_path),
            "media_type": "text/tab-separated-values",
            "role": f"dataset_view_selection:{VIEW_ID}",
        },
        "ontology_artifact": {
            **_artifact(ontology_path),
            "media_type": "application/json",
            "role": f"dataset_view_ontology:{VIEW_ID}",
        },
        "derivation_artifact": {
            **_artifact(derivation_path),
            "media_type": "application/json",
            "role": f"dataset_view_derivation:{VIEW_ID}",
        },
        "authority_artifacts": [
            {
                **_artifact(authority_paths[authority]),
                "media_type": "application/json",
                "role": f"dataset_view_authority:{VIEW_ID}:{authority}",
            }
            for authority in (
                "rights",
                "topology",
                "donor_join",
                "labels",
                "qc",
                "reference",
                "exposure_audit",
            )
        ],
        "ready": True,
    }
    generated_contract.parent.mkdir(parents=True, exist_ok=True)
    _write_json_exclusive(generated_contract, contract)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    subparsers = value.add_subparsers(dest="command", required=True)
    select = subparsers.add_parser("select")
    select.add_argument("--mapping-metadata", required=True, type=Path)
    select.add_argument("--author-labels", required=True, type=Path)
    select.add_argument("--output", required=True, type=Path)
    assemble_parser = subparsers.add_parser("assemble")
    assemble_parser.add_argument("--selection", required=True, type=Path)
    assemble_parser.add_argument("--source-lock", required=True, type=Path)
    assemble_parser.add_argument("--atac-export", required=True, type=Path)
    assemble_parser.add_argument("--output", required=True, type=Path)
    assemble_parser.add_argument("--authority-root", required=True, type=Path)
    assemble_parser.add_argument("--generated-contract", required=True, type=Path)
    return value


def main() -> int:
    args = parser().parse_args()
    if args.command == "select":
        select_cells(
            mapping_metadata=args.mapping_metadata,
            author_labels=args.author_labels,
            output=args.output,
        )
        print(args.output.resolve().as_posix())
    else:
        assemble(
            selection=args.selection,
            source_lock_path=args.source_lock,
            atac_export=args.atac_export,
            output=args.output,
            authority_root=args.authority_root,
            generated_contract=args.generated_contract,
        )
        print(args.generated_contract.resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
