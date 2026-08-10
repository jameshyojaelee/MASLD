#!/usr/bin/env python3
"""Audit the Yakubovsky MATLAB source join without reading program outcomes.

This script is deliberately narrower than ``audit_yakubovsky_schema.py``.  It
dereferences only identifiers and covariates needed to establish the source
chain

    donor -> barcode -> expression axis -> coordinates -> zonation -> lipid

and never reads the expression matrix values or a program registry.  A shared
row count is recorded as structural alignment, not as an authoritative
continuous-lipid join: the pinned public code must also document how the lipid
measurement was attached to the barcode.  The separately documented Loupe
barcode path is audited as a binary lipid-zone fallback and never interpreted
as a gradient.  These distinctions implement the source gate fixed in
docs/plans/2026-08-07_paper_program/10_SPATIAL_ACQUISITION_AND_SOURCE_GATES.md.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
import pandas as pd


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
PAPER_DOI = "10.1038/s41586-026-10377-y"

# Nature 2026, Methods ("Lipid droplet quantification") and Extended Data
# Fig. 10a.  These are source-paper inclusion counts, not values inferred from
# the deposited object.
PAPER_LIPID_SPOTS = {"M1": 3992, "M2": 1806, "M3": 2721, "P6": 3361}

# Donors named by the pinned source for its source-defined Loupe category
# import (a2 lines 174--196).  A deposited category file exists for only a
# subset; missing files are untestable, never reconstructed from geometry.
CODED_ORDINAL_DONORS = ["P6", "P7", "P14", "P18", "P17", "M1"]

MATLAB_IMPORT_ORDER = [
    "P2", "P3", "P6", "P7", "P14", "P17", "P18", "P21",
    "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8",
]

# Explicit assignments in a2_visium_12_patients_initial_QCs_for_github.m,
# lines 82--97.  Despite its historical filename and README wording, that
# script declares and reorders all 16 human samples.
MATLAB_FINAL_ORDER = [
    "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M1",
    "P7", "P6", "P14", "P18", "P17", "P2", "P21", "P3",
]

EXPECTED_SOURCE_FILES = {
    "human_samples_metadata.xlsx": {
        "bytes": 21_271,
        "md5": "ff13c3f429d1a3f84a6035a03e72639f",
    },
    "v.mat": {
        "bytes": 3_996_202_856,
        "md5": "73f2ae74d2984363511d063af2873b0a",
    },
    "zon_struct_all_full.mat": {
        "bytes": 203_729_828,
        "md5": "8dc2e38d58c84e16a5143c9a60d02146",
    },
}


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def shape_string(obj: h5py.Dataset | h5py.Group) -> str:
    shape = getattr(obj, "shape", None)
    return "x".join(str(value) for value in shape) if shape is not None else ""


def decode_matlab_char(dataset: h5py.Dataset) -> str:
    values = np.asarray(dataset[...]).reshape(-1, order="F")
    return "".join(chr(int(value)) for value in values if int(value) != 0)


def decode_matlab_cellstr(dataset: h5py.Dataset, handle: h5py.File) -> list[str]:
    strings: list[str] = []
    for reference in np.asarray(dataset[...]).reshape(-1, order="F"):
        if not reference:
            strings.append("")
            continue
        target = handle[reference]
        if not isinstance(target, h5py.Dataset):
            raise TypeError(f"Expected MATLAB char dataset, found {target.name}")
        strings.append(decode_matlab_char(target))
    return strings


def raw_object_reference_addresses(dataset: h5py.Dataset) -> np.ndarray:
    """Read MATLAB cell references as object-header addresses without dereferencing.

    h5py's high-level object-reference iterator resolves each target through the
    HDF5 metadata index.  MATLAB v7.3 files can contain more than a million
    referenced character objects, making that random lookup path extremely
    expensive.  HDF5's reference memory type exposes the exact on-disk object
    addresses directly; a subsequent sequential object traversal can then
    decode only the requested targets.
    """
    if h5py.check_dtype(ref=dataset.dtype) is not h5py.Reference:
        raise TypeError(f"Expected object-reference dataset, found {dataset.dtype}")
    addresses = np.empty(dataset.shape, dtype=np.uint64)
    dataset.id.read(
        h5py.h5s.ALL,
        h5py.h5s.ALL,
        addresses,
        mtype=h5py.h5t.STD_REF_OBJ,
    )
    return addresses.reshape(-1, order="F")


def decode_matlab_cellstr_bulk(
    datasets: dict[int, h5py.Dataset], handle: h5py.File
) -> dict[int, list[str]]:
    """Decode several MATLAB cell-string arrays in one sequential traversal."""
    address_vectors = {
        key: raw_object_reference_addresses(dataset)
        for key, dataset in datasets.items()
    }
    wanted = {
        int(address)
        for addresses in address_vectors.values()
        for address in addresses
        if int(address) != 0
    }
    resolved: dict[int, str] = {}

    def visitor(_name: str, obj: h5py.Dataset | h5py.Group) -> bool | None:
        if not isinstance(obj, h5py.Dataset):
            return None
        address = int(h5py.h5o.get_info(obj.id).addr)
        if address not in wanted:
            return None
        resolved[address] = decode_matlab_char(obj)
        return True if len(resolved) == len(wanted) else None

    if wanted:
        handle.visititems(visitor)
    missing = sorted(wanted.difference(resolved))
    if missing:
        examples = ";".join(str(value) for value in missing[:5])
        raise KeyError(
            f"Could not resolve {len(missing)} MATLAB cell-string references; "
            f"example object addresses: {examples}"
        )
    return {
        key: [resolved[int(address)] if int(address) else "" for address in addresses]
        for key, addresses in address_vectors.items()
    }


def numeric_vector(dataset: h5py.Dataset) -> np.ndarray:
    return np.asarray(dataset[...], dtype=float).reshape(-1, order="F")


def read_integrity_manifest(path: Path, data_dir: Path) -> tuple[bool, list[dict[str, Any]]]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows: list[dict[str, Any]] = []
    all_pass = True
    for filename, expected in EXPECTED_SOURCE_FILES.items():
        matches = frame.loc[frame["file"] == filename]
        manifest_present = len(matches) == 1
        observed = matches.iloc[0].to_dict() if manifest_present else {}
        source_path = data_dir / filename
        size_match = source_path.is_file() and source_path.stat().st_size == expected["bytes"]
        checksum_match = (
            manifest_present
            and observed.get("expected_md5") == expected["md5"]
            and observed.get("observed_md5") == expected["md5"]
        )
        passed = bool(manifest_present and size_match and checksum_match)
        all_pass = all_pass and passed
        rows.append(
            {
                "filename": filename,
                "manifest_present": manifest_present,
                "observed_bytes_now": source_path.stat().st_size if source_path.exists() else "",
                "expected_bytes": expected["bytes"],
                "expected_md5": expected["md5"],
                "manifest_observed_md5": observed.get("observed_md5", ""),
                "manifest_sha256": observed.get("sha256", ""),
                "integrity_pass": passed,
                "verification_scope": "current_byte_size_plus_upstream_checksum_verified_manifest",
            }
        )
    return all_pass, rows


def source_code_audit(source_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    import_path = source_dir / "Matlab_functions/import_visium_data_funcion_for_github.m"
    a1_path = source_dir / "Matlab_scripts/a1_import_data_raw_and_filtered_for_github.m"
    a2_path = source_dir / "Matlab_scripts/a2_visium_12_patients_initial_QCs_for_github.m"
    a3_path = source_dir / "Matlab_scripts/a3_visium_12_patients_zonation_reconstruction_for_github.m"
    source_files = list(source_dir.rglob("*.m")) + list(source_dir.rglob("*.groovy"))
    import_text = import_path.read_text(errors="replace")
    a1_text = a1_path.read_text(errors="replace")
    a2_text = a2_path.read_text(errors="replace")
    a3_text = a3_path.read_text(errors="replace")
    all_source_text = "\n".join(path.read_text(errors="replace") for path in source_files)
    ordinal_join_tokens = (
        "lipid_zones.csv",
        "strcmpi(t.spot_name,lipid_spots(j))",
        "v{i}.ind_lipid_spots=ind_lipid_spots",
    )
    authoritative_ordinal_join = all(token in a2_text for token in ordinal_join_tokens)

    checks = [
        {
            "concept": "donor_identity",
            "source_path": str(a1_path.relative_to(source_dir)),
            "source_lines": "21,27-33",
            "deposited_field": "v[i].patient",
            "code_support": "pati list assigns the patient field explicitly",
            "authoritative_for_spot_join": True,
        },
        {
            "concept": "barcode_expression_axis",
            "source_path": str(import_path.relative_to(source_dir)),
            "source_lines": "31-36",
            "deposited_field": "v[i].spot_name; v[i].mat",
            "code_support": "count-table column names become spot_name and the same columns become mat",
            "authoritative_for_spot_join": all(
                token in import_text
                for token in ("t.spot_name=T.Properties.VariableNames(2:end)", "t.mat=table2array(T(:,2:end))")
            ),
        },
        {
            "concept": "barcode_coordinates",
            "source_path": str(import_path.relative_to(source_dir)),
            "source_lines": "21-25,47-53",
            "deposited_field": "v[i].spot_name; v[i].coor",
            "code_support": "stored spot names are normalized and matched to tissue-position barcodes before coordinate assignment",
            "authoritative_for_spot_join": "strcmpi(barcode_space,str)" in import_text,
        },
        {
            "concept": "synchronized_qc_filter",
            "source_path": str(a2_path.relative_to(source_dir)),
            "source_lines": "61-73",
            "deposited_field": "v[i].mat; v[i].spot_name; v[i].coor",
            "code_support": "the identical ind_out vector deletes expression columns, spot names, and coordinate rows",
            "authoritative_for_spot_join": all(
                token in a2_text
                for token in ("t.mat(:,ind_out)=[]", "t.spot_name(ind_out)=[]", "t.coor(ind_out,:)=[]")
            ),
        },
        {
            "concept": "continuous_zonation",
            "source_path": str(a3_path.relative_to(source_dir)),
            "source_lines": "69-91",
            "deposited_field": "v[i].eta",
            "code_support": "eta is computed from the expression columns already aligned to spot_name",
            "authoritative_for_spot_join": "v{i}.eta=sum_pp./(sum_pp+sum_pc)" in a3_text,
        },
        {
            "concept": "categorical_zonation",
            "source_path": str(a3_path.relative_to(source_dir)),
            "source_lines": "159-176",
            "deposited_field": "v[i].zon_struct.zone_index",
            "code_support": "zone_index is derived from eta and median-filtered without reordering spots",
            "authoritative_for_spot_join": all(
                token in a3_text
                for token in ("extract_zonation_for_github", "median_zone_filter_for_github", "zone_index_med")
            ),
        },
        {
            "concept": "quantitative_lipid",
            "source_path": "NOT_PRESENT_IN_PINNED_SOURCE",
            "source_lines": "",
            "deposited_field": "v[i].lipid_percentage; v[i].tissue_percentage",
            "code_support": "paper Methods describe QuPath measurement per spot, but the pinned code contains no construction or barcode-keyed import of these arrays",
            "authoritative_for_spot_join": "lipid_percentage" in all_source_text,
        },
        {
            "concept": "binary_lipid_loupe_category",
            "source_path": str(a2_path.relative_to(source_dir)),
            "source_lines": "174-196",
            "deposited_field": "v[i].ind_lipid_spots when present",
            "code_support": "source-defined lipid-zone barcodes are matched explicitly to spot_name before ind_lipid_spots is stored; this supports only a binary lipid-zone versus non-lipid-zone label",
            "authoritative_for_spot_join": authoritative_ordinal_join,
        },
    ]

    declaration = import_text.splitlines()[0]
    declared_name = declaration.split("=", 1)[1].split("(", 1)[0].strip()
    call_name = "import_visium_data_function_for_github"
    metadata_assignment_is_positional = "v{i}.metadata = metadata(i,:)" in a1_text
    loupe_helper_present = any(path.stem == "read_loupe_barcode" for path in source_files)
    diagnostics = [
        {
            "diagnostic_id": "YAK-CODE-01",
            "severity": "reproducibility_blocker",
            "observed": declared_name,
            "expected": call_name,
            "pass": declared_name == call_name,
            "detail": "Pinned filename/declaration spell 'funcion'; a1 calls 'function'. The source scripts are not runnable as named.",
        },
        {
            "diagnostic_id": "YAK-META-01",
            "severity": "metadata_warning",
            "observed": "positional metadata(i,:) assignment" if metadata_assignment_is_positional else "not_found",
            "expected": "join metadata by Patient Code",
            "pass": not metadata_assignment_is_positional,
            "detail": "The deposited workbook contains 17 rows including P4 while the import list has 16 samples; external metadata must be rejoined by Patient Code.",
        },
        {
            "diagnostic_id": "YAK-CODE-02",
            "severity": "reproducibility_warning",
            "observed": "present" if loupe_helper_present else "read_loupe_barcode helper not deposited",
            "expected": "deposited helper used by a2 lines 180-181",
            "pass": loupe_helper_present,
            "detail": "The audit does not execute the missing helper: it independently parses the public CSV and verifies the stored positive-index barcodes, so this warning does not by itself invalidate a matching stored ordinal join.",
        },
        {
            "diagnostic_id": "YAK-LIPID-01",
            "severity": "continuous_lipid_gate_blocker",
            "observed": "no lipid_percentage construction/import in pinned source",
            "expected": "barcode-keyed QuPath export or deposited code proving row alignment",
            "pass": "lipid_percentage" in all_source_text,
            "detail": "Equal vector lengths alone are a positional reconstruction and do not satisfy the authoritative lipid join gate.",
        },
        {
            "diagnostic_id": "YAK-ORD-01",
            "severity": "ordinal_source_check",
            "observed": "explicit barcode matching into ind_lipid_spots" if authoritative_ordinal_join else "required barcode-matching tokens absent",
            "expected": "lipid_zones.csv -> spot_name exact match -> ind_lipid_spots",
            "pass": authoritative_ordinal_join,
            "detail": "This source path authorizes a binary category only; it does not authorize a continuous gradient or ordered dose-response interpretation.",
        },
    ]
    return checks, diagnostics


def classify_lipid_scale(values: np.ndarray) -> str:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return "no_finite_values"
    minimum = float(np.min(finite))
    maximum = float(np.max(finite))
    if minimum >= 0 and maximum <= 1:
        return "stored_0_to_1_fraction_paper_reports_percent_spot_area"
    if minimum >= 0 and maximum <= 100:
        return "stored_0_to_100_compatible_with_percent_spot_area"
    return "outside_expected_percentage_range"


def inspect_v_object(
    path: Path, authoritative_ordinal_source: bool
) -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, Any]]
]:
    schema_rows: list[dict[str, Any]] = []
    join_rows: list[dict[str, Any]] = []
    donor_summary: dict[str, dict[str, Any]] = {}

    with h5py.File(path, "r") as handle:
        root = handle["v"]
        references = np.asarray(root[...]).reshape(-1, order="F")
        if len(references) != 16:
            raise ValueError(f"Expected 16 v cells, found {len(references)}")

        groups: dict[int, h5py.Group] = {}
        for cell_index, reference in enumerate(references, start=1):
            group = handle[reference]
            if not isinstance(group, h5py.Group):
                raise TypeError(f"v cell {cell_index} does not reference a struct group")
            groups[cell_index] = group
        barcodes_by_cell = decode_matlab_cellstr_bulk(
            {index: group["spot_name"] for index, group in groups.items()},
            handle,
        )

        for cell_index, group in groups.items():
            patient = decode_matlab_char(group["patient"])
            feature = decode_matlab_char(group["main_feature"])
            barcodes = barcodes_by_cell[cell_index]
            n_spots = len(barcodes)
            n_genes = int(np.prod(group["gene_name"].shape))
            coordinates_h5 = np.asarray(group["coor"][...], dtype=float)
            coordinates = coordinates_h5.T
            eta = numeric_vector(group["eta"])
            zone = numeric_vector(group["zon_struct/zone_index"])
            mat = group["mat"]

            lipid = numeric_vector(group["lipid_percentage"]) if "lipid_percentage" in group else None
            tissue = numeric_vector(group["tissue_percentage"]) if "tissue_percentage" in group else None
            ordinal_indices_raw = (
                numeric_vector(group["ind_lipid_spots"])
                if "ind_lipid_spots" in group else None
            )
            morpho_path = (
                decode_matlab_char(group["lipid_morpho_file_path"])
                if "lipid_morpho_file_path" in group
                else ""
            )

            matrix_axis_matches = mat.shape == (n_spots, n_genes)
            coordinate_axis_matches = coordinates.shape == (n_spots, 2)
            eta_axis_matches = eta.size == n_spots
            zone_axis_matches = zone.size == n_spots
            barcode_unique = len(set(barcodes)) == n_spots and all(barcodes)
            final_order_matches = patient == MATLAB_FINAL_ORDER[cell_index - 1]
            ordinal_field_present = ordinal_indices_raw is not None
            ordinal_indices_valid = bool(
                ordinal_field_present
                and np.isfinite(ordinal_indices_raw).all()
                and np.allclose(ordinal_indices_raw, np.round(ordinal_indices_raw))
                and len(np.unique(ordinal_indices_raw)) == len(ordinal_indices_raw)
                and np.all(ordinal_indices_raw >= 1)
                and np.all(ordinal_indices_raw <= n_spots)
            )
            ordinal_indices = (
                {int(value) for value in ordinal_indices_raw}
                if ordinal_indices_valid else set()
            )

            for field_name in (
                "patient", "main_feature", "spot_name", "gene_name", "coor", "mat",
                "eta", "zon_struct/zone_index", "lipid_percentage", "tissue_percentage",
                "lipid_morpho_file_path", "ind_lipid_spots",
            ):
                if field_name == "zon_struct/zone_index":
                    obj = group["zon_struct/zone_index"]
                elif field_name in group:
                    obj = group[field_name]
                else:
                    obj = None
                schema_rows.append(
                    {
                        "file": path.name,
                        "v_cell_index_1based": cell_index,
                        "patient": patient,
                        "struct_hdf5_path": group.name,
                        "field": field_name,
                        "hdf5_path": obj.name if obj is not None else "",
                        "present": obj is not None,
                        "object_type": type(obj).__name__ if obj is not None else "",
                        "shape": shape_string(obj) if obj is not None else "",
                        "dtype": str(getattr(obj, "dtype", "")) if obj is not None else "",
                    }
                )

            structural_complete_count = 0
            structural_lipid_complete_count = 0
            ordinal_join_complete_count = 0
            for spot_index, barcode in enumerate(barcodes, start=1):
                coordinate_complete = (
                    coordinate_axis_matches
                    and np.isfinite(coordinates[spot_index - 1, :]).all()
                )
                eta_complete = eta_axis_matches and np.isfinite(eta[spot_index - 1])
                zone_complete = zone_axis_matches and np.isfinite(zone[spot_index - 1])
                structural_complete = bool(
                    matrix_axis_matches
                    and barcode_unique
                    and coordinate_complete
                    and eta_complete
                    and zone_complete
                )
                structural_complete_count += int(structural_complete)

                lipid_value = (
                    float(lipid[spot_index - 1])
                    if lipid is not None and lipid.size == n_spots
                    else np.nan
                )
                tissue_value = (
                    float(tissue[spot_index - 1])
                    if tissue is not None and tissue.size == n_spots
                    else np.nan
                )
                lipid_structural = bool(structural_complete and np.isfinite(lipid_value))
                if patient in PAPER_LIPID_SPOTS:
                    structural_lipid_complete_count += int(lipid_structural)
                ordinal_label_available = bool(ordinal_field_present and ordinal_indices_valid)
                ordinal_lipid_zone = (
                    spot_index in ordinal_indices if ordinal_label_available else ""
                )
                ordinal_join_complete = bool(
                    ordinal_label_available
                    and structural_complete
                    and authoritative_ordinal_source
                )
                ordinal_join_complete_count += int(ordinal_join_complete)

                join_rows.append(
                    {
                        "dataset": "Yakubovsky_2026",
                        "v_cell_index_1based": cell_index,
                        "patient": patient,
                        "main_feature": feature,
                        "spot_index_1based": spot_index,
                        "barcode": barcode,
                        "composite_spot_id": f"{patient}:{barcode}",
                        "x_coordinate": coordinates[spot_index - 1, 0] if coordinate_axis_matches else "",
                        "y_coordinate": coordinates[spot_index - 1, 1] if coordinate_axis_matches else "",
                        "eta": eta[spot_index - 1] if eta_axis_matches else "",
                        "zone_index": zone[spot_index - 1] if zone_axis_matches else "",
                        "paper_primary_lipid_donor": patient in PAPER_LIPID_SPOTS,
                        "lipid_percentage": lipid_value if np.isfinite(lipid_value) else "",
                        "tissue_percentage": tissue_value if np.isfinite(tissue_value) else "",
                        "expression_axis_present": matrix_axis_matches,
                        "coordinate_complete": coordinate_complete,
                        "zonation_complete": eta_complete and zone_complete,
                        "structural_expression_coordinate_zonation_join": structural_complete,
                        "structural_lipid_alignment": lipid_structural,
                        "authoritative_lipid_join": False,
                        "lipid_join_basis": (
                            "same_struct_row_order_only; quantitative_lipid_generation_code_absent"
                            if lipid is not None
                            else "no_quantitative_lipid_field"
                        ),
                        "ordinal_lipid_label_available": ordinal_label_available,
                        "ordinal_lipid_zone": ordinal_lipid_zone,
                        "source_defined_binary_lipid_class": (
                            "lipid_zone" if ordinal_lipid_zone is True
                            else "non_lipid_zone" if ordinal_lipid_zone is False
                            else ""
                        ),
                        "authoritative_ordinal_lipid_join": ordinal_join_complete,
                    }
                )

            finite_lipid = lipid[np.isfinite(lipid)] if lipid is not None else np.array([])
            lipid_non_degenerate = bool(
                finite_lipid.size >= 2 and np.unique(finite_lipid).size > 1
                and float(np.nanstd(finite_lipid)) > 0
            )
            n_ordinal_positive = len(ordinal_indices) if ordinal_indices_valid else 0
            ordinal_non_degenerate = bool(
                ordinal_indices_valid and 0 < n_ordinal_positive < n_spots
            )
            donor_summary[patient] = {
                "v_cell_index_1based": cell_index,
                "patient": patient,
                "main_feature": feature,
                "struct_hdf5_path": group.name,
                "expected_final_order_patient": MATLAB_FINAL_ORDER[cell_index - 1],
                "final_order_matches_source_code": final_order_matches,
                "n_spots": n_spots,
                "n_genes": n_genes,
                "n_unique_barcodes": len(set(barcodes)),
                "matrix_hdf5_shape": shape_string(mat),
                "coordinate_hdf5_shape": shape_string(group["coor"]),
                "eta_hdf5_shape": shape_string(group["eta"]),
                "zone_hdf5_shape": shape_string(group["zon_struct/zone_index"]),
                "matrix_spot_axis_matches": matrix_axis_matches,
                "coordinate_axis_matches": coordinate_axis_matches,
                "eta_axis_matches": eta_axis_matches,
                "zone_axis_matches": zone_axis_matches,
                "barcode_unique_within_donor": barcode_unique,
                "n_structural_expression_coordinate_zonation_join": structural_complete_count,
                "structural_expression_coordinate_zonation_join_fraction": (
                    structural_complete_count / n_spots if n_spots else np.nan
                ),
                "paper_primary_lipid_donor": patient in PAPER_LIPID_SPOTS,
                "paper_expected_lipid_spots": PAPER_LIPID_SPOTS.get(patient, ""),
                "paper_spot_count_matches": (
                    n_spots == PAPER_LIPID_SPOTS[patient]
                    if patient in PAPER_LIPID_SPOTS else "not_applicable"
                ),
                "lipid_field_present": lipid is not None,
                "tissue_field_present": tissue is not None,
                "lipid_morpho_file_path": morpho_path,
                "n_finite_lipid_values": int(finite_lipid.size),
                "n_structural_lipid_aligned_spots": structural_lipid_complete_count,
                "structural_lipid_join_fraction": (
                    structural_lipid_complete_count / int(finite_lipid.size)
                    if patient in PAPER_LIPID_SPOTS and finite_lipid.size
                    else "not_applicable"
                ),
                "lipid_measurement_coverage_fraction": (
                    int(finite_lipid.size) / n_spots
                    if patient in PAPER_LIPID_SPOTS and n_spots
                    else "not_applicable"
                ),
                "lipid_min": float(np.min(finite_lipid)) if finite_lipid.size else "",
                "lipid_max": float(np.max(finite_lipid)) if finite_lipid.size else "",
                "lipid_n_unique": int(np.unique(finite_lipid).size) if finite_lipid.size else 0,
                "lipid_non_degenerate": lipid_non_degenerate,
                "lipid_scale_audit": classify_lipid_scale(lipid) if lipid is not None else "not_applicable",
                "authoritative_lipid_join_verified": False,
                "ordinal_candidate_donor": patient in CODED_ORDINAL_DONORS,
                "ind_lipid_spots_field_present": ordinal_field_present,
                "ordinal_indices_integer_unique_in_range": ordinal_indices_valid,
                "n_lipid_zone_spots": n_ordinal_positive if ordinal_field_present else "",
                "n_non_lipid_zone_spots": (
                    n_spots - n_ordinal_positive if ordinal_field_present else ""
                ),
                "n_valid_zonated_spots": structural_complete_count,
                "ordinal_category_non_degenerate": ordinal_non_degenerate,
                "ordinal_join_fraction": (
                    ordinal_join_complete_count / n_spots
                    if ordinal_field_present and n_spots else "not_applicable"
                ),
                "authoritative_ordinal_join_structurally_verified": bool(
                    ordinal_field_present
                    and ordinal_indices_valid
                    and ordinal_join_complete_count == n_spots
                    and authoritative_ordinal_source
                ),
            }

    return schema_rows, join_rows, donor_summary


def normalize_loupe_barcode(value: str) -> str:
    """Normalize only the delimiter change made by the pinned Visium import."""
    return str(value).strip().upper().replace("_", "-")


def audit_ordinal_loupe_join(
    source_dir: Path,
    join_rows: list[dict[str, Any]],
    donor_summary: dict[str, dict[str, Any]],
    authoritative_ordinal_source: bool,
) -> list[dict[str, Any]]:
    """Independently reconcile stored positive indices to public Loupe barcodes."""
    rows_by_patient: dict[str, list[dict[str, Any]]] = {}
    for row in join_rows:
        rows_by_patient.setdefault(str(row["patient"]), []).append(row)

    audit_rows: list[dict[str, Any]] = []
    for patient in CODED_ORDINAL_DONORS:
        public_path = source_dir / "Loupe_categories" / patient / "lipid_zones.csv"
        public_barcodes: list[str] = []
        if public_path.is_file():
            with public_path.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle)
                if not reader.fieldnames:
                    raise ValueError(f"No header in {public_path}")
                barcode_field = next(
                    (field for field in reader.fieldnames if field.strip().lower() == "barcode"),
                    reader.fieldnames[0],
                )
                public_barcodes = [
                    normalize_loupe_barcode(row.get(barcode_field, ""))
                    for row in reader
                    if str(row.get(barcode_field, "")).strip()
                ]

        public_set = set(public_barcodes)
        donor_rows = rows_by_patient.get(patient, [])
        stored_positive = {
            normalize_loupe_barcode(row["barcode"])
            for row in donor_rows
            if row.get("ordinal_lipid_zone") is True
        }
        retained_barcodes = {
            normalize_loupe_barcode(row["barcode"]) for row in donor_rows
        }
        stored_in_public = stored_positive & public_set
        public_retained = public_set & retained_barcodes
        stored_to_public_fraction = (
            len(stored_in_public) / len(stored_positive) if stored_positive else 0.0
        )
        public_retention_fraction = (
            len(public_retained) / len(public_set) if public_set else 0.0
        )
        summary = donor_summary.get(patient, {})
        independent_join_pass = bool(
            authoritative_ordinal_source
            and public_path.is_file()
            and bool(summary.get("ind_lipid_spots_field_present", False))
            and bool(summary.get("ordinal_indices_integer_unique_in_range", False))
            and len(stored_positive) > 0
            and stored_to_public_fraction >= 0.90
        )
        audit_rows.append(
            {
                "patient": patient,
                "coded_in_a2": True,
                "public_loupe_file_present": public_path.is_file(),
                "public_file": (
                    str(public_path.relative_to(source_dir)) if public_path.is_file() else ""
                ),
                "n_public_positive_barcodes": len(public_set),
                "public_positive_barcodes_unique": len(public_set) == len(public_barcodes),
                "ind_lipid_spots_field_present": summary.get(
                    "ind_lipid_spots_field_present", False
                ),
                "n_stored_positive_indices": len(stored_positive),
                "n_stored_positive_barcodes_in_public_csv": len(stored_in_public),
                "stored_positive_to_public_csv_fraction": stored_to_public_fraction,
                "n_public_positive_barcodes_retained": len(public_retained),
                "public_positive_retention_fraction": public_retention_fraction,
                "n_stored_positive_not_in_public_csv": len(stored_positive - public_set),
                "stored_positive_not_in_public_csv_examples": ";".join(
                    sorted(stored_positive - public_set)[:5]
                ),
                "n_public_positive_not_retained": len(public_set - retained_barcodes),
                "public_positive_not_retained_examples": ";".join(
                    sorted(public_set - retained_barcodes)[:5]
                ),
                "source_code_lines": "a2_visium_12_patients_initial_QCs_for_github.m:174-196",
                "comparison_normalization": "uppercase_and_underscore_to_hyphen_only",
                "independent_ordinal_barcode_join_pass": independent_join_pass,
            }
        )
    return audit_rows


def inspect_aggregate_zonation(path: Path, final_order: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with h5py.File(path, "r") as handle:
        root = handle["zon_struct_all_full"]
        references = np.asarray(root[...]).reshape(-1, order="F")
        for cell_index, reference in enumerate(references, start=1):
            group = handle[reference]
            indices = numeric_vector(group["orig_ind_in_v"]).astype(int).tolist()
            patients = [
                final_order[index - 1] for index in indices
                if 1 <= index <= len(final_order)
            ]
            rows.append(
                {
                    "cell_index_1based": cell_index,
                    "struct_hdf5_path": group.name,
                    "category": decode_matlab_char(group["category"]),
                    "normalization_method": decode_matlab_char(group["normalization_method"]),
                    "orig_ind_in_v_1based": ";".join(str(value) for value in indices),
                    "resolved_patients": ";".join(patients),
                    "n_genes": int(np.prod(group["gene_name"].shape)),
                    "mn_shape": shape_string(group["mn"]),
                    "has_spot_barcodes": "spot_name" in group,
                    "has_coordinates": "coor" in group,
                    "has_lipid": "lipid_percentage" in group,
                    "join_role": "aggregate_zonation_only_not_a_spot_level_join_bridge",
                }
            )
    return rows


def git_value(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], text=True, stderr=subprocess.STDOUT
    ).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--source-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    source_dir = args.source_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    started = utc_now()

    manifest_path = output_dir / "lean_download_manifest.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing checksum-verified lean manifest: {manifest_path}")
    integrity_pass, integrity_rows = read_integrity_manifest(manifest_path, data_dir)

    code_map, code_diagnostics = source_code_audit(source_dir)
    ordinal_code_map = next(
        row for row in code_map if row["concept"] == "binary_lipid_loupe_category"
    )
    authoritative_ordinal_source = bool(
        ordinal_code_map["authoritative_for_spot_join"]
    )
    schema_rows, join_rows, donor_summary = inspect_v_object(
        data_dir / "v.mat", authoritative_ordinal_source
    )
    ordinal_audit_rows = audit_ordinal_loupe_join(
        source_dir, join_rows, donor_summary, authoritative_ordinal_source
    )
    ordinal_audit_by_patient = {
        str(row["patient"]): row for row in ordinal_audit_rows
    }
    for patient, summary in donor_summary.items():
        audit = ordinal_audit_by_patient.get(patient, {})
        summary.update(
            {
                "ordinal_public_loupe_file_present": audit.get(
                    "public_loupe_file_present", False
                ),
                "ordinal_n_public_positive_barcodes": audit.get(
                    "n_public_positive_barcodes", ""
                ),
                "ordinal_n_stored_positive_barcodes_in_public_csv": audit.get(
                    "n_stored_positive_barcodes_in_public_csv", ""
                ),
                "ordinal_stored_positive_to_public_csv_fraction": audit.get(
                    "stored_positive_to_public_csv_fraction", "not_applicable"
                ),
                "ordinal_public_positive_retention_fraction": audit.get(
                    "public_positive_retention_fraction", "not_applicable"
                ),
                "ordinal_independent_barcode_join_pass": audit.get(
                    "independent_ordinal_barcode_join_pass", False
                ),
            }
        )
    aggregate_rows = inspect_aggregate_zonation(
        data_dir / "zon_struct_all_full.mat", MATLAB_FINAL_ORDER
    )

    clinical = pd.read_excel(data_dir / "human_samples_metadata.xlsx", sheet_name="Clinical Data")
    metadata_codes = clinical["Patient Code"].astype(str).tolist()
    metadata_positions = {patient: index for index, patient in enumerate(metadata_codes, start=1)}
    object_codes = list(donor_summary)
    object_positions = {
        patient: int(summary["v_cell_index_1based"])
        for patient, summary in donor_summary.items()
    }
    all_patients = sorted(set(metadata_codes) | set(object_codes) | set(MATLAB_IMPORT_ORDER))

    cohort_rows: list[dict[str, Any]] = []
    for patient in all_patients:
        summary = donor_summary.get(patient, {})
        import_position = (
            MATLAB_IMPORT_ORDER.index(patient) + 1 if patient in MATLAB_IMPORT_ORDER else ""
        )
        metadata_position = metadata_positions.get(patient, "")
        cohort_rows.append(
            {
                "patient": patient,
                "donor_class": "live_donor" if patient.startswith("M") else "adjacent_normal",
                "in_clinical_metadata": patient in metadata_positions,
                "clinical_metadata_row_1based": metadata_position,
                "in_matlab_import_list": patient in MATLAB_IMPORT_ORDER,
                "matlab_import_position_1based": import_position,
                "positional_metadata_assignment_would_match": (
                    metadata_position == import_position
                    if metadata_position != "" and import_position != "" else False
                ),
                "in_v_object": patient in donor_summary,
                "v_cell_index_1based": object_positions.get(patient, ""),
                "v_order_matches_a2_reorder": summary.get("final_order_matches_source_code", ""),
                "paper_primary_lipid_donor": patient in PAPER_LIPID_SPOTS,
                "has_quantitative_lipid_field": summary.get("lipid_field_present", False),
                "observed_spots": summary.get("n_spots", ""),
                "paper_expected_lipid_spots": PAPER_LIPID_SPOTS.get(patient, ""),
                "source_role": (
                    "paper_primary_lipid_and_visium"
                    if patient in PAPER_LIPID_SPOTS
                    else "visium_expression_sample" if patient in donor_summary
                    else "metadata_only"
                ),
                "resolution": (
                    "P4_metadata_only_not_in_16_sample_import"
                    if patient == "P4" else "included_in_v_object"
                ),
            }
        )

    primary_summaries = [
        donor_summary[patient]
        for patient in PAPER_LIPID_SPOTS
        if patient in donor_summary
    ]
    primary_donors_present = all(patient in donor_summary for patient in PAPER_LIPID_SPOTS)
    primary_counts_match = primary_donors_present and all(
        int(summary["n_spots"]) == int(summary["paper_expected_lipid_spots"])
        for summary in primary_summaries
    )
    n_primary_counts_matching = sum(
        int(summary["n_spots"]) == int(summary["paper_expected_lipid_spots"])
        for summary in primary_summaries
    )
    primary_finite_lipid_total = sum(
        int(summary["n_finite_lipid_values"]) for summary in primary_summaries
    )
    primary_lipid_spot_total = sum(
        int(summary["n_spots"]) for summary in primary_summaries
    )
    lipid_measurement_coverage_fraction = (
        primary_finite_lipid_total / primary_lipid_spot_total
        if primary_lipid_spot_total > 0 else 0.0
    )
    structural_join_fraction = (
        sum(int(summary["n_structural_lipid_aligned_spots"]) for summary in primary_summaries)
        / primary_finite_lipid_total
        if primary_finite_lipid_total > 0 else 0.0
    )
    continuous_passing_donors: list[str] = []
    for summary in primary_summaries:
        donor_pass = bool(
            int(summary["n_spots"]) == int(summary["paper_expected_lipid_spots"])
            and int(summary["n_structural_lipid_aligned_spots"]) >= 100
            and bool(summary["lipid_non_degenerate"])
            and float(summary["structural_lipid_join_fraction"]) >= 0.90
        )
        summary["continuous_source_gate_donor_pass"] = donor_pass
        summary["continuous_source_gate_disposition"] = (
            "pass_structural_only_authoritative_join_still_required"
            if donor_pass else "fail_count_join_spots_or_variation"
        )
        if donor_pass:
            continuous_passing_donors.append(str(summary["patient"]))
    donor_thresholds_pass = len(continuous_passing_donors) >= 3
    # The prespecified source gate requires at least 90% authoritative joining
    # among lipid-tested barcodes; it does not require every spot in every
    # non-lipid reference donor to have finite zonation.  Preserve the complete
    # 16-donor object and its aligned axes, while allowing a donor's explicitly
    # invalid spots to be excluded when at least 90% and at least 100 valid
    # spots remain.  The ordinal donor gate below independently applies the
    # same threshold to each actually tested donor.
    expression_zonation_pass = (
        len(donor_summary) == 16
        and object_codes == MATLAB_FINAL_ORDER
        and all(
            bool(summary["matrix_spot_axis_matches"])
            and bool(summary["coordinate_axis_matches"])
            and bool(summary["eta_axis_matches"])
            and bool(summary["zone_axis_matches"])
            and bool(summary["barcode_unique_within_donor"])
            and int(summary["n_structural_expression_coordinate_zonation_join"])
            >= 100
            and float(
                summary[
                    "structural_expression_coordinate_zonation_join_fraction"
                ]
            )
            >= 0.90
            for summary in donor_summary.values()
        )
    )
    minimum_expression_zonation_join_fraction = min(
        float(
            summary["structural_expression_coordinate_zonation_join_fraction"]
        )
        for summary in donor_summary.values()
    )
    quantitative_lipid_code_map = next(
        row for row in code_map if row["concept"] == "quantitative_lipid"
    )
    authoritative_lipid_join = bool(
        quantitative_lipid_code_map["authoritative_for_spot_join"]
    )
    continuous_lipid_gate_pass = bool(
        integrity_pass
        and expression_zonation_pass
        and donor_thresholds_pass
        and authoritative_lipid_join
    )
    ordinal_passing_donors: list[str] = []
    for patient in CODED_ORDINAL_DONORS:
        summary = donor_summary.get(patient, {})
        audit = ordinal_audit_by_patient.get(patient, {})
        donor_pass = bool(
            summary
            and audit.get("public_loupe_file_present", False)
            and summary.get("ind_lipid_spots_field_present", False)
            and summary.get("ordinal_indices_integer_unique_in_range", False)
            and int(summary.get("n_spots", 0)) >= 100
            and int(summary.get("n_valid_zonated_spots", 0)) >= 100
            and summary.get("ordinal_category_non_degenerate", False)
            and float(summary.get("ordinal_join_fraction", 0.0)) >= 0.90
            and float(audit.get("stored_positive_to_public_csv_fraction", 0.0)) >= 0.90
            and audit.get("independent_ordinal_barcode_join_pass", False)
        )
        summary["ordinal_source_gate_donor_pass"] = donor_pass
        if donor_pass:
            summary["ordinal_source_gate_disposition"] = "pass"
            ordinal_passing_donors.append(patient)
        elif not audit.get("public_loupe_file_present", False):
            summary["ordinal_source_gate_disposition"] = (
                "untestable_no_public_loupe_category_file"
            )
        else:
            summary["ordinal_source_gate_disposition"] = "fail_source_or_join_criterion"

    for patient, summary in donor_summary.items():
        if patient not in CODED_ORDINAL_DONORS:
            summary["ordinal_source_gate_donor_pass"] = False
            summary["ordinal_source_gate_disposition"] = "not_source_coded_for_ordinal_lipid"

    ordinal_lipid_gate_pass = bool(
        integrity_pass
        and expression_zonation_pass
        and authoritative_ordinal_source
        and len(ordinal_passing_donors) >= 3
    )
    terminal_verdict = (
        "pass_continuous_lipid" if continuous_lipid_gate_pass
        else "pass_ordinal_lipid" if ordinal_lipid_gate_pass
        else "pass_zonation_only" if integrity_pass and expression_zonation_pass
        else "fail_unjoinable" if integrity_pass
        else "fail_integrity"
    )

    criteria = [
        {
            "criterion_id": "integrity",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": True,
            "observed": "3/3 lean files match current bytes and upstream checksum manifest",
            "pass": integrity_pass,
            "evidence": "lean_download_manifest.tsv plus current byte sizes",
        },
        {
            "criterion_id": "human_expression_donors",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": True,
            "observed": len(donor_summary),
            "pass": len(donor_summary) == 16 and object_codes == MATLAB_FINAL_ORDER,
            "evidence": "v.mat root v cell array and explicit a2 lines 82-97",
        },
        {
            "criterion_id": "expression_barcode_coordinate_zonation_join",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": True,
            "observed": (
                f"{len(donor_summary)}/16 donor structs have aligned axes; "
                "each retains at least 100 valid spots and at least 90% "
                "expression-barcode-coordinate-zonation joining; minimum "
                f"fraction={minimum_expression_zonation_join_fraction:.6f}"
            ),
            "pass": expression_zonation_pass,
            "evidence": (
                "yakubovsky_join_audit.tsv, donor schema summary, pinned source "
                "field map, and the prespecified >=90% source-gate threshold"
            ),
        },
        {
            "criterion_id": "paper_lipid_donors_and_counts",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": False,
            "observed": ";".join(
                f"{patient}:{donor_summary[patient]['n_spots']}"
                for patient in PAPER_LIPID_SPOTS if patient in donor_summary
            ),
            "pass": n_primary_counts_matching >= 3,
            "evidence": "Nature Methods/Extended Data Fig. 10a versus v.mat",
        },
        {
            "criterion_id": "structural_lipid_join_at_least_90_percent",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": False,
            "observed": structural_join_fraction,
            "pass": structural_join_fraction >= 0.90,
            "evidence": (
                "finite lipid observations with valid expression/coordinates/"
                "zonation divided by all finite lipid observations; finite/all-"
                "tissue coverage is reported separately"
            ),
        },
        {
            "criterion_id": "minimum_donors_spots_and_variation",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": False,
            "observed": (
                f"{len(continuous_passing_donors)}/{len(PAPER_LIPID_SPOTS)} "
                f"paper donors pass structurally: "
                f"{';'.join(continuous_passing_donors) or 'none'}"
            ),
            "pass": donor_thresholds_pass,
            "evidence": "yakubovsky_donor_schema_summary.tsv",
        },
        {
            "criterion_id": "lipid_unit_documented",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": False,
            "observed": "percent of Visium spot area from QuPath fat-droplet pixels",
            "pass": True,
            "evidence": f"Nature Methods, {PAPER_DOI}; deposited field lipid_percentage",
        },
        {
            "criterion_id": "authoritative_barcode_to_lipid_join",
            "required_for_continuous_lipid": True,
            "required_for_ordinal_lipid": False,
            "observed": "lipid vectors are positionally aligned in v.mat, but their construction/import is absent from the pinned source",
            "pass": authoritative_lipid_join,
            "evidence": "pinned-source search plus yakubovsky_source_field_map.tsv",
        },
        {
            "criterion_id": "source_defined_ordinal_barcode_join",
            "required_for_continuous_lipid": False,
            "required_for_ordinal_lipid": True,
            "observed": "lipid_zones.csv is matched to spot_name and stored as ind_lipid_spots",
            "pass": authoritative_ordinal_source,
            "evidence": "pinned a2 source lines 174-196",
        },
        {
            "criterion_id": "at_least_three_public_ordinal_donors",
            "required_for_continuous_lipid": False,
            "required_for_ordinal_lipid": True,
            "observed": f"{len(ordinal_passing_donors)}/6 coded donors pass: {';'.join(ordinal_passing_donors) or 'none'}",
            "pass": len(ordinal_passing_donors) >= 3,
            "evidence": "yakubovsky_ordinal_loupe_join_audit.tsv",
        },
        {
            "criterion_id": "ordinal_spot_count_and_nondegenerate_category",
            "required_for_continuous_lipid": False,
            "required_for_ordinal_lipid": True,
            "observed": ";".join(
                f"{patient}:{donor_summary[patient]['n_lipid_zone_spots']}/"
                f"{donor_summary[patient]['n_non_lipid_zone_spots']}/"
                f"{donor_summary[patient]['n_spots']}"
                for patient in ordinal_passing_donors
            ),
            "pass": len(ordinal_passing_donors) >= 3,
            "evidence": "stored ind_lipid_spots, valid 1-based indices, and retained spot universe",
        },
        {
            "criterion_id": "ordinal_expression_zonation_and_public_barcode_join_at_least_90_percent",
            "required_for_continuous_lipid": False,
            "required_for_ordinal_lipid": True,
            "observed": ";".join(
                f"{patient}:retained={float(donor_summary[patient]['ordinal_join_fraction']):.6f},"
                f"stored_to_public={float(ordinal_audit_by_patient[patient]['stored_positive_to_public_csv_fraction']):.6f}"
                for patient in ordinal_passing_donors
            ),
            "pass": len(ordinal_passing_donors) >= 3,
            "evidence": "yakubovsky_join_audit.tsv and independent normalized barcode comparison",
        },
        {
            "criterion_id": "ordinal_semantic_scope",
            "required_for_continuous_lipid": False,
            "required_for_ordinal_lipid": True,
            "observed": "binary source-defined lipid-zone versus non-lipid-zone category",
            "pass": True,
            "evidence": "no category severity or ordered exposure is deposited",
        },
    ]

    extra_lipid_donors = sorted(
        patient for patient, summary in donor_summary.items()
        if bool(summary["lipid_field_present"]) and patient not in PAPER_LIPID_SPOTS
    )
    missing_ordinal_files = [
        patient for patient in CODED_ORDINAL_DONORS
        if not ordinal_audit_by_patient[patient]["public_loupe_file_present"]
    ]
    blocker_rows = code_diagnostics + [
        {
            "diagnostic_id": "YAK-LIPID-02",
            "severity": "scope_guard",
            "observed": ";".join(extra_lipid_donors) if extra_lipid_donors else "none",
            "expected": "primary lipid donors fixed to M1;M2;M3;P6",
            "pass": True,
            "detail": "Any extra quantitative field is excluded from the primary test because the paper's image-quality gate names exactly four donors.",
        },
        {
            "diagnostic_id": "YAK-ZON-01",
            "severity": "information",
            "observed": f"{len(aggregate_rows)} aggregate zonation cells without barcodes/coordinates/lipid",
            "expected": "spot-level join bridge",
            "pass": False,
            "detail": "zon_struct_all_full.mat summarizes zonation categories and cannot establish the missing barcode-to-lipid mapping.",
        },
        {
            "diagnostic_id": "YAK-ORD-02",
            "severity": "ordinal_coverage_limit",
            "observed": ";".join(missing_ordinal_files) if missing_ordinal_files else "none",
            "expected": "public lipid_zones.csv for every source-coded donor",
            "pass": not missing_ordinal_files,
            "detail": "P14, P18, and P17 are untestable for the binary category because their source-coded Loupe files are not deposited; labels must not be reconstructed.",
        },
        {
            "diagnostic_id": "YAK-ORD-03",
            "severity": "terminology_guard",
            "observed": "one binary membership label per retained spot",
            "expected": "binary lipid-zone versus non-lipid-zone interpretation only",
            "pass": True,
            "detail": "The Loupe category has no deposited severity or ordering. Never call it a continuous lipid gradient or an ordered dose-response.",
        },
        {
            "diagnostic_id": "YAK-ORD-04",
            "severity": "interpretation_limit",
            "observed": "non_lipid_zone is the complement of Loupe lipid-zone membership among retained spots",
            "expected": "category contrast, not proof of zero lipid",
            "pass": True,
            "detail": "The complement must not be called lipid-free; it is only the source-defined non-lipid-zone comparison class.",
        },
    ]

    field_map_fields = [
        "concept", "source_path", "source_lines", "deposited_field",
        "code_support", "authoritative_for_spot_join",
    ]
    schema_fields = list(schema_rows[0])
    join_fields = list(join_rows[0])
    donor_fields = list(next(iter(donor_summary.values())))
    cohort_fields = list(cohort_rows[0])
    aggregate_fields = list(aggregate_rows[0])
    criteria_fields = list(criteria[0])
    blocker_fields = list(blocker_rows[0])
    ordinal_audit_fields = list(ordinal_audit_rows[0])

    atomic_write_tsv(output_dir / "yakubovsky_source_field_map.tsv", code_map, field_map_fields)
    atomic_write_tsv(output_dir / "yakubovsky_hdf5_join_schema.tsv", schema_rows, schema_fields)
    atomic_write_tsv(output_dir / "yakubovsky_join_audit.tsv", join_rows, join_fields)
    atomic_write_tsv(
        output_dir / "yakubovsky_donor_schema_summary.tsv",
        list(donor_summary.values()), donor_fields,
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_ordinal_loupe_join_audit.tsv",
        ordinal_audit_rows, ordinal_audit_fields,
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_object_cohort_reconciliation.tsv",
        cohort_rows, cohort_fields,
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_aggregate_zonation_inventory.tsv",
        aggregate_rows, aggregate_fields,
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_source_gate_criteria.tsv", criteria, criteria_fields
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_unresolved_blockers.tsv", blocker_rows, blocker_fields
    )
    atomic_write_tsv(
        output_dir / "yakubovsky_integrity_recheck.tsv",
        integrity_rows, list(integrity_rows[0]),
    )

    ordinal_counts = ";".join(
        f"{patient}:{donor_summary[patient]['n_lipid_zone_spots']}/"
        f"{donor_summary[patient]['n_non_lipid_zone_spots']}/"
        f"{donor_summary[patient]['n_spots']}"
        for patient in ordinal_passing_donors
    )
    if continuous_lipid_gate_pass:
        gate_reason = "The authoritative continuous lipid source gate passed."
        next_action = "Proceed with the prespecified continuous lipid model."
    elif ordinal_lipid_gate_pass:
        gate_reason = (
            "The quantitative lipid arrays pass structural checks but remain unauthorized "
            "because their barcode-keyed construction/import is absent. The independent "
            "source-defined Loupe barcode audit passes for at least three donors, authorizing "
            "only a binary lipid-zone versus non-lipid-zone analysis."
        )
        next_action = (
            "Proceed only with the source-defined binary lipid-zone versus non-lipid-zone "
            "analysis in the passing donors; do not describe this exposure as continuous, "
            "graded, ordered, or dose-responsive."
        )
    elif integrity_pass and expression_zonation_pass:
        gate_reason = (
            "Expression and zonation join, but neither an authoritative continuous lipid "
            "join nor the minimum three-donor ordinal fallback passed."
        )
        next_action = "Use the deposit as a zonation reference only."
    else:
        gate_reason = "A required integrity or expression-zonation join criterion failed."
        next_action = "Stop dataset-native inference and resolve the failing source criterion."

    source_gate = {
        "release_id": RELEASE_ID,
        "dataset": "yakubovsky2026",
        "status": terminal_verdict,
        "terminal_verdict": terminal_verdict,
        "source_gate_pass": continuous_lipid_gate_pass or ordinal_lipid_gate_pass,
        "integrity_gate_pass": integrity_pass,
        "zonation_gate_pass": expression_zonation_pass,
        "continuous_lipid_gate_pass": continuous_lipid_gate_pass,
        "ordinal_lipid_gate_pass": ordinal_lipid_gate_pass,
        "authoritative_barcode_join_verified": expression_zonation_pass,
        "authoritative_lipid_join_verified": authoritative_lipid_join,
        "authoritative_continuous_lipid_join_verified": authoritative_lipid_join,
        "authoritative_ordinal_lipid_join_verified": bool(
            authoritative_ordinal_source and ordinal_lipid_gate_pass
        ),
        "n_expression_donors": len(donor_summary),
        "n_paper_lipid_donors": len(PAPER_LIPID_SPOTS),
        "n_paper_lipid_donors_present": sum(
            patient in donor_summary for patient in PAPER_LIPID_SPOTS
        ),
        "paper_lipid_spot_counts_match": primary_counts_match,
        "n_paper_lipid_donors_with_matching_spot_counts": n_primary_counts_matching,
        "structural_lipid_join_fraction": structural_join_fraction,
        "lipid_measurement_coverage_fraction": lipid_measurement_coverage_fraction,
        "minimum_donor_thresholds_pass": donor_thresholds_pass,
        "n_continuous_structurally_passing_donors": len(continuous_passing_donors),
        "continuous_structurally_passing_donors": ";".join(
            continuous_passing_donors
        ),
        "n_coded_ordinal_donors": len(CODED_ORDINAL_DONORS),
        "n_public_ordinal_donors": sum(
            bool(row["public_loupe_file_present"]) for row in ordinal_audit_rows
        ),
        "n_passing_ordinal_donors": len(ordinal_passing_donors),
        "passing_ordinal_donors": ";".join(ordinal_passing_donors),
        "ordinal_donor_spot_counts_lipid_nonlipid_total": ordinal_counts,
        "unavailable_ordinal_donors": ";".join(missing_ordinal_files),
        "ordinal_effect_interpretation": "binary_lipid_zone_vs_non_lipid_zone",
        "reason": gate_reason,
        "next_action": next_action,
        "v2_registry_read": False,
    }
    gate_fields = list(source_gate)
    atomic_write_tsv(output_dir / "yakubovsky_source_gate.tsv", [source_gate], gate_fields)
    atomic_write_tsv(output_dir / "gate_status.tsv", [source_gate], gate_fields)
    atomic_write_json(output_dir / "yakubovsky_source_gate.json", source_gate)

    repo_commit = git_value(source_dir, "rev-parse", "HEAD")
    repo_origin = git_value(source_dir, "remote", "get-url", "origin")
    repo_dirty = bool(git_value(source_dir, "status", "--porcelain"))
    execution = {
        "release_id": RELEASE_ID,
        "dataset": "yakubovsky2026",
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256_file(Path(__file__).resolve()),
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "superseded_slurm_job_ids": "19631161;19631452;19635568",
        "superseded_job_states": (
            "CANCELLED;COMPLETED_BUT_GATE_LOGIC_SUPERSEDED;"
            "CANCELLED_BEFORE_OUTPUT"
        ),
        "supersession_reason": (
            "19631161 cancelled after a micromamba shared-process lock stall; "
            "19631452 preserved in the acquisition archive after independent "
            "rederivation found an incorrect 100-percent all-reference-donor "
            "join requirement instead of the prespecified 90-percent lipid-"
            "tested-barcode threshold; 19635568 cancelled before output when a "
            "second audit found the continuous fallback must require at least "
            "three usable donors rather than all four paper donors"
        ),
        "gate_logic_version": "prespecified_90pct_lipid_tested_v2",
        "cellstr_resolution": "raw_object_reference_address_plus_sequential_hdf5_visit",
        "started_utc": started,
        "completed_utc": utc_now(),
        "github_origin": repo_origin,
        "github_commit": repo_commit,
        "github_dirty": repo_dirty,
        "expression_matrix_values_read": False,
        "v2_registry_read": False,
        "exit_state": "pass",
    }
    atomic_write_tsv(
        output_dir / "execution_manifest_source_join.tsv", [execution], list(execution)
    )

    print(json.dumps(source_gate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
