"""Build and validate the immutable MASLD atlas contract.

This module intentionally reads HDF5 directly. Contract construction must not
materialize the 1.23-million-cell expression matrix or depend on AnnData's
in-memory behavior.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Iterable, Iterator

import h5py
import numpy as np

from .config import canonical_json_bytes, repo_path, write_json_exclusive


class ContractError(RuntimeError):
    """Raised when atlas metadata or counts violate the frozen contract."""


PAIRING_FILES = {
    "GSE136103": "data/GSE136103/metadata/donor_pairing.csv",
    "GSE185477": "data/GSE185477/metadata/donor_pairing.csv",
    "GSE202379": "data/GSE202379/metadata/donor_pairing.csv",
    "GSE244832": "data/GSE244832/metadata/donor_pairing.csv",
}

CONTROL_CONDITIONS = {"Healthy", "Healthy control", "NORMAL"}


@dataclass(frozen=True)
class LibraryRecord:
    dataset: str
    library_id: str
    assay_id: str
    donor_id: str
    n_cells: int
    native_condition: str
    harmonized_stage: str
    preparation: str
    technical_batch: str
    strict_reference: bool
    primary_query: bool
    query_control: bool
    analysis_eligible: bool
    exclusion_reason: str
    provenance: str


def _decode_array(dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[:]
    if values.dtype.kind in {"S", "O", "U"}:
        return dataset.asstr()[:]
    return values


def read_h5ad_column(group: h5py.Group, key: str) -> tuple[np.ndarray, np.ndarray]:
    """Return categorical values and integer codes for an obs/var column."""
    node = group[key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        return _decode_array(node["categories"]), node["codes"][:].astype(np.int64)
    values = _decode_array(node)
    categories, codes = np.unique(values, return_inverse=True)
    return categories, codes.astype(np.int64)


def _read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _pairing_maps(config: dict[str, Any]) -> tuple[dict[str, str], dict[str, str]]:
    library_to_donor: dict[str, str] = {}
    library_to_condition: dict[str, str] = {}
    for dataset, relpath in PAIRING_FILES.items():
        path = repo_path(config, relpath)
        if not path.exists():
            raise ContractError(f"missing donor-pairing authority: {path}")
        for row in _read_csv(path):
            donor = f"{dataset}_{row['donor_id']}"
            for library in row["rna_srrs"].split(";"):
                if library in library_to_donor:
                    raise ContractError(f"library appears twice in donor pairing: {library}")
                library_to_donor[library] = donor
                library_to_condition[library] = row["condition"]
    return library_to_donor, library_to_condition


def _liver_maps(config: dict[str, Any]) -> tuple[dict[str, str], dict[str, str]]:
    run_table = repo_path(config, "Liver_Atlas/metadata/SraRunTable.csv")
    sample_info = repo_path(config, "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv")
    run_rows = {row["Run"]: row for row in _read_csv(run_table)}
    info_rows = {
        row["characteristics: shortFileName"]: row
        for row in _read_csv(sample_info, delimiter="\t")
    }
    run_to_donor: dict[str, str] = {}
    run_to_assay: dict[str, str] = {}
    for run, row in run_rows.items():
        short = row["shortfilename"]
        if short not in info_rows:
            continue
        match = re.search(r"H\d+", info_rows[short]["title"])
        if match is None:
            # The run table also contains non-human or non-scRNA assays. Only
            # canonical human atlas runs are required to resolve below.
            continue
        run_to_donor[run] = f"Liver_Atlas_{match.group()}"
        run_to_assay[run] = row["Experiment"]
    return run_to_donor, run_to_assay


def _liver_histology(config: dict[str, Any]) -> dict[str, dict[str, str]]:
    path = Path(config["_config_path"]).parent / "reference/liver_atlas_histology_v1.tsv"
    rows = _read_csv(path, delimiter="\t")
    if len(rows) != 19:
        raise ContractError(f"expected 19 Liver Atlas histology rows, found {len(rows)}")
    clean = {row["donor_id"] for row in rows if row["strict_clean_reference"] == "true"}
    if clean != {"H06", "H07", "H10", "H22"}:
        raise ContractError(f"unexpected clean Liver Atlas donor set: {sorted(clean)}")
    return {row["donor_id"]: row for row in rows}


def harmonize_stage(dataset: str, condition: str, liver_clean: bool = False) -> str:
    if dataset == "Liver_Atlas":
        return "Healthy" if liver_clean else "Non-clean/unknown liver"
    if condition in CONTROL_CONDITIONS:
        return "Healthy"
    if condition in {"NAFLD", "MASL"}:
        return "Steatosis"
    if condition in {"NASH", "MASH", "NASH w/o cirrhosis", "NASH with cirrhosis"}:
        return "Steatohepatitis"
    if condition == "end stage":
        return "End-stage liver disease"
    if condition == "Cirrhotic":
        return "Cirrhosis (mixed etiology)"
    return condition or "Unknown"


def _library_metadata(atlas: Path) -> tuple[list[str], dict[str, dict[str, str]], Counter[str]]:
    with h5py.File(atlas, "r") as handle:
        obs = handle["obs"]
        required = {"sample", "dataset", "condition", "preparation_method", "cell_type"}
        missing = sorted(required - set(obs.keys()))
        if missing:
            raise ContractError(f"canonical H5AD missing obs columns: {missing}")
        sample_cat, sample_codes = read_h5ad_column(obs, "sample")
        dataset_cat, dataset_codes = read_h5ad_column(obs, "dataset")
        condition_cat, condition_codes = read_h5ad_column(obs, "condition")
        prep_cat, prep_codes = read_h5ad_column(obs, "preparation_method")
        n_cells = len(sample_codes)
        if not all(len(x) == n_cells for x in (dataset_codes, condition_codes, prep_codes)):
            raise ContractError("obs columns have inconsistent lengths")
        counts = Counter(str(sample_cat[i]) for i in sample_codes)
        metadata: dict[str, dict[str, str]] = {}
        for i in range(n_cells):
            library = str(sample_cat[sample_codes[i]])
            values = {
                "dataset": str(dataset_cat[dataset_codes[i]]),
                "condition": str(condition_cat[condition_codes[i]]),
                "preparation": str(prep_cat[prep_codes[i]]),
            }
            prior = metadata.setdefault(library, values)
            if prior != values:
                raise ContractError(f"library metadata varies between cells: {library}")
    return sorted(metadata), metadata, counts


def build_library_records(config: dict[str, Any]) -> list[LibraryRecord]:
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    libraries, observed, cell_counts = _library_metadata(atlas)
    paired_donor, paired_condition = _pairing_maps(config)
    liver_donor, liver_assay = _liver_maps(config)
    histology = _liver_histology(config)
    exclusions = set(config["input"]["excluded_analysis_libraries"])
    strict_config = {
        (dataset, donor)
        for dataset, donors in config["roles"]["strict_reference"].items()
        for donor in donors
    }
    primary_datasets = set(config["roles"]["primary_query_datasets"])
    records: list[LibraryRecord] = []
    for library in libraries:
        meta = observed[library]
        dataset = meta["dataset"]
        assay = library
        provenance = "canonical_h5ad_obs"
        if dataset == "Liver_Atlas":
            if library not in liver_donor or library not in liver_assay:
                raise ContractError(f"unmapped canonical Liver Atlas run: {library}")
            donor = liver_donor[library]
            donor_short = donor.split("_", 2)[-1]
            assay = liver_assay[library]
            row = histology[donor_short]
            liver_clean = row["strict_clean_reference"] == "true"
            condition = "Histology-clean" if liver_clean else "Non-clean/unknown liver"
            provenance = "SraRunTable+GSE192740_sampleInfo+Guilliams_Table_S6"
        elif library in paired_donor:
            donor = paired_donor[library]
            donor_short = donor.split("_", 1)[1]
            condition = paired_condition[library]
            liver_clean = False
            provenance = f"{dataset}/metadata/donor_pairing.csv"
        else:
            donor = f"{dataset}_{library}"
            donor_short = library
            condition = meta["condition"]
            liver_clean = False
        strict_reference = (dataset, donor_short) in strict_config
        primary_query = dataset in primary_datasets
        query_control = primary_query and condition in CONTROL_CONDITIONS
        analysis_eligible = library not in exclusions
        exclusion_reason = "" if analysis_eligible else "prespecified_low_cell_library"
        stage = harmonize_stage(dataset, condition, liver_clean)
        technical_batch = f"{dataset}|{meta['preparation']}"
        records.append(
            LibraryRecord(
                dataset=dataset,
                library_id=library,
                assay_id=assay,
                donor_id=donor,
                n_cells=int(cell_counts[library]),
                native_condition=condition,
                harmonized_stage=stage,
                preparation=meta["preparation"],
                technical_batch=technical_batch,
                strict_reference=strict_reference,
                primary_query=primary_query,
                query_control=query_control,
                analysis_eligible=analysis_eligible,
                exclusion_reason=exclusion_reason,
                provenance=provenance,
            )
        )
    return records


def summarize_contract(records: Iterable[LibraryRecord]) -> dict[str, dict[str, int]]:
    rows = list(records)

    def summary(selected: list[LibraryRecord]) -> dict[str, int]:
        return {
            "libraries": len(selected),
            "donors": len({row.donor_id for row in selected}),
            "cells": sum(row.n_cells for row in selected),
        }

    analyzed = [row for row in rows if row.analysis_eligible]
    reference = [row for row in analyzed if row.strict_reference]
    query = [row for row in analyzed if row.primary_query]
    query_donors = {row.donor_id for row in query}
    control_donors = {row.donor_id for row in query if row.query_control}
    return {
        "descriptive": summary(rows),
        "analyzed": summary(analyzed),
        "strict_reference": summary(reference),
        "primary_query": {
            **summary(query),
            "controls": len(control_donors),
            "cases": len(query_donors - control_donors),
        },
    }


def validate_expected_contract(config: dict[str, Any], observed: dict[str, Any]) -> None:
    expected = config["expected_contract"]
    failures: list[str] = []
    for roster, wanted in expected.items():
        found = observed.get(roster, {})
        for key, value in wanted.items():
            if found.get(key) != value:
                failures.append(f"{roster}.{key}: expected {value}, observed {found.get(key)}")
    if failures:
        raise ContractError("atlas contract failed:\n  " + "\n  ".join(failures))


def hash_gene_order(atlas: Path) -> str:
    with h5py.File(atlas, "r") as handle:
        genes = _decode_array(handle["raw/var/_index"])
    digest = hashlib.sha256()
    for gene in genes:
        digest.update(str(gene).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def sha256_path(path: str | Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_tree(path: str | Path) -> tuple[str, list[dict[str, Any]]]:
    """Hash a directory by ordered relative paths and file contents."""
    root = Path(path).resolve()
    if not root.is_dir():
        raise ContractError(f"model bundle is not a directory: {root}")
    entries: list[dict[str, Any]] = []
    digest = hashlib.sha256()
    for current in sorted(root.rglob("*"), key=lambda value: value.relative_to(root).as_posix()):
        if current.is_symlink():
            raise ContractError(f"symlink is forbidden in immutable bundle: {current}")
        if not current.is_file():
            continue
        relative = current.relative_to(root).as_posix()
        file_hash = sha256_path(current)
        size = current.stat().st_size
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(size).encode("ascii"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\0")
        entries.append({"path": relative, "size_bytes": size, "sha256": file_hash})
    if not entries:
        raise ContractError(f"immutable bundle contains no files: {root}")
    return digest.hexdigest(), entries


def hash_raw_counts(atlas: Path, chunk_size: int = 4_000_000) -> tuple[str, dict[str, Any]]:
    """Stream and validate canonical raw CSR components.

    The digest includes shape, indptr, indices, and data. Float-backed counts are
    accepted only when every value is finite, non-negative, and exactly integral.
    """
    digest = hashlib.sha256()
    with h5py.File(atlas, "r") as handle:
        raw = handle["raw/X"]
        shape = tuple(int(x) for x in raw.attrs["shape"])
        digest.update(json.dumps(shape).encode("ascii"))
        for name in ("indptr", "indices", "data"):
            data = raw[name]
            digest.update(name.encode("ascii"))
            digest.update(str(data.dtype).encode("ascii"))
            for start in range(0, len(data), chunk_size):
                values = data[start : start + chunk_size]
                if name == "data":
                    if not np.all(np.isfinite(values)):
                        raise ContractError("raw/X contains non-finite counts")
                    if np.any(values < 0) or not np.all(values == np.floor(values)):
                        raise ContractError("raw/X is not non-negative integer count data")
                digest.update(np.ascontiguousarray(values).view(np.uint8))
        stats = {"shape": list(shape), "nnz": int(len(raw["data"])), "encoding": "csr_matrix"}
    return digest.hexdigest(), stats


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


class DeterministicGzipTextWriter:
    """Text writer with mtime=0 so identical manifests are byte-identical."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.raw = None
        self.compressed = None
        self.text = None

    def __enter__(self):
        self.raw = self.path.open("xb")
        self.compressed = gzip.GzipFile(
            filename="", mode="wb", fileobj=self.raw, mtime=0
        )
        self.text = io.TextIOWrapper(self.compressed, encoding="utf-8", newline="")
        return self.text

    def __exit__(self, exc_type, exc, traceback):
        try:
            if self.text is not None:
                self.text.close()
        finally:
            # GzipFile deliberately leaves a caller-owned file object open.
            if self.raw is not None and not self.raw.closed:
                self.raw.close()
        return False


def write_contract_outputs(
    config: dict[str, Any], output: Path, *, hash_counts: bool = True, write_cells: bool = True
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False)
    records = build_library_records(config)
    observed = summarize_contract(records)
    validate_expected_contract(config, observed)
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    gene_hash = hash_gene_order(atlas)
    raw_hash, raw_stats = hash_raw_counts(atlas) if hash_counts else (None, {})

    library_rows = [asdict(row) for row in records]
    fields = list(library_rows[0])
    write_tsv(output / "library_manifest.tsv", library_rows, fields)

    assay_groups: dict[str, list[LibraryRecord]] = defaultdict(list)
    donor_groups: dict[str, list[LibraryRecord]] = defaultdict(list)
    for row in records:
        assay_groups[row.assay_id].append(row)
        donor_groups[row.donor_id].append(row)
    assay_rows = []
    for assay, group in sorted(assay_groups.items()):
        assay_rows.append({
            "assay_id": assay,
            "dataset": group[0].dataset,
            "donor_id": group[0].donor_id,
            "libraries": ";".join(sorted(x.library_id for x in group)),
            "n_libraries": len(group),
            "n_cells": sum(x.n_cells for x in group),
        })
    write_tsv(output / "assay_manifest.tsv", assay_rows, list(assay_rows[0]))
    donor_rows = []
    for donor, group in sorted(donor_groups.items()):
        conditions = sorted({x.native_condition for x in group})
        stages = sorted({x.harmonized_stage for x in group})
        if len(conditions) != 1 or len(stages) != 1:
            raise ContractError(f"donor metadata conflict: {donor}: {conditions}, {stages}")
        donor_rows.append({
            "donor_id": donor,
            "dataset": group[0].dataset,
            "native_condition": conditions[0],
            "harmonized_stage": stages[0],
            "strict_reference": any(x.strict_reference for x in group),
            "primary_query": any(x.primary_query for x in group),
            "query_control": any(x.query_control for x in group),
            "analysis_eligible": any(x.analysis_eligible for x in group),
            "n_libraries_descriptive": len(group),
            "n_libraries_analyzed": sum(x.analysis_eligible for x in group),
            "n_cells_descriptive": sum(x.n_cells for x in group),
            "n_cells_analyzed": sum(x.n_cells for x in group if x.analysis_eligible),
        })
    write_tsv(output / "donor_manifest.tsv", donor_rows, list(donor_rows[0]))

    cell_order_hash = None
    if write_cells:
        cell_order_hash = _write_cell_manifest(
            atlas, output / "cell_manifest.tsv.gz", records
        )

    manifest_names = ["library_manifest.tsv", "assay_manifest.tsv", "donor_manifest.tsv"]
    if write_cells:
        manifest_names.append("cell_manifest.tsv.gz")
    manifest_hashes = {name: sha256_path(output / name) for name in manifest_names}

    lock = {
        "schema_version": "masld-cl-contract-v1",
        "config_sha256": config["_config_sha256"],
        "atlas_realpath": os.path.realpath(atlas),
        "atlas_size_bytes": atlas.stat().st_size,
        "observed_contract": observed,
        "gene_order_sha256": gene_hash,
        "raw_counts_sha256": raw_hash,
        "raw_counts": raw_stats,
        "cell_order_sha256": cell_order_hash,
        "manifest_sha256": manifest_hashes,
        "production_ready": raw_hash is not None and write_cells,
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output / "contract_lock.json", lock)
    return lock


def _write_cell_manifest(atlas: Path, path: Path, records: list[LibraryRecord]) -> str:
    by_library = {row.library_id: row for row in records}
    digest = hashlib.sha256()
    with h5py.File(atlas, "r") as handle, DeterministicGzipTextWriter(path) as out:
        obs = handle["obs"]
        sample_cat, sample_codes = read_h5ad_column(obs, "sample")
        cell_cat, cell_codes = read_h5ad_column(obs, "_index")
        type_cat, type_codes = read_h5ad_column(obs, "cell_type")
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "cell_id", "source_cell_id", "library_id", "assay_id", "donor_id",
            "dataset", "cell_type", "strict_reference", "primary_query",
            "query_control", "analysis_eligible",
        ])
        seen: set[str] = set()
        for i in range(len(sample_codes)):
            library = str(sample_cat[sample_codes[i]])
            row = by_library[library]
            source_id = str(cell_cat[cell_codes[i]])
            cell_id = f"{library}|{source_id}"
            if cell_id in seen:
                raise ContractError(f"duplicate compound cell ID: {cell_id}")
            seen.add(cell_id)
            digest.update(cell_id.encode("utf-8"))
            digest.update(b"\0")
            writer.writerow([
                cell_id, source_id, library, row.assay_id, row.donor_id, row.dataset,
                str(type_cat[type_codes[i]]), row.strict_reference, row.primary_query,
                row.query_control, row.analysis_eligible,
            ])
    if len(seen) != sum(row.n_cells for row in records):
        raise ContractError("cell manifest row count changed during construction")
    return digest.hexdigest()


def hash_atlas_cell_order(atlas: Path) -> str:
    digest = hashlib.sha256()
    with h5py.File(atlas, "r") as handle:
        obs = handle["obs"]
        sample_cat, sample_codes = read_h5ad_column(obs, "sample")
        cell_cat, cell_codes = read_h5ad_column(obs, "_index")
        for sample_code, cell_code in zip(sample_codes, cell_codes):
            cell_id = f"{sample_cat[sample_code]}|{cell_cat[cell_code]}"
            digest.update(str(cell_id).encode("utf-8"))
            digest.update(b"\0")
    return digest.hexdigest()


def verify_contract_lock(config: dict[str, Any], lock_path: str | Path, *, full_hash: bool) -> dict[str, Any]:
    with Path(lock_path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-contract-v1":
        raise ContractError("unsupported contract lock schema")
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if hashlib.sha256(canonical_json_bytes(payload)).hexdigest() != lock.get("lock_sha256"):
        raise ContractError("contract lock content hash mismatch")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise ContractError("contract lock was built with a different config")
    validate_expected_contract(config, lock.get("observed_contract", {}))
    atlas = repo_path(config, config["input"]["atlas_h5ad"])
    if os.path.realpath(atlas) != lock.get("atlas_realpath"):
        raise ContractError("canonical atlas realpath changed")
    if atlas.stat().st_size != lock.get("atlas_size_bytes"):
        raise ContractError("canonical atlas size changed")
    if hash_gene_order(atlas) != lock.get("gene_order_sha256"):
        raise ContractError("canonical gene order changed")
    if not lock.get("production_ready") or not lock.get("raw_counts_sha256"):
        raise ContractError("contract lock lacks a complete raw-count hash")
    lock_dir = Path(lock_path).resolve().parent
    expected_manifests = lock.get("manifest_sha256", {})
    required_manifests = {
        "library_manifest.tsv", "assay_manifest.tsv", "donor_manifest.tsv",
        "cell_manifest.tsv.gz",
    }
    if set(expected_manifests) != required_manifests:
        raise ContractError("contract lock does not name all four immutable manifests")
    for name, expected_hash in expected_manifests.items():
        path = lock_dir / name
        if not path.is_file() or sha256_path(path) != expected_hash:
            raise ContractError(f"contract manifest changed or is missing: {name}")
    if full_hash:
        observed_hash, _ = hash_raw_counts(atlas)
        if observed_hash != lock["raw_counts_sha256"]:
            raise ContractError("canonical raw counts changed")
        if hash_atlas_cell_order(atlas) != lock.get("cell_order_sha256"):
            raise ContractError("canonical cell identity or order changed")
    return lock
