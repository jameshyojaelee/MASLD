#!/usr/bin/env python3
"""Independently verify target-value-free GSE296875 model inputs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import h5py
import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-observed-multiome-model-input-verification-v1"
VERIFICATION_ID = "gse296875_observed_multiome_model_input_verification_20260825"
HEX64 = re.compile(r"^[0-9a-f]{64}$")


class ModelInputVerificationError(ValueError):
    """Raised when a model input is not isolated from held target values."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ModelInputVerificationError("JSON object required")
    return value


def _strings(values: Sequence[Any]) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _salted_hash(namespace: str, *parts: str) -> str:
    return sha256("\0".join((namespace, *parts)).encode("utf-8")).hexdigest()


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise ModelInputVerificationError(f"{label} drifted")
    return path


def _tsv(path: Path, fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(fields):
            raise ModelInputVerificationError(f"{path.name} schema differs")
        return [dict(row) for row in reader]


def _csr(group: Any, shape: tuple[int, int]) -> dict[str, Any]:
    if set(group) != {"data", "indices", "indptr", "shape"}:
        raise ModelInputVerificationError("CSR schema differs")
    observed_shape = tuple(int(value) for value in group["shape"][:])
    data = np.asarray(group["data"][:])
    indices = np.asarray(group["indices"][:], dtype=np.int64)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    if (
        observed_shape != shape
        or data.ndim != 1
        or indices.shape != data.shape
        or indptr.shape != (shape[0] + 1,)
        or int(indptr[0]) != 0
        or int(indptr[-1]) != data.size
        or np.any(indptr[1:] < indptr[:-1])
        or np.any(indices < 0)
        or np.any(indices >= shape[1])
        or np.any(data < 0)
        or not np.equal(data, np.floor(data)).all()
    ):
        raise ModelInputVerificationError("CSR values or axes differ")
    return {
        "shape": list(shape),
        "nnz": int(data.size),
        "sum": int(data.sum()),
        "data_sha256": sha256(data.tobytes()).hexdigest(),
        "indices_sha256": sha256(indices.tobytes()).hexdigest(),
        "indptr_sha256": sha256(indptr.tobytes()).hexdigest(),
    }


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("verification_id") != VERIFICATION_ID:
        raise ModelInputVerificationError("verification identity differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise ModelInputVerificationError("verification firewall is open")
    model = _tree(root, config["model_input_artifact"], "model input artifact")
    masks = _tree(root, config["mask_plan_artifact"], "mask plan artifact")
    raw = config.get("raw_donor_identifier_reference")
    if not isinstance(raw, dict) or raw.get("allowed_dataset") != "obs/donor_id" or raw.get("biological_count_values_read") is not False:
        raise ModelInputVerificationError("raw donor reference contract differs")
    lexical = reject_symlink_components(root / str(raw.get("path", "")), label="raw donor reference")
    raw_path = lexical.resolve(strict=True)
    raw_path.relative_to(root)
    if not raw_path.is_file() or _digest(raw_path) != raw.get("sha256"):
        raise ModelInputVerificationError("raw donor reference drifted")
    return {"model": model, "masks": masks, "raw": raw_path}


def verify(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    resolved = validate_config(root, config)
    expected = config["expected"]
    parent = _json(resolved["model"] / "receipt.json")
    children = parent.get("child_artifacts")
    if not isinstance(children, list) or len(children) != expected["children"]:
        raise ModelInputVerificationError("child roster differs")
    with h5py.File(resolved["raw"], "r") as handle:
        raw_donors = set(_strings(handle["obs/donor_id"][:]))
    if len(raw_donors) != expected["donors"]:
        raise ModelInputVerificationError("raw donor census differs")
    input_masks = _tsv(
        resolved["masks"] / "observed_atac_input_peaks.tsv",
        ("input_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "source_genomic_fold", "held_genomic_fold"),
    )
    target_masks = _tsv(
        resolved["masks"] / "target_peaks.tsv",
        ("target_hash", "chromosome", "bed_start_0based", "bed_end_half_open", "genomic_fold"),
    )
    planned_units = _tsv(
        resolved["masks"] / "donor_lineage_units.tsv",
        ("donor_hash", "lineage", "donor_fold", "nuclei"),
    )
    fold_receipts: list[dict[str, Any]] = []
    reference_rows: tuple[str, ...] | None = None
    reference_rna: dict[str, Any] | None = None
    observed_folds: set[int] = set()
    for child_record in children:
        fold = int(child_record["held_genomic_fold"])
        if fold in observed_folds:
            raise ModelInputVerificationError("duplicate held genomic fold")
        observed_folds.add(fold)
        child = resolved["model"] / str(child_record["path"])
        verify_frozen_tree(child)
        if _digest(child / "ARTIFACTS.json") != child_record.get("artifacts_sha256"):
            raise ModelInputVerificationError("child artifact drifted")
        receipt = _json(child / "receipt.json")
        with h5py.File(child / "model_input.h5", "r") as handle:
            if (
                set(handle) != {"rows", "rna", "observed_atac_input", "targets_without_values"}
                or handle.attrs.get("schema_version") != "masld-bench-observed-multiome-model-input-h5-v1"
                or int(handle.attrs.get("held_genomic_fold", -1)) != fold
                or bool(handle.attrs.get("target_atac_values_present", True))
            ):
                raise ModelInputVerificationError("model HDF5 root differs")
            rows = handle["rows"]
            if set(rows) != {"row_hash", "unit_hash", "lineage", "donor_fold", "nuclei", "rna_observed_mask", "atac_observed_mask"}:
                raise ModelInputVerificationError("row schema differs")
            row_hashes = tuple(_strings(rows["row_hash"][:]))
            unit_hashes = _strings(rows["unit_hash"][:])
            lineages = _strings(rows["lineage"][:])
            donor_folds = [int(value) for value in rows["donor_fold"][:]]
            nuclei = np.asarray(rows["nuclei"][:], dtype=np.int64)
            if (
                len(row_hashes) != expected["donor_lineage_units"]
                or len(set(row_hashes)) != len(row_hashes)
                or any(not HEX64.fullmatch(value) for value in row_hashes + tuple(unit_hashes))
                or len(set(unit_hashes)) != expected["donors"]
                or len(set(lineages)) != expected["lineages"]
                or int(nuclei.sum()) != expected["nuclei"]
                or not np.all(rows["rna_observed_mask"][:] == 1)
                or not np.all(rows["atac_observed_mask"][:] == 1)
            ):
                raise ModelInputVerificationError("row census or masks differ")
            expected_row_hashes = [
                _salted_hash(
                    "gse296875:observed_multiome:row:v1",
                    row["donor_hash"],
                    row["lineage"],
                )
                for row in planned_units
            ]
            if (
                list(row_hashes) != expected_row_hashes
                or unit_hashes != [row["donor_hash"] for row in planned_units]
                or lineages != [row["lineage"] for row in planned_units]
                or donor_folds != [int(row["donor_fold"]) for row in planned_units]
                or list(nuclei) != [int(row["nuclei"]) for row in planned_units]
            ):
                raise ModelInputVerificationError("row identity differs from mask plan")
            unit_lineages: dict[str, set[str]] = defaultdict(set)
            unit_folds: dict[str, set[int]] = defaultdict(set)
            for unit, lineage, donor_fold in zip(unit_hashes, lineages, donor_folds, strict=True):
                unit_lineages[unit].add(lineage)
                unit_folds[unit].add(donor_fold)
            fold_census = Counter(next(iter(values)) for values in unit_folds.values())
            if (
                any(len(values) != expected["lineages"] for values in unit_lineages.values())
                or any(len(values) != 1 for values in unit_folds.values())
                or [fold_census[index] for index in range(5)] != expected["donor_fold_census"]
            ):
                raise ModelInputVerificationError("donor grouping differs")
            if raw_donors.intersection(row_hashes) or raw_donors.intersection(unit_hashes):
                raise ModelInputVerificationError("raw donor ID entered model artifact")
            if reference_rows is None:
                reference_rows = row_hashes
            elif row_hashes != reference_rows:
                raise ModelInputVerificationError("row order differs across genomic folds")
            rna = handle["rna"]
            ensembl_ids = _strings(rna["ensembl_id"][:]) if "ensembl_id" in rna else []
            gene_names = _strings(rna["gene_name"][:]) if "gene_name" in rna else []
            if (
                set(rna) != {"counts_csr", "ensembl_id", "gene_name"}
                or rna.attrs.get("scale") != "raw_integer_donor_by_lineage_sum"
                or len(ensembl_ids) != expected["rna_features"]
                or len(gene_names) != expected["rna_features"]
            ):
                raise ModelInputVerificationError("RNA schema differs")
            rna_identity = _csr(rna["counts_csr"], (expected["donor_lineage_units"], expected["rna_features"]))
            if reference_rna is None:
                reference_rna = rna_identity
            elif rna_identity != reference_rna:
                raise ModelInputVerificationError("RNA values differ across genomic folds")
            observed = handle["observed_atac_input"]
            input_rows = [row for row in input_masks if int(row["held_genomic_fold"]) == fold]
            if (
                set(observed) != {"counts_csr", "input_hash", "chromosome", "bed_start_0based", "bed_end_half_open"}
                or observed.attrs.get("scale") != "raw_integer_donor_by_lineage_sum"
            ):
                raise ModelInputVerificationError("observed ATAC input schema differs")
            observed_identity = _csr(observed["counts_csr"], (expected["donor_lineage_units"], expected["input_peaks_per_fold"]))
            if (
                _strings(observed["input_hash"][:]) != [row["input_hash"] for row in input_rows]
                or _strings(observed["chromosome"][:]) != [row["chromosome"] for row in input_rows]
                or [int(value) for value in observed["bed_start_0based"][:]] != [int(row["bed_start_0based"]) for row in input_rows]
                or [int(value) for value in observed["bed_end_half_open"][:]] != [int(row["bed_end_half_open"]) for row in input_rows]
                or any(int(row["source_genomic_fold"]) == fold for row in input_rows)
            ):
                raise ModelInputVerificationError("observed ATAC mask differs")
            targets = handle["targets_without_values"]
            target_rows = [row for row in target_masks if int(row["genomic_fold"]) == fold]
            if (
                set(targets) != {"target_hash", "chromosome", "bed_start_0based", "bed_end_half_open"}
                or bool(targets.attrs.get("atac_values_present", True))
                or _strings(targets["target_hash"][:]) != [row["target_hash"] for row in target_rows]
                or _strings(targets["chromosome"][:]) != [row["chromosome"] for row in target_rows]
                or [int(value) for value in targets["bed_start_0based"][:]] != [int(row["bed_start_0based"]) for row in target_rows]
                or [int(value) for value in targets["bed_end_half_open"][:]] != [int(row["bed_end_half_open"]) for row in target_rows]
                or set(_strings(targets["target_hash"][:])).intersection(_strings(observed["input_hash"][:]))
            ):
                raise ModelInputVerificationError("target-without-values schema differs")
        if (
            receipt.get("held_genomic_fold") != fold
            or receipt.get("donor_lineage_units") != expected["donor_lineage_units"]
            or receipt.get("target_atac_values_read") is not False
            or receipt.get("target_atac_values_written") is not False
            or receipt.get("rna_identity") != rna_identity
            or receipt.get("observed_atac_input_identity") != observed_identity
        ):
            raise ModelInputVerificationError("child receipt target firewall differs")
        fold_receipts.append({"held_genomic_fold": fold, "rna_identity": rna_identity, "observed_atac_identity": observed_identity})
    if observed_folds != set(range(expected["children"])):
        raise ModelInputVerificationError("held genomic folds differ")
    return {
        "schema_version": "masld-bench-observed-multiome-model-input-verification-receipt-v1",
        "verification_id": VERIFICATION_ID,
        "model_input_artifacts_sha256": config["model_input_artifact"]["artifacts_sha256"],
        "folds_verified": len(fold_receipts),
        "fold_receipts": sorted(fold_receipts, key=lambda row: row["held_genomic_fold"]),
        "donors": expected["donors"],
        "donor_lineage_units": expected["donor_lineage_units"],
        "raw_donor_identifier_reference_read": True,
        "biological_count_values_read_from_raw_source": False,
        "target_atac_values_present_in_model_artifacts": False,
        "development_metric_calculated": False,
        "model_fit": False,
        "sealed_outcomes_read": False,
        "promotion_gate_passed": True
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = verify(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_model_inputs.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_model_input_verification", "promotion_gate_passed": True, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
