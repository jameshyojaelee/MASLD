#!/usr/bin/env python3
"""Materialize ARC-adjusted GSE296875 bigWig insertions on the exchange axis.

GSE296875 Cell Ranger ARC endpoints are already Tn5-adjusted; this script
counts ``start`` and ``end-1`` without another shift.  Explicit +4/-5 applies
only to construction of the custom-Bowtie GSE244832 query authority.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from hashlib import sha256
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import (
    freeze_tree,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-gse244832-task-native-blind-baselines-v1"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
ROLES = ("valid", "test")
_WORKER_WINDOWS: dict[str, tuple[tuple[str, int, int], ...]] = {}


class SourceWindowMaterializationError(ValueError):
    """Raised when source tracks or the fixed-window/Tn5 requirement differs."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SourceWindowMaterializationError(f"JSON object required: {path}")
    return value


def _resolve_tree(
    root: Path,
    record: Mapping[str, Any],
    artifact_class: str,
    *,
    verify_payload: bool,
) -> Path:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    manifest_path = path / "ARTIFACTS.json"
    manifest = (
        verify_frozen_tree(path)
        if verify_payload
        else read_json(manifest_path)
    )
    if (
        path.is_symlink()
        or digest(manifest_path) != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise SourceWindowMaterializationError(f"{artifact_class} authority differs")
    return path


def validate_config(root: Path, config_path: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    config = read_json(config_path.resolve(strict=True))
    if (
        config.get("schema_version") != SCHEMA
        or config.get("source_dataset_id") != "gse296875"
        or config.get("target_dataset_id") != "gse244832"
        or config.get("models") != list(("lsi", "observed_atac_glm", "observed_atac_only"))
    ):
        raise SourceWindowMaterializationError("campaign identity differs")
    for field in (
        "condition_values_read",
        "phenotype_values_read",
        "rna_assay_read",
        "development_outcomes_read",
        "sealed_data_read",
        "metrics_calculated",
        "target_role_atac_consumed",
        "target_dataset_fit_or_adaptation",
        "missing_as_zero",
    ):
        if config.get("firewall", {}).get(field) is not False:
            raise SourceWindowMaterializationError(f"campaign firewall opened: {field}")
    paths = {
        "registration": _resolve_tree(
            root,
            config["execution_registration"],
            "gse244832_observed_atac_execution_registration",
            verify_payload=True,
        ),
        "axis": _resolve_tree(
            root,
            config["exchange_axis"],
            "gse244832_reference_guarded_atac_exchange_axis",
            verify_payload=True,
        ),
        "source": _resolve_tree(
            root,
            config["source_bigwigs"],
            "gse296875_deduplicated_tn5_bigwigs",
            verify_payload=False,
        ),
        "source_fragments": _resolve_tree(
            root,
            config["source_fragments"],
            "gse296875_donor_lineage_fragments",
            verify_payload=False,
        ),
    }
    registration = read_json(paths["registration"] / "execution_contract.json")
    source_contract = read_json(paths["source"] / "contract.json")
    fragment_contract = read_json(paths["source_fragments"] / "contract.json")
    if (
        registration.get("axis_artifacts_sha256")
        != config["exchange_axis"]["artifacts_sha256"]
        or registration.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or registration.get("sequence_execution_authorized") is not False
        or source_contract.get("source_coordinates_are_already_tn5_adjusted") is not True
        or source_contract.get("source_coordinate_transform")
        != "cellranger_arc_alignment_plus_4_minus_5"
        or source_contract.get("additional_shift_applied") is not False
        or source_contract.get("left_insertion_position") != "fragment_start"
        or source_contract.get("right_insertion_position")
        != "fragment_end_minus_1_BED_exclusive"
        or source_contract.get("output_signal_contract")
        != "two_insertions_per_unique_fragment_record"
        or source_contract.get("fragments_artifacts_sha256")
        != config["source_fragments"]["artifacts_sha256"]
        or fragment_contract.get("source_coordinates")
        != "cellranger_arc_2_0_0_tn5_adjusted"
        or fragment_contract.get("assay_signal_unit")
        != "one_unique_fragment_record_not_readSupport"
    ):
        raise SourceWindowMaterializationError("source/target Tn5 geometry differs")
    return config, paths


def read_windows(path: Path) -> dict[str, tuple[tuple[str, int, int], ...]]:
    by_role: dict[str, list[tuple[str, int, int]]] = {role: [] for role in ROLES}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "window_index",
            "window_id",
            "role",
            "contig",
            "start",
            "end",
        ):
            raise SourceWindowMaterializationError("fixed-window fields differ")
        seen: set[str] = set()
        expected_index = 0
        contigs: dict[str, set[str]] = {role: set() for role in ROLES}
        for row in reader:
            role = row["role"]
            start, end = int(row["start"]), int(row["end"])
            if (
                role not in by_role
                or int(row["window_index"]) != expected_index
                or row["window_id"] in seen
                or start < 0
                or end - start != 1000
            ):
                raise SourceWindowMaterializationError("fixed-window identity differs")
            expected_index += 1
            seen.add(row["window_id"])
            contigs[role].add(row["contig"])
            by_role[role].append((row["contig"], start, end))
    if (
        any(len(rows) != 16_000 for rows in by_role.values())
        or contigs["valid"] & contigs["test"]
    ):
        raise SourceWindowMaterializationError("whole-contig rotation geometry differs")
    return {role: tuple(rows) for role, rows in by_role.items()}


def read_source_units(source: Path) -> list[dict[str, Any]]:
    records = {
        row["path"]: row for row in read_json(source / "ARTIFACTS.json")["artifacts"]
    }
    units: list[dict[str, Any]] = []
    with (source / "bigwig_manifest.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "donor_id",
            "lineage_id",
            "analysis_role",
            "path",
            "size_bytes",
            "sha256",
        }
        if not required.issubset(reader.fieldnames or ()):
            raise SourceWindowMaterializationError("source bigWig manifest fields differ")
        for row in reader:
            if row["analysis_role"] == "primary" and row["lineage_id"] in LINEAGES:
                artifact = records.get(row["path"])
                if (
                    artifact is None
                    or int(row["size_bytes"]) != int(artifact["size_bytes"])
                    or row["sha256"] != artifact["sha256"]
                ):
                    raise SourceWindowMaterializationError("source bigWig binding differs")
                units.append(
                    {
                        "donor_id": row["donor_id"],
                        "lineage_id": row["lineage_id"],
                        "path": row["path"],
                        "size_bytes": int(row["size_bytes"]),
                        "sha256": row["sha256"],
                    }
                )
    donors = sorted({row["donor_id"] for row in units}, key=int)
    if (
        len(donors) != 39
        or len(units) != 39 * 4
        or {
            (row["donor_id"], row["lineage_id"]) for row in units
        }
        != {(donor, lineage) for donor in donors for lineage in LINEAGES}
    ):
        raise SourceWindowMaterializationError("source donor-lineage rectangle differs")
    return sorted(units, key=lambda row: (int(row["donor_id"]), LINEAGES.index(row["lineage_id"])))


def _worker_init(windows: Mapping[str, Sequence[tuple[str, int, int]]]) -> None:
    global _WORKER_WINDOWS
    _WORKER_WINDOWS = {role: tuple(values) for role, values in windows.items()}


def _read_one_bigwig(source_root: str, unit: Mapping[str, Any]) -> tuple[str, str, Any]:
    import numpy as np
    import pyBigWig

    path = (Path(source_root) / str(unit["path"])).resolve(strict=True)
    if path.is_symlink() or path.stat().st_size != int(unit["size_bytes"]):
        raise SourceWindowMaterializationError("source bigWig path or size differs")
    values = np.zeros((len(ROLES), 16_000), dtype=np.uint32)
    handle = pyBigWig.open(str(path))
    try:
        for role_index, role in enumerate(ROLES):
            for window_index, (contig, start, end) in enumerate(_WORKER_WINDOWS[role]):
                observed = handle.stats(contig, start, end, type="sum", exact=True)[0]
                total = 0.0 if observed is None else float(observed)
                rounded = round(total)
                if total < 0 or abs(total - rounded) > 1.0e-5 or rounded > np.iinfo(np.uint32).max:
                    raise SourceWindowMaterializationError("source bigWig count is not uint32-exact")
                values[role_index, window_index] = rounded
    finally:
        handle.close()
    if np.any(values.sum(axis=1, dtype=np.uint64) == 0):
        raise SourceWindowMaterializationError("source donor-lineage role has zero support")
    return str(unit["donor_id"]), str(unit["lineage_id"]), values


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_runtime_lock(output: Path) -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"],
        check=True,
        capture_output=True,
        text=True,
    )
    write_text_exclusive(output / "pip_freeze.txt", completed.stdout)
    write_text_exclusive(
        output / "python_version.txt",
        f"{platform.python_version()}\n{sys.executable}\n",
    )


def materialize(
    root: Path,
    config_path: Path,
    output: Path,
    workers: int,
) -> dict[str, Any]:
    import numpy as np

    if output.exists() or not 1 <= workers <= 16:
        raise SourceWindowMaterializationError("output or worker contract differs")
    root = root.resolve(strict=True)
    config, paths = validate_config(root, config_path)
    windows = read_windows(paths["axis"] / "windows.tsv")
    units = read_source_units(paths["source"])
    output.mkdir(parents=True, mode=0o750)
    counts = np.zeros((39, 4, 2, 16_000), dtype=np.uint32)
    donor_order = sorted({row["donor_id"] for row in units}, key=int)
    donor_index = {donor: index for index, donor in enumerate(donor_order)}
    completed: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_worker_init,
        initargs=(windows,),
    ) as executor:
        futures = {
            executor.submit(_read_one_bigwig, str(paths["source"]), unit): unit
            for unit in units
        }
        for future in as_completed(futures):
            donor, lineage, values = future.result()
            counts[donor_index[donor], LINEAGES.index(lineage)] = values
            completed.append(
                {
                    "donor_id": donor,
                    "lineage_id": lineage,
                    "valid_total": int(values[0].sum(dtype=np.uint64)),
                    "test_total": int(values[1].sum(dtype=np.uint64)),
                }
            )
    if len(completed) != 156 or np.any(counts.sum(axis=3, dtype=np.uint64) == 0):
        raise SourceWindowMaterializationError("source fixed-window count rectangle differs")
    np.save(output / "source_counts.uint32.npy", counts, allow_pickle=False)
    write_tsv(
        output / "source_axis.tsv",
        ("donor_index", "donor_id", "lineage_index", "lineage_id"),
        [
            {
                "donor_index": donor_idx,
                "donor_id": donor,
                "lineage_index": lineage_idx,
                "lineage_id": lineage,
            }
            for donor_idx, donor in enumerate(donor_order)
            for lineage_idx, lineage in enumerate(LINEAGES)
        ],
    )
    write_tsv(
        output / "materialization_qc.tsv",
        ("donor_id", "lineage_id", "valid_total", "test_total"),
        sorted(completed, key=lambda row: (int(row["donor_id"]), LINEAGES.index(row["lineage_id"]))),
    )
    write_json_exclusive(
        output / "contract.json",
        {
            "schema_version": "masld-bench-gse296875-gse244832-source-window-counts-v1",
            "dataset_id": "gse296875",
            "donors": 39,
            "lineages": list(LINEAGES),
            "roles": list(ROLES),
            "windows_per_role": 16_000,
            "shape": list(counts.shape),
            "dtype": "uint32",
            "biological_unit": "donor_lineage_pseudobulk",
            "signal_unit": "deduplicated_tn5_insertion_count",
            "source_fragment_transform": "cellranger_arc_plus_4_minus_5_already_applied",
            "additional_shift_applied": False,
            "exchange_axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
            "source_bigwigs_artifacts_sha256": config["source_bigwigs"]["artifacts_sha256"],
            "source_fragments_artifacts_sha256": config["source_fragments"]["artifacts_sha256"],
            "execution_registration_artifacts_sha256": config["execution_registration"]["artifacts_sha256"],
            "source_bigwig_payload_hashes_reused_from_frozen_manifest": True,
            "condition_or_phenotype_read": False,
            "target_dataset_values_read": False,
            "development_outcomes_read": False,
            "metrics_calculated": False,
        },
    )
    write_runtime_lock(output)
    artifact_sha = freeze_tree(
        output,
        {
            "artifact_class": "gse296875_gse244832_source_window_counts",
            "dataset_id": "gse296875",
            "donor_lineage_units": 156,
            "roles": list(ROLES),
            "exchange_axis_artifacts_sha256": config["exchange_axis"]["artifacts_sha256"],
            "source_bigwigs_artifacts_sha256": config["source_bigwigs"]["artifacts_sha256"],
            "source_fragments_artifacts_sha256": config["source_fragments"]["artifacts_sha256"],
            "target_dataset_values_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return {"output": str(output), "artifacts_sha256": artifact_sha, "units": 156}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", dest="config_path", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=16)
    arguments = parser.parse_args()
    print(json.dumps(materialize(**vars(arguments)), sort_keys=True))


if __name__ == "__main__":
    main()
