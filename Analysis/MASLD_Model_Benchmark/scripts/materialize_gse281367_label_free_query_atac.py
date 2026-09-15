#!/usr/bin/env python3
"""Materialize rotation-specific GSE281367 query ATAC without scored contigs."""

from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import csv
from dataclasses import dataclass
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-gse281367-label-free-query-atac-materialization-v1"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
ROLES = ("valid", "test")


class QueryATACMaterializationError(ValueError):
    """Raised when query ATAC could expose target contigs or label values."""


@dataclass(frozen=True)
class WindowIndex:
    rows: tuple[tuple[str, int, int, str], ...]
    by_contig: Mapping[str, tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]]


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise QueryATACMaterializationError(f"JSON object required: {path}")
    return value


def read_reference_lengths(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] in result or int(fields[1]) <= 0:
                raise QueryATACMaterializationError("reference FAI differs")
            result[fields[0]] = int(fields[1])
    if not all(f"chr{value}" in result for value in [*range(1, 23), "X", "Y"]):
        raise QueryATACMaterializationError("reference lacks primary contigs")
    return result


def read_windows(path: Path) -> dict[str, WindowIndex]:
    by_role: dict[str, list[tuple[str, int, int, str]]] = {role: [] for role in ROLES}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("window_index", "window_id", "role", "contig", "start", "end"):
            raise QueryATACMaterializationError("window fields differ")
        for row in reader:
            role = row["role"]
            if role not in by_role:
                raise QueryATACMaterializationError("window role differs")
            start, end = int(row["start"]), int(row["end"])
            if start < 0 or end - start != 1000:
                raise QueryATACMaterializationError("window geometry differs")
            by_role[role].append((row["contig"], start, end, row["window_id"]))
    if any(len(rows) != 16000 for rows in by_role.values()):
        raise QueryATACMaterializationError("role-specific window count differs")
    if {row[0] for row in by_role["valid"]} & {row[0] for row in by_role["test"]}:
        raise QueryATACMaterializationError("query and scored-target contigs overlap")
    result: dict[str, WindowIndex] = {}
    for role, rows in by_role.items():
        contigs: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
        for index, (contig, start, end, _window) in enumerate(rows):
            contigs[contig].append((start, end, index))
        packed = {}
        for contig, values in contigs.items():
            values.sort()
            packed[contig] = tuple(tuple(item[offset] for item in values) for offset in range(3))
        result[role] = WindowIndex(tuple(rows), packed)
    return result


def windows_at(index: WindowIndex, contig: str, position: int) -> Iterable[int]:
    values = index.by_contig.get(contig)
    if values is None:
        return ()
    starts, ends, indexes = values
    offset = bisect_right(starts, position) - 1
    hits: list[int] = []
    while offset >= 0 and ends[offset] > position:
        hits.append(indexes[offset])
        offset -= 1
    return hits


def read_membership(path: Path) -> tuple[dict[str, dict[str, str]], dict[str, set[str]]]:
    primary: dict[str, dict[str, str]] = defaultdict(dict)
    all_barcodes: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "dataset_id", "donor_id", "raw_barcode", "lineage_id", "analysis_role", "outer_fold"
        ):
            raise QueryATACMaterializationError("label-free membership fields differ")
        for row in reader:
            if row["dataset_id"] != "gse281367":
                raise QueryATACMaterializationError("membership dataset differs")
            donor, barcode = row["donor_id"], row["raw_barcode"]
            if barcode in all_barcodes[donor]:
                raise QueryATACMaterializationError("donor-local barcode duplicated")
            all_barcodes[donor].add(barcode)
            if row["analysis_role"] == "primary":
                if row["lineage_id"] not in LINEAGES:
                    raise QueryATACMaterializationError("primary lineage differs")
                primary[donor][barcode] = row["lineage_id"]
    if sorted(all_barcodes) != [f"Z{index:02d}" for index in range(1, 13)]:
        raise QueryATACMaterializationError("membership donor roster differs")
    return dict(primary), dict(all_barcodes)


def process_donor(
    donor_id: str,
    fragment_path: Path,
    expected_size: int,
    primary: Mapping[str, str],
    all_barcodes: set[str],
    windows: Mapping[str, WindowIndex],
    reference_lengths: Mapping[str, int],
) -> dict[str, Any]:
    import numpy as np

    if fragment_path.stat().st_size != expected_size:
        raise QueryATACMaterializationError(f"fragment size differs: {donor_id}")
    profiles = {role: np.zeros((4, 16000, 20), dtype=np.uint32) for role in ROLES}
    counts = {role: np.zeros((4, 16000), dtype=np.uint32) for role in ROLES}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    seen_primary: set[str] = set()
    prior_contig: str | None = None
    prior_start = -1
    closed: set[str] = set()
    selected_rows = 0
    with gzip.open(fragment_path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise QueryATACMaterializationError(f"fragment width differs: {donor_id}:{line_number}")
            contig, start_text, end_text, barcode, multiplicity_text = fields[:5]
            try:
                start, end, multiplicity = int(start_text), int(end_text), int(multiplicity_text)
            except ValueError as error:
                raise QueryATACMaterializationError("fragment numeric field differs") from error
            if (
                contig not in reference_lengths or start < 0 or end <= start
                or end > reference_lengths[contig] or multiplicity < 1
            ):
                raise QueryATACMaterializationError("fragment geometry differs")
            if contig != prior_contig:
                if contig in closed:
                    raise QueryATACMaterializationError("fragment contig reappears")
                if prior_contig is not None:
                    closed.add(prior_contig)
                prior_contig, prior_start = contig, -1
            if start < prior_start:
                raise QueryATACMaterializationError("fragments are unsorted")
            prior_start = start
            lineage = primary.get(barcode)
            if lineage is None:
                continue
            seen_primary.add(barcode)
            selected_rows += 1
            lineage_offset = lineage_index[lineage]
            for position in (start, end - 1):
                for role in ROLES:
                    for window_offset in windows_at(windows[role], contig, position):
                        window_start = windows[role].rows[window_offset][1]
                        bin_offset = (position - window_start) // 50
                        if not 0 <= bin_offset < 20:
                            raise QueryATACMaterializationError("profile bin differs")
                        profiles[role][lineage_offset, window_offset, bin_offset] += 1
                        counts[role][lineage_offset, window_offset] += 1
    if set(primary) != seen_primary:
        raise QueryATACMaterializationError(f"primary barcodes absent from fragments: {donor_id}")
    for role in ROLES:
        if not np.array_equal(profiles[role].sum(axis=2, dtype=np.uint64), counts[role]):
            raise QueryATACMaterializationError("profile/count conservation differs")
    return {
        "donor_id": donor_id,
        "profiles": profiles,
        "counts": counts,
        "selected_fragment_rows": selected_rows,
        "primary_barcodes": len(primary),
        "all_membership_barcodes": len(all_barcodes),
        "gzip_crc_verified_to_eof": True,
    }


def validate_config(root: Path, config: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    if config.get("schema_version") != SCHEMA or config.get("dataset_id") != "gse281367":
        raise QueryATACMaterializationError("query materialization identity differs")
    for field in (
        "condition_values_read", "phenotype_values_read", "evaluator_outcome_artifact_read",
        "metrics_calculated", "raw_multiplicity_exported", "genome_wide_library_total_exported",
        "missing_as_zero",
    ):
        if config.get(field) is not False:
            raise QueryATACMaterializationError(f"query firewall opened: {field}")
    axis_spec = config["axis_authority"]
    axis_root = (root / axis_spec["path"]).resolve(strict=True)
    axis = verify_frozen_tree(axis_root)
    if (
        _digest(axis_root / "ARTIFACTS.json") != axis_spec["artifacts_sha256"]
        or axis.get("metadata", {}).get("artifact_class") != axis_spec["artifact_class"]
        or _digest(axis_root / "label_free_cell_membership.tsv.gz") != axis_spec["membership_sha256"]
        or _digest(axis_root / "windows.tsv") != axis_spec["windows_sha256"]
    ):
        raise QueryATACMaterializationError("axis authority differs")
    source_spec = config["source_authority"]
    source_root = (root / source_spec["path"]).resolve(strict=True)
    source_manifest = json.loads((source_root / "ARTIFACTS.json").read_text(encoding="utf-8"))
    inventory = source_root / source_spec["inventory_path"]
    if (
        _digest(source_root / "ARTIFACTS.json") != source_spec["artifacts_sha256"]
        or source_manifest.get("metadata", {}).get("artifact_class") != source_spec["artifact_class"]
        or _digest(inventory) != source_spec["inventory_sha256"]
    ):
        raise QueryATACMaterializationError("source authority differs")
    reference = Path(config["reference"]["fai_path"]).resolve(strict=True)
    if _digest(reference) != config["reference"]["fai_sha256"]:
        raise QueryATACMaterializationError("reference FAI differs")
    implementation = config["implementation"]
    for prefix in ("builder", "test"):
        path = (root / implementation[f"{prefix}_path"]).resolve(strict=True)
        if _digest(path) != implementation[f"{prefix}_sha256"]:
            raise QueryATACMaterializationError(f"implementation differs: {prefix}")
    if config.get("firewall", {}).get("scored_target_contigs_exported_in_query_values") is not False:
        raise QueryATACMaterializationError("scored-target firewall differs")
    return axis_root, inventory, reference


def _read_source_inventory(path: Path, fragment_root: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("path", "size_bytes", "sha256"):
            raise QueryATACMaterializationError("source inventory fields differ")
        for row in reader:
            path_value = Path(row["path"])
            try:
                relative = path_value.relative_to(fragment_root)
            except ValueError:
                continue
            if relative.parts[-2:] != ("outs", "fragments.tsv.gz") or len(relative.parts) != 3:
                continue
            donor = relative.parts[0]
            result[donor] = row
    if sorted(result) != [f"Z{index:02d}" for index in range(1, 13)]:
        raise QueryATACMaterializationError("fragment inventory donor roster differs")
    return result


def _verify_fragment_source(record: Mapping[str, str]) -> str:
    path = Path(record["path"]).resolve(strict=True)
    if path.is_symlink() or path.stat().st_size != int(record["size_bytes"]):
        raise QueryATACMaterializationError(f"fragment path or size differs: {path}")
    observed = _digest(path)
    if observed != record["sha256"]:
        raise QueryATACMaterializationError(f"fragment hash differs: {path}")
    return observed


def _write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_role_output(
    output: Path,
    role: str,
    windows: WindowIndex,
    profiles: Any,
    counts: Any,
    qc_rows: Sequence[Mapping[str, Any]],
    authorities: Mapping[str, Any],
) -> None:
    import numpy as np

    target_role = "test" if role == "valid" else "valid"
    output.mkdir(parents=True, mode=0o750)
    np.save(output / "query_unit_fragment_profile.uint32.npy", profiles)
    np.save(output / "query_unit_fragment_total.uint32.npy", counts)
    _write_tsv(
        output / "query_windows.tsv",
        ("context_index", "window_id", "contig", "start", "end"),
        (
            {"context_index": index, "window_id": row[3], "contig": row[0], "start": row[1], "end": row[2]}
            for index, row in enumerate(windows.rows)
        ),
    )
    _write_tsv(
        output / "query_axis.tsv",
        ("donor_index", "donor_id", "lineage_index", "lineage_id"),
        (
            {
                "donor_index": donor_index,
                "donor_id": f"Z{donor_index + 1:02d}",
                "lineage_index": lineage_index,
                "lineage_id": lineage,
            }
            for donor_index in range(12)
            for lineage_index, lineage in enumerate(LINEAGES)
        ),
    )
    _write_tsv(
        output / "processing_qc.tsv",
        ("donor_id", "selected_fragment_rows", "primary_barcodes", "all_membership_barcodes", "gzip_crc_verified_to_eof"),
        qc_rows,
    )
    write_json_exclusive(
        output / "query_contract.json",
        {
            "schema_version": "masld-bench-gse281367-label-free-query-atac-v1",
            "context_role": role,
            "scored_target_role": target_role,
            "donors": 12,
            "lineages": list(LINEAGES),
            "windows": 16000,
            "profile_bins": 20,
            "scored_target_contigs_exported_in_query_values": False,
            "query_transform_fit": False,
            "condition_values_read": False,
            "phenotype_values_read": False,
            "evaluator_outcome_artifact_read": False,
            "metrics_calculated": False,
            "rna_state": "structurally_missing",
            "model_input_files": [
                "query_unit_fragment_profile.uint32.npy",
                "query_unit_fragment_total.uint32.npy",
                "query_windows.tsv",
                "query_axis.tsv",
                "query_contract.json",
                "input_authorities.json",
                "ARTIFACTS.json",
            ],
            "audit_only_files": ["processing_qc.tsv"],
            "champion_claim_allowed": False,
        },
    )
    write_json_exclusive(output / "input_authorities.json", authorities)
    freeze_tree(
        output,
        {
            "artifact_class": "gse281367_label_free_query_atac",
            "dataset_id": "gse281367",
            "context_role": role,
            "scored_target_role": target_role,
            "condition_values_read": False,
            "evaluator_outcome_artifact_read": False,
            "metrics_calculated": False,
            "champion_claim_allowed": False,
            "status": "passed",
        },
    )


def materialize(root: Path, config_path: Path, valid_output: Path, test_output: Path, workers: int) -> None:
    import numpy as np

    if valid_output.exists() or test_output.exists() or not 1 <= workers <= 6:
        raise QueryATACMaterializationError("output or worker contract differs")
    root = root.resolve(strict=True)
    config = _load_json(config_path.resolve(strict=True))
    axis_root, inventory_path, reference_path = validate_config(root, config)
    windows = read_windows(axis_root / "windows.tsv")
    primary, all_barcodes = read_membership(axis_root / "label_free_cell_membership.tsv.gz")
    reference = read_reference_lengths(reference_path)
    fragment_root = Path(config["fragment_root"]).resolve(strict=True)
    inventory = _read_source_inventory(inventory_path, fragment_root)
    donors = sorted(inventory)
    with ThreadPoolExecutor(max_workers=min(workers, 4)) as executor:
        verified_hashes = list(executor.map(_verify_fragment_source, [inventory[donor] for donor in donors]))
    if len(set(verified_hashes)) != len(donors):
        raise QueryATACMaterializationError("fragment source hashes are duplicated")
    profiles = {role: np.zeros((12, 4, 16000, 20), dtype=np.uint32) for role in ROLES}
    counts = {role: np.zeros((12, 4, 16000), dtype=np.uint32) for role in ROLES}
    qc: dict[str, dict[str, Any]] = {}
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {}
        for donor in donors:
            record = inventory[donor]
            path = Path(record["path"]).resolve(strict=True)
            future = executor.submit(
                process_donor,
                donor,
                path,
                int(record["size_bytes"]),
                primary[donor],
                all_barcodes[donor],
                windows,
                reference,
            )
            futures[future] = donor
        for future in as_completed(futures):
            result = future.result()
            donor = result["donor_id"]
            offset = donors.index(donor)
            for role in ROLES:
                profiles[role][offset] = result["profiles"][role]
                counts[role][offset] = result["counts"][role]
            qc[donor] = {key: value for key, value in result.items() if key not in {"profiles", "counts"}}
    authorities = {
        "campaign_id": config["campaign_id"],
        "campaign_sha256": _digest(config_path.resolve(strict=True)),
        "axis_artifacts_sha256": config["axis_authority"]["artifacts_sha256"],
        "source_artifacts_sha256": config["source_authority"]["artifacts_sha256"],
        "source_inventory_sha256": config["source_authority"]["inventory_sha256"],
        "reference_fai_sha256": config["reference"]["fai_sha256"],
        "builder_sha256": config["implementation"]["builder_sha256"],
        "fragment_sha256_by_donor": dict(zip(donors, verified_hashes, strict=True)),
        "condition_or_phenotype_authority_included": False,
        "evaluator_outcome_authority_included": False,
    }
    write_role_output(
        valid_output,
        "valid",
        windows["valid"],
        profiles["valid"],
        counts["valid"],
        [qc[donor] for donor in donors],
        authorities,
    )
    write_role_output(
        test_output,
        "test",
        windows["test"],
        profiles["test"],
        counts["test"],
        [qc[donor] for donor in donors],
        authorities,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--valid-output", required=True, type=Path)
    parser.add_argument("--test-output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    arguments = parser.parse_args()
    materialize(arguments.root, arguments.config, arguments.valid_output, arguments.test_output, arguments.workers)


if __name__ == "__main__":
    main()
