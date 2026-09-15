#!/usr/bin/env python3
"""Build evaluator-only fixed-window ATAC outcomes from raw fragments.

The primary unit is one deduplicated fragment row.  The fifth fragment column
is retained only in a secondary multiplicity-weighted total because its PCR
semantics differ across sources.  GSE281367 Cell Ranger ATAC fragments are
already Tn5-adjusted and use start/end-1 cut sites.  GSE244832 was generated
directly from Bowtie2 template coordinates, so its cut sites require +4/-5.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage")
MIN_CELLS = 50
PROFILE_BINS = 20
EXPECTED_WINDOW_WIDTH = 1000
EXPECTED_WINDOWS_PER_ROLE = 16_000
MEMBERSHIP_ARTIFACT_CLASS = "atac_transport_cell_membership"
SOURCE_AUDIT_ARTIFACT_CLASS = "atac_transport_source_admission_audit"
CELLRANGER_FRAGMENT_FORMAT_URL = (
    "https://www.10xgenomics.com/support/software/cell-ranger-atac/"
    "latest/analysis/fragments-file"
)


class ATACOutcomeError(RuntimeError):
    """Raised when a raw-fragment outcome does not meet the frozen requirements."""


@dataclass(frozen=True)
class WindowIndex:
    rows: tuple[tuple[str, int, int, str, str], ...]
    by_contig: Mapping[str, tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]]


@dataclass
class DonorResult:
    dataset_id: str
    donor_id: str
    unit_profile: Any
    unit_total: Any
    multiplicity_total: Any
    genome_unit_fragments: Any
    genome_multiplicity_fragments: Any
    qc: dict[str, Any]


def _sha256_file(path: Path) -> str:
    from hashlib import sha256

    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_artifact(
    root: Path,
    *,
    expected_sha256: str,
    artifact_class: str,
    cohort_ids: Sequence[str],
) -> dict[str, Any]:
    manifest_path = root / "ARTIFACTS.json"
    if _sha256_file(manifest_path) != expected_sha256:
        raise ATACOutcomeError(f"ARTIFACTS SHA-256 differs: {root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    metadata = manifest.get("metadata", {})
    if metadata.get("artifact_class") != artifact_class:
        raise ATACOutcomeError(f"artifact class differs: {root}")
    if metadata.get("cohort_family_ids") != list(cohort_ids):
        raise ATACOutcomeError(f"cohort family IDs differ: {root}")
    return manifest


def read_reference_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ATACOutcomeError(f"malformed FAI row {line_number}")
            contig, length_text = fields[:2]
            length = int(length_text)
            if contig in lengths or length <= 0:
                raise ATACOutcomeError(f"invalid FAI contig at row {line_number}")
            lengths[contig] = length
    if not all(f"chr{value}" in lengths for value in [*range(1, 23), "X", "Y"]):
        raise ATACOutcomeError("reference lacks a primary GRCh38 contig")
    return lengths


def read_fold_roles(path: Path) -> dict[str, str]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if set(payload) != {"train", "valid", "test"}:
        raise ATACOutcomeError("genomic fold roles differ")
    role_by_contig: dict[str, str] = {}
    for role, contigs in payload.items():
        for contig in contigs:
            if contig in role_by_contig:
                raise ATACOutcomeError(f"contig occurs in two genomic roles: {contig}")
            role_by_contig[contig] = role
    return role_by_contig


def read_windows(
    path: Path,
    *,
    role_by_contig: Mapping[str, str],
    reference_lengths: Mapping[str, int],
    expected_width: int = EXPECTED_WINDOW_WIDTH,
    expected_per_role: int | None = EXPECTED_WINDOWS_PER_ROLE,
) -> WindowIndex:
    rows: list[tuple[str, int, int, str, str]] = []
    seen_ids: set[str] = set()
    role_counts: Counter[str] = Counter()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 10:
                raise ATACOutcomeError(f"region row {line_number} is not narrowPeak-like")
            contig, start_text, end_text, name = fields[:4]
            prefix = name.split("|", 1)[0]
            if prefix not in {"fixed_valid", "fixed_test"}:
                continue
            role = prefix.removeprefix("fixed_")
            start, end = int(start_text), int(end_text)
            if role_by_contig.get(contig) != role:
                raise ATACOutcomeError(f"window/contig role differs at row {line_number}")
            if (
                contig not in reference_lengths
                or start < 0
                or end - start != expected_width
                or end > reference_lengths[contig]
            ):
                raise ATACOutcomeError(f"window geometry differs at row {line_number}")
            pieces = name.split("|")
            if len(pieces) != 3 or not pieces[1] or not pieces[2]:
                raise ATACOutcomeError(f"window identifier differs at row {line_number}")
            window_id = pieces[1]
            if window_id in seen_ids:
                raise ATACOutcomeError(f"duplicate window ID: {window_id}")
            seen_ids.add(window_id)
            role_counts[role] += 1
            rows.append((contig, start, end, window_id, role))
    if expected_per_role is not None and role_counts != {
        "test": expected_per_role,
        "valid": expected_per_role,
    }:
        raise ATACOutcomeError(f"fixed-window role counts differ: {dict(role_counts)}")
    if not rows:
        raise ATACOutcomeError("no fixed validation/test windows")
    indexed: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for index, (contig, start, end, _window_id, _role) in enumerate(rows):
        indexed[contig].append((start, end, index))
    by_contig: dict[str, tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]] = {}
    for contig, entries in indexed.items():
        entries.sort()
        starts = tuple(entry[0] for entry in entries)
        ends = tuple(entry[1] for entry in entries)
        indexes = tuple(entry[2] for entry in entries)
        if any(left >= right for left, right in zip(starts, starts[1:])):
            raise ATACOutcomeError(f"window starts are not unique on {contig}")
        if any(left >= right for left, right in zip(ends, ends[1:])):
            raise ATACOutcomeError(f"window ends are not increasing on {contig}")
        by_contig[contig] = (starts, ends, indexes)
    return WindowIndex(tuple(rows), by_contig)


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


def tn5_positions(dataset_id: str, start: int, end: int) -> tuple[int, int]:
    """Return assay-native cut sites from a 0-based half-open fragment."""

    if dataset_id == "gse281367":
        return start, end - 1
    if dataset_id == "gse244832":
        if end - start < 10:
            raise ATACOutcomeError("unadjusted custom fragment is shorter than 10 bp")
        return start + 4, end - 5
    raise ATACOutcomeError(f"unsupported fragment coordinate contract: {dataset_id}")


def read_membership(
    path: Path,
) -> tuple[
    dict[tuple[str, str], dict[str, str]],
    dict[tuple[str, str], set[str]],
    dict[tuple[str, str], dict[str, int]],
    dict[tuple[str, str], dict[str, str]],
]:
    primary: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
    all_barcodes: dict[tuple[str, str], set[str]] = defaultdict(set)
    cell_counts: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {lineage: 0 for lineage in LINEAGES}
    )
    attributes: dict[tuple[str, str], dict[str, str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "dataset_id",
            "donor_id",
            "raw_barcode",
            "lineage_id",
            "analysis_role",
            "condition",
            "outer_fold",
        }
        if not required <= set(reader.fieldnames or ()):
            raise ATACOutcomeError("membership fields differ")
        for row in reader:
            key = (row["dataset_id"], row["donor_id"])
            barcode = row["raw_barcode"]
            if barcode in all_barcodes[key]:
                raise ATACOutcomeError(f"duplicate donor barcode: {key}/{barcode}")
            all_barcodes[key].add(barcode)
            observed = {
                "condition": row["condition"],
                "outer_fold": row["outer_fold"],
            }
            if key in attributes and attributes[key] != observed:
                raise ATACOutcomeError(f"donor attributes differ: {key}")
            attributes[key] = observed
            if row["analysis_role"] == "primary":
                lineage = row["lineage_id"]
                if lineage not in LINEAGES:
                    raise ATACOutcomeError(f"unregistered primary lineage: {lineage}")
                primary[key][barcode] = lineage
                cell_counts[key][lineage] += 1
    expected = {"gse244832": 18, "gse281367": 12}
    observed = Counter(dataset for dataset, _donor in attributes)
    if observed != expected:
        raise ATACOutcomeError(f"membership donor counts differ: {dict(observed)}")
    return dict(primary), dict(all_barcodes), dict(cell_counts), attributes


def read_source_inventory(path: Path) -> dict[Path, dict[str, str]]:
    result: dict[Path, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or ()) != {"path", "size_bytes", "sha256"}:
            raise ATACOutcomeError("source inventory fields differ")
        for row in reader:
            source = Path(row["path"]).resolve(strict=True)
            if source in result:
                raise ATACOutcomeError(f"duplicate source inventory path: {source}")
            result[source] = row
    return result


def _fragment_path(dataset_id: str, donor_id: str, root: Path) -> Path:
    if dataset_id == "gse281367":
        return (root / donor_id / "outs" / "fragments.tsv.gz").resolve(strict=True)
    if dataset_id == "gse244832":
        return (root / f"{donor_id}_fragments.tsv.gz").resolve(strict=True)
    raise ATACOutcomeError(f"unsupported fragment dataset: {dataset_id}")


def process_donor(
    *,
    dataset_id: str,
    donor_id: str,
    fragment_path: Path,
    expected_size: int,
    primary_by_barcode: Mapping[str, str],
    all_barcodes: set[str],
    windows: WindowIndex,
    reference_lengths: Mapping[str, int],
    profile_bins: int = PROFILE_BINS,
) -> DonorResult:
    import numpy as np

    if fragment_path.stat().st_size != expected_size:
        raise ATACOutcomeError(f"fragment size drift: {fragment_path}")
    n_windows = len(windows.rows)
    profile = np.zeros((len(LINEAGES), n_windows, profile_bins), dtype=np.uint32)
    unit_total = np.zeros((len(LINEAGES), n_windows), dtype=np.uint32)
    multiplicity_total = np.zeros((len(LINEAGES), n_windows), dtype=np.uint64)
    genome_unit = np.zeros(len(LINEAGES), dtype=np.uint64)
    genome_multiplicity = np.zeros(len(LINEAGES), dtype=np.uint64)
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}
    seen_primary: set[str] = set()
    seen_all: set[str] = set()
    fragment_rows = 0
    selected_rows = 0
    ignored_unregistered_rows = 0
    comment_rows = 0
    prior_contig: str | None = None
    prior_start = -1
    closed_contigs: set[str] = set()
    with gzip.open(fragment_path, "rt", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                comment_rows += 1
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ATACOutcomeError(f"fragment width differs: {fragment_path}:{line_number}")
            contig, start_text, end_text, barcode, multiplicity_text = fields[:5]
            try:
                start, end, multiplicity = (
                    int(start_text),
                    int(end_text),
                    int(multiplicity_text),
                )
            except ValueError as error:
                raise ATACOutcomeError(
                    f"non-integer fragment field: {fragment_path}:{line_number}"
                ) from error
            if (
                contig not in reference_lengths
                or start < 0
                or end <= start
                or end > reference_lengths[contig]
                or multiplicity < 1
            ):
                raise ATACOutcomeError(
                    f"invalid fragment geometry/count: {fragment_path}:{line_number}"
                )
            if contig != prior_contig:
                if contig in closed_contigs:
                    raise ATACOutcomeError(f"fragment contig reappears: {fragment_path}:{contig}")
                if prior_contig is not None:
                    closed_contigs.add(prior_contig)
                prior_contig, prior_start = contig, -1
            if start < prior_start:
                raise ATACOutcomeError(f"fragments are unsorted: {fragment_path}:{line_number}")
            prior_start = start
            fragment_rows += 1
            if barcode in all_barcodes:
                seen_all.add(barcode)
            lineage = primary_by_barcode.get(barcode)
            if lineage is None:
                if barcode not in all_barcodes:
                    ignored_unregistered_rows += 1
                continue
            seen_primary.add(barcode)
            selected_rows += 1
            lineage_offset = lineage_index[lineage]
            genome_unit[lineage_offset] += 1
            genome_multiplicity[lineage_offset] += multiplicity
            for position in tn5_positions(dataset_id, start, end):
                for window_offset in windows_at(windows, contig, position):
                    window_start = windows.rows[window_offset][1]
                    window_end = windows.rows[window_offset][2]
                    profile_offset = (
                        (position - window_start) * profile_bins // (window_end - window_start)
                    )
                    if not 0 <= profile_offset < profile_bins:
                        raise ATACOutcomeError("internal profile-bin calculation failed")
                    profile[lineage_offset, window_offset, profile_offset] += 1
                    unit_total[lineage_offset, window_offset] += 1
                    multiplicity_total[lineage_offset, window_offset] += multiplicity
    unseen_primary = sorted(set(primary_by_barcode) - seen_primary)
    if unseen_primary:
        raise ATACOutcomeError(
            f"primary membership barcodes absent from fragments: {dataset_id}/{donor_id} "
            f"n={len(unseen_primary)} first={unseen_primary[:3]}"
        )
    if not np.array_equal(profile.sum(axis=2, dtype=np.uint64), unit_total):
        raise ATACOutcomeError(f"profile/count conservation failed: {dataset_id}/{donor_id}")
    return DonorResult(
        dataset_id=dataset_id,
        donor_id=donor_id,
        unit_profile=profile,
        unit_total=unit_total,
        multiplicity_total=multiplicity_total,
        genome_unit_fragments=genome_unit,
        genome_multiplicity_fragments=genome_multiplicity,
        qc={
            "dataset_id": dataset_id,
            "donor_id": donor_id,
            "fragment_path": str(fragment_path),
            "fragment_rows": fragment_rows,
            "comment_rows": comment_rows,
            "selected_primary_fragment_rows": selected_rows,
            "ignored_unregistered_fragment_rows": ignored_unregistered_rows,
            "membership_barcodes": len(all_barcodes),
            "membership_barcodes_seen": len(seen_all),
            "primary_barcodes": len(primary_by_barcode),
            "primary_barcodes_seen": len(seen_primary),
            "gzip_crc_verified_to_eof": True,
            "coordinates_0_based_half_open": True,
        },
    )


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _write_windows(path: Path, windows: WindowIndex) -> None:
    rows = [
        {
            "window_index": index,
            "window_id": window_id,
            "role": role,
            "contig": contig,
            "start": start,
            "end": end,
        }
        for index, (contig, start, end, window_id, role) in enumerate(windows.rows)
    ]
    _write_tsv(
        path,
        ("window_index", "window_id", "role", "contig", "start", "end"),
        rows,
    )


def build_outputs(args: argparse.Namespace) -> dict[str, Any]:
    import numpy as np

    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    evaluator = output / "evaluator_outcomes"
    evaluator.mkdir()
    verify_artifact(
        args.source_audit,
        expected_sha256=args.source_artifacts_sha256,
        artifact_class=SOURCE_AUDIT_ARTIFACT_CLASS,
        cohort_ids=("gse244832", "gse281367"),
    )
    verify_artifact(
        args.membership,
        expected_sha256=args.membership_artifacts_sha256,
        artifact_class=MEMBERSHIP_ARTIFACT_CLASS,
        cohort_ids=("gse244832", "gse281367"),
    )
    if _sha256_file(args.regions) != args.regions_sha256:
        raise ATACOutcomeError("fixed-region source SHA-256 differs")
    if _sha256_file(args.fold_roles) != args.fold_roles_sha256:
        raise ATACOutcomeError("genomic-fold source SHA-256 differs")
    if _sha256_file(args.reference_fai) != args.reference_fai_sha256:
        raise ATACOutcomeError("reference FAI SHA-256 differs")
    if (
        _sha256_file(args.gse281_fragment_builder)
        != args.gse281_fragment_builder_sha256
        or _sha256_file(args.gse244_fragment_builder)
        != args.gse244_fragment_builder_sha256
    ):
        raise ATACOutcomeError("fragment-builder provenance SHA-256 differs")
    reference_lengths = read_reference_lengths(args.reference_fai)
    windows = read_windows(
        args.regions,
        role_by_contig=read_fold_roles(args.fold_roles),
        reference_lengths=reference_lengths,
    )
    _write_windows(evaluator / "windows.tsv", windows)
    primary, all_barcodes, cell_counts, attributes = read_membership(
        args.membership / "cell_membership.tsv.gz"
    )
    inventory = read_source_inventory(args.source_audit / "source_inventory.tsv")
    roots = {
        "gse281367": args.gse281_fragments,
        "gse244832": args.gse244_fragments,
    }
    keys = sorted(attributes)
    donor_order = {
        dataset: [donor for observed, donor in keys if observed == dataset]
        for dataset in ("gse244832", "gse281367")
    }
    datasets: dict[str, dict[str, Any]] = {}
    for dataset, donors in donor_order.items():
        dataset_root = evaluator / dataset
        dataset_root.mkdir()
        shape = (len(donors), len(LINEAGES), len(windows.rows))
        datasets[dataset] = {
            "donors": donors,
            "profile": np.lib.format.open_memmap(
                dataset_root / "unit_fragment_profile.uint32.npy",
                mode="w+",
                dtype=np.uint32,
                shape=(*shape, args.profile_bins),
            ),
            "unit": np.lib.format.open_memmap(
                dataset_root / "unit_fragment_total.uint32.npy",
                mode="w+",
                dtype=np.uint32,
                shape=shape,
            ),
            "multiplicity": np.lib.format.open_memmap(
                dataset_root / "raw_multiplicity_total.uint64.npy",
                mode="w+",
                dtype=np.uint64,
                shape=shape,
            ),
            "genome_unit": np.lib.format.open_memmap(
                dataset_root / "genome_unit_fragments.uint64.npy",
                mode="w+",
                dtype=np.uint64,
                shape=(len(donors), len(LINEAGES)),
            ),
            "genome_multiplicity": np.lib.format.open_memmap(
                dataset_root / "genome_multiplicity_fragments.uint64.npy",
                mode="w+",
                dtype=np.uint64,
                shape=(len(donors), len(LINEAGES)),
            ),
        }
    tasks: list[dict[str, Any]] = []
    for dataset, donor in keys:
        fragment_path = _fragment_path(dataset, donor, roots[dataset])
        if fragment_path not in inventory:
            raise ATACOutcomeError(f"fragment absent from frozen inventory: {fragment_path}")
        tasks.append(
            {
                "dataset_id": dataset,
                "donor_id": donor,
                "fragment_path": fragment_path,
                "expected_size": int(inventory[fragment_path]["size_bytes"]),
                "primary_by_barcode": primary[(dataset, donor)],
                "all_barcodes": all_barcodes[(dataset, donor)],
                "windows": windows,
                "reference_lengths": reference_lengths,
                "profile_bins": args.profile_bins,
            }
        )
    qc_rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_key = {
            executor.submit(process_donor, **task): (task["dataset_id"], task["donor_id"])
            for task in tasks
        }
        for future in as_completed(future_to_key):
            dataset, donor = future_to_key[future]
            result = future.result()
            offset = donor_order[dataset].index(donor)
            destination = datasets[dataset]
            destination["profile"][offset] = result.unit_profile
            destination["unit"][offset] = result.unit_total
            destination["multiplicity"][offset] = result.multiplicity_total
            destination["genome_unit"][offset] = result.genome_unit_fragments
            destination["genome_multiplicity"][offset] = (
                result.genome_multiplicity_fragments
            )
            qc_rows.append(result.qc)
    axis_rows: list[dict[str, Any]] = []
    eligible_rows = 0
    for dataset, payload in datasets.items():
        for value in payload.values():
            if hasattr(value, "flush"):
                value.flush()
        dataset_root = evaluator / dataset
        eligible = np.zeros((len(payload["donors"]), len(LINEAGES)), dtype=np.bool_)
        for donor_offset, donor in enumerate(payload["donors"]):
            key = (dataset, donor)
            for lineage_offset, lineage in enumerate(LINEAGES):
                cells = cell_counts[key][lineage]
                is_eligible = cells >= MIN_CELLS
                eligible[donor_offset, lineage_offset] = is_eligible
                eligible_rows += int(is_eligible)
                axis_rows.append(
                    {
                        "dataset_id": dataset,
                        "donor_index": donor_offset,
                        "donor_id": donor,
                        "lineage_index": lineage_offset,
                        "lineage_id": lineage,
                        "condition": attributes[key]["condition"],
                        "outer_fold": attributes[key]["outer_fold"],
                        "cells": cells,
                        "eligible_min_50_cells": str(is_eligible).lower(),
                    }
                )
        np.save(dataset_root / "eligible_min_50_cells.bool.npy", eligible)
    _write_tsv(
        evaluator / "donor_lineage_axis.tsv",
        (
            "dataset_id",
            "donor_index",
            "donor_id",
            "lineage_index",
            "lineage_id",
            "condition",
            "outer_fold",
            "cells",
            "eligible_min_50_cells",
        ),
        sorted(axis_rows, key=lambda row: (row["dataset_id"], row["donor_index"], row["lineage_index"])),
    )
    qc_fields = (
        "dataset_id",
        "donor_id",
        "fragment_path",
        "fragment_rows",
        "comment_rows",
        "selected_primary_fragment_rows",
        "ignored_unregistered_fragment_rows",
        "membership_barcodes",
        "membership_barcodes_seen",
        "primary_barcodes",
        "primary_barcodes_seen",
        "gzip_crc_verified_to_eof",
        "coordinates_0_based_half_open",
    )
    _write_tsv(
        output / "fragment_qc.tsv",
        qc_fields,
        sorted(qc_rows, key=lambda row: (row["dataset_id"], row["donor_id"])),
    )
    summary = {
        "schema_version": "masld-bench-atac-transport-outcomes-v1",
        "artifact_role": "evaluator_only_development_outcomes",
        "cohort_family_ids": ["gse244832", "gse281367"],
        "donors": {dataset: len(donors) for dataset, donors in donor_order.items()},
        "windows": len(windows.rows),
        "window_roles": dict(sorted(Counter(row[4] for row in windows.rows).items())),
        "window_width": EXPECTED_WINDOW_WIDTH,
        "profile_bins": args.profile_bins,
        "profile_bin_width": EXPECTED_WINDOW_WIDTH // args.profile_bins,
        "primary_count_unit": "deduplicated_fragment_row_tn5_insertions",
        "secondary_count_unit": "raw_fifth_column_multiplicity_tn5_insertions",
        "coordinate_contract_by_dataset": {
            "gse281367": {
                "source": "Cell_Ranger_ATAC_2.1.0_already_Tn5_adjusted",
                "cut_sites": ["fragment_start", "fragment_end_minus_1"],
                "local_builder_sha256": args.gse281_fragment_builder_sha256,
                "official_format": CELLRANGER_FRAGMENT_FORMAT_URL,
            },
            "gse244832": {
                "source": "custom_Bowtie2_raw_template_geometry",
                "cut_sites": ["fragment_start_plus_4", "fragment_end_minus_5"],
                "local_builder_sha256": args.gse244_fragment_builder_sha256,
            },
        },
        "lineages": list(LINEAGES),
        "eligibility_rule": f"cells >= {MIN_CELLS}",
        "eligible_donor_lineages": eligible_rows,
        "fragment_files_read_to_eof": len(qc_rows),
        "source_hash_authority": str(args.source_audit / "source_inventory.tsv"),
        "model_inputs_may_read_evaluator_outcomes": False,
        "dataset_activated_for_development_transport": True,
        "same_cell_pairing_claimed": False,
    }
    (output / "outcome_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-audit", type=Path, required=True)
    parser.add_argument("--source-artifacts-sha256", required=True)
    parser.add_argument("--membership", type=Path, required=True)
    parser.add_argument("--membership-artifacts-sha256", required=True)
    parser.add_argument("--regions", type=Path, required=True)
    parser.add_argument("--regions-sha256", required=True)
    parser.add_argument("--fold-roles", type=Path, required=True)
    parser.add_argument("--fold-roles-sha256", required=True)
    parser.add_argument("--reference-fai", type=Path, required=True)
    parser.add_argument("--reference-fai-sha256", required=True)
    parser.add_argument("--gse281-fragment-builder", type=Path, required=True)
    parser.add_argument("--gse281-fragment-builder-sha256", required=True)
    parser.add_argument("--gse244-fragment-builder", type=Path, required=True)
    parser.add_argument("--gse244-fragment-builder-sha256", required=True)
    parser.add_argument("--gse281-fragments", type=Path, required=True)
    parser.add_argument("--gse244-fragments", type=Path, required=True)
    parser.add_argument("--profile-bins", type=int, default=PROFILE_BINS)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1 or args.profile_bins < 1:
        parser.error("workers and profile bins must be positive")
    if EXPECTED_WINDOW_WIDTH % args.profile_bins:
        parser.error("profile bins must divide the 1-kb fixed window")
    return args


def main() -> int:
    summary = build_outputs(parse_args())
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
