#!/usr/bin/env python3
"""Build raw-count, deposited-sample pseudobulk without source-derived outcomes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import tarfile
from typing import BinaryIO, Iterable, TextIO


SAMPLES = {
    "GSM7660623": ("SAMN36706212", "DIFF_HRG_P_Seq2"),
    "GSM7660624": ("SAMN36706211", "P_Seq2"),
    "GSM7660625": ("SAMN36706210", "PSG_1"),
    "GSM7660626": ("SAMN36706209", "PSG_2"),
    "GSM7660627": ("SAMN36706208", "PSG_3"),
}


class MaterializationError(RuntimeError):
    """Raised when source topology differs from the frozen contract."""


def file_digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def bytes_digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def safe_name(value: str) -> None:
    pure = PurePosixPath(value)
    if pure.is_absolute() or ".." in pure.parts:
        raise MaterializationError(f"unsafe archive member: {value}")


def find_outer_member(archive: tarfile.TarFile, sample: str, suffix: str) -> tarfile.TarInfo:
    matches = [
        member
        for member in archive.getmembers()
        if member.isfile() and member.name.startswith(sample) and member.name.endswith(suffix)
    ]
    if len(matches) != 1:
        raise MaterializationError(
            f"expected one {sample} member ending {suffix}, observed {len(matches)}"
        )
    safe_name(matches[0].name)
    return matches[0]


def read_small_nested_members(
    outer: tarfile.TarFile,
    outer_member: tarfile.TarInfo,
    suffixes: Iterable[str],
) -> dict[str, bytes]:
    wanted = set(suffixes)
    found: dict[str, bytes] = {}
    stream = outer.extractfile(outer_member)
    if stream is None:
        raise MaterializationError(f"cannot open {outer_member.name}")
    with tarfile.open(fileobj=stream, mode="r|gz") as nested:
        for member in nested:
            safe_name(member.name)
            if not member.isfile():
                continue
            match = next((suffix for suffix in wanted if member.name.endswith(suffix)), None)
            if match is None:
                continue
            handle = nested.extractfile(member)
            if handle is None:
                raise MaterializationError(f"cannot read {member.name}")
            found[match] = handle.read()
            if set(found) == wanted:
                break
    if set(found) != wanted:
        raise MaterializationError(f"missing nested members: {sorted(wanted - set(found))}")
    return found


def parse_guide_reference(payload: bytes) -> tuple[list[dict[str, str]], dict[str, str]]:
    text = io.StringIO(payload.decode("utf-8-sig"), newline="")
    rows = list(csv.DictReader(text))
    expected = {
        "id", "name", "read", "pattern", "sequence", "feature_type",
        "target_gene_id", "target_gene_name",
    }
    if not rows or set(rows[0]) != expected:
        raise MaterializationError("guide reference schema differs")
    guide_to_target: dict[str, str] = {}
    for row in rows:
        for key in (row["id"], row["name"]):
            prior = guide_to_target.setdefault(key, row["target_gene_name"])
            if prior != row["target_gene_name"]:
                raise MaterializationError(f"guide maps to multiple targets: {key}")
    return rows, guide_to_target


def classify_call(features: list[str], guide_to_target: dict[str, str]) -> tuple[str, str]:
    try:
        targets = [guide_to_target[feature] for feature in features]
    except KeyError as error:
        raise MaterializationError(f"unmapped guide call: {error.args[0]}") from error
    noncontrols = sorted({target for target in targets if target != "Non-Targeting"})
    has_control = any(target == "Non-Targeting" for target in targets)
    if not noncontrols:
        return "control", "Non-Targeting"
    if has_control:
        return "mixed_control_target", "|".join(noncontrols)
    if len(noncontrols) == 1:
        return "single_gene", noncontrols[0]
    return "multi_gene", "|".join(noncontrols)


def parse_calls(
    payload: bytes,
    guide_to_target: dict[str, str],
) -> tuple[dict[str, dict[str, str | int]], Counter[tuple[str, str, str]]]:
    text = io.StringIO(payload.decode("utf-8-sig"), newline="")
    calls: dict[str, dict[str, str | int]] = {}
    families: Counter[tuple[str, str, str]] = Counter()
    for row in csv.DictReader(text):
        barcode = row["cell_barcode"]
        features = row["feature_call"].split("|")
        if int(row["num_features"]) != len(features):
            raise MaterializationError(f"guide multiplicity differs for {barcode}")
        if barcode in calls:
            raise MaterializationError(f"duplicate guide call for {barcode}")
        classification, target = classify_call(features, guide_to_target)
        family = "|".join(sorted(features))
        calls[barcode] = {
            "classification": classification,
            "target_gene": target,
            "guide_family": family,
            "num_guides": len(features),
        }
        families[(classification, target, family)] += 1
    if not calls:
        raise MaterializationError("guide-call table is empty")
    return calls, families


def decoded_gzip_text(stream: BinaryIO) -> TextIO:
    return io.TextIOWrapper(gzip.GzipFile(fileobj=stream), encoding="utf-8-sig", newline="")


def read_lines(stream: BinaryIO) -> list[str]:
    with decoded_gzip_text(stream) as text:
        return [line.rstrip("\r\n") for line in text]


def process_matrix_archive(
    outer: tarfile.TarFile,
    outer_member: tarfile.TarInfo,
    calls: dict[str, dict[str, str | int]],
) -> tuple[list[tuple[str, str, str]], list[str], dict[str, list[int]], dict[str, int], int]:
    barcodes: list[str] | None = None
    features: list[tuple[str, str, str]] | None = None
    groups_by_column: list[str | None] | None = None
    sums: dict[str, list[int]] | None = None
    classification_counts: Counter[str] = Counter()
    observed_nnz = 0
    stream = outer.extractfile(outer_member)
    if stream is None:
        raise MaterializationError(f"cannot open {outer_member.name}")
    with tarfile.open(fileobj=stream, mode="r|gz") as nested:
        for member in nested:
            safe_name(member.name)
            if not member.isfile():
                continue
            handle = nested.extractfile(member)
            if handle is None:
                raise MaterializationError(f"cannot read {member.name}")
            if member.name.endswith("/barcodes.tsv.gz"):
                barcodes = read_lines(handle)
                if len(barcodes) != len(set(barcodes)):
                    raise MaterializationError("duplicate expression barcode")
            elif member.name.endswith("/features.tsv.gz"):
                rows = [line.split("\t") for line in read_lines(handle)]
                if any(len(row) != 3 for row in rows):
                    raise MaterializationError("expression feature schema differs")
                features = [(row[0], row[1], row[2]) for row in rows]
            elif member.name.endswith("/matrix.mtx.gz"):
                if barcodes is None or features is None:
                    raise MaterializationError("matrix precedes axes")
                groups_by_column = []
                for barcode in barcodes:
                    call = calls.get(barcode)
                    if call is None:
                        classification_counts["unassigned"] += 1
                        groups_by_column.append(None)
                        continue
                    classification = str(call["classification"])
                    classification_counts[classification] += 1
                    if classification == "control":
                        groups_by_column.append("control")
                    elif classification == "single_gene":
                        groups_by_column.append(str(call["target_gene"]))
                    else:
                        groups_by_column.append(None)
                groups = sorted({group for group in groups_by_column if group is not None})
                if "control" not in groups:
                    raise MaterializationError("no non-targeting control cells")
                sums = {group: [0] * len(features) for group in groups}
                with decoded_gzip_text(handle) as matrix:
                    shape: tuple[int, int, int] | None = None
                    for line in matrix:
                        if not line.strip() or line.startswith("%"):
                            continue
                        fields = line.split()
                        if shape is None:
                            if len(fields) != 3:
                                raise MaterializationError("matrix shape differs")
                            shape = tuple(int(value) for value in fields)
                            if shape[:2] != (len(features), len(barcodes)):
                                raise MaterializationError("matrix axes differ")
                            continue
                        if len(fields) != 3:
                            raise MaterializationError("matrix entry differs")
                        row_index, column_index, value = (int(field) for field in fields)
                        observed_nnz += 1
                        group = groups_by_column[column_index - 1]
                        if group is not None:
                            sums[group][row_index - 1] += value
                    if shape is None or observed_nnz != shape[2]:
                        raise MaterializationError("matrix nonzero count differs")
                break
    if barcodes is None or features is None or groups_by_column is None or sums is None:
        raise MaterializationError("matrix archive is incomplete")
    return features, barcodes, sums, dict(classification_counts), observed_nnz


def write_outputs(
    output: Path,
    sample: str,
    biosample: str,
    deposited_series_label: str,
    source: Path,
    guide_rows: list[dict[str, str]],
    guide_reference_payload: bytes,
    calls: dict[str, dict[str, str | int]],
    families: Counter[tuple[str, str, str]],
    features: list[tuple[str, str, str]],
    barcodes: list[str],
    sums: dict[str, list[int]],
    classification_counts: dict[str, int],
    observed_nnz: int,
) -> None:
    output.mkdir(parents=True, exist_ok=False)
    guide_path = output / "guide_crosswalk.csv"
    with guide_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(guide_rows[0]))
        writer.writeheader()
        writer.writerows(guide_rows)
    with (output / "guide_family_summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["classification", "target_gene", "guide_family", "num_cells"])
        for (classification, target, family), count in sorted(families.items()):
            writer.writerow([classification, target, family, count])
    with gzip.open(output / "cell_assignments.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["cell_barcode", "classification", "target_gene", "guide_family", "num_guides"])
        for barcode in barcodes:
            call = calls.get(barcode)
            if call is None:
                writer.writerow([barcode, "unassigned", "", "", 0])
            else:
                writer.writerow([
                    barcode,
                    call["classification"],
                    call["target_gene"],
                    call["guide_family"],
                    call["num_guides"],
                ])
    groups = sorted(sums, key=lambda value: (value != "control", value))
    with gzip.open(output / "raw_target_pseudobulk.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_index", "gene_id", "gene_name", "feature_type", *groups])
        for index, (gene_id, gene_name, feature_type) in enumerate(features):
            writer.writerow([
                index + 1,
                gene_id,
                gene_name,
                feature_type,
                *(sums[group][index] for group in groups),
            ])
    axis_payload = "".join(
        f"{index + 1}\t{gene_id}\t{gene_name}\t{feature_type}\n"
        for index, (gene_id, gene_name, feature_type) in enumerate(features)
    ).encode("utf-8")
    summary = {
        "schema_version": "masld-bench-gse238219-raw-pseudobulk-v1",
        "sample": sample,
        "biosample": biosample,
        "deposited_series_label": deposited_series_label,
        "source_path": source.as_posix(),
        "source_sha256": file_digest(source),
        "gene_axis_sha256": bytes_digest(axis_payload),
        "guide_reference_source_sha256": bytes_digest(guide_reference_payload),
        "guide_crosswalk_sha256": file_digest(guide_path),
        "num_expression_features": len(features),
        "num_expression_cells": len(barcodes),
        "num_matrix_nonzero_entries": observed_nnz,
        "num_guide_called_cells": len(calls),
        "cell_assignment_counts": classification_counts,
        "pseudobulk_columns": groups,
        "primary_single_gene_rule": "one unique non-control target gene and no non-targeting guide; multiple guides to that same gene are retained and their exact guide family is recorded",
        "excluded_from_single_gene_pseudobulk": ["unassigned", "mixed_control_target", "multi_gene"],
        "normalization_performed": False,
        "feature_selection_performed": False,
        "source_perturbation_effect_tables_read": False,
        "biological_unit_status": "deposited_sample_label_not_independence_proof",
        "cells_guides_files_as_replicates": False,
        "gse313774_accessed": False,
    }
    (output / "materialization_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def materialize(source: Path, sample: str, output: Path) -> None:
    if sample not in SAMPLES or not source.is_file() or output.exists():
        raise MaterializationError("materialization request differs")
    biosample, deposited_series_label = SAMPLES[sample]
    with tarfile.open(source, mode="r:") as outer:
        guide_member = find_outer_member(outer, sample, "_crispr_analysis.tar.gz")
        matrix_member = find_outer_member(outer, sample, "_filtered_feature_bc_matrix.tar.gz")
        small = read_small_nested_members(
            outer,
            guide_member,
            ("/feature_reference.csv", "/protospacer_calls_per_cell.csv"),
        )
        guide_rows, guide_to_target = parse_guide_reference(small["/feature_reference.csv"])
        calls, families = parse_calls(
            small["/protospacer_calls_per_cell.csv"], guide_to_target
        )
        features, barcodes, sums, classification_counts, observed_nnz = process_matrix_archive(
            outer, matrix_member, calls
        )
    write_outputs(
        output,
        sample,
        biosample,
        deposited_series_label,
        source,
        guide_rows,
        small["/feature_reference.csv"],
        calls,
        families,
        features,
        barcodes,
        sums,
        classification_counts,
        observed_nnz,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--sample", choices=tuple(SAMPLES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    materialize(args.source, args.sample, args.output)


if __name__ == "__main__":
    main()
