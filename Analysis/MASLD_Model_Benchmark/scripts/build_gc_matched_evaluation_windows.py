#!/usr/bin/env python3
"""Add deterministic cCRE-free, GC-matched genomic evaluation windows."""

from __future__ import annotations

import argparse
from bisect import bisect_left
import csv
import gzip
from hashlib import sha256
import heapq
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


SEED = 20260824
INPUT_LENGTH = 2114
OUTPUT_LENGTH = 1000
GRID_STEP = 2500
GC_CALIPER = 0.04
PRIMARY_CONTIGS = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
CCRE_WINDOW_FIELDS = (
    "contig",
    "output_start",
    "output_end",
    "window_id",
    "genomic_fold",
    "window_class",
    "ccre_class",
    "ccre_id",
    "ccre_start",
    "ccre_end",
    "input_start",
    "input_end",
    "selection_hash",
)


class BackgroundWindowError(ValueError):
    """Raised when outcome-independent background matching does not meet its requirements."""


def stable_integer(*values: object, seed: int = SEED) -> int:
    payload = "\0".join((str(seed), *(str(value) for value in values)))
    return int.from_bytes(sha256(payload.encode("utf-8")).digest(), "big")


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise BackgroundWindowError("forbidden interval is invalid")
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return tuple((start, end) for start, end in merged)


def overlaps(
    merged: Sequence[tuple[int, int]], starts: Sequence[int], start: int, end: int
) -> bool:
    if start < 0 or end <= start or len(merged) != len(starts):
        raise BackgroundWindowError("overlap query is invalid")
    index = bisect_left(starts, end)
    return index > 0 and merged[index - 1][1] > start


def gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    if len(sequence) != INPUT_LENGTH or set(sequence).difference("ACGT"):
        raise BackgroundWindowError("evaluation sequence is not unambiguous A/C/G/T")
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def match_gc_targets(
    targets: Sequence[tuple[float, str]],
    candidates: Sequence[tuple[float, int]],
) -> dict[str, tuple[float, int]]:
    """Greedily match ascending targets to a unique nearest candidate."""

    ordered_targets = sorted(targets, key=lambda row: (row[0], row[1]))
    ordered_candidates = sorted(candidates, key=lambda row: (row[0], row[1]))
    if len(ordered_candidates) < len(ordered_targets):
        raise BackgroundWindowError("fewer background candidates than targets")
    lower: list[tuple[float, int]] = []
    candidate_index = 0
    matched: dict[str, tuple[float, int]] = {}
    for target_gc, target_id in ordered_targets:
        while (
            candidate_index < len(ordered_candidates)
            and ordered_candidates[candidate_index][0] <= target_gc
        ):
            candidate_gc, center = ordered_candidates[candidate_index]
            heapq.heappush(lower, (-candidate_gc, center))
            candidate_index += 1
        lower_difference = (
            target_gc + lower[0][0] if lower else math.inf
        )
        upper_difference = (
            ordered_candidates[candidate_index][0] - target_gc
            if candidate_index < len(ordered_candidates)
            else math.inf
        )
        if lower_difference <= upper_difference:
            negative_gc, center = heapq.heappop(lower)
            matched[target_id] = (-negative_gc, center)
        elif candidate_index < len(ordered_candidates):
            matched[target_id] = ordered_candidates[candidate_index]
            candidate_index += 1
        else:
            raise BackgroundWindowError("GC matcher exhausted its candidates")
    if len(matched) != len(targets) or len({value[1] for value in matched.values()}) != len(matched):
        raise BackgroundWindowError("GC matching is not one-to-one")
    return matched


def read_tsv(path: Path, expected_fields: Sequence[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(expected_fields):
            raise BackgroundWindowError(f"{path.name} fields differ")
        return [dict(row) for row in reader]


def write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: Sequence[float], fraction: float) -> float:
    if not values or not 0 <= fraction <= 1:
        raise BackgroundWindowError("percentile input is invalid")
    ordered = sorted(values)
    index = min(len(ordered) - 1, math.ceil(fraction * len(ordered)) - 1)
    return ordered[max(index, 0)]


def build_windows(
    *,
    split_contract: Path,
    ccre_bed: Path,
    fasta: Path,
    output: Path,
    diagnostic_artifacts_sha256: str,
    seed: int = SEED,
    grid_step: int = GRID_STEP,
    gc_caliper: float = GC_CALIPER,
) -> None:
    import pyfaidx

    split_contract = split_contract.resolve(strict=True)
    ccre_bed = ccre_bed.resolve(strict=True)
    fasta = fasta.resolve(strict=True)
    if output.exists() or grid_step < INPUT_LENGTH:
        raise BackgroundWindowError("invalid output or grid-step contract")
    windows = read_tsv(
        split_contract / "ccre_evaluation_windows.tsv", CCRE_WINDOW_FIELDS
    )
    if len(windows) != 80_000:
        raise BackgroundWindowError("cCRE evaluation window count differs")
    folds = read_tsv(
        split_contract / "genomic_folds.tsv",
        ("genomic_fold", "contigs", "total_bp"),
    )
    fold_by_contig = {}
    for row in folds:
        for contig in row["contigs"].split(","):
            fold_by_contig[contig] = int(row["genomic_fold"])
    if tuple(sorted(fold_by_contig, key=lambda value: PRIMARY_CONTIGS.index(value))) != PRIMARY_CONTIGS:
        raise BackgroundWindowError("genomic fold contig universe differs")

    forbidden_raw: dict[str, list[tuple[int, int]]] = {
        contig: [] for contig in PRIMARY_CONTIGS
    }
    with gzip.open(ccre_bed, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise BackgroundWindowError("ENCODE cCRE BED row differs")
            contig = fields[0]
            if contig in forbidden_raw:
                forbidden_raw[contig].append((int(fields[1]), int(fields[2])))
    forbidden = {
        contig: merge_intervals(values) for contig, values in forbidden_raw.items()
    }
    forbidden_starts = {
        contig: tuple(start for start, _end in values)
        for contig, values in forbidden.items()
    }

    reference = pyfaidx.Fasta(
        str(fasta),
        as_raw=True,
        sequence_always_upper=True,
        rebuild=False,
        read_ahead=1_000_000,
    )
    try:
        chrom_sizes = {contig: len(reference[contig]) for contig in PRIMARY_CONTIGS}
        targets_by_fold: dict[int, list[tuple[float, str]]] = {
            fold: [] for fold in range(5)
        }
        target_by_id = {}
        for row in sorted(
            windows, key=lambda value: (PRIMARY_CONTIGS.index(value["contig"]), int(value["input_start"]))
        ):
            contig = row["contig"]
            start, end = int(row["input_start"]), int(row["input_end"])
            sequence = str(reference[contig][start:end])
            fraction = gc_fraction(sequence)
            targets_by_fold[fold_by_contig[contig]].append(
                (fraction, row["window_id"])
            )
            target_by_id[row["window_id"]] = {**row, "gc_fraction": fraction}

        candidates_by_fold: dict[int, list[tuple[float, int]]] = {
            fold: [] for fold in range(5)
        }
        candidate_location: dict[int, tuple[str, int]] = {}
        half = INPUT_LENGTH // 2
        for contig in PRIMARY_CONTIGS:
            size = chrom_sizes[contig]
            offset = stable_integer("background_grid", contig, seed=seed) % grid_step
            first_center = half + offset
            candidates = []
            for center in range(first_center, size - half, grid_step):
                start = center - half
                end = start + INPUT_LENGTH
                if overlaps(
                    forbidden[contig], forbidden_starts[contig], start, end
                ):
                    continue
                sequence = str(reference[contig][start:end])
                try:
                    fraction = gc_fraction(sequence)
                except BackgroundWindowError:
                    continue
                token = PRIMARY_CONTIGS.index(contig) * 1_000_000_000 + center
                if token in candidate_location:
                    raise BackgroundWindowError("background token is duplicated")
                candidate_location[token] = (contig, center)
                candidates.append((fraction, token))
            candidates_by_fold[fold_by_contig[contig]].extend(candidates)

        matches = {}
        for fold in range(5):
            if len(candidates_by_fold[fold]) < len(targets_by_fold[fold]):
                raise BackgroundWindowError(
                    f"fold {fold} has fewer cCRE-free windows than targets"
                )
            values = match_gc_targets(
                targets_by_fold[fold], candidates_by_fold[fold]
            )
            for target_id, (background_gc, token) in values.items():
                contig, center = candidate_location[token]
                matches[target_id] = (contig, background_gc, center)
    finally:
        reference.close()

    background_rows = []
    attempted_differences = []
    retained_differences = []
    class_attempted = {value: 0 for value in sorted({row["ccre_class"] for row in target_by_id.values()})}
    class_retained = {value: 0 for value in class_attempted}
    for target_id in sorted(target_by_id):
        target = target_by_id[target_id]
        contig, background_gc, center = matches[target_id]
        if fold_by_contig[contig] != int(target["genomic_fold"]):
            raise BackgroundWindowError("background genomic fold differs from target")
        input_start = center - INPUT_LENGTH // 2
        input_end = input_start + INPUT_LENGTH
        output_start = center - OUTPUT_LENGTH // 2
        output_end = output_start + OUTPUT_LENGTH
        difference = abs(float(target["gc_fraction"]) - background_gc)
        attempted_differences.append(difference)
        class_attempted[target["ccre_class"]] += 1
        if difference > gc_caliper:
            continue
        retained_differences.append(difference)
        class_retained[target["ccre_class"]] += 1
        background_rows.append(
            {
                "contig": contig,
                "output_start": output_start,
                "output_end": output_end,
                "window_id": f"bg_{target_id}",
                "genomic_fold": target["genomic_fold"],
                "window_class": "ccre_free_background",
                "ccre_class": target["ccre_class"],
                "matched_window_id": target_id,
                "ccre_id": "",
                "input_start": input_start,
                "input_end": input_end,
                "gc_fraction": f"{background_gc:.8f}",
                "gc_absolute_difference": f"{difference:.8f}",
            }
        )
    summary = {
        "attempted_pairs": len(attempted_differences),
        "retained_pairs": len(retained_differences),
        "retained_fraction": len(retained_differences) / len(attempted_differences),
        "gc_caliper": gc_caliper,
        "attempted_mean_absolute_gc_difference": sum(attempted_differences)
        / len(attempted_differences),
        "attempted_p95_absolute_gc_difference": percentile(
            attempted_differences, 0.95
        ),
        "retained_mean_absolute_gc_difference": sum(retained_differences)
        / len(retained_differences),
        "retained_median_absolute_gc_difference": percentile(
            retained_differences, 0.5
        ),
        "retained_p95_absolute_gc_difference": percentile(
            retained_differences, 0.95
        ),
        "retained_maximum_absolute_gc_difference": max(retained_differences),
        "per_class_coverage": {
            value: {
                "attempted": class_attempted[value],
                "retained": class_retained[value],
                "fraction": class_retained[value] / class_attempted[value],
            }
            for value in class_attempted
        },
    }
    if (
        summary["retained_fraction"] < 0.75
        or summary["retained_maximum_absolute_gc_difference"] > gc_caliper
    ):
        raise BackgroundWindowError(f"GC matching quality gate failed: {summary}")

    output.mkdir(mode=0o750)
    fields = (
        "contig",
        "output_start",
        "output_end",
        "window_id",
        "genomic_fold",
        "window_class",
        "ccre_class",
        "matched_window_id",
        "ccre_id",
        "input_start",
        "input_end",
        "gc_fraction",
        "gc_absolute_difference",
    )
    write_tsv(output / "background_windows.tsv", fields, background_rows)
    combined_rows = []
    matched_target_ids = {row["matched_window_id"] for row in background_rows}
    for window_id in sorted(target_by_id):
        target = target_by_id[window_id]
        matched_window_id = (
            f"bg_{window_id}" if window_id in matched_target_ids else ""
        )
        combined_rows.append(
            {
                "contig": target["contig"],
                "output_start": target["output_start"],
                "output_end": target["output_end"],
                "window_id": window_id,
                "genomic_fold": target["genomic_fold"],
                "window_class": "encode_ccre",
                "ccre_class": target["ccre_class"],
                "matched_window_id": matched_window_id,
                "ccre_id": target["ccre_id"],
                "input_start": target["input_start"],
                "input_end": target["input_end"],
                "gc_fraction": f"{float(target['gc_fraction']):.8f}",
                "gc_absolute_difference": "0.00000000",
            }
        )
    combined_rows.extend(background_rows)
    combined_rows.sort(
        key=lambda row: (
            int(row["genomic_fold"]),
            PRIMARY_CONTIGS.index(str(row["contig"])),
            int(row["input_start"]),
            str(row["window_id"]),
        )
    )
    write_tsv(output / "evaluation_windows.tsv", fields, combined_rows)
    (output / "gc_match_summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    contract = {
        "schema_version": "masld-bench-genomic-evaluation-windows-v2",
        "source_positive_windows": "ENCODE_v4_ENCFF420VPZ",
        "positive_windows": 80_000,
        "background_windows": len(background_rows),
        "background_source": "deterministic_whole_genome_grid",
        "background_grid_step": grid_step,
        "background_excludes_all_encode_ccres": True,
        "background_inputs_nonoverlapping": True,
        "matching": "same_genomic_fold_nearest_input_sequence_GC_without_replacement",
        "background_role": "secondary_common_support_only",
        "all_positive_ccre_windows_remain_primary": True,
        "full_positive_set_background_AUPRC_prohibited": True,
        "gc_caliper": gc_caliper,
        "input_length": INPUT_LENGTH,
        "output_length": OUTPUT_LENGTH,
        "seed": seed,
        "outcome_independent": True,
        "held_atac_used": False,
        "gc_quality_gate": {
            "retained_fraction_min": 0.75,
            "retained_maximum_absolute_difference": gc_caliper,
        },
        "observed_gc_summary": summary,
        "design_diagnostic_artifacts_sha256": diagnostic_artifacts_sha256,
    }
    (output / "contract.json").write_text(
        json.dumps(contract, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-contract", type=Path, required=True)
    parser.add_argument("--ccre-bed", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostic-artifacts-sha256", required=True)
    arguments = parser.parse_args()
    build_windows(
        split_contract=arguments.split_contract,
        ccre_bed=arguments.ccre_bed,
        fasta=arguments.fasta,
        output=arguments.output,
        diagnostic_artifacts_sha256=arguments.diagnostic_artifacts_sha256,
    )
    print(arguments.output.resolve(strict=True))


if __name__ == "__main__":
    main()
